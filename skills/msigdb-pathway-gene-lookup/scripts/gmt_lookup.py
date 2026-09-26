#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MSigDB GMT 通路 <-> 基因 双向快速检索工具（纯标准库，无第三方依赖）

用法概览：
    python gmt_lookup.py guide "人肝癌 bulk rna-seq"    # 向导：按组织/实验推荐该用哪几个集合
    python gmt_lookup.py list                          # 列出所有集合与统计
    python gmt_lookup.py find "TNF signaling"          # 正向：通路名 -> 基因
    python gmt_lookup.py show HALLMARK_APOPTOSIS       # 正向：精确集合名 -> 全部基因
    python gmt_lookup.py gene TP53,EGFR                # 反向：基因 -> 参与的通路
    python gmt_lookup.py build --force                 # 强制重建索引

数据目录解析优先级： --data-dir > 环境变量 MSIGDB_DIR > <skill>/config.json > 内置默认
索引文件默认落在数据目录下 .msigdb_index.sqlite3；GMT 文件指纹（文件名+大小+mtime）变化时自动重建。
"""

from __future__ import annotations

import argparse
import array
import csv
import hashlib
import json
import os
import re
import sqlite3
import sys
import time

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUNDLED_DATA_DIR = os.path.join(SKILL_DIR, "data")
DEFAULT_DB_NAME = ".msigdb_index.sqlite3"
BATCH = 20000

# ---------------------------------------------------------------- 集合元信息

COLL_LABEL = {
    "H": "Hallmark 标志性基因集",
    "C1": "位置型基因集（染色体臂/细胞遗传带）",
    "C2": "专家审编基因集（KEGG / REACTOME / WP / BIOCARTA / PID 等）",
    "C3": "调控靶标基因集（转录因子靶标 TFT / miRNA 靶标 MIR）",
    "C4": "计算型基因集（共表达模块 / 邻域）",
    "C5": "本体基因集（GO BP / CC / MF、HPO）",
    "C6": "致癌签名（扰动实验）",
    "C7": "免疫签名（免疫相关扰动实验）",
    "C8": "细胞类型签名（单细胞标记）",
    "C9": "计算扰动签名（DepMap CRISPR / CCLE）",
}

# set_id 前缀 -> 来源数据库
SOURCE_PREFIX = [
    ("HALLMARK_", "HALLMARK"),
    ("KEGG_", "KEGG"),
    ("REACTOME_", "REACTOME"),
    ("WP_", "WIKIPATHWAYS"),
    ("BIOCARTA_", "BIOCARTA"),
    ("PID_", "PID"),
    ("GOBP_", "GO:BP"),
    ("GOCC_", "GO:CC"),
    ("GOMF_", "GO:MF"),
    ("HP_", "HPO"),
    ("MIR_", "MIR"),
    ("TFT_", "TFT"),
]

# 经典“信号通路”来源在排序中加权；扰动/签名/靶标类集合降权（生物学家检索时通常不想要它们）
SOURCE_BONUS = {
    "HALLMARK": 60,
    "KEGG": 70,
    "REACTOME": 65,
    "WIKIPATHWAYS": 60,
    "BIOCARTA": 55,
    "PID": 55,
    "GO:BP": 20,
    "GO:CC": 10,
    "GO:MF": 10,
    "HPO": 0,
    "MIR": -40,
    "TFT": -40,
}
COLL_PENALTY = {"C1": -60, "C3": -40, "C4": -60, "C6": -50, "C7": -70, "C8": -70, "C9": -70}

CURATED_COND = "collection='C2' AND source_db IN ('KEGG','REACTOME','WIKIPATHWAYS','BIOCARTA','PID')"

# --collection 别名 -> SQL 条件
ALIAS_COND = {
    "H": "collection='H'",
    "HALLMARK": "collection='H'",
    "C1": "collection='C1'",
    "C2": "collection='C2'",
    "C3": "collection='C3'",
    "C4": "collection='C4'",
    "C5": "collection='C5'",
    "C6": "collection='C6'",
    "C7": "collection='C7'",
    "C8": "collection='C8'",
    "C9": "collection='C9'",
    "KEGG": "source_db='KEGG'",
    "REACTOME": "source_db='REACTOME'",
    "WP": "source_db='WIKIPATHWAYS'",
    "WIKIPATHWAYS": "source_db='WIKIPATHWAYS'",
    "BIOCARTA": "source_db='BIOCARTA'",
    "PID": "source_db='PID'",
    "GO": "source_db IN ('GO:BP','GO:CC','GO:MF')",
    "GOBP": "source_db='GO:BP'",
    "GOCC": "source_db='GO:CC'",
    "GOMF": "source_db='GO:MF'",
    "HP": "source_db='HPO'",
    "HPO": "source_db='HPO'",
    "MIR": "source_db='MIR'",
    "MIRNA": "source_db='MIR'",
    "TFT": "source_db='TFT'",
    "TPATHWAY": CURATED_COND,
    "PATHWAY": "(collection='H' OR (" + CURATED_COND + "))",
    "CURATED": "(collection='H' OR (" + CURATED_COND + "))",
    # 语义别名：按「研究方面」选集合时用，与 C1..C9 完全等价
    "IMMUNE": "collection='C7'",
    "IMMUNESIGDB": "collection='C7'",
    "CELLTYPE": "collection='C8'",
    "CELLTYPES": "collection='C8'",
    "MARKER": "collection='C8'",
    "PERTURB": "collection='C9'",
    "PERTURBATION": "collection='C9'",
    "DEPMAP": "collection='C9'",
    "ONCOGENIC": "collection='C6'",
    "CANCER": "collection IN ('C6','C9')",
    "REGULATORY": "collection='C3'",
    "TARGETS": "collection='C3'",
    "POSITIONAL": "collection='C1'",
    "LOCATION": "collection='C1'",
    "ONTOLOGY": "collection='C5'",
    "COMPUTATIONAL": "collection='C4'",
    "MODULE": "collection='C4'",
    "CGP": "collection='C2' AND source_db=''",
    "CP": CURATED_COND,
    "CLASSIC": CURATED_COND,
}


# ---------------------------------------------------------------- 工具函数


def eprint(*a):
    print(*a, file=sys.stderr)


def norm_key(s: str) -> str:
    """大写 + 非字母数字转下划线：'PI3K-Akt signaling pathway' -> 'PI3K_AKT_SIGNALING_PATHWAY'"""
    return re.sub(r"[^A-Z0-9]+", "_", s.upper()).strip("_")


def squash(s: str) -> str:
    """只留字母数字，规避 NF-kB -> NFKB 这类分隔符差异"""
    return re.sub(r"[^A-Z0-9]", "", s.upper())


def source_db_of(set_id: str) -> str:
    for pre, name in SOURCE_PREFIX:
        if set_id.startswith(pre):
            return name
    return ""


def collection_of(fname: str) -> str:
    """h.all.v2026.1... -> H ; c2.all.v... -> C2 ; 兜底用文件名主干大写"""
    m = re.match(r"^([A-Za-z]+\d*)\.", os.path.basename(fname))
    return m.group(1).upper() if m else os.path.splitext(os.path.basename(fname))[0].upper()


def version_of(fname: str) -> str:
    m = re.search(r"v(\d+(?:\.\d+)*)", os.path.basename(fname))
    return ("v" + m.group(1)) if m else ""


def fmt_int(n) -> str:
    return f"{int(n):,}"


# ---------------------------------------------------------------- 配置/路径


def resolve_config() -> dict:
    cfg_path = os.path.join(SKILL_DIR, "config.json")
    if os.path.exists(cfg_path):
        try:
            with open(cfg_path, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception as exc:  # 配置损坏不应致命
            eprint(f"[warn] 读取 config.json 失败，忽略：{exc}")
    return {}


def resolve_data_dir(cli_value) -> str:
    """GMT 目录解析优先级：

    1. 命令行 --data-dir
    2. 环境变量 MSIGDB_DIR
    3. config.json 里的 data_dir
    4. skill 自带的 <SKILL_DIR>/data/   ← 兜底，保证装完即用

    前三者属于「用户明确指定」，路径不存在就报错（避免掩盖拼写错误）；
    都没指定时才用自带数据。所以本 skill 不依赖任何外部目录。
    """
    cfg = resolve_config()
    explicit = cli_value or os.environ.get("MSIGDB_DIR") or cfg.get("data_dir")
    if explicit:
        p = os.path.expanduser(str(explicit)).replace("\\", "/")
        if not os.path.isdir(p):
            raise SystemExit(f"指定的 GMT 数据目录不存在：{p}")
        return p
    if os.path.isdir(BUNDLED_DATA_DIR):
        return BUNDLED_DATA_DIR
    raise SystemExit(
        f"没找到 GMT 数据目录。请确认 {BUNDLED_DATA_DIR} 存在，"
        "或用 --data-dir / 环境变量 MSIGDB_DIR / config.json 的 data_dir 指定自己的目录。"
    )


def find_gmt_files(data_dir: str) -> list:
    out = []
    for name in sorted(os.listdir(data_dir)):
        if name.lower().endswith(".gmt"):
            p = os.path.join(data_dir, name)
            if os.path.isfile(p):
                out.append(p)
    return out


def fingerprint(files: list) -> str:
    rec = {}
    for p in files:
        st = os.stat(p)
        rec[os.path.basename(p)] = [st.st_size, int(st.st_mtime)]
    return json.dumps(rec, sort_keys=True)


# ---------------------------------------------------------------- 索引构建


def parse_gmt(path: str):
    """逐行产出 (set_id, url, genes[])，自动跳过空行/列数不足的行并集内去重"""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        first = True
        for raw in fh:
            line = raw.rstrip("\r\n")
            if first:
                line = line.lstrip("\ufeff")
                first = False
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            set_id = parts[0].strip()
            if not set_id:
                continue
            url = parts[1].strip()
            seen = set()
            genes = []
            for g in parts[2:]:
                g = g.strip().upper()
                if g and g not in seen:
                    seen.add(g)
                    genes.append(g)
            if genes:
                yield set_id, url, genes


def build_index(data_dir: str, db_path: str, quiet=False) -> dict:
    files = find_gmt_files(data_dir)
    if not files:
        raise SystemExit(f"在 {data_dir} 下没找到任何 .gmt 文件")

    if os.path.exists(db_path):
        try:
            os.remove(db_path)
        except OSError as exc:
            raise SystemExit(f"无法覆盖旧索引 {db_path}：{exc}")

    t0 = time.time()
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    conn.executescript(
        """
        CREATE TABLE sets(
            set_pk      INTEGER PRIMARY KEY,
            set_id      TEXT NOT NULL,
            collection  TEXT NOT NULL,
            version     TEXT,
            source_file TEXT NOT NULL,
            source_db   TEXT,
            n_genes     INTEGER NOT NULL,
            url         TEXT,
            genes       TEXT NOT NULL,
            is_primary  INTEGER NOT NULL DEFAULT 1
        );
        CREATE INDEX idx_sets_id   ON sets(set_id);
        CREATE INDEX idx_sets_prim ON sets(is_primary, set_id);
        CREATE INDEX idx_sets_coll ON sets(collection);
        CREATE INDEX idx_sets_db   ON sets(source_db);
        CREATE TABLE gene_index(
            gene    TEXT PRIMARY KEY,
            set_pks BLOB NOT NULL
        ) WITHOUT ROWID;
        CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);
        """
    )

    set_rows = []
    gene_map = {}  # gene -> array('I') of set_pk
    pk = 0
    per_file = []
    seen_id = {}  # set_id -> 首次出现的文件名
    dups = {}  # set_id -> [文件名...]
    n_primary = 0

    for path in files:
        fname = os.path.basename(path)
        coll = collection_of(path)
        ver = version_of(path)
        n_file = 0
        for set_id, url, genes in parse_gmt(path):
            if set_id in seen_id:
                # 同名 set_id 出现在多个文件中：保留全部行，但只有首次出现的标记为主版本
                dups.setdefault(set_id, [seen_id[set_id]]).append(fname)
                is_primary = 0
            else:
                seen_id[set_id] = fname
                is_primary = 1
                n_primary += 1
            pk += 1
            n_file += 1
            set_rows.append(
                (
                    pk,
                    set_id,
                    coll,
                    ver,
                    fname,
                    source_db_of(set_id),
                    len(genes),
                    url,
                    "\t".join(genes),
                    is_primary,
                )
            )
            for g in genes:
                arr = gene_map.get(g)
                if arr is None:
                    gene_map[g] = array.array("I", [pk])
                else:
                    arr.append(pk)
        per_file.append((fname, coll, n_file))
        if not quiet:
            eprint(f"  解析 {fname:<34} 集合 {n_file:>6}")

    conn.executemany(
        "INSERT INTO sets(set_pk,set_id,collection,version,source_file,source_db,n_genes,url,genes,is_primary)"
        " VALUES(?,?,?,?,?,?,?,?,?,?)",
        set_rows,
    )

    gi_rows = [(g, arr.tobytes()) for g, arr in gene_map.items()]
    for i in range(0, len(gi_rows), BATCH):
        conn.executemany("INSERT INTO gene_index(gene,set_pks) VALUES(?,?)", gi_rows[i : i + BATCH])

    conn.executemany(
        "INSERT INTO meta(key,value) VALUES(?,?)",
        [
            ("fingerprint", fingerprint(files)),
            ("built_at", time.strftime("%Y-%m-%d %H:%M:%S")),
            ("data_dir", data_dir),
            ("n_files", str(len(files))),
            ("n_sets", str(n_primary)),
            ("n_rows", str(len(set_rows))),
            ("n_genes", str(len(gene_map))),
            ("n_entries", str(sum(r[6] for r in set_rows))),
            ("n_dups", str(len(dups))),
        ],
    )
    conn.commit()
    conn.execute("PRAGMA optimize")
    conn.commit()
    conn.close()

    return {
        "n_files": len(files),
        "n_sets": n_primary,
        "n_rows": len(set_rows),
        "n_genes": len(gene_map),
        "n_entries": sum(r[6] for r in set_rows),
        "seconds": round(time.time() - t0, 1),
        "db_path": db_path,
        "per_file": per_file,
        "dups": dups,
    }


def _needs_build(data_dir, db_path):
    if not os.path.exists(db_path):
        return True
    files = find_gmt_files(data_dir)
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
        row = conn.execute("SELECT value FROM meta WHERE key='fingerprint'").fetchone()
        conn.close()
    except sqlite3.DatabaseError:
        return True
    return (not row) or row[0] != fingerprint(files)


def _connect_ro(db_path):
    """以只读方式打开索引，并调优 I/O（mmap + 大页缓存）。"""
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        conn.execute("PRAGMA mmap_size=268435456")  # 256MB 内存映射
        conn.execute("PRAGMA cache_size=-65536")  # 64MB 页缓存
        conn.execute("PRAGMA temp_store=MEMORY")
        conn.execute("PRAGMA query_only=ON")
    except sqlite3.DatabaseError:
        pass
    return conn


def _writable(path: str) -> bool:
    try:
        os.makedirs(path, exist_ok=True)
        probe = os.path.join(path, ".write_probe")
        with open(probe, "w", encoding="utf-8") as fh:
            fh.write("")
        os.remove(probe)
        return True
    except OSError:
        return False


def _user_cache_dir() -> str:
    """跨平台用户缓存目录：Linux/macOS 用 XDG_CACHE_HOME，Windows 用 LOCALAPPDATA。"""
    if os.name == "nt":
        base = os.environ.get("LOCALAPPDATA") or os.path.expanduser("~")
        return os.path.join(base, "msigdb-gmt-lookup", "cache")
    base = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    return os.path.join(base, "msigdb-gmt-lookup")


def default_index_path(data_dir: str) -> str:
    """索引默认放数据目录；若数据目录只读（例如 skill 装在只读位置），
    退回用户缓存目录（Windows: %LOCALAPPDATA%/msigdb-gmt-lookup/cache；
    Linux/macOS: $XDG_CACHE_HOME/msigdb-gmt-lookup），保证仍然可用。"""
    cand = os.path.join(data_dir, DEFAULT_DB_NAME)
    if os.path.exists(cand) or _writable(data_dir):
        return cand
    cache = _user_cache_dir()
    os.makedirs(cache, exist_ok=True)
    tag = hashlib.md5(os.path.abspath(data_dir).encode("utf-8")).hexdigest()[:12]
    eprint(f"[info] 数据目录不可写，索引改放到 {cache}")
    return os.path.join(cache, f"index-{tag}.sqlite3")


def open_index(data_dir, cli_db=None, auto_build=True, quiet=False):
    db_path = (
        cli_db or os.environ.get("MSIGDB_INDEX") or default_index_path(data_dir)
    )
    if _needs_build(data_dir, db_path):
        if not auto_build:
            raise SystemExit(f"索引缺失或已过期：{db_path}（去掉 --no-build 以自动重建）")
        eprint("[info] 首次使用或 GMT 有变化，正在建立索引（约 10 秒，只需一次）…")
        stats = build_index(data_dir, db_path, quiet=quiet)
        eprint(
            f"[info] 索引完成：{fmt_int(stats['n_sets'])} 个基因集 / {fmt_int(stats['n_genes'])} 个基因 / "
            f"{fmt_int(stats['n_entries'])} 条基因条目，耗时 {stats['seconds']}s -> {db_path}"
        )
        warn_dups(stats)
    return _connect_ro(db_path), db_path


def warn_dups(stats):
    """同名 set_id 出现在多个 GMT 文件时的提示（默认只查主版本）。"""
    dups = stats.get("dups") or {}
    if not dups:
        return
    eprint(
        f"[warn] 有 {len(dups)} 个 set_id 出现在多个 GMT 文件中；索引已保留全部行，"
        "默认检索只返回首次出现的主版本。"
    )
    eprint("       要查某个非主版本：加 --keep-dups 检索，再用 --file 过滤来源文件。")
    for k in list(dups)[:5]:
        eprint(f"       例：{k}  <- {', '.join(dups[k])}")


# ---------------------------------------------------------------- 过滤条件


def build_where(collections, files=None, primary_only=False):
    """编译 --collection / --file / 主版本 过滤条件，返回 (sql, params)"""
    groups, params = [], []
    if primary_only:
        groups.append("is_primary=1")
    if collections:
        conds = []
        for raw in collections:
            key = str(raw).strip().upper()
            if not key:
                continue
            cond = ALIAS_COND.get(key)
            if cond is None:
                raise SystemExit(
                    f"未知的集合/来源别名：{raw}\n可用别名：{', '.join(sorted(ALIAS_COND))}"
                )
            conds.append(f"({cond})")
        if conds:
            groups.append("(" + " OR ".join(conds) + ")")
    if files:
        subs = []
        for f in files:
            subs.append("source_file LIKE ?")
            params.append(f"%{f}%")
        groups.append("(" + " OR ".join(subs) + ")")
    return (" AND ".join(groups), params)


# ---------------------------------------------------------------- 正向检索


def fetch_genes_map(conn, pks):
    """批量取多个 set_pk 的基因列表，避免逐个查询（分批规避 SQLite 变量上限）。"""
    out = {}
    pks = list(dict.fromkeys(pks))
    for i in range(0, len(pks), 500):
        chunk = pks[i : i + 500]
        for row in conn.execute(
            f"SELECT set_pk,genes FROM sets WHERE set_pk IN ({','.join('?'*len(chunk))})", chunk
        ):
            out[row["set_pk"]] = row["genes"].split("\t")
    return out


# 集合名里高频出现的英文功能词。模糊匹配（近似层）计算命中率时忽略它们，
# 否则 "CORRELATED_WITH_ABL2" 会因为所有 C9 集合都含 CORRELATED/WITH 而全量命中。
STOPWORDS = {
    "WITH", "WITHOUT", "AND", "OF", "TO", "BY", "VIA", "IN", "FOR", "THE", "FROM",
    "ON", "OR", "NOT", "AT", "AS", "IS", "ARE", "ITS", "A", "AN", "UP", "DN",
}


def _is_subseq(needle, hay):
    n = len(needle)
    if n == 0 or n > len(hay):
        return False
    for i in range(len(hay) - n + 1):
        if hay[i : i + n] == needle:
            return True
    return False


def rank_sets(query, rows):
    """分层打分：精确 > 前缀 > 词组 > 全词 > 包含 > 近似，再叠加来源权重。

    关键点：连字符/空格会被归一化成下划线，'NF-kB' -> ['NF','KB']，而被检索的集合
    里常写作 NFKB（单个 token）。所以同时用「归一化分词」和「压平整体」两套 token 去匹配，
    取较高分，避免 NF-kB / PI3K-Akt 这类写法因分隔符差异漏掉或降级。
    """
    qn = norm_key(query)
    qs = squash(query)
    qtok = [t for t in qn.split("_") if t]
    if not qtok:
        return []
    variants = [qtok]
    if qs and qs not in qtok:
        variants.append([qs])

    out = []
    for r in rows:
        sid = r["set_id"]
        stok = sid.split("_")
        ssq = squash(sid)
        best, how = 0, ""
        for vt in variants:
            if sid == qn:
                s, h = 1000, "精确"
            elif ssq == qs:
                s, h = 950, "精确(忽略分隔符)"
            elif sid.startswith(qn + "_"):
                s, h = 880, "前缀"
            elif qs and ssq.startswith(qs):
                s, h = 800, "前缀(忽略分隔符)"
            elif _is_subseq(vt, stok):
                s, h = 760, "词组"
            elif all(t in stok for t in vt):
                s, h = 700, "全词命中"
            elif qn in sid:
                s, h = 640, "包含"
            elif qs and qs in ssq:
                s, h = 600, "包含(忽略分隔符)"
            else:
                core = [t for t in vt if t not in STOPWORDS] or vt
                hit = sum(1 for t in core if any(t in st for st in stok))
                ratio = hit / len(core)
                if ratio >= 0.6:
                    s, h = 420 + int(ratio * 100), "近似"
                else:
                    continue
            if s > best:
                best, how = s, h
        if not best:
            continue
        best += SOURCE_BONUS.get(r["source_db"] or "", 0) + COLL_PENALTY.get(r["collection"], 0)
        out.append((best, len(stok), len(sid), how, r))
    out.sort(key=lambda x: (-x[0], x[1], x[2], x[4]["set_id"]))
    return out


def cmd_find(conn, args):
    """支持一次传入多个关键词，共用一次进程启动（本机 python 启动约 1.7s，多查几次更划算）。"""
    queries = args.query if isinstance(args.query, list) else [args.query]
    multi = len(queries) > 1
    all_res, all_pks = [], []
    for qi, q in enumerate(queries, 1):
        if multi and args.format == "table":
            print(f"\n{'#'*100}\n### 关键词 {qi}/{len(queries)}\n{'#'*100}")
        res = _find_one(conn, args, q, emit=(args.format == "table"))
        all_res.extend(res)
        all_pks.extend(r["set_pk"] for r in res)

    if args.csv and all_res:
        write_find_csv(args.csv, all_res)
        print(f"[已导出] {args.csv}（长格式：set_id / 集合 / 来源 / 基因数 / 匹配 / 基因）")
    if args.gmt and all_pks:
        write_gmt(args.gmt, conn, list(dict.fromkeys(all_pks)))
        print(f"[已导出] {args.gmt}（标准 GMT 格式，可直接 read.gmt 读入）")
    return all_res


def _find_one(conn, args, query, emit=True):
    where, params = build_where(args.collection, args.file, primary_only=not args.keep_dups)
    sql = "SELECT set_pk,set_id,collection,source_file,source_db,n_genes FROM sets"
    if where:
        sql += " WHERE " + where
    rows = conn.execute(sql, params).fetchall()
    if not rows:
        print("索引为空，或给定的过滤条件没有命中任何基因集。")
        if args.collection or args.file:
            print("提示：去掉 --collection / --file 可扩大到全部集合（但优先确认是不是该换个集合）。")
        return []

    ranked = rank_sets(query, rows)
    if args.exact:
        ranked = [x for x in ranked if x[3] in ("精确", "精确(忽略分隔符)")]
    limit = len(ranked) if args.all else max(1, args.limit)
    picked = ranked[:limit]

    suggestion = ""
    if not ranked:
        toks = [t for t in norm_key(query).split("_") if len(t) >= 3]
        bits = []
        if args.collection or args.file:
            bits.append(
                "先确认这个方向是不是不在当前集合里——对照 guide 里各集合的「适合的场景」，"
                "该换集合就换集合（例如找细胞类型标记要用 C8，找免疫状态要用 C7），"
                "别直接放弃限定"
            )
        if len(toks) > 1:
            bits.append("拆成更短的关键词单独检索：" + "、".join(f'find "{t}"' for t in toks[:4]))
        bits.append("换同义词（如 WNT / CATENIN、NFKB / RELA、PI3K / AKT）")
        suggestion = "未命中。" + "；".join(bits) + "。"
    elif len(ranked) <= 2 and (args.collection or args.file) and not args.exact:
        suggestion = "（命中较少：先想想是不是该换个集合，确要扩大范围再去掉 --collection / --file）"

    return render_forward(conn, args, query, picked, len(ranked), suggestion, emit=emit)


def render_forward(conn, args, query, picked, total_hits, suggestion="", emit=True):
    scope = ",".join(args.collection) if args.collection else "全部集合"
    if args.file:
        scope += " 文件≈" + ",".join(args.file)
    out = []

    def say(*a):
        if emit:
            print(*a)

    say(f'MSigDB 正向检索: "{query}"    范围={scope}')
    say(
        f"命中 {fmt_int(total_hits)} 个基因集，显示 {len(picked)} 个"
        + ("（--all 可全部显示）" if total_hits > len(picked) else "")
    )
    say("=" * 100)
    if suggestion:
        say(suggestion)

    preview = args.preview
    if len(picked) == 1 and preview == 12:
        preview = 0  # 单个结果默认给全量基因

    gene_maps = fetch_genes_map(conn, [r["set_pk"] for _s, _a, _b, _h, r in picked])
    for i, (score, _nt, _nl, how, r) in enumerate(picked, 1):
        genes = gene_maps.get(r["set_pk"], [])
        out.append(
            {
                "set_pk": r["set_pk"],
                "set_id": r["set_id"],
                "collection": r["collection"],
                "source_db": r["source_db"] or "",
                "source_file": r["source_file"],
                "n_genes": r["n_genes"],
                "match": how,
                "genes": genes,
            }
        )
        say(f"[{i}] {r['set_id']}")
        say(
            f"    集合={r['collection']}  来源={(r['source_db'] or '-')}  基因数={fmt_int(r['n_genes'])}  "
            f"匹配={how}  文件={r['source_file']}"
        )
        if preview and preview > 0 and len(genes) > preview:
            say(f"    基因(前{preview}): " + ", ".join(genes[:preview]) + f"  …(共 {len(genes)})")
        else:
            say("    基因: " + ", ".join(genes))
        say()
    return out


def write_find_csv(path, results):
    with open(path, "w", encoding="utf-8-sig", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["set_id", "collection", "source_db", "source_file", "n_genes", "match", "gene"])
        for r in results:
            for g in r["genes"]:
                w.writerow(
                    [r["set_id"], r["collection"], r["source_db"], r["source_file"], r["n_genes"], r["match"], g]
                )


def write_gmt(path, conn, pks):
    """导出标准 GMT：set_id<TAB>url<TAB>gene1<TAB>gene2..."""
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        for pk in pks:
            row = conn.execute("SELECT set_id,url,genes FROM sets WHERE set_pk=?", (pk,)).fetchone()
            fh.write("\t".join([row["set_id"], row["url"] or row["set_id"]]) + "\t" + row["genes"] + "\n")


def cmd_show(conn, args):
    # 指定了 --file 就按文件精确取版本；否则默认只取主版本（除非 --keep-dups）
    file_cond, fparams = build_where(None, args.file, primary_only=(not args.file and not args.keep_dups))
    target = norm_key(args.set_id)
    sql = "SELECT set_pk,set_id,collection,source_file,source_db,n_genes FROM sets WHERE set_id=?"
    params = [target]
    if file_cond:
        sql += " AND " + file_cond
        params += fparams
    rows = conn.execute(sql, params).fetchall()

    if not rows:  # 忽略分隔符的宽松匹配
        sql2 = "SELECT set_pk,set_id,collection,source_file,source_db,n_genes FROM sets WHERE REPLACE(set_id,'_','')=?"
        p2 = [squash(args.set_id)]
        if file_cond:
            sql2 += " AND " + file_cond
            p2 += fparams
        loose = conn.execute(sql2, p2).fetchall()
        if len(loose) == 1:
            rows = loose
        elif len(loose) > 1:
            print(f"'{args.set_id}' 模糊匹配到 {len(loose)} 个集合，请用完整名称重试：")
            for r in loose[:40]:
                print(f"  {r['set_id']}  ({r['collection']}, {r['n_genes']} genes, {r['source_file']})")
            return []
        else:
            print(f"没找到集合：{args.set_id}")
            print('提示：先用 find "关键词" 检索，再用返回的完整 set_id 调用 show。')
            return []

    if len(rows) > 1:
        print(f"'{args.set_id}' 在多个 GMT 文件中存在 {len(rows)} 个版本，请用 --file 指定：")
        for r in rows:
            print(f"  {r['set_id']}  ({r['collection']}, {r['n_genes']} genes, {r['source_file']})")
        return []

    row = rows[0]
    genes = conn.execute("SELECT genes FROM sets WHERE set_pk=?", (row["set_pk"],)).fetchone()[0].split("\t")
    if args.format == "genes":
        print("\n".join(genes))
    else:
        print(f"集合: {row['set_id']}")
        print(
            f"集合族={row['collection']}  来源={(row['source_db'] or '-')}  文件={row['source_file']}  "
            f"基因数={fmt_int(row['n_genes'])}"
        )
        print("=" * 100)
        print(", ".join(genes))
    if args.csv:
        with open(args.csv, "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["set_id", "collection", "source_db", "n_genes", "gene"])
            for g in genes:
                w.writerow([row["set_id"], row["collection"], row["source_db"] or "", row["n_genes"], g])
        print(f"[已导出] {args.csv}")
    if args.gmt:
        write_gmt(args.gmt, conn, [row["set_pk"]])
        print(f"[已导出] {args.gmt}")
    return [{"set_id": row["set_id"], "genes": genes}]


# ---------------------------------------------------------------- 反向检索


def cmd_gene(conn, args):
    genes = [g.strip().upper() for g in re.split(r"[,\s;]+", args.genes) if g.strip()]
    if not genes:
        raise SystemExit('请提供至少一个基因 symbol，例如：gene TP53 或 gene "TP53,EGFR"')

    hit_sets, missing = {}, []
    for g in genes:
        row = conn.execute("SELECT set_pks FROM gene_index WHERE gene=?", (g,)).fetchone()
        if row is None:
            missing.append(g)
            continue
        pks = array.array("I")
        pks.frombytes(row[0])
        hit_sets[g] = set(pks)

    if not hit_sets:
        print(f"索引中没有这些基因：{', '.join(missing)}")
        print("提示：请用 MSigDB 官方 symbol（TP53 / NFKB1 / MTOR），不要用蛋白名或别名。")
        return []

    if len(hit_sets) == 1:
        mode, merged = "single", list(hit_sets.values())[0]
    elif args.mode == "intersect":
        mode, merged = "intersect", set.intersection(*hit_sets.values())
    else:
        mode, merged = "union", set.union(*hit_sets.values())

    if not merged:
        print("给定基因没有共同出现的基因集（多基因时请去掉 --mode intersect）。")
        return []

    where, wparams = build_where(args.collection, args.file, primary_only=not args.keep_dups)
    sql = (
        "SELECT set_pk,set_id,collection,source_file,source_db,n_genes FROM sets"
        f" WHERE set_pk IN ({','.join('?'*len(merged))})"
    )
    params = list(merged)
    if where:
        sql += " AND " + where
        params += wparams
    rows = conn.execute(sql, params).fetchall()
    rows = sorted(rows, key=lambda r: (r["collection"], -(r["n_genes"] or 0), r["set_id"]))

    scope = ",".join(args.collection) if args.collection else "全部集合"
    mode_txt = {"single": "单基因", "union": "并集", "intersect": "交集"}[mode]
    print(f"MSigDB 反向检索: {', '.join(genes)}    范围={scope}    模式={mode_txt}")
    if missing:
        print(f"注意：未找到 {', '.join(missing)}（可能是别名或非官方 symbol）")
    print(f"命中 {fmt_int(len(rows))} 个基因集")
    print("=" * 100)

    bycoll = {}
    for r in rows:
        bycoll.setdefault(r["collection"], []).append(r)
    print("按集合族分布：")
    for coll in sorted(bycoll, key=lambda c: -len(bycoll[c])):
        dbs = {}
        for r in bycoll[coll]:
            dbs[r["source_db"] or "-"] = dbs.get(r["source_db"] or "-", 0) + 1
        detail = "  ".join(f"{k}:{v}" for k, v in sorted(dbs.items(), key=lambda x: -x[1])[:6])
        print(f"  {coll:<3} {fmt_int(len(bycoll[coll])):>6} 个    {detail}")
    if not rows:
        print("（当前过滤条件下没有命中）")
        return []
    print()

    limit = len(rows) if args.all else max(1, args.limit)
    for i, r in enumerate(rows[:limit], 1):
        print(f"[{i}] {r['set_id']}")
        print(
            f"    集合={r['collection']}  来源={(r['source_db'] or '-')}  "
            f"基因数={fmt_int(r['n_genes'])}  文件={r['source_file']}"
        )
    if len(rows) > limit:
        print(f"\n… 另有 {fmt_int(len(rows)-limit)} 个未显示（--all 显示全部，-n 调整）")

    gene_maps = fetch_genes_map(conn, [r["set_pk"] for r in rows])
    out = []
    for r in rows:
        allg = gene_maps.get(r["set_pk"], [])
        out.append(
            {
                "set_id": r["set_id"],
                "collection": r["collection"],
                "source_db": r["source_db"] or "",
                "source_file": r["source_file"],
                "n_genes": r["n_genes"],
                "matched_genes": ",".join(g for g in genes if g in allg),
                "genes": allg,
            }
        )
    if args.csv:
        with open(args.csv, "w", encoding="utf-8-sig", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["set_id", "collection", "source_db", "source_file", "n_genes", "matched_genes", "genes"])
            for r in out:
                w.writerow(
                    [
                        r["set_id"],
                        r["collection"],
                        r["source_db"],
                        r["source_file"],
                        r["n_genes"],
                        r["matched_genes"],
                        " ".join(r["genes"]),
                    ]
                )
        print(f"[已导出] {args.csv}")
    if args.gmt:
        write_gmt(args.gmt, conn, [r["set_pk"] for r in rows])
        print(f"[已导出] {args.gmt}")
    return out


# ---------------------------------------------------------------- list / build


def cmd_list(conn, args):
    rows = conn.execute(
        "SELECT collection, COUNT(*) n, SUM(n_genes) g, MIN(version) v"
        " FROM sets WHERE is_primary=1 GROUP BY collection ORDER BY collection"
    ).fetchall()
    meta = dict(conn.execute("SELECT key,value FROM meta").fetchall())
    total_sets = sum(r["n"] for r in rows)
    total_entries = sum(r["g"] or 0 for r in rows)
    n_rows = int(meta.get("n_rows", total_sets))
    print(f"数据目录 : {meta.get('data_dir','')}")
    print(f"索引构建 : {meta.get('built_at','')}   索引文件 : {args._db}")
    print(
        f"规模     : {fmt_int(total_sets)} 个基因集 / {fmt_int(meta.get('n_genes',0))} 个唯一基因 / "
        f"{fmt_int(total_entries)} 条基因条目"
    )
    if n_rows != total_sets:
        print(
            f"重复处理 : 共 {fmt_int(n_rows)} 行，其中 {fmt_int(n_rows-total_sets)} 行为同名 set_id 的非主版本"
            "（默认不返回，加 --keep-dups 可查）"
        )
    print("=" * 100)
    print(f"{'集合':<6}{'版本':<10}{'基因集数':>10}{'基因条目':>12}   说明")
    for r in rows:
        print(
            f"{r['collection']:<6}{(r['v'] or ''):<10}{fmt_int(r['n']):>10}{fmt_int(r['g'] or 0):>12}   "
            f"{COLL_LABEL.get(r['collection'],'')}"
        )
    print()
    print("来源库分布（可用于 --collection 过滤）：")
    dbs = conn.execute(
        "SELECT source_db db, COUNT(*) n FROM sets WHERE source_db<>'' AND is_primary=1"
        " GROUP BY source_db ORDER BY n DESC"
    ).fetchall()
    print("  " + "   ".join(f"{d['db']}:{fmt_int(d['n'])}" for d in dbs))
    print()
    print("GMT 文件：")
    for r in conn.execute(
        "SELECT source_file f, SUM(is_primary) p, COUNT(*) n FROM sets"
        " GROUP BY source_file ORDER BY source_file"
    ):
        extra = "" if r["n"] == r["p"] else f"（含 {fmt_int(r['n']-r['p'])} 个非主版本）"
        print(f"  {r['f']:<38}{fmt_int(r['p']):>8}{extra}")
    print()
    print("可选 --collection 别名：")
    print("  " + ", ".join(sorted(ALIAS_COND)))
    print()
    print("提示：PATHWAY = 经典信号通路（H + C2 的 KEGG/REACTOME/WP/BIOCARTA/PID）")
    print("      不知道该用哪几个集合？先跑 guide（见 SKILL.md 第 0 节的选择流程）")


# ---------------------------------------------------------------- 集合选择向导

# 十大集合速查：适合什么组织、什么实验、能回答什么（guide 子命令打印）
GUIDE_SETS = [
    dict(
        coll="H",
        alias="H / HALLMARK",
        n="50 集",
        pathway="是（通路级，已收敛）",
        content="50 条被专家统一收敛的标志性过程：凋亡、缺氧、EMT、炎症、增殖(E2F/MYC)、"
        "PI3K-AKT-mTOR、代谢、干性等",
        tissue="任何组织、任何样本（人类 symbols 版）",
        exp="bulk RNA-seq / 芯片的差异表达、临床队列、多组学打分、GSEA 富集",
        use_when="想知道「这批基因整体落在哪些核心生物学过程」——最稳的首选，结果干净好解释",
        limit="粒度粗，问不到具体分子机制通路；不反映细胞类型构成",
        scenes=[
            "常规 bulk 差异基因的第一层富集：先看整体落在哪些核心过程",
            "临床队列 / 多组学分组比较（高低风险、应答 vs 无应答）",
            "跨数据集、跨平台比较：50 条口径统一，不同研究能对上号",
            "用 GSVA / ssGSEA 给每个样本打 50 个过程的分，做连续变量分析",
            "论文里需要一张「条数少、干净、好解释」的富集图",
        ],
        output="GSEA / ORA 富集图；GSVA 打分成「样本 × 过程」矩阵",
    ),
    dict(
        coll="C1",
        alias="C1 / POSITIONAL / LOCATION",
        n="302 集",
        pathway="不是（染色体位置，不是功能通路）",
        content="按染色体臂 / 细胞遗传带（cytoband）划分的基因位置集合",
        tissue="任何组织；实际主要用在肿瘤遗传学、核型与拷贝数研究",
        exp="CNA / CNV / 核型分析、染色体臂扩增缺失、把变异区间映射到基因",
        use_when="要看「这条染色体臂 / 这个区域上的基因整体是否异常」时",
        limit="做常规表达富集会有位置偏倚（脚本默认已对其降权 -60），别当功能通路用",
        scenes=[
            "肿瘤 CNA / CNV：把扩增或缺失区段映射成基因清单",
            "染色体臂级研究（1q 扩增、9p 缺失这类经典事件）",
            "白血病等血液肿瘤的细胞遗传学分组",
            "位置偏倚对照：判断某个基因是否只是「跟着染色体臂一起变」的假阳性",
            "全基因组区间 / 区域内基因的整体行为",
        ],
        output="CNA 区段 → 基因注释；做富集时充当位置偏倚对照",
    ),
    dict(
        coll="C2",
        alias="C2 / PATHWAY / CP / CGP",
        n="7,670 集",
        pathway="是（真正意义上的信号通路几乎都在这）",
        content="专家审编：CP 经典通路 4,115（KEGG 844 + REACTOME 1,839 + WikiPathways 925 + "
        "BioCarta 292 + PID 196）+ CGP 扰动签名 3,574",
        tissue="任何组织、任何实验体系",
        exp="bulk 差异表达富集、通路打分、多组学；CGP 子库对应「处理 / 过表达 / 敲低」类细胞实验",
        use_when="用户说「我要看信号通路」——就是这里；别名 --collection PATHWAY = H + 这部分经典通路",
        limit="KEGG 子库是 LEGACY 186 集，不含 PI3K-Akt / TNF 整通路（见 references/msigdb-collections.md）",
        scenes=[
            "用户说「找信号通路」时的主战场（KEGG / REACTOME / WP / BioCarta / PID）",
            "要机制级细节：具体到受体、激酶、级联步骤——REACTOME 做得最细",
            "画通路图 / 机制示意图，需要这条通路的完整基因清单（WikiPathways 有配图）",
            "多数据库交叉验证：同一条通路在 KEGG 与 REACTOME 都命中，结论更稳",
            "扰动类实验：CGP 子库是药物 / 过表达 / 敲低的表达签名集",
        ],
        output="GSEA / ORA 富集；导出 gmt 供 clusterProfiler 或 GSEA 软件直接读入",
    ),
    dict(
        coll="C3",
        alias="C3 / REGULATORY / TFT / MIR",
        n="3,714 集",
        pathway="半是（调控靶标，不是信号级联）",
        content="调控靶标：TFT（转录因子靶基因，来自 ChIP / 预测）+ MIR（miRNA 靶基因）",
        tissue="任何组织；肿瘤、发育、免疫方向都用得多",
        exp="ChIP-seq / CUT&Tag / ATAC-seq、miRNA 测序或 mimics / inhibitor、转录调控研究",
        use_when="做过或打算做「某个 TF / 某个 miRNA 的下游靶基因是什么」",
        limit="默认检索对 MIR/TFT 降权（生物学家检索时通常不想要），选 C3 时不要再叠 PATHWAY",
        scenes=[
            "ChIP-seq / CUT&Tag / ATAC 拿到 peak 后，问「这是哪个转录因子的靶标」",
            "miRNA 测序或 mimics / inhibitor 实验，找下游靶基因",
            "单基因通路分析：只关心某一个 TF 的下游会牵出什么",
            "上游调控反推：把差异基因对应回可能的调控 TF",
            "用 TFT 集合做转录因子活性打分（ssGSEA 给 TF 排序）",
        ],
        output="TF / miRNA 靶基因清单；转录因子活性打分表",
    ),
    dict(
        coll="C4",
        alias="C4 / COMPUTATIONAL / MODULE",
        n="1,006 集",
        pathway="不是（计算出来的共表达模块）",
        content="计算型：3CA 癌症共表达模块、CGN 共表达邻域、CM 模块",
        tissue="任何组织；肿瘤与复杂疾病的网络分析最常见",
        exp="没有先验通路假设的探索性分析、WGCNA 类共表达网络",
        use_when="想数据驱动地先看「数据自己聚成了哪些模块」",
        limit="集合名是抽象模块编号，名字本身没有生物学含义；解释要靠模块内的基因",
        scenes=[
            "完全没有先验假设的探索：先看数据自己聚成哪些模块，再逐个注释",
            "共表达网络 / WGCNA 的模块与已知模块做比对",
            "肿瘤分子亚型分型（3CA 癌症共表达模块）",
            "找「这套数据特有」的基因群，而不是套公共通路",
        ],
        output="模块 → 基因群清单；再拿模块去做注释或分型",
    ),
    dict(
        coll="C5",
        alias="C5 / GO / GOBP / GOCC / GOMF / HPO",
        n="16,283 集",
        pathway="GO:BP 算通路级，CC/MF 是功能与组分，HPO 是表型",
        content="本体：GO:BP 7,538 + GO:CC 1,080 + GO:MF 1,872 + HPO 5,793",
        tissue="GO 部分任何组织；HPO 部分面向人类临床与遗传病表型",
        exp="标准 GO 富集（任何转录组 / 蛋白组）；HPO 用于患者表型、变异解读、罕见病",
        use_when="要穷尽式功能注释，或做的是临床遗传 / 表型方向",
        limit="GO 层级冗余严重，结果动辄上千条，必须按层级或去冗余后再挑",
        scenes=[
            "标准 GO 富集，BP / CC / MF 分开报——审稿人最认的口径",
            "穷尽式注释，不想漏掉细小的功能项",
            "GO:CC 做定位研究（膜 / 胞质 / 核 / 细胞器）",
            "GO:MF 做分子功能（酶活、结合、转运）",
            "HPO 做临床：患者表型、罕见病、变异致病性解读、表型-基因关联",
        ],
        output="GO 富集三件套（BP / CC / MF）；HPO 表型富集",
    ),
    dict(
        coll="C6",
        alias="C6 / ONCOGENIC",
        n="189 集",
        pathway="是（致癌通路激活状态签名）",
        content="致癌签名：癌基因激活 / 抑癌基因失活 / 突变型肿瘤的扰动签名",
        tissue="肿瘤组织、癌细胞系",
        exp="肿瘤 bulk RNA-seq、细胞系 panel、药物敏感性关联",
        use_when="想判断「样本里哪条致癌通路被激活了」（KRAS 激活、E2F 靶标、抑癌失活等）",
        limit="只有肿瘤方向适用；非肿瘤样本别用",
        scenes=[
            "判断样本里哪条致癌通路被激活（KRAS、MYC、E2F、β-catenin 等）",
            "肿瘤分子分型与亚型间比较",
            "抑癌基因失活的间接推断",
            "癌旁 vs 肿瘤：看哪种致癌程序被打开",
            "细胞系 panel 里把签名与药物敏感性关联",
        ],
        output="致癌通路激活状态打分；肿瘤分型依据",
    ),
    dict(
        coll="C7",
        alias="C7 / IMMUNE",
        n="5,219 集",
        pathway="否（免疫细胞状态签名）",
        content="免疫签名：ImmuneSigDB（免疫细胞类型 / 状态 / 扰动）+ VAX（疫苗应答）",
        tissue="血液、脾、淋巴结、肿瘤免疫微环境等任何含免疫细胞的组织",
        exp="免疫浸润与炎症研究、感染与自身免疫、肿瘤免疫与免疫治疗应答、bulk 免疫组成推断",
        use_when="研究里只要有「免疫 / 炎症」这条线，就该把它拉进来",
        limit="区分度极高，非免疫样本也会命中；结论必须结合细胞组成一起解释",
        scenes=[
            "肿瘤免疫微环境：免疫浸润程度、热肿瘤 vs 冷肿瘤",
            "bulk 数据推断免疫细胞组成（免疫签名打分，或做解卷积的参考）",
            "感染与疫苗研究：抗病毒应答、VAX 疫苗应答签名",
            "自身免疫与炎症性疾病：类风湿、IBD、银屑病等",
            "免疫治疗 / 免疫检查点研究的应答分层",
            "单细胞里给免疫亚群做状态注释（耗竭、记忆、活化）",
        ],
        output="免疫浸润打分；免疫细胞状态标签",
    ),
    dict(
        coll="C8",
        alias="C8 / CELLTYPE / MARKER",
        n="866 集",
        pathway="否（细胞类型标记）",
        content="细胞类型签名：单细胞研究里积累的各组织细胞类型标记集",
        tissue="任何有细胞类型异质性的组织（脑、肝、肾、肿瘤、外周全血…）",
        exp="scRNA-seq / snRNA-seq / 空间转录组的细胞类型注释、bulk 组成解卷积（CIBERSORTx / MuSiC）",
        use_when="要做细胞类型注释，或想知道 bulk 样本里各细胞类型的占比",
        limit="只有 866 集、覆盖有限；注释不到就退回「marker 基因 + H 打分」的路子",
        scenes=[
            "scRNA-seq / snRNA-seq 的细胞类型注释（配合 marker 基因一起看）",
            "空间转录组 spot / 区域的细胞构成",
            "bulk 解卷积（CIBERSORTx / MuSiC）的参考集",
            "组织构成变化：肿瘤 vs 正常、纤维化、脂肪肝这类",
            "类器官 / 共培养体系确认细胞身份",
        ],
        output="细胞类型标签；bulk 样本各细胞类型占比",
    ),
    dict(
        coll="C9",
        alias="C9 / PERTURB / DEPMAP",
        n="62 集",
        pathway="否（扰动响应签名）",
        content="计算扰动签名：DepMap CRISPR 敲除 + CCLE 细胞系表达 / 药物数据推导的响应签名",
        tissue="细胞系为主（肿瘤细胞系最典型）",
        exp="CRISPR screen、siRNA / shRNA 敲低、药物处理转录组、靶点验证",
        use_when="想回答「我敲了这个基因 / 用了这个药，转录组变化像不像文献里某个已知扰动」",
        limit="只有 62 集、覆盖面窄；问机制通路别用它",
        scenes=[
            "CRISPR screen / 敲低实验：敲掉之后转录组像哪个已知扰动",
            "药物机制推断（MoA）：药物响应像哪个基因扰动",
            "靶点验证与脱靶评估",
            "细胞系 panel 里把药物敏感性与 DepMap 基因依赖关联",
        ],
        output="扰动相似度匹配结果；靶点 / 药物机制的佐证",
    ),
]

# 场景 -> 推荐集合（keys 全部小写，对用户描述做子串匹配）
SCENARIOS = [
    dict(
        name="肿瘤 / 癌症 bulk 转录组（差异表达 → 通路）",
        keys=[
            "肿瘤", "癌", "cancer", "tumor", "tumour", "oncolog", "腺癌", "癌旁", "活检",
            "实体瘤", "肝癌", "肺癌", "胃癌", "乳腺癌", "结直肠癌", "白血病",
        ],
        main=["H", "C2"],
        extra=["C6", "C5"],
        why="H 给核心生物学过程、C2 给具体信号通路（KEGG/REACTOME/WP/PID）；再看 C6 有没有致癌通路被激活",
        note="C2 的 KEGG 子库是 LEGACY 186 集，不含 PI3K-Akt / TNF 整通路，用 REACTOME + H 替代",
    ),
    dict(
        name="免疫 / 炎症 / 感染 / 自身免疫 / 肿瘤免疫微环境",
        keys=[
            "免疫", "炎症", "感染", "自身免疫", "淋巴", "巨噬", "t细胞", "b细胞", "干扰素",
            "immune", "inflamm", "infiltrat", "vaccine", "疫苗", "免疫治疗", "pd-1", "pd1",
        ],
        main=["C7", "H"],
        extra=["C2", "C8"],
        why="C7 是 ImmuneSigDB + VAX 免疫签名，专治免疫/炎症细胞状态；H 补整体过程；要拆细胞组成加 C8",
        note="非免疫组织也可能命中免疫签名，结论要结合细胞组成一起解释",
    ),
    dict(
        name="单细胞 / 空间转录组 / 细胞类型注释 / 组成解卷积",
        keys=[
            "单细胞", "单核", "scrna", "snrna", "single cell", "single-cell", "细胞类型", "cell type",
            "marker", "标记基因", "空间转录", "spatial", "细胞组成", "解卷积", "deconvol",
            "cibersort", "musiq", "亚群", "cluster 注释",
        ],
        main=["C8", "H"],
        extra=["C2", "C7"],
        why="C8 是单细胞来源的细胞类型签名，做注释 / 占比 / 解卷积最直接；H 用来给 cluster 打基因集分数",
        note="C8 只有 866 集，注释不到就退回「marker 基因 + H 打分」的路子",
    ),
    dict(
        name="CRISPR / 敲低 / 过表达 / 药物处理（功能基因组学与扰动）",
        keys=[
            "crispr", "敲除", "敲低", "敲减", "sirna", "shrna", "过表达", "扰动", "perturb",
            "depmap", "药物处理", "药物筛选", "耐药", "靶点验证", "screen", "抑制剂", "treatment",
        ],
        main=["C9", "C2"],
        extra=["H", "C6"],
        why="C9 是 DepMap CRISPR / CCLE 计算扰动签名，把「敲了谁 / 用了什么药」的响应映射到已知扰动；"
        "C2 的 CGP 子库是海量扰动签名",
        note="CGP 子库用 --collection CGP（= C2 中 3,574 个非经典通路来源的集）",
    ),
    dict(
        name="ChIP-seq / ATAC-seq / 转录因子 / miRNA 靶标",
        keys=[
            "chip", "cut&tag", "cut and tag", "atac", "转录因子", "mirna", "microrna", "靶基因",
            "靶标", "结合位点", "motif", "启动子", "regulon", "tft",
        ],
        main=["C3"],
        extra=["C2", "H"],
        why="C3 存 TFT（TF 靶标）与 MIR（miRNA 靶标），正好回答「这个 TF / miRNA 的下游是什么」",
        note="选 C3 时不要再叠 PATHWAY，否则默认降权会让 TFT/MIR 沉底",
    ),
    dict(
        name="基础生物学过程（代谢 / 缺氧 / 凋亡 / 增殖 / 应激），非肿瘤体系",
        keys=[
            "代谢", "metabol", "缺氧", "hypoxia", "凋亡", "apoptos", "增殖", "proliferat",
            "细胞周期", "cell cycle", "应激", "stress", "自噬", "autophag", "衰老", "senescence",
            "发育", "develop", "干细胞", "stem", "分化", "differentiat", "小鼠", "大鼠", "斑马鱼",
        ],
        main=["H", "C2"],
        extra=["C5"],
        why="H 的 50 条正好就是这些核心过程；C2 补具体机制通路；要穷尽式功能注释再叠 C5(GO:BP)",
        note="非人类物种注意：本库是 Hs symbols，做小鼠请先转成大写人类 symbol 写法",
    ),
    dict(
        name="临床 / 遗传病 / 罕见病 / 患者表型",
        keys=[
            "临床", "患者", "遗传病", "罕见病", "表型", "phenotype", "hpo", "变异", "突变",
            "孟德尔", "家系", "疾病", "诊断",
        ],
        main=["C5", "H"],
        extra=["C2"],
        why="C5 含 HPO 5,793 个表型 / 疾病集，是临床与遗传学表型分析的唯一来源；表达侧仍用 H / C2",
        note="只用 HPO 时写 --collection HPO（或 HP）",
    ),
    dict(
        name="拷贝数 / 染色体 / 核型 / 染色体臂",
        keys=[
            "拷贝数", "cna", "cnv", "染色体", "chromosom", "核型", "karyotyp", "扩增", "amplif",
            "缺失", "deletion", "染色体臂", "cytoband", "位置",
        ],
        main=["C1"],
        extra=["C6", "C2"],
        why="C1 是染色体细胞遗传带位置集，把 CNA 区域映射到基因；功能后果再看 C6 / C2",
        note="只做 CNA 分析时才主用 C1；常规表达富集不要选它",
    ),
    dict(
        name="没有明确通路假设，想数据驱动地找模块",
        keys=[
            "没有假设", "数据驱动", "模块", "module", "共表达", "co-express", "wgcna", "聚类",
            "无假设", "探索", "未知机制",
        ],
        main=["C4", "C2"],
        extra=["H"],
        why="C4 是计算型共表达模块（3CA / CGN / CM），适合先看数据自己聚出的模块；C2 / H 给可读的通路名",
        note="C4 的集合名是抽象模块号，名字本身没有生物学含义",
    ),
    dict(
        name="只是查一条已知通路包含哪些基因 / 导出 gmt 给 GSEA",
        keys=[
            "查通路", "某条通路", "包含哪些基因", "基因列表", "导出", "gmt", "gsea", "富集",
            "enrichment", "通路基因", "信号通路",
        ],
        main=["H", "C2"],
        extra=["C5"],
        why="经典信号通路就是 H + C2 的 KEGG/REACTOME/WP/BIOCARTA/PID，别名 --collection PATHWAY 一步到位",
        note="免疫通路请叠 C7；细胞类型标记请叠 C8；具体机制通路用 REACTOME / WP 更全",
    ),
]

GUIDE_SHORT = {
    "H": "H Hallmark 标志性过程",
    "C1": "C1 位置型（染色体带）",
    "C2": "C2 专家审编经典通路（KEGG/REACTOME/WP/PID/BIOCARTA）",
    "C3": "C3 调控靶标（TF / miRNA）",
    "C4": "C4 计算型共表达模块",
    "C5": "C5 本体（GO / HPO）",
    "C6": "C6 致癌签名",
    "C7": "C7 免疫签名",
    "C8": "C8 细胞类型签名",
    "C9": "C9 计算扰动签名（DepMap / CCLE）",
}

GUIDE_DEFAULT = [
    ("-c PATHWAY", "经典信号通路（H + KEGG/REACTOME/WP/BIOCARTA/PID，约 3,900 集）：通用默认，覆盖大多数需求"),
    ("-c PATHWAY,GOBP", "再加 C5 的 GO:BP，做穷尽式功能注释时用"),
    ("-c H,C2,C5,C7,C8", "确实想「都看看」时的五集合折中组合（不含 C1/C3/C4/C9 这些专用签名）"),
]


def _coll_text(colls):
    return " + ".join(GUIDE_SHORT.get(c, c) for c in colls)


def _print_guide_set(g):
    print(f"【{g['coll']}】{g['alias']}   （{g['n']}）")
    print(f"  是信号通路吗 : {g['pathway']}")
    print(f"  内容         : {g['content']}")
    print(f"  适合组织     : {g['tissue']}")
    print(f"  适合实验     : {g['exp']}")
    print("  适合的场景   :")
    for s in g.get("scenes") or []:
        print(f"    · {s}")
    if g.get("output"):
        print(f"  典型产出     : {g['output']}")
    print(f"  什么时候选它 : {g['use_when']}")
    print(f"  不适合 / 坑  : {g['limit']}")
    print()


def cmd_guide(args):
    q = " ".join(args.keyword or []).strip().lower()
    hits = []
    if q:
        for sc in SCENARIOS:
            matched = [k for k in sc["keys"] if k in q]
            if matched:
                # 打分 = 命中关键词总长度 + 命中个数：越长的词越具体，
                # 这样「单细胞 / 细胞类型」能压过同样命中的宽泛词（如「免疫」）。
                hits.append((sum(len(k) for k in matched) + len(matched), sc))
        hits.sort(key=lambda x: -x[0])

    print("MSigDB 十大集合选择向导（MSigDB v2026.1.Hs，10 个集合 / 35,361 个基因集）")
    print("=" * 92)

    if q and not hits:
        print(f"查询：{q}")
        print("没有匹配到预设场景 —— 下面是全部 10 个集合的适用面对照，按用户实际的组织 + 实验类型挑。")
        print()

    want = None
    if hits:
        print(f"查询：{q}")
        print()
        for rank, (score, sc) in enumerate(hits[:3]):
            tag = "主推" if rank == 0 else "备选"
            print(f"[{tag}] {sc['name']}")
            print(f"  建议集合 : {_coll_text(sc['main'])}")
            print(f"             --collection {','.join(sc['main'])}")
            if sc.get("extra"):
                print(f"  可选叠加 : {_coll_text(sc['extra'])}")
                print(f"             --collection {','.join(sc['extra'])}")
            print(f"  理由     : {sc['why']}")
            if sc.get("note"):
                print(f"  注意     : {sc['note']}")
            print()
        want = set(hits[0][1]["main"]) | set(hits[0][1].get("extra") or [])
        print("-" * 92)
        if args.all:
            print("下面展开全部 10 个集合；每个集合的「适合的场景」清单可直接挑 2~3 条讲给用户。")
        else:
            print(f"下面只展开主推场景涉及的集合（{'/'.join(sorted(want))}）。")
            print("每个集合都带了「适合的场景」清单，转述给用户时挑 2~3 条最贴合他体系的即可；")
            print("想看全部 10 个加 --all（或不带关键词再跑一次 guide）。")
        print()

    for g in GUIDE_SETS:
        if want is not None and g["coll"] not in want and not args.all:
            continue
        _print_guide_set(g)

    print("-" * 92)
    print("默认建议（用户说「都看看 / 你决定」时按这个顺序）：")
    for cmd, why in GUIDE_DEFAULT:
        print(f"  {cmd:<20} {why}")
    print()
    print("选完集合后，再按用户要的方向检索（都带 --collection）：")
    print('  gmt_lookup.py find "TNF signaling" --collection H,C2 -n 10')
    print("  gmt_lookup.py gene FOXO1 --collection PATHWAY")
    print("提示：具体操作流程与「先问组织 + 实验再选集合」的强制步骤见 SKILL.md 第 0 节。")


def cmd_build(data_dir, args):
    stats = build_index(data_dir, args._db, quiet=args.quiet)
    print("索引构建完成")
    print(f"  数据目录 : {data_dir}")
    print(f"  索引文件 : {stats['db_path']}")
    print(
        f"  规模     : {fmt_int(stats['n_sets'])} 个基因集 / {fmt_int(stats['n_genes'])} 个唯一基因 / "
        f"{fmt_int(stats['n_entries'])} 条基因条目"
    )
    print(f"  耗时     : {stats['seconds']}s")
    for f, coll, n in stats["per_file"]:
        print(f"    {coll:<4}{f:<38}{fmt_int(n):>8}")
    warn_dups(stats)
    return stats


# ---------------------------------------------------------------- CLI


def _common_parser():
    """所有子命令共享的选项。默认值用 SUPPRESS，这样放在子命令前后都能生效。"""
    c = argparse.ArgumentParser(add_help=False)
    c.add_argument(
        "--data-dir",
        default=argparse.SUPPRESS,
        help="GMT 目录（默认用 skill 自带的 data/，可用本项指向自己的目录）",
    )
    c.add_argument(
        "--db",
        default=argparse.SUPPRESS,
        help="索引文件路径（默认 <data-dir>/.msigdb_index.sqlite3；数据目录只读时自动放到用户缓存）",
    )
    c.add_argument("--no-build", action="store_true", default=argparse.SUPPRESS, help="索引缺失时不要自动重建")
    c.add_argument("--file", action="append", default=argparse.SUPPRESS, help="按 GMT 文件名（子串）过滤，可重复")
    c.add_argument("--gmt", default=argparse.SUPPRESS, help="把命中的基因集导出为标准 GMT 文件")
    c.add_argument(
        "--keep-dups",
        action="store_true",
        default=argparse.SUPPRESS,
        help="同时返回同名 set_id 的非主版本（另有旧版 GMT 文件时才需要）",
    )
    return c


def main(argv=None):
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except Exception:
        pass

    cfg = resolve_config()
    common = _common_parser()
    p = argparse.ArgumentParser(
        prog="gmt_lookup.py",
        parents=[common],
        description="MSigDB GMT 通路 <-> 基因 双向快速检索",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "示例：\n"
            "  gmt_lookup.py guide '人肝癌 bulk rna-seq 差异表达'   # 先定范围：推荐该查哪几个集合\n"
            "  gmt_lookup.py list\n"
            '  gmt_lookup.py find "TNF signaling"\n'
            '  gmt_lookup.py find "TNF signaling" "WNT signaling" apoptosis   # 一次查多个关键词\n'
            "  gmt_lookup.py find apoptosis --collection H --csv out.csv\n"
            '  gmt_lookup.py find "PI3K-Akt" --preview 0\n'
            "  gmt_lookup.py show HALLMARK_ADIPOGENESIS --gmt subset.gmt\n"
            "  gmt_lookup.py gene TP53 --collection PATHWAY\n"
            '  gmt_lookup.py gene "TP53,EGFR" --mode intersect\n'
            "\n"
            "范围别名 --collection：H/HALLMARK, C1..C9, KEGG, REACTOME, WP, BIOCARTA, PID,\n"
            "                      GO/GOBP/GOCC/GOMF, HP/HPO, MIR, TFT, PATHWAY(经典通路), CURATED,\n"
            "                      语义别名 IMMUNE(C7), CELLTYPE(C8), PERTURB(C9), ONCOGENIC(C6),\n"
            "                      REGULATORY(C3), POSITIONAL(C1), ONTOLOGY(C5), COMPUTATIONAL(C4),\n"
            "                      CANCER(C6+C9), CP/CLASSIC(C2 经典通路), CGP(C2 扰动签名)\n"
        ),
    )

    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("list", parents=[common], help="列出集合与统计").set_defaults(func="list")

    s = sub.add_parser(
        "guide", parents=[common], help="选择向导：按组织/实验类型推荐该用哪几个集合"
    )
    s.add_argument(
        "keyword",
        nargs="*",
        help="可选：用户的组织/样本/实验类型描述，如 '人肝癌 bulk rna-seq 差异表达'",
    )
    s.add_argument(
        "--all",
        "-a",
        action="store_true",
        help="即使命中了场景，也把全部 10 个集合的适合场景一并展开",
    )
    s.set_defaults(func="guide")

    s = sub.add_parser("build", parents=[common], help="重建索引")
    s.add_argument("--quiet", action="store_true", help="不打印逐文件进度")
    s.set_defaults(func="build")

    s = sub.add_parser("find", parents=[common], help="正向：通路名/关键词 -> 基因")
    s.add_argument("query", nargs="+", help="一个或多个通路名/关键词，如 'TNF signaling' MAPK apoptosis")
    s.add_argument("--collection", "-c", action="append", default=[], help="限定集合/来源，可重复")
    s.add_argument("--limit", "-n", type=int, default=int(cfg.get("default_limit", 20)), help="最多显示几个基因集")
    s.add_argument("--all", action="store_true", help="显示全部命中")
    s.add_argument("--exact", action="store_true", help="只要名称精确匹配")
    s.add_argument("--preview", type=int, default=12, help="每个基因集预览前 N 个基因；0=全部")
    s.add_argument("--format", choices=["table", "json", "genes"], default="table")
    s.add_argument("--csv", default=None, help="导出长格式 CSV")
    s.set_defaults(func="find")

    s = sub.add_parser("show", parents=[common], help="正向：精确集合名 -> 全部基因")
    s.add_argument("set_id", help="完整集合名，如 HALLMARK_APOPTOSIS")
    s.add_argument("--format", choices=["inline", "genes"], default="inline")
    s.add_argument("--csv", default=None, help="导出 CSV")
    s.set_defaults(func="show")

    s = sub.add_parser("gene", parents=[common], help="反向：基因 -> 参与的通路")
    s.add_argument("genes", help="一个或多个基因 symbol，逗号/空格分隔")
    s.add_argument("--collection", "-c", action="append", default=[], help="限定集合/来源")
    s.add_argument("--mode", choices=["union", "intersect"], default="union", help="多基因：并集或交集")
    s.add_argument("--limit", "-n", type=int, default=40, help="最多显示几个基因集")
    s.add_argument("--all", action="store_true", help="显示全部命中")
    s.add_argument("--format", choices=["table", "json"], default="table")
    s.add_argument("--csv", default=None, help="导出 CSV")
    s.set_defaults(func="gene")

    args = p.parse_args(argv)

    # 共享选项在子命令前后都可能出现，统一用 getattr 兜底
    collapsed = []
    for c in getattr(args, "collection", None) or []:
        collapsed.extend([x for x in re.split(r"[,\s]+", str(c)) if x])
    args.collection = collapsed
    args.file = [x for x in (getattr(args, "file", None) or []) if x]
    args.data_dir = getattr(args, "data_dir", None)
    args.db = getattr(args, "db", None)
    args.no_build = getattr(args, "no_build", False)
    args.keep_dups = getattr(args, "keep_dups", False)
    args.gmt = getattr(args, "gmt", None)

    # guide 是纯静态向导，不需要索引，先处理掉，避免无谓的索引路径解析与重建
    if args.func == "guide":
        cmd_guide(args)
        return

    data_dir = resolve_data_dir(args.data_dir)
    # 索引路径只解析一次，避免 default_index_path 的只读回退提示被打印两次
    args._db = args.db or os.environ.get("MSIGDB_INDEX") or default_index_path(data_dir)

    if args.func == "build":
        cmd_build(data_dir, args)
        return

    conn, db_path = open_index(data_dir, cli_db=args._db, auto_build=not args.no_build)
    args._db = db_path
    try:
        if args.func == "list":
            cmd_list(conn, args)
        elif args.func == "find":
            res = cmd_find(conn, args)
            if args.format == "json":
                print(json.dumps(res, ensure_ascii=False, indent=1))
            elif args.format == "genes":
                seen, uniq = set(), []
                for r in res:
                    for g in r["genes"]:
                        if g not in seen:
                            seen.add(g)
                            uniq.append(g)
                print(" ".join(uniq))
        elif args.func == "show":
            cmd_show(conn, args)
        elif args.func == "gene":
            res = cmd_gene(conn, args)
            if args.format == "json":
                print(json.dumps(res, ensure_ascii=False, indent=1))
    finally:
        conn.close()


if __name__ == "__main__":
    main()
