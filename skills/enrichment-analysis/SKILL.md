---
name: enrichment-analysis
description: 基因列表的富集分析与美化出图，四步式：① ORA（GO BP/CC/MF、KEGG、MSigDB 自定义 GMT），输入是给定基因列表；② GSEA，输入是全部基因按 logFC 排序、不做筛选（结果单独用 NES 条形图展示，不进美化主图）；③ 把全部 p<0.05 的通路统计归纳成功能大类（输出排名前 20 大类的表格 <前缀>_func_top20_<方向>.csv 供用户查阅，并做基因重叠去冗余挑代表通路）后，停下让用户选出图方式；④ 出图：默认均衡配额模式（--mode=even）固定 20 条通路、BP/CC/MF/KEGG 平均各分 5 条，某一板块通路数量不够时富余名额自动匀给其他板块；也可按全库排名取前 N 条或只画某个功能方向（这两种旧模式仍受 --max_terms=15 上限约束）。出图用 gground + ggprism 版式（左侧分类圆角色块 + 通路名 + 斜体基因名 + 基因数圆点，右侧 -log10(p.adjust) 圆角柱）。支持人类与小鼠（org.Hs.eg.db / org.Mm.eg.db，KEGG hsa / mmu），可按全部/上调/下调分别出三套结果。当用户提到 GO 富集、KEGG 富集、通路富集、富集分析、富集图、美化富集图、气泡图/条形图富集、GSEA、基因集富集、Hallmark 富集、MSigDB 富集、通路功能分类、通路归纳、通路统计、功能方向、前五类通路、top10 通路、top20 大类、enrichGO、enrichKEGG、clusterProfiler、gground、富集结果表 等需求时使用。
agent_created: true
---

# 富集分析：ORA（GO / KEGG / MSigDB）+ GSEA + 归纳 + 美化主图

## 流程总览（★ 第 3 步必须停下问用户）

```
上游（可选但最常用）：geo-microarray-analysis 02_deg_plots.R
        └─ 产出 <前缀>_DEG.csv / <前缀>_all.csv / <前缀>_for_enrichment.txt
           ↓ 富集侧一条命令接手：--geo_dir=<geo 产物目录>
0. 先分清 ORA 还是 GSEA（见下表）
1. 01_enrich_ora.R   出全部被检验条目 + 显著标记        -> 9 张表
2. 02_enrich_gsea.R  （可选支线）全部基因排序的 GSEA    -> 表 + NES 条形图
3. 04_enrich_summary.R  把全部 p<0.05 的通路统计归纳成功能方向
        ★★ 跑完必须停：把报告交给用户，让他选 A 还是 B ★★
4. 03_enrich_plot.R  按用户的选择出美化主图
        A) --mode=even（默认）           固定 20 条：BP/CC/MF/KEGG 各 5，不足匀给其他板块
        B) --func_csv=... --category=... 某个功能方向
        C) --mode=global --top=10        全库排名前 10（旧挑法，仍可用）
```

## ★ 出图条数规则（03 步，2026-09-25 起）

1. **美化主图固定 20 条、均衡分配（默认 `--mode=even`）**：整张图条数由
   `--even_total`（默认 20）决定，在实际出现的分类间平均分配——BP/CC/MF/KEGG
   四类齐全时**各 5 条**；某板块通过筛选的通路不足配额时，富余名额按
   **BP→CC→MF→KEGG 顺序轮转**匀给还有存货的板块，总条数仍尽量凑满 20
   （全部库存都不足 20 时有多少画多少，日志会告警）。even 模式下去重先做、
   配额在去重后的池子上算，`--top` 与 `--max_terms` 都不生效（给了会打印提示）。
2. **旧挑法仍可用**：`--mode=global`（全库按 p.adjust 排名取前 N）与
   `--mode=per_ontology`（每类各取前 N，参考代码行为）保留，且仍受
   `--max_terms=15`（默认）硬上限约束。
3. **GSEA 结果不进美化主图**：`03_enrich_plot.R` 读表时检测 GSEA 专有列
   （NES / setSize），发现即报错退出。GSEA 结果**单独展示**，用 02 步自带的
   NES 条形图 `<前缀>_GSEA_top_NES.pdf/png`（按 |NES| 排序、红=激活/蓝=抑制，
   `--plot_top` 控制条数）和可选的富集曲线 `--curve=N`。**不要**把 GSEA 表复制成
   ORA 命名喂给 03，也不要把 GSEA 的显著通路手动混进 ORA 主图 —— 两种口径
   （排序富集 vs 超代表）画在一张图里审稿人会挑。

**为什么第 3 步必须停**：上百条 p<0.05 的通路里，光看排名前 10 会丢掉"这个方向整体
在讲什么"；而按功能方向出图又需要用户先知道有哪些方向、每个方向多少条。归纳结果
是用户做选择的依据，**不要让脚本或 agent 替他选方向**。

## ★ 从 GEO 差异分析直接接手（最常用的一条路）

上一步是 `geo-microarray-analysis` 技能的第二步（`02_deg_plots.R`）。它跑完会在产物
目录写一份交接清单 **`<前缀>_for_enrichment.txt`**（key=value 文本，可人工核对），
里面写清了对比式、分组与例数、上下调条数、以及 `gene_col`/`logfc_col`/`p_col`/
`change_col` 这些列名。

于是富集侧**只要一个 `--geo_dir`** 就够，不用手写八个参数：

```bash
SK="<TOOLKIT>/skills/enrichment-analysis/scripts"

# ① ORA（吃 geo 的 <前缀>_DEG.csv，默认不再二次筛选）
Rscript --vanilla "$SK/01_enrich_ora.R" \
  --geo_dir="<geo 的产物目录>" --gmt_sets=H --direction=all,up,down

# ② GSEA（吃 geo 的 <前缀>_all.csv，即全部基因）
Rscript --vanilla "$SK/02_enrich_gsea.R" \
  --geo_dir="<geo 的产物目录>" --gmt_sets=H --plot=TRUE

# ③ 归纳（★ 跑完把 *_func_report.txt 交给用户，等他选 A / B）
Rscript --vanilla "$SK/04_enrich_summary.R" \
  --in_dir="<geo 的产物目录>/enrich" --prefix=<前缀> --direction=all,up,down

# ④ 按用户的选择出图（A 或 B，见第三步/第四步）
```

`--geo_dir` 会自动做这些事（都在日志的「[geo 交接]」块里逐条打出来，可核对）：

| 自动设定 | 值 | 来源 |
|---|---|---|
| ORA 输入 | `<前缀>_DEG.csv` | 清单 `deg_file`（找不到就退到 `_all.csv` 并告警） |
| GSEA 输入 | `<前缀>_all.csv` | 清单 `all_file` |
| `--ora_p` / `--ora_logfc` | `1` / `0`（**不再二次筛选**） | geo 已经筛过了；让 all/up/down 与 geo 口径逐字一致 |
| `--split_by` | `change` | 用 geo 的 `change` 列（UP/DOWN）拆上下调，**不靠 logFC 正负** |
| `--gene_col` / `--logfc_col` / `--p_col` | 清单里的列名 | `rownames` / `logFC` / `P.Value`（geo 默认） |
| `--species` | `human` | 清单 `species`（**geo 不解析芯片物种，小鼠芯片请显式加 `--species=mouse`**） |
| `--outdir` | `<geo 目录>/enrich` | 固定约定，后续 04/03 直接读这个目录 |
| `--prefix` | `<前缀>` | 清单 `prefix` |

几个要点：

- **`--geo_dir` 只要指向 `02_deg_plots.R` 的产物目录**。⚠️ geo 的 02 步**没有 `--outdir`**，
  产物落在它的**运行目录**（进程 CWD），所以别去猜 `<prefix>_DEG.csv` 在哪，
  直接看清单里的 `deg_file` 或让 `--geo_dir` 指对。
- 读不到清单时会退化为**按文件名约定**找 `<前缀>_DEG.csv` / `<前缀>_all.csv`；
  目录里有多个数据集时要用 `--geo_prefix=` 指认（只有一个时自动认）。
- **想让 ORA 用 geo 的 `_all.csv` 重新筛**（比如换个阈值）：加 `--geo_ora_source=all
  --ora_p_type=adj.P.Val --ora_logfc=1`，或干脆显式给 `--deg=` 覆盖自动识别。
- 用户已经明确说"用哪批基因做富集"时，**以用户为准**，别被清单里的默认值带跑。
- 交接清单只带 `species=human`（geo 流程不解析平台物种）。**见 `mouse` 字样或用户提到
  小鼠芯片，必须显式加 `--species=mouse`。**

## 先分清两条路线（这一步不做对，后面全错）

| | **ORA**（`01_enrich_ora.R`） | **GSEA**（`02_enrich_gsea.R`） |
|---|---|---|
| 输入 | **已经筛好的基因列表**（显著差异基因） | **全部基因**，按 logFC 排序，**不做任何筛选** |
| 问的问题 | 这一小撮基因里，哪条通路被**过度代表**？ | 某条通路的基因是否**系统性堆在排序列表两端**？ |
| 阈值 | 必须给 `--ora_p` / `--ora_logfc` | **没有阈值**，给了也不生效 |

用户在同一个任务里常常两个都要。**默认按这个分岔问清**，不要自己替他把"全部基因"塞进 ORA。

## 前置：环境

- **R**：任意 R ≥ 4.2（包装器自动探测；也可用环境变量 `RSCRIPT` 指定）
- **已装的包**：clusterProfiler 4.18.2、org.Hs.eg.db、enrichplot、DOSE、
  ggprism 1.0.7、**gground 1.0.1**、tidyverse、fgsea。
  `gground` 是 GitHub 包（`devtools::install_github("dxsbiocc/gground")`），
  装好后无需重装。小鼠的 `org.Mm.eg.db` 见 `references/pitfalls.md` 5.2。

### 怎么启动（★ 首选：直接调 Rscript）

某些沙箱环境会拦掉 bash 嵌套调用，此时**直接用绝对路径调 `Rscript`**，在同一条命令里先把环境变量设好：

```bash
export TMP=/tmp/rtmp TEMP=/tmp/rtmp TMPDIR=/tmp/rtmp   # Windows: 指向纯 ASCII 目录如 C:/Rtmp
export COMSPEC='C:\WINDOWS\system32\cmd.exe'
unset LC_ALL LANG LC_CTYPE LC_COLLATE LC_TIME
SK="<TOOLKIT>/skills/enrichment-analysis/scripts"
export ENRICH_SKILL_SCRIPTS="$SK"

Rscript --vanilla "$SK/01_enrich_ora.R" <参数...>
```

这五件事**必须在 R 启动前办好**，R 起来之后在脚本里补救不了：

| # | 变量 | 为什么 |
|---|---|---|
| 1 | `TMP/TEMP/TMPDIR` → 纯 ASCII 目录 | `tempdir()` 在 R 启动时固定；中文用户名会让它变成坏路径 |
| 2 | `COMSPEC` | 不设则 R 的 `shell()` 报 `'/c' not found` |
| 3 | **`unset LC_ALL LANG LC_*`** | 环境里的 `LC_ALL=C.UTF-8` 会让 Windows 版 R 退化成 `C` locale，此时**含非 ASCII 字符的路径完全无法寻址**，`file.exists` 静默失败 |
| 4 | `ENRICH_SKILL_SCRIPTS` → scripts 目录 | 三个脚本靠它找 `lib_enrich_common.R` |
| 5 | 路径写 `D:/...`，别写 `/d/...` | POSIX 路径交给 `Rscript.exe` 会段错误（139、无产物无报错） |

**bash 可用时的便捷入口**是 `run_enrich.sh`，它把这五件事都包好了
（还会自愈 PATH，防精简 Git 环境缺 coreutils）：

```bash
bash "$SK/run_enrich.sh" "$SK/01_enrich_ora.R" --deg=... --species=human
```

- 已知无害现象：R 退出时报 `Segmentation fault`（139），发生在工作全部完成之后。
  **判断成败看产物文件，不要看退出码。**

## 第一步：ORA（GO / KEGG / MSigDB）

```bash
SK="<TOOLKIT>/skills/enrichment-analysis/scripts"
bash "$SK/run_enrich.sh" "$SK/01_enrich_ora.R" \
  --deg="<差异分析表>" --ora_p_type=P.Value --ora_logfc=1 \
  --species=human --outdir=./enrich --prefix=<前缀> --gmt_sets=H
```

### 输入怎么给（两种，可以同时给）

| 给法 | 含义 | 什么时候用 |
|---|---|---|
| `--gene_list=<文件>` | 直接给基因列表（每行一个，逗号/空白/分号分隔都行） | 用户手上就是一张基因清单；结果那一套叫 `list`（可用 `--list_name` 改名） |
| `--deg=<差异分析表>` | 按 `--ora_p` / `--ora_p_type` / `--ora_logfc` 筛出显著基因，再拆成 `all` / `up` / `down` | 输入是 GEO 差异分析的输出表 |

> ⚠️ **最容易犯的错：把"全部基因"的差异表直接丢给 ORA。**
> 实例：`GSE62452_T_all.csv` 有 23307 行，`adj.P.Val<0.05` 就有 **10960** 个基因；
> 拿 1 万个基因做 ORA，几乎每条通路都显著，结果没有信息量。
> 脚本会在基因集 > 2000 时打印醒目告警。正确做法是给**已经筛过的**基因：
>
> - 直接用 GEO 技能的 `<前缀>_DEG.csv`（那已经是 `|logFC|>1 & P<0.05` 的 292 个基因），或
> - 明确加 `--ora_logfc=1`（该数据集口径：`|logFC|>1 & P.Value<0.05` → 上调 177 / 下调 115）
>
> **`--ora_logfc` 对 `all` 也生效**，所以 `all` 恒等于 `up ∪ down`，三套结果的基因数对得上。

### 常用参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `--species=human\|mouse` | human | 也接受 `hsa`/`mmu`/`人`/`小鼠` |
| `--geo_dir=<目录>` | 无 | **★ 从 geo 差异分析直接接手**（详见上文「从 GEO 差异分析直接接手」） |
| `--geo_prefix=<前缀>` | 自动 | geo 的 `--prefix`；目录里有多个数据集时必须给 |
| `--geo_manifest=<文件>` | 自动 | 直接给 `<前缀>_for_enrichment.txt` |
| `--geo_ora_source=deg\|all` | deg | ORA 吃 `_DEG.csv` 还是 `_all.csv`（后者会按阈值重新筛） |
| `--split_by=logfc\|change` | geo 模式为 change | 上下调按 logFC 拆还是按 geo 的 `change` 列拆 |
| `--orgdb_sqlite=<路径>` | 无 | 没装注释包时直接给 OrgDb 的 `.sqlite` |
| `--id_type=SYMBOL` | SYMBOL | 表里是 ENTREZID / ENSEMBL 就改这里（`--id_type=ENTREZID` 时跳过 ID 转换） |
| `--gene_col` / `--logfc_col` / `--p_col` | 自动识别 | GEO 技能产出的表列名是 `logFC` / `P.Value` / `adj.P.Val`，行名在**无名第一列**（脚本已能识别为 `rownames1`） |
| `--ora_p_type=auto` | auto | auto = `adj.P.Val` 优先，其次 `P.Value`；也可写 `P.Value` |
| `--ont=ALL` | ALL | GO 本体：`ALL` / `BP` / `CC` / `MF` |
| `--kegg=auto` | auto | `auto`（在线失败退 GMT）/ `online` / `gmt` / `none` |
| `--gmt_sets=H,C2` | 无 | 要跑的 MSigDB 集合，见下表 |
| `--universe=none` | none | `none` = clusterProfiler 默认（与参考代码一致）；`deg` = 用输入表里**全部**基因作背景（更严格） |
| `--direction=all,up,down` | 三套 | 要跑哪几套 |
| `--pvalue/--padj/--qvalue` | 0.05 三个 | 只决定结果表里 `sig` 列的真假，**不影响计算**（计算恒用 cutoff=1，全部条目都留下） |

> ⚠️ **小鼠（`--species=mouse`）先用前自检**：`requireNamespace("org.Mm.eg.db")`。
> Windows 上 Bioconductor 注释包可能没有二进制、源码安装要 Rtools，装不上时按
> `references/pitfalls.md` 5.2 直接 `loadDb()` 一个现成的 `org.Mm.eg.sqlite`。
> 用小鼠前先 `requireNamespace("org.Mm.eg.db")`，或按
> `references/pitfalls.md` 5.2 把 `org.Mm.eg.sqlite` 放到 `<orgdb sqlite 目录>/`。
> **别把"注释包没装"报成"没富集结果"。**

`--gmt_sets` 可用的名字（GMT 数据自动从 msigdb 技能目录找）：

| 名字 | 是什么 | 名字 | 是什么 |
|---|---|---|---|
| `H`/`HALLMARK` | 50 条标志性过程 | `GOBP`/`GOCC`/`GOMF` | C5 的 GO 三支 |
| `C1`…`C9` | 整个集合 | `IMMUNE`（C7） | 免疫签名 |
| `KEGG`/`REACTOME`/`WP`/`PID`/`BIOCARTA` | C2 里的对应来源 | `CELLTYPE`（C8） | 细胞类型标记 |
| `PATHWAY` 无此名，用 `H,C2` 组合 | | `PERTURB`（C9） | 扰动签名 |

> ⚠️ **MSigDB 的 GMT 只有人类 symbol 版**。小鼠样本要用 MSigDB，得自己先做同源转换
> （小鼠用 GO/KEGG 不受影响）。

### 产出

| 文件 | 内容 |
|---|---|
| `<前缀>_GO_<方向>.csv` | GO 富集**完整**表（所有被检验条目），列含 `ONTOLOGY/ID/Description/GeneRatio/pvalue/p.adjust/qvalue/geneID/geneName/Count/sig` |
| `<前缀>_KEGG_<方向>.csv` | KEGG 同上 |
| `<前缀>_GMT_<方向>.csv` | MSigDB 同上，多一列 `collection` 标出集合名 |
| `<前缀>_<库>_<方向>_sig.csv` | 只含 `sig = TRUE` 的表（对齐原来的使用习惯） |
| `<前缀>_ORA_<方向>_input_genes.txt` / `_input_entrez.txt` | **本次实际用的基因清单**，可复现 |
| `<前缀>_id_map_<方向>.csv` | SYMBOL ↔ ENTREZID 转换明细（看丢了多少） |
| `<前缀>_ORA_summary.csv` / `.txt` | 各方向、各库的被检验条目数与显著条目数 |

> 表里 `geneID` 是 **ENTREZID**、`geneName` 是 **SYMBOL**（参考代码里 `setReadable`
> 之后只剩 symbol 一列，这里两列都留着）。

## 第二步：GSEA（全部基因，不筛）

```bash
bash "$SK/run_enrich.sh" "$SK/02_enrich_gsea.R" \
  --deg="<含全部基因的差异表>" --species=human \
  --outdir=./enrich --prefix=<前缀> --gmt_sets=H --plot=TRUE
```

- 排序列默认自动找 `logFC`；想改用 t 值排序就 `--stat_col=t`。
- **`--padj` 的 GSEA 惯例是 0.25**（不是 0.05），脚本默认就是 0.25。
- 重复基因默认保留 `|logFC|` 最大的一条。
- `--seed=1234` 已固定（fgsea 的 P 值来自随机置换，不固定就不可复现）；
  `--nperm=10000` 控制简化置换次数。
- 产出：`<前缀>_GSEA_{GO,KEGG,GMT}.csv` + `_sig.csv`、`_GSEA_ranked_list.txt`
  （**实际用的排序向量**，一行 `基因<TAB>排序值`）、`_GSEA_summary.csv`、
  `_GSEA_top_NES.pdf/png`（显著性最强的基因集 NES 条形图，红=激活、蓝=抑制）。
- `--curve=5` 可额外给前 5 个基因集画 GSEA 富集曲线（`--curve_which=pos|neg|both`）。
- **★ GSEA 的展示就到 02 步为止**：NES 条形图 + 富集曲线就是它的发表级出图，
  **不要**把 GSEA 结果塞进 03 的美化主图（03 检测到 GSEA 表会直接报错拒绝）。
  GSEA 显著通路动辄几百上千条，也没有"挑 15 条画主图"的需求；
  要引用具体通路时引用 `_GSEA_*_sig.csv` 表格即可。

## 第三步：归纳（统计 + 功能归类）★ 跑完必须停下问用户

```bash
SK="<TOOLKIT>/skills/enrichment-analysis/scripts"
bash "$SK/run_enrich.sh" "$SK/04_enrich_summary.R" \
  --in_dir=./enrich --prefix=<前缀> --direction=all,up,down --p_type=padj --p=0.05
```

它把**全部**通过阈值的通路（不是前 10 条）做四件事：

1. **筛选**：默认 `p.adjust ≤ 0.05`；想按你原话的原始 p 就用 `--p_type=p`
   （脚本会同时把两种口径的条数都打出来，方便对照）
2. **归纳**：按 `references/functional-rules.tsv` 把每条通路归到一个功能大类
   （**84 条规则 / 20 个功能大类**，按顺序先匹配到的优先）
3. **去冗余**：同一个功能大类内部，按基因重叠（Jaccard ≥ `--jaccard`，默认 0.5）
   折叠高度重复的通路，留下**代表通路**（例：42 条 ECM 相关通路 → 24 条代表）
4. **按通路数量排名**，报告前 `--top_cats`（默认 **20**，正好覆盖全部定义大类）个功能方向，
   并把这份排名**单独落成 CSV 表格**（`<前缀>_func_top20_<方向>.csv`）供用户查阅

产出：

| 文件 | 用途 |
|---|---|
| **`<前缀>_func_report.txt`** | **交给用户的报告**：每个方向前 20 个功能方向 + 每个方向的代表通路 + 未归类清单 + 下面两个选项 |
| **`<前缀>_func_top20_<方向>.csv`** | **排名前 20 大类清单（供用户查阅）**：rank/通路数/代表数/最小 p/库构成；「未归类」不参与排名；实际命中不足 20 个时全列 |
| `<前缀>_func_pathways_<方向>.csv` | 逐条通路的归类结果（含 `func_cat`/`func_sub`/`representative`），**第四步 B 方案的输入** |
| `<前缀>_func_counts_<方向>.csv` | 功能方向汇总（通路数、代表数、最小 p、库构成） |
| `<前缀>_func_counts_by_direction.csv` | 所有方向合并的汇总表 |
| `<前缀>_func_rules_used.tsv` | 本次实际生效的规则表（便于审计与复现） |

> ⚠️ 合并汇总表叫 `_func_counts_by_direction.csv`，**不是** `_func_counts_all.csv` ——
> 后者是「方向 = all」那张逐方向表，两者同名会互相覆盖（2026-09-18 修）。

**跑完把报告摆给用户，明确问他：**

> 全部 p<0.05 的通路归纳下来是这几个方向（前 20 大类清单见 `*_func_top20_*.csv`，
> 报告里带通路数/代表数，这里摘前几个）：
> 1. 细胞外基质与黏附 44 条（代表 26）
> 2. 代谢与能量 33 条（代表 17）…
>
> 你想怎么出图？
> **A) 均衡配额（默认）**：固定 20 条，BP/CC/MF/KEGG 各 5 条，某板块不够时
>    富余名额自动匀给其他板块
> **B) 只看某个功能方向**（比如「细胞外基质与黏附」，我会自动只画该方向的代表通路，
>    不重复；受 `--max_terms=15` 上限约束）

**在用户回答之前不要出图。**

- 「未归类」**不参与排名**（它不是功能方向），脚本会单列出来。实测 GSE62452
  三个方向补完规则后**未归类已清零**（0 / 0 / 0 条）。
- **规则表按顺序匹配**，所以"细胞外基质"必须排在"信号转导"前面
  （否则 `ECM-receptor interaction` 会被 `receptor` 抢走）；
  "骨骼、牙齿与矿化"要排在"发育与分化"前面（否则 `bone/cartilage development`
  会被 `development` 抢走）。

### 20 个功能大类（要跟用户介绍可选方向时照这个念）

细胞死亡与自噬 / 细胞周期与增殖 / 免疫与炎症 / 血液与凝血 / 细胞外基质与黏附 /
细胞骨架与运动 / 代谢与能量 / 酶活性与催化 / 翻译与蛋白稳态 / 转录与表观遗传 /
信号转导 / 物质运输与离子 / 应激与损伤修复 / **骨骼、牙齿与矿化** / 发育与分化 /
神经与感觉 / 内分泌与激素 / 细胞器与亚细胞定位 / **稳态与内环境** / 疾病与感染
（另有一个 `未归类` 兜底，不参与排名）

加粗的两个是 2026-09-18 为把未归类压到 0 而新增的：
- **骨骼、牙齿与矿化**：bone / ossification / osteo / chondro / cartilage / skelet /
  odontogenesis / tooth / mineralization / bone resorption
- **稳态与内环境**：兜住所有前面分类没接住的 `homeostasis` 泛指项
  （`tissue homeostasis`、`anatomical structure homeostasis` 等）

## 第四步：按用户的选择出美化主图

### A) 均衡配额（默认；用户说"就按默认出"或直接要图）

```bash
bash "$SK/run_enrich.sh" "$SK/03_enrich_plot.R" \
  --in_dir=./enrich --prefix=<前缀> --direction=up --db=GO,KEGG --mode=even
```

`--mode=even`（2026-09-25 起为默认，可省略）= 固定 `--even_total=20` 条，
BP/CC/MF/KEGG 各 5 条；某板块不足配额时富余名额按 BP→CC→MF→KEGG 顺序轮转
匀给还有存货的板块（日志会打印配额/库存/实际三行，逐项核对）。配额在去重后
的池子上算，`--top`、`--max_terms` 在此模式不生效。

> 旧的按排名挑法仍在：`--mode=global --top=10`（全库排名前 10，默认的
> 撞车去重可能把条数压下去，要足额加 `--dedup_padj=FALSE`）。

### B) 只看某个功能方向（用户说"看 ECM 那一块"）

```bash
bash "$SK/run_enrich.sh" "$SK/03_enrich_plot.R" \
  --in_dir=./enrich --prefix=<前缀> --direction=up --db=GO,KEGG,GMT \
  --func_csv=./enrich/<前缀>_func_pathways_up.csv \
  --category="细胞外基质与黏附" --top=0
```

- `--func_csv` + `--category=<功能大类>`：只画该方向的通路
- ⚠ **2026-09-25 起 even 是默认模式，B 方案同样固定 20 条均衡配额**（在该方向
  的代表通路里按 BP/CC/MF/KEGG 各 5 条挑，不足匀给；`--top`/`--max_terms` 不生效）。
  想保持旧行为（该方向代表通路全画、`--max_terms=15` 截断）就显式加
  `--mode=per_ontology --top=0`；确认要全画才给 `--max_terms=0`
- `--use_representative=TRUE`（默认）只画 04 步挑出的代表通路
- 用了 `--func_csv` 时默认**关掉** `--dedup_padj`（04 已按基因重叠挑过代表，
  再按 p 值去重会误删不同的通路；even 模式下同理，配额在未去重池子上算）
- 想不出现在可用的大类名，脚本在筛空时会直接列出全部可选值

### 版式参数（两种方式通用）

| 参数 | 默认 | 说明 |
|---|---|---|
| `--direction` | all | 画哪一套：`all`/`up`/`down`/`list` |
| `--db=GO,KEGG` | GO,KEGG | 想叠 MSigDB 就写 `--db=GO,KEGG,GMT`（分类会多出 H/C2… 各自一块） |
| **`--mode=even`** | even | **`even` 均衡配额（默认）**：总条数固定 `--even_total`，各分类平均分、不足时轮转匀给；`per_ontology` 每个分类各取前 N 条（参考代码行为）；`global` 全库取前 N 条 |
| **`--even_total=20`** | 20 | even 模式整张图的总条数；某分类不足配额时富余名额自动匀给其他分类 |
| `--top=5` | 5 | 挑几条（仅 per_ontology/global 生效）；**`0` = 不限制**（配合 `--category` 用）；even 模式下忽略 |
| **`--max_terms=15`** | 15 | **仅 per_ontology/global 生效的硬上限**：挑条 + 去重后按 p.adjust 截断并打印告警；`0` = 关闭上限。even 模式下忽略（总条数由 `--even_total` 决定） |
| `--func_csv` `--category` | 无 | 按 04 步的功能归类表筛选；`--category` 给功能大类名 |
| `--use_representative=TRUE` | TRUE | 用了 `--func_csv` 时只画代表通路（04 已按基因重叠去冗余） |
| `--fig_filter=padj` `--fig_thr=0.05` | padj ≤ 0.05 | 过滤口径；**某类被筛空时会自动退化为"用该类全部条目"并打印告警**，不会静默出空图 |
| `--pal_idx=3` | 3 | 三套内置配色（3 = 参考代码最终使用的那套） |
| `--dedup_padj=TRUE` | TRUE | 同一 p.adjust 只留 Count 最多的一条（参考代码的行为）；用了 `--func_csv` 时默认关闭 |
| **`--gene_pos=below`** | below | 基因名放彩色柱**下方**（默认，不压柱子）；`inside` = 参考代码原样（压在柱子里） |
| **`--bar_height=0.6`** | 0.6 | 彩色柱在 y 方向的高度（参考代码写死 0.6）。基因名的位置**跟着它走** |
| **`--gene_dy`** | 自动 | 基因名的 y 偏移（行单位）；不给就按柱高自动算。想手动微调再给 |
| **`--cat_label_chars=0`** | 自动 | 分类标签（BP/CC/MF/KEGG）每行最多几个字符。自动时 "KEGG" 会折成两行，色块往左加宽到刚好包住 |
| `--cat_size=3.88` | 3.88 | 分类标签字号 |
| `--max_genes=0` | 0 | 每条通路最多画几个基因名，0 = 全列（基因名很长时可设 20 防止冲出画布） |
| `--auto_gene_size=TRUE` | TRUE | 按最长基因名字符数自动缩字号（默认）；想改加宽画布用 `--auto_width=TRUE` |
| `--legend=right\|bottom\|none` | right | 图例位置 |
| `--width/--height/--format/--dpi` | 12/8/pdf,png/600 | **高度双向自适应**：条目多自动加高、条目少自动收矮（4 行 → 4.6 英寸）。显式写 `--height=` 就完全按你给的来 |

> **基因名会很长**（一条通路富集到 40 个基因时字符串可能 260+ 字符）。
> 默认策略是**自动缩基因名字号**把它塞进画布（实测 261 字符时 3.5 → 2.27），
> 好处是不裁切、不丢内容。另外两个选项：
> `--auto_width=TRUE`（加宽画布、保持字号，可能到 17 英寸宽）、
> `--max_genes=20`（截断并标注 `... (+n)`）。
>
> **基因名位置**：默认贴在彩色柱**下沿**（`--gene_pos=below`），位置由
> `--bar_height` 推出来，所以调柱高时文字会跟着柱子动；参考代码是把基因名
> 压在柱子**里面**（`vjust = 2.6`，位置只跟字号有关），要那个效果写
> `--gene_pos=inside`。
>
> **左侧分类标签**：色块宽度是"数据单位"，x 轴越长色块越窄，`KEGG` 会被裁成
> `EGG`。脚本会按色块能放几个字符**把标签折行**（KEGG → `KE`/`GG` 两行），
> 并把色块**往左**加宽到刚好包住文字（右边界不动，所以基因数圆点和 x=0 都不移位）。

产出：`<前缀>_enrich_<方向>.pdf` / `.png`、`<前缀>_enrich_<方向>_selected_terms.csv`
（**入选条目清单，先核这张表**）、`_plot.log`。
用 `--category` 时文件名会自动带上功能方向（`<前缀>_enrich_<方向>_细胞外基质与黏附.pdf`），
用 `--mode=global` 时带 `top<N>`。

**配色与参考代码的差别（重要）**：参考代码用裸向量 `pal` 配合
`levels = rev(c('BP','CC','MF','KEGG'))`，四类齐全时
`KEGG=pal[1], MF=pal[2], CC=pal[3], BP=pal[4]`；**但一旦某一类没有条目，
颜色会整体平移一格**（BP 会拿到 KEGG 的颜色）。本技能把颜色**按分类名绑定**，
所以 `all`/`up`/`down` 三张图的同一分类永远是同一个颜色，也便于并排比较。
四类齐全时与参考代码完全一致。

## 关键决策点（不要替用户默默决定）

0. **★ 先看上游有没有 GEO 差异分析产物** —— 用户说"对 GSE62452 做富集"时，
   **先去 geo 的产物目录找 `<前缀>_for_enrichment.txt`**，用 `--geo_dir` 接手，
   不要手工拼 `--deg`/列名/阈值这一串（容易把 geo 的 `P.Value` 口径换成
   `adj.P.Val`，基因数直接换一套）。找不到清单再退回手工路径。
1. **这个任务是 ORA 还是 GSEA** —— 看输入是"筛过的基因"还是"全部基因"。
2. **ORA 用哪一批基因** —— 必须让用户确认显著性口径（`p` / `adj.P.Val`）和
   要不要 `|logFC|` 阈值。**默认 `--ora_logfc=0` 会得到上千个基因，不是他要的。**
3. **物种 / ID 类型** —— SYMBOL 还是 ENTREZID，人还是鼠。
4. **要不要 MSigDB 集合** —— 要的话先问是哪几个（H 最通用，C2 看具体通路，
   C7/C8 看免疫/细胞类型）。**不要一次全选 10 个集合**，C5 一个就有 1.6 万个集。
5. **GSEA 的排序依据** —— 默认 logFC；用 t 值排要显式指定。
6. **★ 出图方式（归纳之后必问）** —— 把 04 步报告和前 20 大类表格（`*_func_top20_*.csv`）
   摆给用户，让他二选一：**A) 均衡配额 20 条**（默认，BP/CC/MF/KEGG 各 5、不足匀给）
   还是 **B) 只看某个功能方向**（受 `--max_terms=15` 上限约束）。
   若选 B，**还要等他指定是哪一个功能方向**。这一问不能省，也不要替他选。
7. **GSEA 结果怎么展示** —— 单独用 02 步的 NES 条形图（`_GSEA_top_NES.pdf/png`），
   **不进 03 的美化主图**。用户要求"GSEA 也画成美化图"时，先解释两种口径
   不能混图，再给 NES 条形图方案。
7. **画哪几个方向** —— `all`/`up`/`down` 各一张（三张），还是只画一张。
8. **功能归类规则表要不要改** —— 04 报告里「未归类」条数偏多时，问用户是补规则
   （编辑 `references/functional-rules.tsv`）还是先照现有归类出图。

## 坑（详见 `references/pitfalls.md`）

- **`adj.P.Val` 与 `P.Value` 差一个数量级**：GSE62452 上 `P.Value<0.05` 有 12473 个基因，
  `adj.P.Val<0.05` 有 10960 个，而 `|logFC|>1 & P<0.05` 只有 292 个。
  口径写错，富集结果完全换一套。
- **功能归类规则表里的短词必须加词边界**：`actin` 会命中 `acting`（把
  "oxidoreductase activity, acting on NAD(P)H" 归成细胞骨架）、`ral` 会命中
  `endochondral`。有风险的短词一律写成 `\bactin\b` / `\bral[ab]?\b`（R 端用
  `perl = TRUE` 的 `grepl`）。
- **规则表按顺序先匹配到的优先**，所以具体规则要排在宽泛规则前面
  （`细胞外基质` 要在 `信号转导` 之前，否则 `ECM-receptor interaction` 被 `receptor` 抢走）。
- **`enrichKEGG` 需要联网**（KEGG REST）。失败时 `--kegg=auto` 会自动退到本地
  MSigDB 的 KEGG 子集，但那是 `KEGG_LEGACY`（186 条），
  **不含 `KEGG_PI3K_AKT_SIGNALING_PATHWAY`、`KEGG_TNF_SIGNALING_PATHWAY`**，
  通路 ID 也不是 `hsaXXXXX`。退化了要在回复里说明。
- **KEGG 每跑一个方向都要下载一次**，三个方向 ~4 分钟。已经跑过就别重跑，
  改阈值用 `03_enrich_plot.R` 从 CSV 重画即可（富集值不随阈值改变）。
- **`factor(Description, levels = Description)` 会崩**：同一 Description 在 BP/CC/MF
  重复出现时 R 4.x 直接报 `duplicated levels`。脚本已改成 `unique()`。
- **`--mode=global --top=10` 不一定真出 10 条**：默认的撞车去重会把 p.adjust 完全
  相同的通路合成一条（实测 10 → 6）。要足额就加 `--dedup_padj=FALSE`。
  （even 模式没这个问题：去重先做，配额在去重后的池子上算，总数稳定。）
- **条数上限现在分两套（2026-09-25 起）**：even 模式（默认）总条数由
  `--even_total=20` 决定、各分类均衡分配，`--top`/`--max_terms` 不生效；
  per_ontology/global 模式仍是 `--top`（每类/全库挑几条）+ `--max_terms=15`
  （整图硬上限）两级限制。用户说"条数多点没关系"时才允许 `--max_terms=0`，
  且要在回复里说明 SCI 投稿不推荐。
- **GSEA 表不能喂给 03**：03 只吃 ORA 产物；喂进 GSEA 表（含 NES/setSize 列）会
  报错退出而不是画出错误口径的图。GSEA 展示走 02 的 NES 条形图。
- 图内文字一律 ASCII（Windows 下中文字形缺失会渲染成乱码 ASCII）。

## 回归基线（数字对不上说明中间出了错）

**均衡配额 + 前 20 大类表格**（`<回归测试目录>/enrich` → 产物落
`<回归测试目录>/_regress_even/`，GSE62452，2026-09-25 实测）：

- `03 --direction=up --db=GO,KEGG`（even 默认）：库存 BP 92 / CC 24 / MF 32 / KEGG 10，
  配额 5/5/5/5 → **正好 20 条，各分类 5 条** ✓
- `--even_total=60`（每类配额 15，KEGG 库存 10 不足）：告警"KEGG 只有 10 条"→
  富余 5 个名额轮转匀给 BP+2 / CC+2 / MF+1 → **BP 17 / CC 17 / MF 16 / KEGG 10，合计 60** ✓
- 极端 `--fig_thr=1e-10`（全库存只剩 5 条）：打印"不足 --even_total=20，有多少画多少" → 5 条 ✓
- `--mode=global --top=10`（旧挑法，up）：去重后 6 条，行为不变 ✓
- `04`（all,up,down）：新写出 `GSE62452_T_func_top20_{all,up,down}.csv`
  （all/up 各 17 行 + 表头、down 11 行 + 表头——实际命中不足 20 个就全列；
  含 rank 列，「未归类」不进表），报告排名段完整列出全部大类 ✓
- 归纳数字与 2026-09-18 基线一致：all 228→132、up 214→121、down 85→40 ✓

**旧行为存档（even 模式之前，2026-09-24 实测）**：

- `--mode=global --top=10`：去重后 8 条 ≤ 15 ✓（与改动前一致，行为不变）
- ECM 方向 `--top=0`（30 条代表通路）：打印"超过上限"告警 → 截到 **15 条**
  （KEGG 1 / MF 4 / CC 4 / BP 6），产物 `GSE62452_T_enrich_up_细胞外基质与黏附.*`
- `--max_terms=0`：显式关闭上限，发育与分化方向出 18 条（行为与改动前一致）✓
- GSEA 表复制成 ORA 命名喂 03：**报错退出**，提示"GSEA 单独用 02 的 NES 条形图" ✓

**与 geo 技能联通**（`--geo_dir`，GSE62452，2026-09-18 实测）：
geo 02 跑出 `177 up / 115 down`（与 geo 基线一致）→ 富集侧只给 `--geo_dir`，
结论与手工给 `--deg=…_all.csv --ora_p_type=P.Value --ora_logfc=1` **完全一致**：

- ORA（吃 `_DEG.csv`，`split_by=change`）：`all 292 = up 177 + down 115` ✓；
  GO 212/196/76、KEGG 14/12/7、Hallmark 2/6/2
- GSEA（吃 `_all.csv`）：GO 8228/2567、KEGG 352/157、H 50/39
- 04 归纳：与手工路径同样得到 all 228→133、up 214→126、down 85→41

**04 步的功能归纳**（`--p_type=padj --p=0.05`，84 条规则 / 20 个大类；实测值）：

| 方向 | 显著通路 | 去冗余代表 | 前 8 个功能方向（通路数 / 代表数） |
|---|---|---|---|
| all | 228（GO 212 / KEGG 14 / GMT 2） | 132 | 细胞外基质与黏附 49/31、代谢与能量 41/20、酶活性与催化 25/12、发育与分化 17/12、应激与损伤修复 16/6、物质运输与离子 14/7、血液与凝血 14/5、信号转导 11/7 |
| up | 214（GO 196 / KEGG 12 / GMT 6） | 121 | 细胞外基质与黏附 51/30、发育与分化 35/18、信号转导 17/10、应激与损伤修复 16/7、酶活性与催化 16/6、细胞骨架与运动 13/8、血液与凝血 13/3、骨骼、牙齿与矿化 11/6 |
| down | 85（GO 76 / KEGG 7 / GMT 2） | 40 | 代谢与能量 35/14、酶活性与催化 16/6、物质运输与离子 9/4、细胞器与亚细胞定位 6/2、内分泌与激素 4/4、应激与损伤修复 4/3、血液与凝血 4/2、免疫与炎症 3/1 |

- **未归类 0 / 0 / 0 条**（补规则前是 19 / 32 / 4）
- 命中大类数：all 16 个 / up 17 个 / down 11 个（总定义 20 个）
- 同阈值下按**原始 pvalue** 算的条数：all 728 / up 655 / down 302（`--p_type=p` 可切过去）
- 耗时 < 1 秒（纯统计，不联网）
- 生物学上说得通：up 是 ECM/黏附 + 发育（EMT）+ 信号转导 + 骨骼矿化，
  down 是代谢/消化吸收 —— 胰腺癌 vs 癌旁，肿瘤侧基质重塑与增殖，
  正常侧腺泡消化功能

**GSE62452_T**（人，GPL6244，130 样本，case 69 / control 61；
输入 `<本地数据目录>/GSE62452_T_all.csv`，23307 基因；
`--ora_logfc=1 --ora_p_type=P.Value` → 292 个基因，与 GEO 技能的
`GSE62452_T_DEG.csv` 完全一致）。实测 2026-09-18：

| 方向 | 输入基因 | 转成 ENTREZID | GO 被检验 / 显著 | KEGG 被检验 / 显著 | MSigDB H 被检验 / 显著 |
|---|---|---|---|---|---|
| all | **292** | 285（丢 7，2.4%） | 3881 / **212** | 226 / **14** | 45 / **2** |
| up | **177** | 173 | 3080 / **196** | 174 / **12** | 39 / **6** |
| down | **115** | 112 | 2056 / **76** | 154 / **7** | 28 / **2** |

- `all = up + down = 292` ✓（`--ora_logfc` 对 all 也生效，脚本会自己核对并打印 ✓/⚠）
- 耗时约 137 秒（三个方向 × GO/KEGG/GMT，KEGG 联网下载占大头）
- 美化图入选条目数（旧行为存档）：all 16 条、up 15 条、down 16 条（去重后每类不足
  5 条是正常的；2026-09-24 起 `--max_terms=15` 生效后 all/down 会截到 15）。
  **2026-09-25 起 even 模式为默认：固定 20 条、各分类均衡分配**（见上方回归基线）

**GSEA**（同一份全表，23307 个基因按 logFC 排序，**不筛选**）：

| 库 | 被检验 | 显著（p<0.05 且 padj<0.25） | 激活 / 抑制 |
|---|---|---|---|
| GO（BP+CC+MF） | 8228 | 2567 | 2170 / 397 |
| KEGG | 352 | 157 | 120 / 37 |
| MSigDB H | 50 | 39 | 31 / 8 |

|NES| 最大的前几条（生物学上说得通：胰腺癌 EMT + 增殖 + 干扰素）：
`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`、`INTERFERON_GAMMA_RESPONSE`、
`INTERFERON_ALPHA_RESPONSE`、`APICAL_JUNCTION`、`E2F_TARGETS`、
`KRAS_SIGNALING_UP`、`MITOTIC_SPINDLE`、`G2M_CHECKPOINT`（均激活）；
`PANCREAS_BETA_CELLS`（抑制）。

耗时约 32 秒（`--go=FALSE --kegg=none` 只跑 MSigDB 时）。

## 资源

- `scripts/run_enrich.sh` —— 环境包装器（bash 可用时用它，被沙箱拦就直接调 Rscript）
- `scripts/01_enrich_ora.R` —— ORA：GO / KEGG / MSigDB，all+up+down 三套
- `scripts/02_enrich_gsea.R` —— GSEA：GO / KEGG / MSigDB，全部基因排序
- `scripts/04_enrich_summary.R` —— **归纳**：全部显著通路 → 功能大类 + 去冗余代表 + 排名报告
- `scripts/03_enrich_plot.R` —— 美化主图（gground + ggprism 版式；`--mode` / `--category`）
- `scripts/lib_enrich_common.R` —— 公共库（参数、物种、读表、GMT、配色、归类规则、Jaccard 去冗余）
- `references/functional-rules.tsv` —— **功能分类规则表（可编辑）**：84 条规则 / **20 个功能大类**，
  按顺序先匹配到的优先。归错或未归类多就直接改这个文件
- `references/figure-style.md` —— 版式与配色的逐项对照（含三套配色的色值、与参考代码的差异清单）
- `references/pitfalls.md` —— 踩坑手册。**改动脚本前先读**
