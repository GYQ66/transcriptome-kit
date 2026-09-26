#!/usr/bin/env Rscript
# ===========================================================================
# 01_string_network.R —— STRING 蛋白互作（PPI）网络：查询 -> 建图 -> 可视化
#
# 功能等价于参考代码（STRINGdb$map -> get_interactions -> igraph -> ggraph），
# 取数方式做了参数化，三选一：
#   --engine=stringdb  用 Bioconductor 的 STRINGdb 包（与参考代码完全一致）
#   --engine=api       直接打 STRING REST API（不需要 STRINGdb 包，推荐）
#   --engine=local     读本地 STRING 数据文件（完全离线）
#   --engine=auto      先 stringdb，失败则 api（默认）
#
# 用法见 --help。图内文字一律 ASCII（Windows 缺中文字形会渲染成乱码）。
# ===========================================================================

SCRIPT_VERSION <- "1.2"

# --- 定位公共库 -----------------------------------------------------------
skill_dir <- Sys.getenv("STRING_SKILL_DIR")
if (!nzchar(skill_dir)) {
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  this <- if (length(a)) tryCatch(normalizePath(sub("^--file=", "", a[1])), error = function(e) "") else ""
  skill_dir <- if (nzchar(this)) dirname(dirname(this)) else ""
}
lib <- file.path(skill_dir, "scripts", "lib_string_common.R")
if (!file.exists(lib)) {
  alt <- file.path(dirname(normalizePath(sub("^--file=", "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))),
        "lib_string_common.R")
  if (file.exists(alt)) lib <- alt else stop("找不到 lib_string_common.R: ", lib)
}
source(lib)

# ===========================================================================
# 参数
# ===========================================================================
DEFAULTS <- list(
  # 输入（三选一）
  genes = "", genes_file = "", table = "",
  gene_col = "", logfc_col = "",
  logfc = "", logfc_file = "",
  # STRING 查询
  species = "9606", version = "11.5", score_threshold = "400",
  network_type = "functional",          # functional | physical
  add_nodes = "0",                      # STRING 自动补充的邻居数（0 = 只查输入基因之间）
  engine = "auto",                      # auto | stringdb | api | local
  api_base = "",                        # 留空则按 version 自动拼
  info_file = "", links_file = "",      # engine=local 用的本地 STRING 文件
  cache_dir = "",                       # 下载缓存目录
  timeout = "120", retries = "2",
  # 基因名处理
  alias_fix = "auto",                   # auto | all | none
  # 图的构建
  prune = "logfc",                      # logfc | none
  keep_isolated = "1",
  # 绘图
  layout = "stress", seed = "42",
  node_size = "6", edge_width_range = "0.1,1.2",
  edge_colour = "grey70", edge_alpha = "0.9",
  palette = "nature", colors = "",
  color_limits = "auto",                # auto 或 "min,max"（如 -2,0）
  color_label = "",                     # 图例标题，留空用 log2 fold-change
  node_colour = "#2171B5",              # 无 logFC 时的统一节点色
  label_size = "3", repel = "1",
  edge_label = "STRINGdb\nconfidence",
  legend = "bottom",                    # bottom | right | left | top | none
  # 输出
  outdir = "./string_out", prefix = "ppi",
  format = "pdf,png,svg",               # 白名单过滤后三格式始终同时输出
  width = "7", height = "7", dpi = "300",
  graphml = "0", log_file = "",
  # 交互式布局微调
  interactive = "1",                    # 默认 1 = 每次都生成 HTML 布局编辑器（0 关闭）
  positions = ""                        # 手动布局坐标 csv（编辑器导出的文件），给了就走 manual 布局
)

USAGE <- "
01_string_network.R —— STRING 蛋白互作网络（查询 / 建图 / 可视化）

【输入基因】三选一
  --genes=A,B,C                 内联基因列表
  --genes_file=genes.csv        基因文件（第一列或 --gene_col 指定列）
  --table=deg.csv               含基因列 + logFC 列的结果表（推荐）
  --gene_col=gene               --table/--genes_file 的基因列名（留空自动猜）
  --logfc_col=logFC             logFC 列名（留空自动猜）

【logFC / 节点着色】可选
  --logfc=-2.0,-1.5,...         与 --genes 一一对应的内联 logFC
  --logfc_file=fc.csv           单独的文件（行数须与基因数一致）
  不给 logFC 时：节点单色（--node_colour），网络结构照常出

【STRING 查询】
  --species=9606                物种 NCBI taxonomy ID（9606人 / 10090小鼠 / 10116大鼠）
  --version=11.5                STRING 版本（11.5 / 12.0）
  --score_threshold=400         置信度阈值 0-1000（400 中高可信；200 低可信）
  --network_type=functional     functional | physical
  --add_nodes=0                 让 STRING 额外补充 N 个邻居（0 = 不补）
  --engine=auto                 auto | stringdb | api | local
  --api_base=                   API 根地址（留空按 version 自动拼）
  --info_file= / --links_file=  engine=local 用的本地 STRING 数据文件(.txt.gz)
  --cache_dir=                  下载缓存目录（留空用 outdir/_cache）
  --timeout=120  --retries=2    单次请求超时(秒) / 失败重试次数

【基因名修正】
  --alias_fix=auto              auto：只给「映射不上」的基因换别名（CoREST->RCOR1、
                                LSD1->KDM1A、MACROH2A2->H2AFY2、H2AX->H2AFX ...）
                                all：无条件全表替换   none：不替换

【建图】
  --prune=logfc                 logfc：只保留有 logFC 的节点（等同参考代码）
                                none：保留全部映射成功的基因
  --keep_isolated=1             把没有任何边的输入基因也加进图
                                （参考代码会静默丢掉孤立节点）

【绘图】
  --layout=stress               stress | fr | kk | dh | circle | graphopt | lgl | mds
  --seed=42                     布局随机种子（保证可复现）
  --node_size=6                 节点大小
  --edge_width_range=0.1,1.2    线宽映射范围
  --edge_colour=grey70  --edge_alpha=0.9
  --palette=nature              nature(参考图配色) | viridis | plasma | rdbu |
                                bluered2 | spectral | bluyl | blues | reds | ...
  --colors=#08306B,#2171B5,...  自定义渐变色（逗号分隔），优先于 --palette
  --color_limits=auto           auto：按 logFC 实际范围；或写 -2,0 固定（参考代码用 -2,0）
  --color_label=                图例标题（默认 log2 fold-change）
  --node_colour=#2171B5         无 logFC 时的节点颜色
  --label_size=3  --repel=1     标签字号 / 是否防重叠
  --legend=bottom               bottom | right | left | top | none

【输出】
  --outdir=./string_out  --prefix=ppi
  --format=pdf,png,svg          三格式始终同时输出（白名单外的值告警并忽略，
                                绝不产出扩展名与内容不符的文件）
  --width=7 --height=7 --dpi=300
  --graphml=0                   1 = 额外导出 Cytoscape 可读的 .graphml
  --log_file=                   把日志同时写进该文件

【交互式布局微调（拖节点摆位）】默认开启
  每次运行都会生成 *_layout_editor.html（自包含、离线可用）+ 最终 pdf/png/svg。
  想按自己摆的布局出图：
  # 1) 浏览器打开 *_layout_editor.html，按住节点拖动，点「导出坐标 CSV」
  # 2) 原命令追加 --positions 重跑
  Rscript 01_string_network.R --table=deg.csv --outdir=./string_out \\
      --positions=./string_out/ppi_positions.csv

  --interactive=1               默认 1；0 = 不生成编辑器（纯一键出图）
  --positions=coords.csv        手动布局坐标（列: symbol 或 STRING_id + x + y）。
                                给了它就走 manual 布局；没覆盖到的节点保留自动布局

【示例】
  # 1) 复刻参考代码的 Nature 图（13 个基因 + 示例 logFC）
  Rscript 01_string_network.R --genes=RCOR3,CoREST,LSD1,RCOR2,MACROH2A2,H2AC7,HDAC2,MDC1,CENPC,UBTF,H2AX,CTCF,MIER3 \\
      --logfc=-2.0,-1.5,-1.4,-1.9,-1.8,-0.6,-0.5,-0.7,-0.4,-0.3,-0.6,-0.8,-0.9 \\
      --color_limits=-2,0 --outdir=./string_out --prefix=kbtbd4

  # 2) 从差异分析结果表出发（自动猜 gene / logFC 列）
  Rscript 01_string_network.R --table=deg.csv --score_threshold=400 \\
      --outdir=./string_out --prefix=deg --format=png

  # 3) 交互式摆位：默认每次都出编辑器；拖完节点导出坐标，--positions 重跑即按摆位出图
  Rscript 01_string_network.R --table=deg.csv --prefix=deg --outdir=./string_out
"

args <- commandArgs(trailingOnly = TRUE)
p <- parse_args(args, DEFAULTS)
if (isTRUE(p$help_flag) || length(args) == 0) {
  cat(USAGE); quit(save = "no", status = 0)
}
given_prune <- any(grepl("^--prune", args))

init_log(p$log_file)

# 拼错的参数名不要静默吞掉（--score_threhold=400 这种最坑）
unknown <- setdiff(names(p), c(names(DEFAULTS), "help_flag"))
if (length(unknown)) {
  msg("[warn] 无法识别的参数（已忽略）：", paste0("--", unknown, collapse = ", "))
  msg("       用 --help 看全部可用参数。")
}

# 布局编辑器（交互模式用）
source(file.path(skill_dir, "scripts", "lib_layout_editor.R"))

# 输出格式归一：pdf+png 始终输出（用户要求）；svg 可选加；非法值告警忽略
fmts <- normalize_formats(p$format)

# ===========================================================================
# 0. 基本检查
# ===========================================================================
rule("=")
msg("STRING PPI 网络  |  01_string_network.R v", SCRIPT_VERSION)
rule("=")

if (!requireNamespace("igraph", quietly = TRUE)) stop("缺少 igraph，请先 install.packages('igraph')")
if (!requireNamespace("ggraph", quietly = TRUE)) stop("缺少 ggraph，请先 install.packages('ggraph')")
if (as_bool(p$repel) && !requireNamespace("ggrepel", quietly = TRUE))
  msg("[warn] 未装 ggrepel，标签防重叠可能失效；install.packages('ggrepel')")

suppressPackageStartupMessages({
  library(igraph); library(ggraph); library(ggplot2)
})

SPECIES_MAP <- c(human = 9606, mouse = 10090, rat = 10116, zebrafish = 7955,
                 fly = 7227, worm = 6239, yeast = 559292)
sp_raw <- trimws(p$species)
species <- if (grepl("^[0-9]+$", sp_raw)) as.integer(sp_raw) else as.integer(SPECIES_MAP[tolower(sp_raw)])
if (is.na(species)) stop("无法识别 --species=", sp_raw, "（用数字 taxid，或 human/mouse/rat/...）")

version <- trimws(p$version)
score_threshold <- max(0, min(1000, as_int(p$score_threshold, 400)))
network_type <- if (tolower(trimws(p$network_type)) %in% c("physical", "p")) "physical" else "functional"
add_nodes <- max(0, as_int(p$add_nodes, 0))
timeout <- as_num(p$timeout, 120); retries <- as_int(p$retries, 2)

outdir <- if (nzchar(p$outdir)) p$outdir else "."
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
prefix <- if (nzchar(p$prefix)) p$prefix else "ppi"
cache_dir <- if (nzchar(p$cache_dir)) p$cache_dir else file.path(outdir, "_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

msg("物种 taxid        : ", species)
msg("STRING 版本       : ", version)
msg("置信度阈值        : ", score_threshold, " (0-1000)")
msg("网络类型          : ", network_type)
msg("输出目录          : ", normalizePath(outdir, mustWork = FALSE))

# ===========================================================================
# 1. 读基因表
# ===========================================================================
gt <- build_gene_table(p)
gt$gene <- trimws(gt$gene)
gt <- gt[!duplicated(gt$gene), , drop = FALSE]
msg("输入基因数        : ", nrow(gt))

alias_mode <- tolower(trimws(p$alias_fix))
if (identical(alias_mode, "all")) {
  a <- apply_alias(gt$gene, mode = "all")
  gt$gene <- a$genes
  if (length(a$changed)) {
    msg("[alias] --alias_fix=all，无条件替换 ", length(a$changed), " 个：")
    for (s in a$changed) msg("        ", s)
  }
}

# ===========================================================================
# 2. 取数引擎（三个）
#   每个引擎返回 list(name=, map=function(genes)->df, edges=function(ids)->raw df)
# ===========================================================================
empty_map <- function() data.frame(gene_symbol = character(0), STRING_id = character(0),
                                   stringsAsFactors = FALSE)

fetch_url <- function(url, dest, timeout, retries) {
  old <- getOption("timeout"); on.exit(options(timeout = old), add = TRUE)
  options(timeout = max(30, timeout))
  methods <- c("wininet", "libcurl", "curl", "wget", "internal", "auto")
  last_err <- ""
  for (i in seq_len(max(1, retries + 1))) {
    for (m in methods) {
      r <- tryCatch({
        utils::download.file(url, dest, quiet = TRUE, mode = "wb", method = m)
        file.exists(dest) && file.size(dest) > 0
      }, error = function(e) { last_err <<- conditionMessage(e); FALSE },
         warning = function(w) { last_err <<- conditionMessage(w); FALSE })
      if (isTRUE(r)) return(TRUE)
    }
    Sys.sleep(1.5 * i)
  }
  msg("[error] 下载失败: ", url, "\n        ", last_err)
  FALSE
}

read_tsv_strict <- function(path)
  utils::read.delim(path, sep = "\t", quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE,
                    fileEncoding = "UTF-8-BOM")

# ---- 引擎 A：STRINGdb（Bioconductor，与参考代码一致）----------------------
make_engine_stringdb <- function() {
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    msg("[engine] STRINGdb 未安装，跳过"); return(NULL)
  }
  msg("[engine] 尝试 STRINGdb ...（首次会下载 aliases/links 数据，几十 MB，慢）")

  # --- 缓存完整性预检 -----------------------------------------------------
  # STRINGdb 的 downloadAbsentFile() 只判断 `file.size > 0`：下载中断留下的残缺
  # .txt.gz 会被当成「已下载」直接使用，映射结果会静默变成垃圾（实测 13 个基因
  # 只映射上 1 个）。这里用 HTTP Content-Length 比对本地长度，发现残缺就拒绝启用。
  sb_url <- function(kind)
    sprintf("https://stringdb-downloads.org/download/protein.%s.v%s/%d.protein.%s.v%s.txt.gz",
            kind, version, species, kind, version)
  head_len <- function(u) {
    if (!requireNamespace("httr", quietly = TRUE)) return(NA_real_)
    h <- tryCatch(httr::HEAD(u, httr::timeout(30)), error = function(e) NULL)
    if (is.null(h)) return(NA_real_)
    cl <- httr::headers(h)[["content-length"]]
    if (is.null(cl)) NA_real_ else suppressWarnings(as.numeric(cl))
  }
  bad <- character(0)
  for (kind in c("aliases", "info", "links")) {
    f <- file.path(cache_dir, basename(sb_url(kind)))
    if (!file.exists(f)) next
    exp <- head_len(sb_url(kind)); act <- file.size(f)
    if (!is.na(exp) && act != exp)
      bad <- c(bad, sprintf("  %s  本地 %d 字节 / 应为 %.0f 字节", basename(f), act, exp))
  }
  if (length(bad)) {
    rule("!")
    msg("[error] STRINGdb 数据缓存残缺，拒绝启用该引擎（残缺文件会让映射静默变成垃圾）：")
    for (b in bad) msg(b)
    msg("  对策：删掉上面这些文件后重跑；或直接用 --engine=api（同一个 STRING ",
        version, " 库，结果同源）。")
    rule("!")
    return(NULL)
  }

  # STRINGdb 内部用 download.file 直接取数据，不吃 fetch_url 的超时设置。
  # 注意：这里不能用 on.exit 还原 —— map()/get_interactions() 的下载发生在本函数返回之后。
  options(timeout = max(timeout, 1800))
  # STRINGdb$new() 会先访问 https://string-db.org/api/tsv-no-header/version 做版本校验，
  # 本机该域名时通时断，所以这里重试几次（包里没有重试）。
  sdb <- NULL
  for (attempt in 1:3) {
    sdb <- tryCatch({
      suppressPackageStartupMessages(library(STRINGdb))
      suppressMessages(suppressWarnings(
        STRINGdb$new(version = version, species = species,
                     score_threshold = score_threshold, input_directory = cache_dir)))
    }, error = function(e) {
      msg("[warn] STRINGdb$new() 第 ", attempt, "/3 次失败：", conditionMessage(e)); NULL
    })
    if (!is.null(sdb)) break
    Sys.sleep(2 * attempt)
  }
  if (is.null(sdb)) {
    msg("[hint] STRINGdb$new() 需要访问 string-db.org 做版本校验，本机该域名常不可达；",
        "改用 --engine=api（同一个 STRING ", version, " 库，结果同源）。")
    return(NULL)
  }
  if (is.null(sdb)) return(NULL)
  list(
    name = "stringdb",
    map = function(genes) {
      if (!length(genes)) return(empty_map())
      df <- data.frame(gene_symbol = as.character(genes), stringsAsFactors = FALSE)
      m <- tryCatch(sdb$map(df, "gene_symbol", removeUnmappedRows = TRUE),
                    error = function(e) { msg("[warn] map() 失败：", conditionMessage(e)); NULL })
      if (is.null(m) || nrow(m) == 0) return(empty_map())
      m <- m[!is.na(m$STRING_id) & nzchar(as.character(m$STRING_id)), , drop = FALSE]
      if (!nrow(m)) return(empty_map())
      data.frame(gene_symbol = as.character(m$gene_symbol),
                 STRING_id = as.character(m$STRING_id), stringsAsFactors = FALSE)
    },
    edges = function(ids) {
      if (length(ids) < 2) return(NULL)
      tryCatch(sdb$get_interactions(ids),
               error = function(e) { msg("[warn] get_interactions() 失败：", conditionMessage(e)); NULL })
    }
  )
}

# ---- 引擎 B：STRING REST API（不依赖 Bioconductor）------------------------
api_base_of <- function(v) {
  if (nzchar(p$api_base)) return(sub("/+$", "", p$api_base))
  vv <- gsub("\\.", "_", trimws(v))
  if (identical(vv, "12") || identical(vv, "12_0")) "https://string-db.org/api" else
    sprintf("https://version%s.string-db.org/api", vv)
}

make_engine_api <- function() {
  base <- api_base_of(version)
  list(
    name = "api",
    map = function(genes) {
      if (!length(genes)) return(empty_map())
      url <- function(gv) sprintf(
        "%s/tsv/get_string_ids?identifiers=%s&species=%d&limit=1&echo_query=1&caller_identity=workbuddy_string_ppi_skill",
        base, paste(utils::URLencode(gv, reserved = TRUE), collapse = "%0d"), species)
      chunks <- split(genes, ceiling(seq_along(genes) / 200))
      out <- list()
      for (k in seq_along(chunks)) {
        gv <- unname(unlist(chunks[[k]]))
        dest <- file.path(cache_dir, sprintf("get_string_ids_%d.tsv", k))
        if (!fetch_url(url(gv), dest, timeout, retries)) next
        d <- tryCatch(read_tsv_strict(dest), error = function(e) NULL)
        if (!is.null(d) && nrow(d)) out[[length(out) + 1]] <- d
        Sys.sleep(1)   # STRING 要求批量调用间隔 >= 1 秒
      }
      if (!length(out)) return(empty_map())
      mid <- do.call(rbind, out)
      if (!all(c("queryItem", "stringId") %in% names(mid))) return(empty_map())
      mid$queryItem <- as.character(mid$queryItem); mid$stringId <- as.character(mid$stringId)
      mid$q_up <- toupper(trimws(mid$queryItem))
      m <- data.frame(gene_symbol = character(0), STRING_id = character(0), stringsAsFactors = FALSE)
      for (g in genes) {
        hit <- which(mid$q_up == toupper(g))
        if (length(hit)) m <- rbind(m, data.frame(gene_symbol = g,
                                                  STRING_id = mid$stringId[hit[1]],
                                                  stringsAsFactors = FALSE))
      }
      m[!duplicated(m$gene_symbol), , drop = FALSE]
    },
    edges = function(ids) {
      if (length(ids) < 2) return(NULL)
      chunks <- split(ids, ceiling(seq_along(ids) / 400))
      el <- list()
      for (k in seq_along(chunks)) {
        iv <- unname(unlist(chunks[[k]]))
        # required_score 在 0-1 口径下才生效；>1 时 STRING 会忽略，故统一传 0-1
        u <- sprintf("%s/tsv/network?identifiers=%s&species=%d&required_score=%.4f&network_type=%s&add_nodes=%d",
                     base, paste(iv, collapse = "%0d"), species,
                     max(0, min(1, score_threshold / 1000)), network_type, add_nodes)
        dest <- file.path(cache_dir, sprintf("network_%d.tsv", k))
        if (!fetch_url(u, dest, timeout, retries)) next
        d <- tryCatch(read_tsv_strict(dest), error = function(e) NULL)
        if (!is.null(d) && nrow(d)) el[[length(el) + 1]] <- d
        Sys.sleep(1)
      }
      if (!length(el)) return(NULL)
      do.call(rbind, el)
    }
  )
}

# ---- 引擎 C：本地 STRING 数据文件（完全离线）------------------------------
make_engine_local <- function() {
  info <- p$info_file; links <- p$links_file
  if (!nzchar(info) || !nzchar(links))
    stop("--engine=local 必须同时给 --info_file（protein.info.*.txt.gz）和 ",
         "--links_file（protein.links.*.txt.gz）")
  msg("[engine] 本地文件: ", info, " | ", links)
  pi <- utils::read.delim(gzfile(info), header = TRUE, sep = "\t", quote = "",
                          comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  nm <- tolower(names(pi))
  idc <- names(pi)[match("protein_external_id", nm)]
  pnc <- names(pi)[match("preferred_name", nm)]
  if (is.na(idc) || is.na(pnc))
    stop("info 文件需要 protein_external_id / preferred_name 两列，现有: ", paste(names(pi), collapse = ", "))
  pi$.up <- toupper(as.character(pi[[pnc]]))
  idx <- split(seq_len(nrow(pi)), pi$.up)
  lk <- NULL
  list(
    name = "local",
    map = function(genes) {
      if (!length(genes)) return(empty_map())
      m <- data.frame(gene_symbol = character(0), STRING_id = character(0), stringsAsFactors = FALSE)
      for (g in genes) {
        i <- idx[[toupper(g)]]
        if (length(i)) m <- rbind(m, data.frame(gene_symbol = g,
                                                STRING_id = as.character(pi[[idc]][i[1]]),
                                                stringsAsFactors = FALSE))
      }
      m[!duplicated(m$gene_symbol), , drop = FALSE]
    },
    edges = function(ids) {
      if (length(ids) < 2) return(NULL)
      if (is.null(lk)) {
        lk <<- utils::read.delim(gzfile(links), header = TRUE, sep = "\t", quote = "",
                                 comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
        lcn <- names(lk)
        cc1 <- intersect(c("protein1", "protein_external_id_a"), lcn)[1]
        cc2 <- intersect(c("protein2", "protein_external_id_b"), lcn)[1]
        ccs <- intersect(c("combined_score", "score"), lcn)[1]
        if (is.na(cc1) || is.na(cc2) || is.na(ccs))
          stop("links 文件列不匹配，现有: ", paste(lcn, collapse = ", "))
        names(lk)[names(lk) == cc1] <- "protein1"
        names(lk)[names(lk) == cc2] <- "protein2"
        names(lk)[names(lk) == ccs] <- "combined_score"
        lk$combined_score <- suppressWarnings(as.numeric(as.character(lk$combined_score)))
        lk$protein1 <- as.character(lk$protein1); lk$protein2 <- as.character(lk$protein2)
      }
      lk[lk$protein1 %in% ids & lk$protein2 %in% ids,
         c("protein1", "protein2", "combined_score"), drop = FALSE]
    }
  )
}

engine <- tolower(trimws(p$engine))
eng <- NULL
if (identical(engine, "stringdb")) {
  eng <- make_engine_stringdb()
  if (is.null(eng)) stop("--engine=stringdb 但 STRINGdb 不可用/初始化失败")
} else if (identical(engine, "api")) {
  eng <- make_engine_api()
} else if (identical(engine, "local")) {
  eng <- make_engine_local()
} else {
  # auto：只有 STRINGdb 的本地缓存已就位时才走 STRINGdb（否则它会现下载几十 MB，很慢且易断）
  have_cache <- any(grepl(sprintf("^%d\\.protein\\.aliases", species), list.files(cache_dir)))
  if (have_cache) {
    eng <- make_engine_stringdb()
  } else {
    msg("[engine] auto：cache_dir 里没有 STRINGdb 数据缓存，直接用 REST API",
        "（同一个 STRING ", version, " 数据库，结果同源）")
  }
  if (is.null(eng)) eng <- make_engine_api()
}
msg("实际使用的引擎    : ", eng$name)

# ===========================================================================
# 3. 映射 + 别名兜底（auto 模式只救「映射不上」的基因）
# ===========================================================================
map_with_alias <- function(genes, mode) {
  m1 <- eng$map(genes)
  got <- m1$gene_symbol
  miss <- setdiff(genes, got)
  subs <- character(0)
  if (length(miss) && !identical(mode, "none")) {
    cand <- unlist(ALIAS_TABLE[miss])
    cand <- cand[nzchar(cand)]
    if (length(cand)) {
      m2 <- eng$map(unname(cand))
      if (nrow(m2)) {
        rev <- setNames(names(cand), unname(cand))
        m2$gene_symbol <- unname(rev[m2$gene_symbol])
        m2$gene_symbol[is.na(m2$gene_symbol)] <- m2$STRING_id[is.na(m2$gene_symbol)]
        m1 <- rbind(m1, m2)
        subs <- paste0(names(cand), " -> ", unname(cand))
      }
    }
  }
  list(mapped = m1[!duplicated(m1$gene_symbol), , drop = FALSE], subs = subs)
}

mr <- map_with_alias(gt$gene, alias_mode)
mapped <- mr$mapped
if (length(mr$subs)) {
  msg("[alias] 以下基因原名映射失败，已用别名重试：")
  for (s in mr$subs) msg("        ", s)
}
unmapped <- setdiff(gt$gene, mapped$gene_symbol)
msg("映射成功          : ", nrow(mapped), " / ", nrow(gt))
# STRINGdb 引擎的数据文件下载不全时会「装作能映射」（返回 1/13 这种），必须报出来
if (identical(eng$name, "stringdb") && nrow(mapped) < 0.5 * nrow(gt)) {
  rule("!")
  msg("[warn] STRINGdb 引擎映射率异常偏低（", nrow(mapped), "/", nrow(gt),
      "），通常是 aliases/links 数据文件没下全。")
  msg("       对策：清掉 ", cache_dir, " 里残缺的 .txt.gz 重跑，或直接用 --engine=api。")
  rule("!")
}
if (length(unmapped)) {
  msg("")
  rule("!")
  msg("[warn] 以下 ", length(unmapped), " 个基因在 STRING ", version, " 中映射不上（不会出现在图里）：")
  for (g in unmapped) msg("        - ", g)
  rule("!")
  msg("")
}
if (nrow(mapped) == 0) stop("没有任何基因映射成功，无法建网络。检查 --species 与基因名。")

unmapped_df <- data.frame(gene = unmapped,
                          suggestion = ifelse(unmapped %in% names(ALIAS_TABLE),
                                              unlist(ALIAS_TABLE[unmapped]), ""),
                          stringsAsFactors = FALSE)
write.csv(unmapped_df, file.path(outdir, paste0(prefix, "_unmapped_genes.csv")),
          row.names = FALSE, fileEncoding = "UTF-8")

# ===========================================================================
# 4. 取边 -> 建图
# ===========================================================================
raw_edges <- eng$edges(mapped$STRING_id)
edges <- normalize_edges(raw_edges)

# STRING 官方 API 各版本 score 标度不一致（见过 0-1 与 0-1000 两种），这里统一到 0-1000
if (nrow(edges)) {
  mx <- max(edges$combined_score, na.rm = TRUE)
  if (is.finite(mx) && mx <= 1.0001) {
    msg("[score] 引擎返回 0-1 标度，已 x1000 归一到 STRING 0-1000 口径")
    edges$combined_score <- round(edges$combined_score * 1000, 0)
  }
}

# --add_nodes>0 时 STRING 会带进输入基因之外的蛋白，顺手记下它们的 symbol
extra_nm <- character(0)
if (!is.null(raw_edges) && nrow(raw_edges) &&
    all(c("stringId_A", "preferredName_A") %in% names(raw_edges))) {
  extra_nm <- c(setNames(as.character(raw_edges$preferredName_A), as.character(raw_edges$stringId_A)),
                setNames(as.character(raw_edges$preferredName_B), as.character(raw_edges$stringId_B)))
}
if (add_nodes == 0) {
  edges <- edges[edges$from %in% mapped$STRING_id & edges$to %in% mapped$STRING_id, , drop = FALSE]
}
msg("原始边数          : ", nrow(edges))
edges <- edges[edges$combined_score >= score_threshold, , drop = FALSE]
msg("过阈值边数        : ", nrow(edges), " (combined_score >= ", score_threshold, ")")

id2symbol <- c(setNames(mapped$gene_symbol, mapped$STRING_id), extra_nm)
sym_of <- function(ids) {
  s <- unname(id2symbol[ids])
  s[is.na(s)] <- ids[is.na(s)]
  s
}

if (nrow(edges) > 0) {
  g <- igraph::graph_from_data_frame(edges[, c("from", "to", "combined_score")], directed = FALSE)
} else {
  g <- igraph::make_empty_graph(n = 0, directed = FALSE)
}
if (as_bool(p$keep_isolated, TRUE)) {
  miss <- setdiff(mapped$STRING_id, igraph::V(g)$name)
  if (length(miss)) {
    g <- igraph::add_vertices(g, nv = length(miss), attr = list(name = miss))
    msg("补回孤立节点      : ", length(miss), " 个（这些基因在当前阈值下没有任何互作）")
  }
}
V(g)$symbol <- sym_of(V(g)$name)
if (igraph::ecount(g) > 0) {
  mx <- max(E(g)$combined_score, na.rm = TRUE)
  E(g)$score_norm <- if (is.finite(mx) && mx > 0) E(g)$combined_score / mx else rep(1, igraph::ecount(g))
} else {
  E(g)$score_norm <- numeric(0)
}

# --- logFC -----------------------------------------------------------------
logFC_named <- stats::setNames(gt$logFC, gt$gene)
V(g)$logFC <- unname(logFC_named[V(g)$symbol])
has_fc <- any(!is.na(V(g)$logFC))
msg("节点数 / 边数     : ", igraph::vcount(g), " / ", igraph::ecount(g))

# --- 剪枝 ------------------------------------------------------------------
prune_mode <- tolower(trimws(p$prune))
if (add_nodes > 0 && !given_prune) {
  prune_mode <- "none"
  msg("[note] --add_nodes>0 且未显式给 --prune，自动改为 prune=none（保留 STRING 补充的邻居）")
}
g2 <- g
if (has_fc && identical(prune_mode, "logfc")) {
  keep <- which(!is.na(V(g)$logFC))
  if (length(keep) < igraph::vcount(g)) {
    msg("剪枝              : 去掉 ", igraph::vcount(g) - length(keep),
        " 个无 logFC 的节点（--prune=logfc，等同参考代码）")
    g2 <- igraph::induced_subgraph(g, vids = keep)
  }
} else if (all(is.na(gt$logFC))) {
  msg("剪枝              : 跳过（输入没有 logFC，--prune 自动视为 none）")
}
msg("最终节点数 / 边数 : ", igraph::vcount(g2), " / ", igraph::ecount(g2))
if (igraph::vcount(g2) == 0) {
  stop("最终图一个节点都没有，无法出图。原因通常是：\n",
       "  - 阈值太高（当前 --score_threshold=", score_threshold, "），一条边都没有；且\n",
       "  - --keep_isolated=0，把没有边的节点也一起去掉了。\n",
       "  对策：下调 --score_threshold（试 200），或设 --keep_isolated=1（默认）。")
}

comp <- igraph::components(g2)
deg <- igraph::degree(g2)
msg("连通分量          : ", comp$no, "（最大分量 ", max(comp$csize), " 个节点）")
if (comp$no > 1) msg("[warn] 网络不连通，有 ", sum(deg == 0),
                     " 个孤立节点；可下调 --score_threshold 或增大 --add_nodes")
if (igraph::ecount(g2) > 0)
  msg("边权范围          : ", round(min(E(g2)$combined_score), 0), " ~ ",
      round(max(E(g2)$combined_score), 0))

# ===========================================================================
# 5. 布局（auto / --positions 手动）+ 出图
# ===========================================================================
set.seed(as_int(p$seed, 42))
layout <- tolower(trimws(p$layout))
if (identical(layout, "stress") && !requireNamespace("graphlayouts", quietly = TRUE)) {
  msg("[warn] 缺 graphlayouts，stress 布局不可用，改用 fr")
  layout <- "fr"
}
interactive_mode <- as_bool(p$interactive, FALSE)
pos_file <- trimws(p$positions)

# --- 自动布局（同时是 manual 坐标未覆盖节点的兜底）--------------------------
lyt <- ggraph::create_layout(g2, layout = layout)
x <- lyt$x; y <- lyt$y

# --- 手动布局：--positions（布局编辑器导出的坐标 csv）-----------------------
manual <- FALSE
if (nzchar(pos_file)) {
  if (!file.exists(pos_file)) stop("找不到 --positions 文件: ", pos_file)
  pf <- utils::read.csv(pos_file, stringsAsFactors = FALSE, check.names = FALSE)
  nm_pf <- tolower(names(pf))
  col_idx <- function(cands) {
    for (cc in cands) { i <- which(nm_pf == cc); if (length(i)) return(i[1]) }
    NA_integer_
  }
  ix <- col_idx(c("x")); iy <- col_idx(c("y"))
  isym <- col_idx(c("symbol", "gene", "name"))
  iid <- col_idx(c("string_id", "id"))
  if (is.na(ix) || is.na(iy))
    stop("--positions 文件缺 x / y 列。需要的列: symbol 或 STRING_id + x + y，现有: ",
         paste(names(pf), collapse = ", "))
  pf$x <- suppressWarnings(as.numeric(as.character(pf[[ix]])))
  pf$y <- suppressWarnings(as.numeric(as.character(pf[[iy]])))
  pf <- pf[is.finite(pf$x) & is.finite(pf$y), , drop = FALSE]
  sym_pf <- if (is.na(isym)) rep("", nrow(pf)) else as.character(pf[[isym]])
  id_pf  <- if (is.na(iid)) rep("", nrow(pf)) else as.character(pf[[iid]])
  hit <- rep(NA_integer_, igraph::vcount(g2))
  for (i in seq_len(igraph::vcount(g2))) {
    j <- which(id_pf == V(g2)$name[i])           # 先按 STRING_id 匹配
    if (!length(j)) j <- which(sym_pf == V(g2)$symbol[i])  # 再按 symbol
    if (length(j)) hit[i] <- j[1]
  }
  got <- sum(!is.na(hit))
  if (got == 0)
    stop("--positions 一个节点都没匹配上（检查 symbol / STRING_id 列的值是否与本网络一致）")
  if (got < igraph::vcount(g2))
    msg("[warn] --positions 只覆盖 ", got, "/", igraph::vcount(g2),
        " 个节点，其余保留自动布局坐标")
  x[!is.na(hit)] <- pf$x[hit[!is.na(hit)]]
  y[!is.na(hit)] <- pf$y[hit[!is.na(hit)]]
  manual <- TRUE
  msg("布局来源          : manual（--positions=", basename(pos_file),
      "，覆盖 ", got, "/", igraph::vcount(g2), " 个节点）")
} else {
  msg("布局来源          : auto（", layout, "，seed=", as_int(p$seed, 42), "）")
}

cols <- build_palette(p$palette, p$colors, n = 7)
lims <- if (!nzchar(p$color_limits) || tolower(trimws(p$color_limits)) %in% c("auto", "na")) {
  if (has_fc) {
    v <- V(g2)$logFC[is.finite(V(g2)$logFC)]
    if (!length(v)) c(0, 1) else { r <- range(v); if (r[1] == r[2]) r + c(-1, 1) else r }
  } else c(0, 1)
} else {
  v <- num_vec(p$color_limits, c(0, 1))
  if (length(v) == 1) c(min(0, v), max(0, v)) else v[1:2]
}
lims <- sort(lims)
msg("配色              : ", p$palette, " (", length(cols), " 段)")
msg("logFC 颜色区间    : ", round(lims[1], 3), " ~ ", round(lims[2], 3))

ewr <- num_vec(p$edge_width_range, c(0.1, 1.2)); if (length(ewr) < 2) ewr <- c(ewr[1], ewr[1])
lab <- if (nzchar(p$color_label)) p$color_label else expression(log[2] * " fold-change")

p0 <- if (manual) ggraph(g2, layout = "manual", x = x, y = y) else ggraph(g2, layout = layout)
if (igraph::ecount(g2) > 0) {
  p0 <- p0 + geom_edge_link(aes(edge_width = score_norm), colour = p$edge_colour,
                            alpha = as_num(p$edge_alpha, 0.9)) +
    scale_edge_width(range = c(ewr[1], ewr[2]), name = p$edge_label)
}
if (has_fc) {
  p0 <- p0 + geom_node_point(aes(colour = logFC), size = as_num(p$node_size, 6)) +
    scale_colour_gradientn(colours = cols, limits = lims, oob = squish_oob, name = lab)
} else {
  p0 <- p0 + geom_node_point(size = as_num(p$node_size, 6), colour = p$node_colour)
}
p0 <- p0 + geom_node_text(aes(label = symbol), repel = as_bool(p$repel, TRUE),
                          size = as_num(p$label_size, 3)) +
  theme_void() + theme(legend.position = tolower(trimws(p$legend)))

figs <- character(0)

# --- 交互模式：生成 HTML 布局编辑器（拖节点 -> 导出坐标 -> --positions 回读）---
if (interactive_mode) {
  node_hex <- rep("#999999", igraph::vcount(g2))
  if (has_fc) {
    tt <- (V(g2)$logFC - lims[1]) / (lims[2] - lims[1])
    tt <- pmin(1, pmax(0, tt))
    pal256 <- grDevices::colorRampPalette(cols)(256)
    fin <- is.finite(tt)
    node_hex[fin] <- pal256[pmin(256, round(tt[fin] * 255) + 1)]
  }
  ed_e <- if (igraph::ecount(g2) > 0) {
    ed0 <- igraph::as_data_frame(g2, what = "edges")
    data.frame(a = match(ed0$from, V(g2)$name), b = match(ed0$to, V(g2)$name),
               w = E(g2)$score_norm, stringsAsFactors = FALSE)
  } else {
    data.frame(a = integer(0), b = integer(0), w = numeric(0))
  }
  ed_e <- ed_e[is.finite(ed_e$a) & is.finite(ed_e$b), , drop = FALSE]
  nd_df <- data.frame(symbol = V(g2)$symbol, STRING_id = V(g2)$name,
                      logFC = V(g2)$logFC, color = node_hex, stringsAsFactors = FALSE)
  meta_txt <- paste0("engine=", eng$name,
                     " | layout=", if (manual) "manual" else layout,
                     " (seed=", as_int(p$seed, 42), ") | ",
                     igraph::vcount(g2), " nodes / ", igraph::ecount(g2),
                     " edges | color=logFC | width=STRING confidence")
  html_f <- file.path(outdir, paste0(prefix, "_layout_editor.html"))
  write_layout_editor_html(html_f, prefix = prefix, nodes = nd_df, edges = ed_e,
                           x = x, y = y, cols = cols, lims = lims,
                           meta = meta_txt, has_fc = has_fc)
  msg("布局编辑器        : ", normalizePath(html_f, mustWork = FALSE))
  if (!manual) {
    pos_csv <- file.path(outdir, paste0(prefix, "_positions_auto.csv"))
    write.csv(data.frame(symbol = V(g2)$symbol, STRING_id = V(g2)$name,
                         x = x, y = y, stringsAsFactors = FALSE),
              pos_csv, row.names = FALSE, fileEncoding = "UTF-8")
    msg("初始坐标快照      : ", normalizePath(pos_csv, mustWork = FALSE))
  }
}

# --- 出图 -------------------------------------------------------------------
# 最终图（pdf+png+svg）总是生成；编辑器是常备的微调入口，不卡出图流程。
figs <- save_plot_formats(p0, outdir, prefix, "ppi_network", fmts,
                          as_num(p$width, 7), as_num(p$height, 7), as_int(p$dpi, 300))
if (interactive_mode) {
  msg("")
  msg("[layout] 要微调布局：浏览器打开 ", prefix, "_layout_editor.html，",
      "按住节点拖动，点「导出坐标 CSV」，")
  msg("         然后原命令追加 --positions=<导出的 csv 路径> 重跑，即按你的摆位重新出图。")
}

# ===========================================================================
# 6. 落盘
# ===========================================================================
nodes_df <- data.frame(
  STRING_id = V(g2)$name,
  symbol = V(g2)$symbol,
  logFC = V(g2)$logFC,
  degree = as.integer(deg[V(g2)$name]),
  stringsAsFactors = FALSE
)
nodes_df <- nodes_df[order(-nodes_df$degree, nodes_df$symbol), , drop = FALSE]
write.csv(nodes_df, file.path(outdir, paste0(prefix, "_nodes.csv")),
          row.names = FALSE, fileEncoding = "UTF-8")

if (igraph::ecount(g2) > 0) {
  ed <- igraph::as_data_frame(g2, what = "edges")
  ed$symbol_from <- sym_of(ed$from); ed$symbol_to <- sym_of(ed$to)
  ed$score_norm <- round(ed$score_norm, 4)
  ed <- ed[, c("from", "to", "symbol_from", "symbol_to", "combined_score", "score_norm")]
  ed <- ed[order(-ed$combined_score), , drop = FALSE]
} else {
  ed <- data.frame(from = character(0), to = character(0), symbol_from = character(0),
                   symbol_to = character(0), combined_score = numeric(0), score_norm = numeric(0))
}
write.csv(ed, file.path(outdir, paste0(prefix, "_edges.csv")),
          row.names = FALSE, fileEncoding = "UTF-8")

if (as_bool(p$graphml, FALSE)) {
  gm <- file.path(outdir, paste0(prefix, "_network.graphml"))
  igraph::write_graph(g2, gm, format = "graphml")
  msg("已导出 graphml    : ", normalizePath(gm))
}
write_run_params(file.path(outdir, paste0(prefix, "_run_params.txt")), p)

hub <- head(nodes_df, 10)
rep_lines <- c(
  "STRING PPI 网络构建报告",
  strrep("=", 72), "",
  paste("生成时间        :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("脚本版本        :", SCRIPT_VERSION),
  paste("引擎            :", eng$name),
  paste("STRING 版本     :", version),
  paste("物种 taxid      :", species),
  paste("置信度阈值      :", score_threshold),
  paste("网络类型        :", network_type),
  paste("add_nodes       :", add_nodes),
  "", strrep("-", 72), "1. 输入与映射", strrep("-", 72),
  paste("输入基因数      :", nrow(gt)),
  paste("别名替换        :", if (length(mr$subs)) paste(mr$subs, collapse = "; ") else "无"),
  paste("映射成功        :", nrow(mapped)),
  paste("未映射          :", length(unmapped),
        if (length(unmapped)) paste0("(", paste(unmapped, collapse = ", "), ")") else ""),
  "", strrep("-", 72), "2. 网络统计", strrep("-", 72),
  paste("剪枝方式        :", if (has_fc && identical(prune_mode, "logfc"))
        "logfc（只保留有 logFC 的节点）" else "none"),
  paste("最终节点数      :", igraph::vcount(g2)),
  paste("最终边数        :", igraph::ecount(g2)),
  paste("连通分量        :", comp$no),
  paste("孤立节点        :", sum(deg == 0)),
  if (igraph::ecount(g2) > 0)
    paste("边权范围        :", round(min(E(g2)$combined_score), 0), "~",
          round(max(E(g2)$combined_score), 0)) else "边权范围        : 无",
  "", strrep("-", 72), "3. 绘图参数", strrep("-", 72),
  paste("布局来源        :", if (manual) paste0("manual (--positions=", basename(pos_file),
        ", seed 仍为 ", as_int(p$seed, 42), ")") else paste0(layout, " (auto, seed=", as_int(p$seed, 42), ")")),
  paste("配色            :", p$palette, "(", paste(cols, collapse = ","), ")"),
  paste("logFC 颜色区间  :", paste(round(lims, 3), collapse = " ~ ")),
  paste("线宽范围        :", paste(ewr, collapse = " ~ ")),
  paste("节点大小/标签   :", as_num(p$node_size, 6), "/", as_num(p$label_size, 3)),
  paste("输出图形        :", paste(basename(figs), collapse = " + "),
        " (", paste(fmts, collapse = "+"), " ", p$width, "x", p$height, "in)"),
  "", strrep("-", 72), "4. 度最高的 10 个节点（hub 参考）", strrep("-", 72)
)
for (i in seq_len(nrow(hub))) {
  rep_lines <- c(rep_lines, sprintf("  %-12s degree=%d  logFC=%s", hub$symbol[i], hub$degree[i],
    ifelse(is.na(hub$logFC[i]), "NA", sprintf("%.3f", hub$logFC[i]))))
}
rep_lines <- c(rep_lines, "", strrep("-", 72), "5. 交付前必须核对", strrep("-", 72),
  "  a. 未映射基因数是否为 0（未映射的基因不会出现在图里，也不会报错）",
  "  b. 别名替换是否符合预期（--alias_fix=none 可关闭）",
  "  c. 置信度阈值：默认 400（中高可信）；网络太稀疏时下调到 200 重跑",
  "  d. 若 --add_nodes>0：STRING 会引入输入基因之外的邻居，节点数会 > 输入基因数",
  "  e. logFC 颜色区间：auto 时取实际范围；要对齐参考原图请显式 --color_limits=-2,0",
  "  f. 布局种子固定 42，换种子图形会变但拓扑不变；--positions 手动布局不受种子影响",
  "     （只对未覆盖的节点用种子兜底）",
  "  g. 网络不连通是真实拓扑，不是脚本出错",
  "  h. 布局来源为 manual 时：核对报告里 --positions 的节点覆盖数是否 = 节点总数",
  "")
writeLines(rep_lines, file.path(outdir, paste0(prefix, "_string_report.txt")))

rule("=")
msg("完成。产物目录: ", normalizePath(outdir))
for (f in figs) msg("  ", basename(f))
if (interactive_mode) {
  msg("  ", prefix, "_layout_editor.html   （浏览器打开，拖节点摆位；--interactive=0 可关闭）")
  if (!manual) msg("  ", prefix, "_positions_auto.csv    （初始坐标快照）")
}
msg("  ", prefix, "_nodes.csv")
msg("  ", prefix, "_edges.csv")
msg("  ", prefix, "_unmapped_genes.csv")
msg("  ", prefix, "_string_report.txt")
msg("  ", prefix, "_run_params.txt")
rule("=")
