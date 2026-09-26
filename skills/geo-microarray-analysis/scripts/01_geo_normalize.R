#!/usr/bin/env Rscript
# ===========================================================================
# 01_geo_normalize.R  ——  【第一步】GEO 芯片标准化 + 分组候选探查
#
# 输入: 用户提供的 GEO 号 + 本地 family.soft.gz + series_matrix.txt.gz
# 输出: 标准化表达矩阵、临床信息、分组候选报告、分组模板
#
# 不联网也能跑。getGEO() 解析上百样本的系列矩阵要 10 分钟以上且内存暴涨，
# 因此本脚本自写流式解析器，只取需要的部分（实测 40MB soft + 14MB matrix
# 约 1 分钟完成，getGEO 同等数据 7 分钟仍未结束）。
#
# 用法:
#   Rscript 01_geo_normalize.R --gse=GSE62452 \
#       --matrix=D:/data/GSE62452_series_matrix.txt.gz \
#       --soft=D:/data/GSE62452_family.soft.gz --dir=./out
# ===========================================================================

options(stringsAsFactors = FALSE, warn = 1)

## ----------------------------- 参数解析 -----------------------------------
parse_args <- function(argv, defaults) {
  for (a in argv) {
    if (!grepl("^--", a)) next
    a <- substring(a, 3L)
    kv <- strsplit(a, "=", fixed = TRUE)[[1]]
    defaults[[kv[1]]] <- if (length(kv) > 1L) paste(kv[-1L], collapse = "=") else TRUE
  }
  defaults
}

OPT <- parse_args(commandArgs(trailingOnly = TRUE), list(
  gse = "", matrix = "", soft = "", dir = ".", prefix = "",
  gpl = "", force_log = "auto", norm = "quantile", dup = "max",
  min_expr = "0", annot_col = "Gene Symbol", sep = "///",
  boxplot = "TRUE", dl_timeout = "900", help = "FALSE"
))

usage <- function() {
  cat("
01_geo_normalize.R —— 【第一步】GEO 芯片标准化 + 分组候选探查

用法: Rscript 01_geo_normalize.R --gse=GSE62452 --matrix=<文件> --soft=<文件> [选项]

输入文件 (用户提供，脚本不联网)
  --matrix=<file>       系列矩阵 GSE*_series_matrix.txt.gz
  --soft=<file>         平台注释 GSE*_family.soft.gz
  --gse=GSEXXXXXX       GEO 号，用于命名与找不到本地文件时的兜底下载

常用
  --dir=.               输出目录，不存在自动创建
  --prefix=GSE62452     输出文件前缀，默认同 --gse
  --gpl=GPLxxxxx        指定注释平台 (GSE 含多平台时用)
  --force_log=auto      auto|TRUE|FALSE，是否强制 log2 转换
  --norm=quantile       归一化方法: quantile|median|cyclicloess|scale|Aquantile|none
  --dup=max             重复基因合并: max(行均值最大, 快) | elementwise(逐元素取最大)
  --min_expr=0          行均值 <= 该值的基因被剔除, 0 表示不过滤
  --annot_col=Gene Symbol  基因列名；找不到时自动探测，再找不到则用
                        gene_assignment 兜底 (Affy ST 类平台没有 Gene Symbol 列)
  --sep=///             一个探针对应多个基因时的分隔符
  --boxplot=TRUE        是否输出归一化前后箱线图

输出 (<dir>/<prefix>*)
  <prefix>.csv / .txt      标准化表达矩阵 (基因 x 样本)
  clinical_<prefix>.csv    样本临床/表型信息
  <prefix>_grouping_candidates.txt   候选分组列报告（候选列+每列取值样本数，供你确定分组）
  <prefix>_group_template.csv        分组模板 (sample,group)，填好后给第二步用
  <prefix>_boxplot.pdf     归一化前后箱线图
  <prefix>_normalized.RData
")
  invisible(NULL)
}

if (isTRUE(OPT$help) || (OPT$gse == "" && OPT$matrix == "")) {
  usage()
  quit(status = if (OPT$gse == "" && OPT$matrix == "") 1L else 0L, save = "no")
}

as_bool <- function(x, default = FALSE) {
  if (length(x) == 0L || is.na(x)) return(default)
  if (is.logical(x)) return(x)
  toupper(as.character(x)) %in% c("TRUE", "T", "YES", "Y", "1")
}

GSE      <- toupper(trimws(OPT$gse))
PREFIX   <- if (nzchar(OPT$prefix)) OPT$prefix else if (nzchar(GSE)) GSE else "GEO"
MIN_EXPR <- as.numeric(OPT$min_expr)
DUP      <- tolower(trimws(OPT$dup))
NORM     <- tolower(trimws(OPT$norm))

## ------------------------------ 环境 --------------------------------------
for (pkg in c("limma", "GEOquery")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("缺少 R 包: ", pkg, "\n请先安装: BiocManager::install(c('GEOquery','limma'))")
  }
}
suppressPackageStartupMessages(library(limma))

# 输入路径必须在 setwd 之前转成绝对路径，否则切目录后相对路径失效
abspath <- function(p) {
  if (!nzchar(p)) return("")
  if (!file.exists(p)) return(p)
  normalizePath(p, winslash = "/", mustWork = FALSE)
}
MATRIX <- abspath(OPT$matrix)
SOFT   <- abspath(OPT$soft)

dir.create(OPT$dir, showWarnings = FALSE, recursive = TRUE)
setwd(OPT$dir)
message("工作目录: ", normalizePath("."))
options(timeout = max(60, as.numeric(OPT$dl_timeout)))

if (!dir.exists(tempdir()) || file.access(tempdir(), 2L) != 0L) {
  stop("R 临时目录不可用: ", tempdir(),
       "\nWindows 中文用户名机器上的常见故障。启动 R 前先设置纯 ASCII 临时目录:",
       "\n  bash:        先把 TMP/TEMP/TMPDIR 指向纯 ASCII 目录（run_geo.sh 已自动处理）",
       "\n  PowerShell:  mkdir D:\\Rtmp; $env:TMP='D:\\Rtmp'; $env:TEMP='D:\\Rtmp'; $env:TMPDIR='D:\\Rtmp'")
}
if (.Platform$OS.type == "windows" && !nzchar(Sys.getenv("COMSPEC"))) {
  cmd <- file.path(Sys.getenv("SystemRoot", "C:\\WINDOWS"), "system32", "cmd.exe")
  if (file.exists(cmd)) {
    Sys.setenv(COMSPEC = cmd)
    message("环境自检: COMSPEC 未设置，已自动设为 ", cmd)
  } else {
    warning("环境自检: COMSPEC 未设置且找不到 cmd.exe，解压步骤可能失败")
  }
}

## ------------------------- 通用文件读取工具 --------------------------------
# 注意: 千万不要给 gzfile()/file() 传 encoding= 参数！
# 实测 gzfile(p, "rt", encoding="latin1") 会让 readLines 只读 51 行就认为 EOF，
# 导致 33364 行的文件只解析到头部。不指定 encoding 才能读全（非法字节仅 warning，已屏蔽）。
open_text <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
}
read_chunk <- function(con, n = 50000L) {
  suppressWarnings(readLines(con, n = n, warn = FALSE))
}

## -------------------- 1. 读取系列矩阵 (本地快速解析) ------------------------
read_series_matrix <- function(path) {
  con <- open_text(path)
  on.exit(close(con), add = TRUE)
  meta <- list(); state <- 0L; body <- character(0)

  repeat {
    ln <- read_chunk(con, 20000L)
    if (!length(ln)) break
    if (state == 0L) {
      i1 <- grep("^!series_matrix_table_begin", ln, useBytes = TRUE)
      if (length(i1)) {
        head_ln <- if (i1[1] > 1L) ln[seq_len(i1[1] - 1L)] else character(0)
        for (l in head_ln) {
          if (!grepl("^!Sample_", l)) next
          parts <- strsplit(l, "\t", fixed = TRUE)[[1]]
          tag <- parts[1]
          val <- if (length(parts) > 1L) parts[-1L] else ""
          meta[[tag]][[length(meta[[tag]]) + 1L]] <- val
        }
        state <- 1L
        ln <- ln[-seq_len(i1[1])]
      } else next
    }
    if (state == 1L) {
      i2 <- grep("^!series_matrix_table_end", ln, useBytes = TRUE)
      if (length(i2)) { body <- c(body, ln[seq_len(i2[1] - 1L)]); break }
      body <- c(body, ln)
    }
  }

  if (length(body) < 2L) stop("系列矩阵里没读到数据表，确认文件是否为 *_series_matrix.txt.gz")
  # 部分 GSE 的表头/首列带双引号 (如 "ID_REF" "GSM1527105" ...)。
  # 因为 read.delim 用了 quote=""，引号会被当字面量留下，必须手动去掉，
  # 否则列名变成 \"GSM1527105\"，第二步按样本名对分组时对不上。
  unq <- function(x) sub('^"(.*)"$', '\\1', x)
  cl <- unq(strsplit(body[1], "\t", fixed = TRUE)[[1]])
  ns <- length(cl) - 1L
  if (ns < 1L) stop("系列矩阵表头只有 1 列，格式异常")

  tf <- tempfile(fileext = ".tsv")
  writeLines(body[-1L], tf)
  d <- read.delim(tf, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE,
                  colClasses = c("character", rep("numeric", ns)))
  unlink(tf)

  m <- as.matrix(d[, -1L, drop = FALSE])
  colnames(m) <- cl[-1L]
  rownames(m) <- make.unique(unq(as.character(d[[1]])))

  gsms <- cl[-1L]
  pdl <- list()
  for (tag in names(meta)) {
    nm0 <- sub("^!Sample_", "", tag)
    for (k in seq_along(meta[[tag]])) {
      v <- as.character(meta[[tag]][[k]])
      # 系列矩阵里的取值一律带双引号: "GPL6244" / "tissue: Pancreatic tumor"
      # 不去掉会让平台号匹配失败、分组列名变成 "tissue:ch1 这种畸形
      v <- sub('^"(.*)"$', '\\1', v)
      if (length(v) == 1L) v <- rep(v, length(gsms))
      if (length(v) != length(gsms)) next
      nm <- if (k > 1L) paste0(nm0, ".", k - 1L) else nm0
      pdl[[nm]] <- v
    }
  }
  pd <- data.frame(pdl, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(pd) <- gsms

  # 把 "age: 52" 这类 characteristics 拆成独立列，GEOquery 也是这么做的
  for (cc in grep("^characteristics_ch", colnames(pd), value = TRUE)) {
    v <- as.character(pd[[cc]])
    if (!any(grepl(":", v, fixed = TRUE))) next
    ch  <- sub("\\.\\d+$", "", sub("^characteristics_", "", cc))
    key <- sub("^\\s*([^:]+)\\s*:.*$", "\\1", v)
    val <- sub("^\\s*[^:]+\\s*:\\s*", "", v)
    ok  <- grepl(":", v, fixed = TRUE)
    for (kk in unique(key[ok])) {
      nmk <- paste0(kk, ":", ch)
      if (nmk %in% colnames(pd)) next
      nv <- rep(NA_character_, nrow(pd))
      nv[ok & key == kk] <- val[ok & key == kk]
      pd[[nmk]] <- nv
    }
  }
  list(expr = m, pd = pd)
}

message("[1/6] 读取系列矩阵 ...")
if (nzchar(MATRIX) && file.exists(MATRIX)) {
  message("  本地文件: ", basename(MATRIX))
  sm <- read_series_matrix(MATRIX)
  dat <- sm$expr
  pd   <- sm$pd
  pd   <- pd[colnames(dat), , drop = FALSE]
} else {
  if (!nzchar(GSE)) stop("既没有 --matrix 也没有 --gse，无法获取表达矩阵")
  message("  本地无系列矩阵，尝试联网下载 ", GSE, " ...")
  suppressMessages(library(GEOquery))
  gset <- GEOquery::getGEO(GSE, destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)
  if (!length(gset)) stop("getGEO 未取到数据，检查 GSE 编号或网络。")
  eset <- gset[[1L]]
  dat  <- Biobase::exprs(eset); dat <- as.matrix(dat); storage.mode(dat) <- "numeric"
  pd   <- Biobase::pData(eset)
}
message("  原始矩阵: ", nrow(dat), " 探针 x ", ncol(dat), " 样本")

## 平台号: 优先 clinical 里的 platform_id，其次 --gpl
if (nzchar(OPT$gpl)) {
  PLAT <- toupper(trimws(OPT$gpl))
} else if ("platform_id" %in% colnames(pd)) {
  pu <- unique(as.character(pd$platform_id))
  pu <- pu[!is.na(pu) & nzchar(pu)]
  PLAT <- if (length(pu)) pu[1L] else ""
  if (length(pu) > 1L) message("  ! 样本涉及多个平台 (", paste(pu, collapse = ", "),
                               ")，默认用 ", PLAT, "；需要别的请用 --gpl 指定")
} else PLAT <- ""

## --------------------------- 2. log2 判断 ----------------------------------
message("[2/6] 判断是否需要 log2 转换 ...")
qx  <- as.numeric(quantile(dat, c(0, .25, .5, .75, .99, 1), na.rm = TRUE))
LogC <- (qx[5] > 100) ||
        (qx[6] - qx[1] > 50 && qx[2] > 0) ||
        (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (tolower(OPT$force_log) == "true")  { LogC <- TRUE;  message("  --force_log=TRUE，强制执行") }
if (tolower(OPT$force_log) == "false") { LogC <- FALSE; message("  --force_log=FALSE，跳过转换") }

logged <- FALSE
if (LogC) {
  n_na <- sum(dat <= 0, na.rm = TRUE)
  ex <- dat; ex[which(ex <= 0)] <- NA
  dat <- log2(ex); logged <- TRUE
  message("  已做 log2 转换", if (n_na > 0) paste0(" (", n_na, " 个 <=0 的值置为 NA)") else "")
} else {
  message("  数据已在 log 尺度，无需转换")
}
cat("  分位数: min=", round(qx[1], 2), " 25%=", round(qx[2], 2), " 中位数=",
    round(qx[3], 2), " 75%=", round(qx[4], 2), " max=", round(qx[6], 2), "\n", sep = "")

## ---------------------- 3. 导出临床信息 + 分组候选 --------------------------
message("[3/6] 导出临床信息并探查分组候选 ...")
write.csv(pd, paste0("clinical_", PREFIX, ".csv"), row.names = TRUE)
message("  表型列数: ", ncol(pd))

report_candidates <- function(pd, prefix) {
  n <- nrow(pd)
  out <- character(0)
  add <- function(...) out <<- c(out, paste0(...))
  add("================ 临床信息表（请据此自行确定分组）================")
  add("样本数: ", n, "  |  表型列数: ", ncol(pd))
  add("")
  add("【候选分组列：取值较少、可作为分组的列，请从中挑选或自行组合】")
  cnt <- 0L
  for (cc in colnames(pd)) {
    v <- trimws(as.character(pd[[cc]]))
    v[is.na(v)] <- ""
    v[v %in% c("NA", "null", "NULL", "---")] <- ""
    # "tissue: xxx" 这类原始列已有更干净的派生列 (tissue:ch1)，跳过避免重复
    if (grepl("^characteristics_ch", cc)) {
      vv <- v[nzchar(v)]
      if (length(vv) && all(grepl("^[^:]{1,40}:\\s*\\S", vv))) next
    }
    tb <- table(v[nzchar(v)])
    if (length(tb) < 2L) next
    if (length(tb) > min(10L, max(2L, floor(n / 2)))) next
    if (min(tb) < 2L || max(tb) > n - 2L) next
    cnt <- cnt + 1L
    add(sprintf("[%d] %s   (%d 个唯一值)", cnt, cc, length(tb)))
    for (nm in names(tb)) add(sprintf("      %-42s %d", nm, tb[[nm]]))
    add("")
  }
  if (cnt == 0L) {
    add("  未找到取值较少的候选列，请直接打开 clinical_", prefix, ".csv 人工判断。")
    add("")
  }
  # 其余高基数列只列名字，提示用户可正则提取
  hi <- character(0)
  for (cc in colnames(pd)) {
    v <- trimws(as.character(pd[[cc]])); v[is.na(v)] <- ""
    tb <- table(v[nzchar(v)])
    if (length(tb) >= 2L && length(tb) <= n &&
        min(tb) >= 2L && max(tb) <= n - 2L) next
    if (length(tb) == 0L) next
    hi <- c(hi, cc)
  }
  if (length(hi)) {
    add("【其他列（取值过多，一般不直接作分组，但可用 --group_regex 从中提取）】")
    add("      ", paste(hi, collapse = ", "))
    add("")
  }
  add("--------------------------------------------------------")
  add("分组是生物学判断，请自行决定，脚本不会替你选：")
  add("  ① 选某列直接用     -> 填好 <prefix>_group_template.csv (sample,group)，第二步 --group_file=...")
  add("  ② 用正则从某列提取 -> 第二步 --group_from=列名 --group_regex='(捕获组)'")
  add("  ③ 直接指名样本     -> 第二步 --group_spec='组A=GSM1,GSM2;组B=GSM3'")
  add("注意: 分组名里的空格/短横会被 make.names 洗成点号 (rheumatoid arthritis -> rheumatoid.arthritis)")
  txt <- paste(out, collapse = "\n")
  writeLines(txt, paste0(prefix, "_grouping_candidates.txt"))
  cat("\n", txt, "\n", sep = "")
}

report_candidates(pd, PREFIX)

tmpl <- data.frame(sample = colnames(dat), group = "", check.names = FALSE)
write.csv(tmpl, paste0(PREFIX, "_group_template.csv"), row.names = FALSE)
message("  分组模板 -> ", PREFIX, "_group_template.csv")

## --------------------------- 4. 平台注释 -----------------------------------
message("[4/6] 读取平台注释 ...")
soft_file <- if (nzchar(SOFT) && file.exists(SOFT)) SOFT else paste0(GSE, "_family.soft.gz")
if (!file.exists(soft_file)) {
  if (!nzchar(GSE)) stop("缺少平台注释文件，请用 --soft 指定 GSE*_family.soft.gz")
  n_num  <- suppressWarnings(as.numeric(sub("^GSE", "", GSE)))
  folder <- if (is.na(n_num)) GSE else paste0("GSE", floor(n_num / 1000), "nnn")
  url    <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/soft/%s_family.soft.gz",
                    folder, GSE, GSE)
  message("  本地无 ", soft_file, "，尝试从 NCBI 下载 (超时 ", OPT$dl_timeout, " 秒) ...")
  ok <- FALSE
  for (attempt in 1:3) {
    tr <- try(utils::download.file(url, soft_file, mode = "wb", quiet = TRUE), silent = TRUE)
    ok <- !inherits(tr, "try-error") && file.exists(soft_file) && file.size(soft_file) > 1000
    if (ok) break
    message("  第 ", attempt, " 次下载失败，重试 ...")
    unlink(soft_file)
  }
  if (!ok) stop("无法获取注释文件。请手动下载并放到工作目录: ", url)
  message("  下载完成: ", round(file.size(soft_file) / 1048576, 2), " MB")
}

# 流式扫描 soft，只取目标平台的 platform_table（40MB 文件约 5 秒）
read_platform_table <- function(path, want_gpl = "") {
  con <- open_text(path)
  on.exit(close(con), add = TRUE)
  state <- 0L; plat <- ""; body <- character(0)
  repeat {
    ln <- read_chunk(con, 50000L)
    if (!length(ln)) break
    if (state == 0L) {
      pos <- 0L
      for (ii in seq_along(ln)) {
        if (grepl("^\\^PLATFORM", ln[ii], useBytes = TRUE)) {
          acc <- trimws(sub("^\\^PLATFORM\\s*=\\s*", "", ln[ii]))
          if (!nzchar(want_gpl) || toupper(acc) == toupper(want_gpl)) {
            plat <- acc; pos <- ii; state <- 1L; break
          }
        }
      }
      if (state == 0L) next
      ln <- ln[-seq_len(pos)]
    }
    if (state == 1L) {
      j <- grep("^!platform_table_begin", ln, useBytes = TRUE)
      if (!length(j)) next
      state <- 2L
      ln <- ln[-seq_len(j[1])]
    }
    if (state == 2L) {
      k <- grep("^!platform_table_end", ln, useBytes = TRUE)
      if (length(k)) { body <- c(body, ln[seq_len(k[1] - 1L)]); break }
      body <- c(body, ln)
    }
  }
  if (length(body) < 2L) return(NULL)
  cl <- strsplit(body[1], "\t", fixed = TRUE)[[1]]
  tf <- tempfile(fileext = ".tsv")
  writeLines(body[-1L], tf)
  d <- read.delim(tf, header = FALSE, sep = "\t", quote = "", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE,
                  colClasses = "character")
  unlink(tf)
  colnames(d) <- cl[seq_len(ncol(d))]
  list(platform = plat, table = d)
}

t0 <- Sys.time()
gpl_res <- read_platform_table(soft_file, PLAT)
if (is.null(gpl_res)) {
  stop("soft 文件中没读到平台表。当前平台: ", if (nzchar(PLAT)) PLAT else "(未指定)",
       "\n试试用 --gpl=GPLxxxx 明确指定。")
}
gpl <- gpl_res$table
message("  平台 ", gpl_res$platform, " (", nrow(gpl), " 行 x ", ncol(gpl),
        " 列), 解析耗时 ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " 秒")

# 基因名提取: 优先 Gene Symbol 类列，否则用 gene_assignment 兜底
extract_symbol <- function(gpl_df, annot_col, sep) {
  cols <- colnames(gpl_df)
  col <- NA_character_; mode <- ""
  cands <- c(annot_col, "Gene Symbol", "Gene_Symbol", "GENE_SYMBOL",
             "Symbol", "SYMBOL", "gene_symbol", "Gene symbol", "GENE")
  for (c0 in cands) if (c0 %in% cols) { col <- c0; mode <- "column"; break }
  if (is.na(col)) {
    hit <- grep("symbol", cols, ignore.case = TRUE, value = TRUE)
    if (length(hit)) { col <- hit[1L]; mode <- "column" }
  }
  if (is.na(col)) {
    for (c0 in c("gene_assignment", "mrna_assignment", "Gene Assignment")) {
      if (c0 %in% cols) { col <- c0; mode <- "assignment"; break }
    }
  }
  if (is.na(col)) {
    stop("平台表里找不到基因名列，也没有 gene_assignment 可兜底。可用列: ",
         paste(cols, collapse = ", "), "\n请用 --annot_col= 指定。")
  }
  v <- as.character(gpl_df[[col]])
  v[is.na(v)] <- ""
  if (mode == "assignment") {
    # 形如 "ENST... // ARF5 // 描述 // ..." 多条用 " /// " 分隔
    v1 <- vapply(strsplit(v, "///", fixed = TRUE),
                 function(x) if (length(x)) x[1] else "", character(1))
    v  <- vapply(strsplit(v1, "//", fixed = TRUE),
                 function(x) if (length(x) > 1L) trimws(x[2]) else "", character(1))
  } else {
    # strsplit("", "///", fixed=TRUE) 返回 character(0)，取 [1] 得 NA —— 先补分隔符
    v <- vapply(strsplit(paste0(v, sep), sep, fixed = TRUE), `[`, character(1), 1L)
    v <- trimws(v)
  }
  v[is.na(v)] <- ""
  list(symbol = v, col = col, mode = mode)
}

SY <- extract_symbol(gpl, OPT$annot_col, OPT$sep)
if (SY$mode == "assignment") {
  message("  平台无 Gene Symbol 列，改用 '", SY$col, "' 字段提取基因名")
} else if (SY$col != OPT$annot_col) {
  message("  基因列 '", OPT$annot_col, "' 不存在，自动改用 '", SY$col, "'")
} else {
  message("  基因列: ", SY$col)
}
if (!"ID" %in% colnames(gpl)) stop("平台表缺少 ID 列。")

ids <- data.frame(probe_id = as.character(gpl[["ID"]]),
                  symbol   = SY$symbol)

## -------------------- 5. 探针匹配 + 重复基因合并 ----------------------------
message("[5/6] 探针匹配与基因合并 ...")
rownames(dat) <- as.character(rownames(dat))
ids <- ids[ids$probe_id %in% rownames(dat), , drop = FALSE]
ids <- ids[!duplicated(ids$probe_id), , drop = FALSE]
if (!nrow(ids)) stop("注释探针与表达矩阵行名无一匹配，检查 --gpl 是否选错。")

dat <- dat[ids$probe_id, , drop = FALSE]
stopifnot(identical(rownames(dat), ids$probe_id))

# 过滤必须对 NA 安全: NA %in% c("","NA") 返回 FALSE，只写 %in% 会把 NA 全留下
bad <- is.na(ids$symbol) | !nzchar(ids$symbol) |
       ids$symbol %in% c("---", "NA", "null", "NULL")
if (any(bad)) message("  剔除无有效基因名的探针 ", sum(bad), " 个")
ids <- ids[!bad, , drop = FALSE]
dat <- dat[ids$probe_id, , drop = FALSE]
message("  注释到 ", nrow(dat), " 个探针 / ", length(unique(ids$symbol)), " 个基因")

if (DUP == "max") {
  rm_vals <- rowMeans(dat, na.rm = TRUE)
  rm_vals[!is.finite(rm_vals)] <- -Inf
  o <- order(ids$symbol, -rm_vals)
  dat <- dat[o, , drop = FALSE]; ids <- ids[o, , drop = FALSE]
  keep <- !duplicated(ids$symbol)
  dat <- dat[keep, , drop = FALSE]; ids <- ids[keep, , drop = FALSE]
} else if (DUP == "elementwise") {
  num <- as.data.frame(lapply(as.data.frame(dat), as.numeric), check.names = FALSE)
  agg <- aggregate(num, by = list(symbol = ids$symbol),
                   FUN = function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
  rownames(agg) <- agg$symbol; agg$symbol <- NULL
  dat <- as.matrix(agg)
  # ⚠ 这里 dat 已经按基因去重完毕（nrow = 唯一基因数），而 ids 仍是**探针级别**
  #   （nrow = 探针数）。必须把 ids 同步到去重后的顺序，否则下面第 508 行的
  #   `rownames(dat) <- ids$symbol` 会因长度不一致报
  #     "Error in dimnames(x) <- dn : length of 'dimnames' [1] not equal to array extent"
  #   （--dup=max 分支里 dat 与 ids 是并行子集的，所以没有这个问题。）
  ids <- ids[match(rownames(dat), ids$symbol), , drop = FALSE]
} else {
  stop("--dup 只支持 max 或 elementwise")
}
rownames(dat) <- ids$symbol
if (anyNA(rownames(dat))) stop("基因名中出现 NA，请检查 --annot_col")
message("  合并后: ", nrow(dat), " 个基因")

if (MIN_EXPR > 0) {
  n0 <- nrow(dat)
  dat <- dat[rowMeans(dat, na.rm = TRUE) > MIN_EXPR, , drop = FALSE]
  message("  低表达过滤 (行均值 > ", MIN_EXPR, "): ", n0, " -> ", nrow(dat))
}

## ----------------------------- 6. 归一化 ------------------------------------
message("[6/6] 归一化 (", NORM, ") ...")
if (NORM != "none") {
  dat_norm <- limma::normalizeBetweenArrays(dat, method = NORM)
  if (as_bool(OPT$boxplot, TRUE)) {
    pdf(paste0(PREFIX, "_boxplot.pdf"), width = max(7, ncol(dat) * 0.35), height = 6)
    boxplot(data.frame(dat), col = "#4DBBD5", main = "Before normalization",
            las = 2, cex.axis = 0.5, outline = FALSE)
    boxplot(data.frame(dat_norm), col = "#4DBBD5",
            main = paste0("After normalization (", NORM, ")"),
            las = 2, cex.axis = 0.5, outline = FALSE)
    invisible(dev.off())
  }
  dat <- dat_norm
} else {
  message("  --norm=none，跳过归一化")
}

## ------------------------------ 导出 ----------------------------------------
write.table(data.frame(ID = rownames(dat), dat, check.names = FALSE),
            file = paste0(PREFIX, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(dat, paste0(PREFIX, ".csv"))
save(dat, pd, file = paste0(PREFIX, "_normalized.RData"))

cat("\n================ 第一步完成 ================\n")
cat("GSE        : ", GSE, "\n")
cat("平台       : ", gpl_res$platform, "\n")
cat("log2       : ", if (logged) "已转换" else "未转换(数据已在 log 尺度)", "\n")
cat("归一化     : ", NORM, "\n")
cat("最终矩阵   : ", nrow(dat), " 基因 x ", ncol(dat), " 样本\n")
cat("输出文件   : ", PREFIX, ".csv, ", PREFIX, ".txt, clinical_", PREFIX, ".csv, ",
    PREFIX, "_grouping_candidates.txt, ", PREFIX, "_group_template.csv\n", sep = "")
cat("\n临床信息表已写入: clinical_", PREFIX, ".csv (", nrow(pd), " 个样本 x ", ncol(pd),
    " 列)\n", sep = "")
cat("请据此（及 _grouping_candidates.txt 候选列）自行确定差异分析的分组。\n")
cat("  确认分组后再跑第二步，脚本不会替你自动决定分组：\n", sep = "")
cat("  Rscript 02_deg_plots.R --expr=", PREFIX, ".csv --group_file=<你填好的分组文件> ",
    "--contrast=<组2>-<组1>\n", sep = "")
