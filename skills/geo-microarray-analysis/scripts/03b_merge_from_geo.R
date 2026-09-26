#!/usr/bin/env Rscript
# ===========================================================================
# 03b_merge_from_geo.R  —— 【编排】从原始 GEO 文件一路做到合并矩阵
#
# 解决的问题: 03_merge_batches.R 要求输入**已标准化**的矩阵，但用户手上通常
# 是 N 组 (GSE*_series_matrix.txt.gz, GSE*_family.soft.gz)。手工做法要跑 N 遍
# 01 步再跑 03 步，且 N 份临床表要人工拼成一张分组表，极易出错。
#
# 本脚本把这条链自动化:
#   [1] 为每个 GSE 调用 01_geo_normalize.R  (逐数据集标准化 + 探针注释 + 去重 + 归一化)
#   [2] 收齐各 GSE 的产物，打印分组候选   <- 在这里 STOP (--stage=normalize)
#   [3] 由各 GSE 的临床表推导统一的样本->分组表
#   [4] 调用 03_merge_batches.R            (取共同基因 + cbind + 批次矫正 + QC)
#
# 关键设计: --stage 把"必须由用户拍板分组"这条铁律做进了命令行。
#   --stage=normalize  只跑到第 [2] 步就停，等用户看完分组候选再继续
#   --stage=merge      跳过标准化，直接用已有的 per_gse 产物做合并
#   --stage=all        一次跑完（适合已经有明确分组方案时）
#
# 用法:
#   # 第一步: 先只做标准化，同时看到每个 GSE 的分组候选
#   Rscript 03b_merge_from_geo.R --gse=GSE55584,GSE57218 --datadir=D:/data \
#           --outdir=./geo_merge --stage=normalize
#
#   # 第二步: 用户指定分组后，跑合并 + 批次矫正
#   Rscript 03b_merge_from_geo.R --gse=GSE55584,GSE57218 --datadir=D:/data \
#           --outdir=./geo_merge --stage=merge \
#           --group_cols='tissue:ch1;tissue:ch1' \
#           --group_maps='RA=case,OA=control;RA=case,Normal=control'
# ===========================================================================

options(stringsAsFactors = FALSE, warn = 1)
msg <- function(...) { base::message(...); flush.console() }

## --------------------------- locale 自愈 ----------------------------------
# 实测（2026-09-17）: 本机环境里带着 LC_ALL=C.UTF-8，而 Windows 版 R 应用不了
# "C.UTF-8"，只会退化到 "C" locale。在 C locale 下 R
#   **完全无法寻址含非 ASCII 字符的路径**：dir.exists/file.exists/dir.create/
#   readLines 一律失败，路径被悄悄改写成乱码
#   （含中文的路径会被探测改写成 8 位十六进制形式）。
# 改成 Windows 原生 UTF-8 locale 后，中文路径读写与中文命令行参数都恢复正常。
#
# 正常入口 run_geo.sh 已经 unset 了 LC_ALL，这里再兜一次：既修自己的文件 API，
# 也把变量从环境里清掉，好让派生出来的 01/03 子进程一启动就是正确的 locale。
if (!grepl("utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
  Sys.unsetenv(c("LC_ALL", "LANG", "LC_CTYPE", "LC_COLLATE", "LC_MONETARY", "LC_TIME"))
  suppressWarnings(tryCatch(Sys.setlocale("LC_CTYPE", ""), error = function(e) NULL))
  if (grepl("utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
    msg("  locale 已修正为 ", Sys.getlocale("LC_CTYPE"))
  } else {
    msg("  ! locale 仍是 ", Sys.getlocale("LC_CTYPE"),
        " —— R 在非 UTF-8 locale 下无法寻址含中文的路径。")
    msg("    请改用 run_geo.sh 启动，或把所有路径改成纯 ASCII。")
  }
}

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
  # --- 输入 ---
  gse = "", datadir = ".", matrices = "", softs = "", gpls = "",
  extra_inputs = "", extra_names = "",
  # --- 转发给 01_geo_normalize.R ---
  norm = "quantile", force_log = "auto", dup = "max", min_expr = "0",
  annot_col = "Gene Symbol", gene_sep = "///", boxplot = "TRUE",
  dl_timeout = "900", allow_download = "FALSE",
  # --- 转发给 03_merge_batches.R ---
  pre = "", file_sep = "auto", gene_case = "asis", toupper_scan = "TRUE",
  match = "intersect", na = "error",
  method = "combat", primary = "combat", combat_prior = "TRUE",
  combat_mean_only = "FALSE", combat_ref = "", combat_mod = "group",
  batch = "", batch_file = "", group_na = "error",
  pca_top = "0", pca_npc = "10", pca_color = "batch", label_samples = "FALSE",
  ellipse = "TRUE", density = "TRUE", density_max = "200",
  cor_heatmap = "FALSE", cor_max = "150",
  fig_style = "classic", fig_title = "",
  fig_width = "7.5", fig_height = "5.5", fig_width_mm = "", fig_height_mm = "",
  pca_width = "", pca_height = "", box_width = "", box_height = "",
  hm_width = "9", hm_height = "8", save_txt = "TRUE", save_raw = "TRUE",
  # --- 分组推导 ---
  group_cols = "", group_maps = "", group_regexes = "", group_file = "", group_col = "",
  clinical_files = "", group_out = "",
  # --- 编排 ---
  outdir = "./geo_merge", prefix = "merged",
  stage = "all", force = "FALSE", stop_on_fail = "TRUE", keep_logs = "TRUE",
  script_dir = "", help = "FALSE"
))

usage <- function() {
  cat("
03b_merge_from_geo.R —— 从原始 GEO 文件一路做到合并矩阵（01 x N + 03 编排）

用法: Rscript 03b_merge_from_geo.R --gse=<GSE1,GSE2,...> --datadir=<目录> [选项]

输入
  --gse=GSE55584,GSE57218,...   要合并的 GEO 号（用 --datadir 自动定位文件）
  --datadir=D:/data             在这些目录下**递归**查找
                                <GSE>_series_matrix.txt.gz 与 <GSE>_family.soft.gz
  --matrices=m1,m2,...          也可直接给出系列矩阵路径（与 --gse 一一对应）
  --softs=s1,s2,...             对应的 family.soft.gz 路径
  --gpls=GPL96,GPL6244,...      可选: 指定注释平台（GSE 含多平台时）
  --extra_inputs=p1,p2,...      可选: 额外并入**已标准化**的矩阵（如 TCGA 表达量）
  --extra_names=n1,n2,...       对应的显示名（默认取文件名）
  --allow_download=FALSE        本地找不到文件时是否允许回退到联网下载

01 步（逐数据集标准化）转发选项
  --norm=quantile  --force_log=auto  --dup=max  --min_expr=0
  --annot_col='Gene Symbol'  --gene_sep=///  --boxplot=TRUE  --dl_timeout=900
  （注意: 这里叫 --gene_sep，对应 01 的 --sep=探针多基因分隔符，与 03 的文件
    分隔符 --file_sep 不是一回事——两个脚本的 --sep 含义不同，故改名避开。）

03 步（合并 + 批次矫正）转发选项
  --pre=none                    03 侧预处理。默认 none: 各数据集已由 01 步标准化过，
                                这里只做取交集 + cbind + 批次矫正，不要再归一化一次。
                                extra_inputs 默认给 auto（未标准化的外部矩阵）。
  --file_sep=auto  --gene_case=asis  --toupper_scan=TRUE  --match=intersect  --na=error
  --method=combat|quantile|both|none   --primary=combat
  --combat_prior=TRUE  --combat_mean_only=FALSE  --combat_ref=  --combat_mod=group|none
  --batch=  --batch_file=  --group_na=error
  --pca_top=0  --pca_npc=10  --pca_color=batch|group  --label_samples=FALSE  --ellipse=TRUE
  --density=TRUE  --density_max=200  --cor_heatmap=FALSE  --cor_max=150
  --fig_style=classic|nature  --fig_title=  --fig_width=7.5  --fig_height=5.5
  --fig_width_mm=  --fig_height_mm=  --pca_width=  --pca_height=  --box_width=  --box_height=
  --hm_width=9  --hm_height=8  --save_txt=TRUE  --save_raw=TRUE

分组（把各 GSE 的临床表自动拼成一张统一的 样本->分组 表）
  --group_cols='tissue:ch1;title'      每个数据集用 clinical_<GSE>.csv 的哪一列
                                       （分号分隔，顺序对应 --gse；留空表示该数据集跳过）
  --group_regexes='re1;re2'            可选: 先从列值里用正则提取（必须一个捕获组）
  --group_maps='RA=case,OA=control;Normal=control,Tumor=case'
                                       列的取值 -> 统一分组名（逗号分隔）
  --group_file=file.csv                已有现成分组表时直接用它（优先于上面的推导）
  --group_col=                         配合 --group_file 指定列名
  --clinical_files=...                 可选: 覆盖自动定位的各 GSE 临床表
                                       (默认 per_gse/clinical_<GSE>.csv，逗号分隔、顺序对齐)
  说明: 先用 --stage=normalize 跑一遍，看每个 GSE 的 <GSE>_grouping_candidates.txt
        里有哪些候选列与取值，再据此写 --group_cols / --group_maps。
        合并完成后，merge/<prefix>_grouping_candidates.txt 会把 N 份临床信息汇总成
        一份（含各数据集候选列、取值计数、跨数据集同名列），**把它交给用户索要分组**。

编排
  --stage=all|normalize|merge|group
                                all       = 一条命令跑完（01 x N -> 03），不中途停
                                normalize = 只跑 01 x N 并打印分组候选后停止（推荐先跑这个）
                                merge     = 跳过标准化，直接合并（可反复重跑调参）
                                group     = **只重建分组表**：用 merge 阶段的合并临床表确定
                                            样本范围，按 --group_cols/--group_maps 重新推导
                                            分组，不重跑 03、不重算批次矫正。
                                            用户改了差异分析要用的分组时用它。
  --outdir=./geo_merge          总输出目录
  --prefix=merged               合并产物的文件名前缀
  --group_out=                  仅 --stage=group: 分组表输出路径
                                (默认 <outdir>/<prefix>_group_derived.csv)
  --force=FALSE                 为 TRUE 时即使某 GSE 已有产物也重跑 01
  --stop_on_fail=TRUE           某个 GSE 失败时是否立即中止
  --keep_logs=TRUE              保留每个 GSE 的 01 步运行日志
  --script_dir=                 01/03 脚本所在目录（一般不用给，会自动定位）

产物
  <outdir>/per_gse/<GSE>.csv  等        每个 GSE 的 01 步全部产物
  <outdir>/logs/01_<GSE>.log            每个 GSE 的标准化日志
  <outdir>/<prefix>_group_derived.csv   编排器推导出的统一分组表（sample,group）
  <outdir>/merge/<prefix>_merged.csv    最终合并矩阵（02 步的输入）
  <outdir>/merge/<prefix>_clinical.csv  合并后的临床信息表（样本名为行名）
  <outdir>/merge/<prefix>_grouping_candidates.txt
                                        各数据集的分组候选汇总 -> 交给用户索要分组
  <outdir>/merge/<prefix>_group_template.csv
                                        分组模板（全部样本已列好）-> 用户填好即可用，
                                        或让用户直接给一张他自己的分组表（--group_file=）
  <outdir>/merge/...                    03 步的全部 QC 产物
", sep = "")
  invisible(NULL)
}

if (isTRUE(OPT$help) || identical(tolower(as.character(OPT$help)), "true")) { usage(); quit(status = 0) }

as_bool <- function(x, default = FALSE) {
  if (length(x) == 0L || is.na(x)) return(default)
  if (is.logical(x)) return(x)
  toupper(as.character(x)) %in% c("TRUE", "T", "YES", "Y", "1")
}
split_list <- function(s, n = 0L) {
  v <- trimws(strsplit(as.character(s), ",", fixed = TRUE)[[1]])
  if (n > 0L && length(v) == 1L && n > 1L) v <- rep(v, n)
  v
}
# 分号分隔的逐数据集列表；长度为 1 时广播到所有数据集
split_per <- function(s, n) {
  v <- trimws(strsplit(as.character(s), ";", fixed = TRUE)[[1]])
  if (length(v) == 1L && n > 1L) v <- rep(v, n)
  if (length(v) != n) stop("这个逐数据集参数需要 ", n, " 段（用分号分隔），实际 ", length(v), " 段: ", s)
  v
}

## --------------------------- 运行环境准备 ---------------------------------
# 子进程（Rscript）的 tempdir() 在它自己启动时就固定了，所以必须在**派生之前**
# 把 TMP/TEMP/TMPDIR 指到纯 ASCII 且可写的目录，否则子进程会撞中文用户名那个坑。
is_ascii <- function(p) !grepl("[^ -~]", p, useBytes = TRUE)
pick_ascii_tmp <- function(extra = NULL) {
  cand <- c(extra, Sys.getenv("TMP"), Sys.getenv("TEMP"), Sys.getenv("TMPDIR"), file.path(tempdir(), "rtmp"), "C:/Rtmp", "D:/Rtmp")
  for (p in cand) {
    if (!nzchar(p) || !is_ascii(p)) next
    d <- gsub("\\\\", "/", p)
    if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(d) && file.access(d, 2L) == 0L) return(d)
  }
  NULL
}

find_rscript <- function() {
  cand <- c(file.path(R.home("bin"), "Rscript.exe"),
            file.path(R.home("bin"), "x64", "Rscript.exe"),
            file.path(R.home(), "bin", "x64", "Rscript.exe"),
            file.path(R.home(), "bin", "Rscript.exe"),
            file.path(R.home("bin"), "Rscript"),
            unname(Sys.which("Rscript")))
  cand <- cand[nzchar(cand)]
  hit <- cand[file.exists(cand)]
  if (!length(hit)) {
    stop("找不到 Rscript.exe。请在 PATH 里加入 R 的 bin 目录，或用 RSCRIPT 环境变量指定。")
  }
  normalizePath(hit[1], winslash = "/", mustWork = FALSE)
}

find_script_dir <- function() {
  d <- Sys.getenv("GEO_SKILL_DIR")
  if (nzchar(d) && dir.exists(d)) return(gsub("\\\\", "/", d))
  if (nzchar(OPT$script_dir) && dir.exists(OPT$script_dir)) return(gsub("\\\\", "/", OPT$script_dir))
  ca <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", ca, value = TRUE)
  if (length(m)) {
    p <- dirname(sub("^--file=", "", m[1]))
    if (dir.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
  }
  stop("无法定位 01/03 脚本所在目录。请用 --script_dir=<scripts 目录> 指定，\n",
       "  或在 run_geo.sh 包装器下运行（它会导出 GEO_SKILL_DIR）。")
}

RSCRIPT <- find_rscript()
SDIR    <- find_script_dir()
S01     <- file.path(SDIR, "01_geo_normalize.R")
S03     <- file.path(SDIR, "03_merge_batches.R")
for (f in c(S01, S03)) {
  if (!file.exists(f)) stop("找不到脚本: ", f, "\n  (--script_dir / GEO_SKILL_DIR 指向的目录不对?)")
}

tmp_ascii <- pick_ascii_tmp(extra = file.path(OPT$outdir, "_rtmp"))
if (is.null(tmp_ascii)) stop("找不到可写的纯 ASCII 临时目录，子进程的 tempdir() 会失效。")
Sys.setenv(TMP = tmp_ascii, TEMP = tmp_ascii, TMPDIR = tmp_ascii)
if (.Platform$OS.type == "windows" && !nzchar(Sys.getenv("COMSPEC"))) {
  Sys.setenv(COMSPEC = file.path(Sys.getenv("SystemRoot", "C:\\WINDOWS"), "system32", "cmd.exe"))
}

## ------------------------------ 工具函数 ----------------------------------

# 子脚本路径若含非 ASCII，R 打开脚本本身可能出问题（pitfalls 1.4.1），与
# run_geo.sh 同样处理: 复制到纯 ASCII 目录再跑。locale 正常后通常用不上，
# 但成本极低，作为兜底保留。
safe_script <- function(p, tag, tmpdir) {
  if (is_ascii(p)) return(gsub("\\\\", "/", p))
  dst <- file.path(tmpdir, paste0("_geo_", tag, ".R"))
  if (!file.copy(p, dst, overwrite = TRUE)) stop("复制脚本失败: ", p, " -> ", dst)
  msg("  ! 脚本路径含非 ASCII，已复制到 ", dst, " 运行 (pitfalls 1.4.1)")
  gsub("\\\\", "/", dst)
}

# 子进程里调 Rscript。
# ⚠ 退出码不可信: 本机 R 退出时固定报 Segmentation fault(139)，与工作是否完成
#   无关（见 pitfalls 1.4），调用方**一律按产物文件判成败**。
run_child <- function(script, args, logfile, tmpdir) {
  t0 <- Sys.time()
  st <- suppressWarnings(system2(RSCRIPT, c("--vanilla", script, args),
                                 stdout = logfile, stderr = logfile, wait = TRUE))
  list(status = st, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}
tail_log <- function(logfile, n = 25L) {
  if (!file.exists(logfile)) return("(无日志)")
  ln <- tryCatch(readLines(logfile, warn = FALSE), error = function(e) character(0))
  ln <- ln[nzchar(ln)]
  paste(utils::tail(ln, n), collapse = "\n")
}

S01 <- safe_script(S01, "01", tmp_ascii)
S03 <- safe_script(S03, "03", tmp_ascii)
msg("  临时目录: ", tmp_ascii, "   locale: ", Sys.getlocale("LC_CTYPE"))

# 只读矩阵 CSV 的表头拿样本名（避免为了取列名去读上百 MB 的文件）
read_samples_of <- function(csv) {
  con <- file(csv, "r"); on.exit(close(con), add = TRUE)
  h <- readLines(con, n = 1L, warn = FALSE)
  if (!length(h)) stop("读不到表头: ", csv)
  v <- strsplit(h, ",", fixed = TRUE)[[1]]
  v <- sub('^"(.*)"$', "\\1", v)
  v[-1]
}

read_clinical <- function(path) {
  d <- read.csv(path, header = TRUE, check.names = FALSE, na.strings = "")
  if (ncol(d) < 2L) stop("临床表列数不足: ", path)
  rn <- as.character(d[[1]])
  d <- d[, -1, drop = FALSE]
  rownames(d) <- rn
  d
}

# 逐数据集把某个临床列的值映射成**统一分组标签**。
#   universe     : 需要分组的样本名（按数据集顺序排列）
#   universe_ds  : 与 universe 等长，每个样本属于哪个数据集
#   cols/maps/regs: 分号分隔的逐数据集参数（顺序与 gses 对应）
# 每个数据集只给自己那些样本打标签，其余保持 NA（交给 03 的 --group_na 处理）。
# 这样 --stage=merge（universe=全部样本）与 --stage=group（universe=合并矩阵里的
# 样本，可能已被 --group_na=drop 裁过）能共用同一段逻辑。
derive_labels <- function(universe, universe_ds, per_dir, gses, cols, maps, regs) {
  lab <- rep(NA_character_, length(universe))
  for (i in seq_along(gses)) {
    g <- gses[i]; col <- cols[i]
    if (!nzchar(col)) { msg("    ", g, ": 未指定列，跳过（该数据集样本分组将为 NA）"); next }
    sel <- which(universe_ds == g)
    if (!length(sel)) { msg("    ", g, ": 没有它的样本，跳过"); next }
    cf <- file.path(per_dir, paste0("clinical_", g, ".csv"))
    if (!file.exists(cf)) {
      stop("找不到临床表: ", cf, "\n  --stage=merge / --stage=group 需要 per_gse 下的完整产物，",
           "请先跑 --stage=normalize 或 --stage=all。")
    }
    cl <- read_clinical(cf)
    if (!col %in% colnames(cl)) {
      stop("数据集 ", g, " 的临床表里没有列 '", col, "'。\n  可选列: ",
           paste(colnames(cl), collapse = ", "))
    }
    v <- trimws(as.character(cl[[col]]))
    if (nzchar(regs[i])) {
      v <- sub(regs[i], "\\1", v)
      msg("    ", g, ": 列 '", col, "' 经正则 '", regs[i], "' 提取")
    }
    if (nzchar(maps[i])) {
      # 映射表 取值=分组名,取值=分组名
      pairs <- trimws(strsplit(maps[i], ",", fixed = TRUE)[[1]])
      pairs <- pairs[nzchar(pairs)]
      keys <- character(0); vals <- character(0)
      for (p in pairs) {
        kv <- strsplit(p, "=", fixed = TRUE)[[1]]
        if (length(kv) < 2L) stop("--group_maps 每段必须是 '取值=分组名'，收到: ", p)
        keys <- c(keys, trimws(kv[1]))
        vals <- c(vals, trimws(paste(kv[-1L], collapse = "=")))
      }
      # 向量化映射。没被映射到的取值（例如要剔除的 rheumatoid arthritis）
      # 保持 NA，交给 03 步的 --group_na 处理。
      # ⚠ 不要写成 `if (!is.na(mp[[x]]))`——键不存在时 mp[[x]] 是 NULL，
      #   `is.na(NULL)` 得 logical(0)，`if` 会报 "argument is of length zero"。
      out <- rep(NA_character_, length(v))
      j <- match(v, keys); ok <- !is.na(j)
      out[ok] <- vals[j[ok]]
      v <- out
    } else {
      v[!nzchar(v) | v %in% c("NA", "null", "NULL")] <- NA_character_
    }
    idx <- match(universe[sel], rownames(cl))
    lab[sel] <- v[idx]
    got <- !is.na(idx)
    # 只统计这个数据集自己那些样本（表里可能出现别的数据集的样本名）
    tb <- table(lab[sel][got], useNA = "no")
    msg(sprintf("    %-12s 列 '%s'  匹配 %d/%d 个样本，分组: %s",
                g, col, sum(got), length(sel),
                paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", ")))
  }
  lab
}

# 在根目录下递归找 <gse>_series_matrix.txt.gz / <gse>_family.soft.gz
find_geo_files <- function(gse, roots) {
  pat_m <- paste0("^", gse, "_series_matrix\\.txt(\\.gz)?$")
  pat_s <- paste0("^", gse, "_family\\.soft(\\.gz)?$")
  r <- character(0)
  for (d in roots) if (dir.exists(d)) r <- c(r, d)
  if (!length(r)) return(list(matrix = "", soft = ""))
  fs <- unlist(lapply(r, function(d) list.files(d, recursive = TRUE, full.names = TRUE)),
               use.names = FALSE)
  if (!length(fs)) return(list(matrix = "", soft = ""))
  b <- basename(fs)
  list(matrix = if (any(grepl(pat_m, b, ignore.case = TRUE)))
         fs[which(grepl(pat_m, b, ignore.case = TRUE))[1]] else "",
       soft   = if (any(grepl(pat_s, b, ignore.case = TRUE)))
         fs[which(grepl(pat_s, b, ignore.case = TRUE))[1]] else "")
}

## -------------------- [1/4] 解析每个数据集的输入 ---------------------------
cat("\n=========== 从原始 GEO 文件到合并矩阵 (01 x N + 03) ===========\n")
stage <- tolower(trimws(as.character(OPT$stage)))
if (!stage %in% c("all", "normalize", "merge", "group")) {
  stop("--stage 只能是 all | normalize | merge | group，收到: ", OPT$stage)
}
do_normalize <- stage %in% c("all", "normalize")   # 是否要跑 01 x N
do_merge     <- stage %in% c("all", "merge")       # 是否要跑 03

outdir  <- gsub("\\\\", "/", OPT$outdir)
per_dir <- file.path(outdir, "per_gse")
mp_dir  <- file.path(outdir, "merge")
log_dir <- file.path(outdir, "logs")
for (d in c(outdir, per_dir, mp_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# 输出目录必须真的能建出来、能寻址。含中文的路径在 locale 不对时 dir.create 会
# 静默失败（见文件开头 locale 自愈那段），这里直接把问题暴露出来而不是让后面崩。
if (!dir.exists(outdir) || file.access(outdir, 2L) != 0L) {
  stop("输出目录不可用: ", outdir, "\n",
       if (!is_ascii(outdir))
         "  该路径含非 ASCII 字符，而当前 locale 下 R 无法寻址它。\n  把 --outdir 换成纯 ASCII 路径（如 D:/geo_merge）即可。\n"
       else "",
       "  当前 locale: ", Sys.getlocale("LC_CTYPE"))
}

gses <- toupper(trimws(split_list(OPT$gse)))
gses <- gses[nzchar(gses)]
mat_in <- trimws(split_list(OPT$matrices))
sof_in <- trimws(split_list(OPT$softs))
gpl_in <- trimws(split_list(OPT$gpls))

if (!length(gses) && !length(mat_in)) {
  stop("必须给 --gse=<GSE1,GSE2,...>（配合 --datadir），或用 --matrices/--softs 直接指定文件。")
}
if (length(mat_in) && length(gses) && length(mat_in) != length(gses)) {
  stop("--matrices 个数 (", length(mat_in), ") 与 --gse 个数 (", length(gses), ") 不一致")
}
n_geo <- if (length(gses)) length(gses) else length(mat_in)
if (!length(gses)) gses <- sub("_series_matrix.*$", "", basename(mat_in))
if (length(gpl_in) == 1L && n_geo > 1L) gpl_in <- rep("", n_geo)
if (!length(gpl_in)) gpl_in <- rep("", n_geo)

roots <- trimws(strsplit(as.character(OPT$datadir), ",", fixed = TRUE)[[1]])
roots <- roots[nzchar(roots)]

dsets <- list()
for (i in seq_len(n_geo)) {
  g <- gses[i]
  m <- if (length(mat_in)) mat_in[i] else ""
  s <- if (length(sof_in)) sof_in[i] else ""
  if (!nzchar(m) || !nzchar(s)) {
    hit <- find_geo_files(g, roots)
    if (!nzchar(m)) m <- hit$matrix
    if (!nzchar(s)) s <- hit$soft
  }
  # 相对路径按第一个 datadir 解析
  fix <- function(p) {
    if (!nzchar(p)) return("")
    if (file.exists(p)) return(normalizePath(p, winslash = "/", mustWork = FALSE))
    if (length(roots) && file.exists(file.path(roots[1], p))) {
      return(normalizePath(file.path(roots[1], p), winslash = "/", mustWork = FALSE))
    }
    p
  }
  dsets[[i]] <- list(gse = g, matrix = fix(m), soft = fix(s), gpl = gpl_in[i])
}

if (do_normalize) {
  msg("[1/4] 检查各数据集输入文件 ...")
  miss <- character(0)
  for (d in dsets) {
    ok_m <- nzchar(d$matrix) && file.exists(d$matrix)
    ok_s <- nzchar(d$soft)   && file.exists(d$soft)
    msg(sprintf("  %-12s matrix=%s  soft=%s", d$gse,
                if (ok_m) basename(d$matrix) else "**未找到**",
                if (ok_s) basename(d$soft)   else "**未找到**"))
    if (!ok_m || !ok_s) miss <- c(miss, d$gse)
  }
  if (length(miss)) {
    msg("")
    msg("以下数据集缺少本地文件: ", paste(miss, collapse = ", "))
    msg("请确认 --datadir 是否正确（会递归查找 <GSE>_series_matrix.txt.gz / <GSE>_family.soft.gz）。")
    if (any(!vapply(roots, is_ascii, logical(1)))) {
      msg("")
      msg("! --datadir 里有含非 ASCII 字符的路径，而当前 locale (", Sys.getlocale("LC_CTYPE"),
          ") 下 R 可能根本无法寻址它 ——")
      msg("  这时 list.files() 会静默返回空，看起来就像\"文件不存在\"。")
      msg("  请用 run_geo.sh 启动（它会清掉 LC_ALL），或把数据放到纯 ASCII 路径下。")
    }
    msg("GEO 下载位置: https://ftp.ncbi.nlm.nih.gov/geo/series/<GSEnnn>/<GSE>/matrix/ 与 .../soft/")
    if (!as_bool(OPT$allow_download)) {
      stop("默认不允许联网下载。确认要联网请加 --allow_download=TRUE。")
    }
    msg("--allow_download=TRUE: 交给 01 步回退到联网下载。")
  }
} else {
  msg("[1/4] --stage=", stage, ": 跳过输入文件检查与标准化，直接使用已有产物")
}

overview <- data.frame(dataset = character(0), n_gene = integer(0), n_sample = integer(0),
                       stringsAsFactors = FALSE)

## -------------------- [2/4] 逐个 GSE 跑 01 步 ------------------------------
if (do_normalize) {
  msg("")
  msg("[2/4] 逐个数据集标准化 (01_geo_normalize.R x ", n_geo, ") ...")
  for (i in seq_len(n_geo)) {
    d <- dsets[[i]]
    csv <- file.path(per_dir, paste0(d$gse, ".csv"))
    logf <- file.path(log_dir, paste0("01_", d$gse, ".log"))
    if (file.exists(csv) && !as_bool(OPT$force)) {
      msg("  [", i, "/", n_geo, "] ", d$gse, ": 已有产物，跳过 (--force=TRUE 可重跑)")
    } else {
      msg("  [", i, "/", n_geo, "] ", d$gse, ": 运行标准化 (日志 -> ", basename(logf), ") ...")
      args <- c(S01,
                paste0("--gse=", d$gse),
                paste0("--dir=", per_dir),
                paste0("--prefix=", d$gse),
                paste0("--norm=", OPT$norm),
                paste0("--force_log=", OPT$force_log),
                paste0("--dup=", OPT$dup),
                paste0("--min_expr=", OPT$min_expr),
                paste0("--annot_col=", OPT$annot_col),
                paste0("--sep=", OPT$gene_sep),
                paste0("--boxplot=", OPT$boxplot),
                paste0("--dl_timeout=", OPT$dl_timeout))
      if (nzchar(d$matrix) && file.exists(d$matrix)) args <- c(args, paste0("--matrix=", d$matrix))
      if (nzchar(d$soft)   && file.exists(d$soft))   args <- c(args, paste0("--soft=", d$soft))
      if (nzchar(d$gpl)) args <- c(args, paste0("--gpl=", d$gpl))
      r <- run_child(S01, args, logf, tmp_ascii)
      # 判成败看产物，不看退出码（本机 R 退出固定 139）
      if (!file.exists(csv)) {
        msg("")
        msg("  !! ", d$gse, " 标准化失败 (退出码 ", r$status, "，用时 ",
            round(r$secs, 1), "s，但未生成 ", basename(csv), ")")
        msg("  ---- 日志尾部 ----")
        msg(paste0("  ", strsplit(tail_log(logf), "\n", fixed = TRUE)[[1]]))
        msg("  ------------------")
        if (as_bool(OPT$stop_on_fail, TRUE)) stop("在 ", d$gse, " 处中止。修正后重跑，已完成的 GSE 会自动跳过。")
        next
      }
      msg("      完成，用时 ", round(r$secs, 1), "s")
    }
    smp <- read_samples_of(csv)
    overview <- rbind(overview, data.frame(dataset = d$gse,
                                           n_gene = NA_integer_, n_sample = length(smp),
                                           stringsAsFactors = FALSE))
  }
} else {
  # merge / group 模式: 只核对已有产物
  msg("[2/4] 核对各数据集已有的标准化产物 ...")
  for (i in seq_len(n_geo)) {
    g <- gses[i]
    csv <- file.path(per_dir, paste0(g, ".csv"))
    if (!file.exists(csv)) {
      stop("--stage=", stage, " 需要已存在的 per_gse/", g, ".csv，但没找到。\n",
           "  请先跑 --stage=normalize（或 --stage=all）完成标准化。")
    }
    if (stage == "group" && !file.exists(file.path(per_dir, paste0("clinical_", g, ".csv")))) {
      stop("--stage=", stage, " 需要 per_gse/clinical_", g, ".csv，但没找到。\n",
           "  请先跑 --stage=normalize（或 --stage=all）完成标准化。")
    }
    overview <- rbind(overview, data.frame(dataset = g, n_gene = NA_integer_,
                                           n_sample = length(read_samples_of(csv)),
                                           stringsAsFactors = FALSE))
  }
}

# 额外并入的已标准化矩阵
extra_in  <- trimws(split_list(OPT$extra_inputs))
extra_in  <- extra_in[nzchar(extra_in)]
extra_nm  <- trimws(split_list(OPT$extra_names))
if (length(extra_in)) {
  if (!length(extra_nm)) extra_nm <- sub("\\.[^.]*$", "", basename(extra_in))
  if (length(extra_nm) != length(extra_in)) stop("--extra_names 个数与 --extra_inputs 不一致")
  for (j in seq_along(extra_in)) {
    p <- extra_in[j]
    if (!file.exists(p)) stop("找不到 --extra_inputs 指定的矩阵: ", p)
    overview <- rbind(overview, data.frame(dataset = extra_nm[j], n_gene = NA_integer_,
                                           n_sample = length(read_samples_of(p)),
                                           stringsAsFactors = FALSE))
  }
  msg("  额外并入 ", length(extra_in), " 个已标准化矩阵: ", paste(extra_nm, collapse = ", "))
}

overview$n_gene <- NULL
msg("")
msg("  各数据集样本数:")
for (k in seq_len(nrow(overview))) {
  msg(sprintf("    %-16s %5d 样本", overview$dataset[k], overview$n_sample[k]))
}

## -------------------- [3/4] 汇总分组候选（STOP 点） ------------------------
all_csv   <- c(file.path(per_dir, paste0(gses, ".csv")), extra_in)
all_names <- c(gses, extra_nm)

msg("")
if (stage != "group") {
  msg("[3/4] 汇报各数据集的分组候选 ...")
  for (g in gses) {
    f <- file.path(per_dir, paste0(g, "_grouping_candidates.txt"))
    if (!file.exists(f)) { msg("  (", g, " 没有候选报告文件，跳过)"); next }
    msg("")
    msg("  ------------- ", g, " -------------")
    for (ln in readLines(f, warn = FALSE)) msg("  ", ln)
  }
}

if (stage == "normalize") {
  msg("")
  msg("===================================================================")
  msg("--stage=normalize 完成: 只做了标准化，尚未合并。")
  msg("")
  msg("!! 下一步必须向用户索要分组（分组是生物学判断，不要自己选）。")
  msg("   把上面每个 GSE 的候选列与取值整理给用户，让他确认:")
  msg("     ① 每个数据集用哪一列  ② 哪些取值算 case / control / 要剔除")
  msg("     ③ 对比方向 (contrast 的分子分母)")
  msg("   或者他直接给一张自己的分组表 (两列: 样本名,分组) -> 合并时 --group_file=")
  msg("")
  msg("拿到答复后写 --group_cols / --group_maps，再跑 --stage=merge:")
  msg("")
  msg("  Rscript 03b_merge_from_geo.R --gse=", paste(gses, collapse = ","),
      " --datadir=<同前>")
  msg("          --outdir=", outdir, " --stage=merge \\")
  msg("          --group_cols='", paste(rep("<列名>", length(gses)), collapse = ";"), "' \\")
  msg("          --group_maps='", paste(rep("<取值A>=case,<取值B>=control", length(gses)), collapse = ";"), "'")
  msg("")
  msg("若各数据集的取值本来就统一，可以省略 --group_maps。")
  msg("合并完成后 merge/<prefix>_grouping_candidates.txt 会把 N 份临床信息汇总成")
  msg("一份（含各数据集候选列、取值计数、跨数据集同名列），那份才是给用户看的主材料。")
  msg("要一次跑完不中途停: --stage=all")
  msg("===================================================================")
  quit(status = 0, save = "no")
}

## -------------------- [4/4] 推导分组表 + 跑 03 步 --------------------------
msg("")
if (stage == "group") {
  msg("[4/4] 重建分组表（不跑 03，不重算批次矫正）...")
} else {
  msg("[4/4] 合并 + 批次矫正 (03_merge_batches.R) ...")
}

# 4a) 分组表: --group_file 优先；否则由各 GSE 的临床表推导
group_file_use <- ""
if (stage == "group") {
  msg("  --stage=group: 跳过合并阶段的分组推导（本次只重建分组表）")
} else if (nzchar(OPT$group_file)) {
  if (!file.exists(OPT$group_file)) stop("找不到 --group_file: ", OPT$group_file)
  group_file_use <- normalizePath(OPT$group_file, winslash = "/", mustWork = FALSE)
  msg("  使用现成的 --group_file（跳过自动推导）")
} else if (nzchar(OPT$group_cols)) {
  cols <- split_per(OPT$group_cols, n_geo)
  maps <- if (nzchar(OPT$group_maps)) split_per(OPT$group_maps, n_geo) else rep("", n_geo)
  regs <- if (nzchar(OPT$group_regexes)) split_per(OPT$group_regexes, n_geo) else rep("", n_geo)

  # 合并后的样本顺序与 03 内部一致: 按数据集顺序 cbind
  all_samples <- unlist(lapply(all_csv, read_samples_of), use.names = FALSE)
  all_ds      <- rep(all_names, times = vapply(all_csv, function(f) length(read_samples_of(f)), 0L))

  msg("  由各数据集临床表推导分组:")
  lab <- derive_labels(all_samples, all_ds, per_dir, gses, cols, maps, regs)
  if (length(extra_in)) {
    msg("    ! 注意: --extra_inputs 的样本无法从 GEO 临床表推导分组，其分组为 NA。")
    msg("      若它们也需要分组，请改用 --group_file 直接给一张覆盖全部样本的分组表。")
  }
  n_na <- sum(is.na(lab))
  if (n_na) {
    msg("  ! 有 ", n_na, " 个样本没有分到组（占 ",
        round(100 * n_na / length(lab), 2), "%）")
    msg("    03 步会按 --group_na=", OPT$group_na, " 处理",
        if (identical(tolower(as.character(OPT$group_na)), "drop"))
          "（这些样本会从合并矩阵中剔除）" else "（默认 error 会直接报错，可改 --group_na=drop）")
    if (n_na == length(lab)) {
      stop("所有样本都没分到组，请检查 --group_cols / --group_maps 是否写对。")
    }
  }
  group_file_use <- file.path(outdir, paste0(OPT$prefix, "_group_derived.csv"))
  write.csv(data.frame(sample = all_samples, group = lab, check.names = FALSE),
            group_file_use, row.names = FALSE)
  msg("  统一分组表 -> ", group_file_use)
} else {
  msg("  未提供 --group_cols/--group_maps/--group_file: 不做分组推导。")
  msg("  ! 此时 ComBat 的 mod=NULL，可能把真实的生物学差异也当成批次效应扣掉。")
  msg("    建议回到 --stage=normalize 看分组候选后再跑。")
  if (identical(tolower(as.character(OPT$combat_mod)), "group")) {
    msg("    (--combat_mod 会自动降级为 none)")
  }
}

## -------------------- --stage=group: 只重建分组表 --------------------------
# 用户看过 merge/<prefix>_grouping_candidates.txt 之后，可能想换一个临床列做差异
# 分析的分组。这时没必要重跑 03（ComBat 很贵），只需按新分组重导分组表。
# 样本范围以**合并临床表**为准（它已对齐最终合并矩阵，若有样本被 --group_na=drop
# 剔除，这里也不会把它们算进来）。
if (stage == "group") {
  cl_merged <- file.path(mp_dir, paste0(OPT$prefix, "_clinical.csv"))
  if (!file.exists(cl_merged)) {
    stop("找不到合并临床表: ", cl_merged, "\n",
         "  --stage=group 需要先跑 --stage=merge（或 all）产出合并矩阵与合并临床表。\n",
         "  若上次合并时没有传临床表，请重跑 --stage=merge。")
  }
  if (!nzchar(OPT$group_cols)) {
    stop("--stage=group 需要 --group_cols（指定每个数据集用 clinical_<GSE>.csv 的哪一列）。\n",
         "  候选列见: ", file.path(mp_dir, paste0(OPT$prefix, "_grouping_candidates.txt")))
  }
  cdf <- read.csv(cl_merged, header = TRUE, check.names = FALSE, row.names = 1, na.strings = "")
  if (!"dataset" %in% colnames(cdf)) {
    stop("合并临床表里没有 dataset 列，无法判断样本属于哪个数据集: ", cl_merged)
  }
  universe    <- rownames(cdf)
  universe_ds <- as.character(cdf$dataset)
  msg("  合并矩阵内的样本: ", length(universe), " 个 (来自 ", length(unique(universe_ds)), " 个数据集)")
  for (g in gses) msg(sprintf("    %-14s %5d 个样本", g, sum(universe_ds == g)))

  cols <- split_per(OPT$group_cols, n_geo)
  maps <- if (nzchar(OPT$group_maps)) split_per(OPT$group_maps, n_geo) else rep("", n_geo)
  regs <- if (nzchar(OPT$group_regexes)) split_per(OPT$group_regexes, n_geo) else rep("", n_geo)
  lab <- derive_labels(universe, universe_ds, per_dir, gses, cols, maps, regs)

  n_na <- sum(is.na(lab))
  if (n_na) {
    msg("  ! ", n_na, " 个样本没有分到组（占 ", round(100 * n_na / length(lab), 2),
        "%），02 步会自动剔除它们")
  }
  if (n_na == length(lab)) stop("所有样本都没分到组，请检查 --group_cols / --group_maps 是否写对。")
  tb  <- table(lab, useNA = "no")
  msg("  分组计数: ", paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", "))

  out <- if (nzchar(OPT$group_out)) OPT$group_out
         else file.path(outdir, paste0(OPT$prefix, "_group_derived.csv"))
  write.csv(data.frame(sample = universe, group = lab, check.names = FALSE), out, row.names = FALSE)

  merged_csv_g <- file.path(mp_dir, paste0(OPT$prefix, "_merged.csv"))
  cat("\n=========== 分组表已重建（未重算批次矫正）===========\n")
  cat("合并矩阵  : ", merged_csv_g, "\n", sep = "")
  cat("分组表    : ", out, "  (", length(lab), " 个样本)\n", sep = "")
  cat("分组计数  : ", paste(sprintf("%s=%d", names(tb), as.integer(tb)), collapse = ", "), "\n", sep = "")
  cat("临床信息  : ", cl_merged, "\n", sep = "")
  cat("候选报告  : ", file.path(mp_dir, paste0(OPT$prefix, "_grouping_candidates.txt")), "\n", sep = "")
  cat("\n下一步 (差异分析):\n")
  cat("  bash run_geo.sh 02_deg_plots.R --expr=", merged_csv_g,
      " --group_file=", out, " --contrast=<组2>-<组1>\n", sep = "")
  cat("\n注 1: 本次只换了差异分析用的分组，合并矩阵的 ComBat 仍是按原分组做的 mod。\n")
  cat("      要让 ComBat 也保护新分组，请带 --group_cols/--group_maps 重跑 --stage=merge。\n")
  cat("注 2: --contrast 用 make.names 之后的分组名，含空格/中文的组名要先确认实际取值。\n")
  quit(status = 0, save = "no")
}

# 4b) 组装 03 的参数
pre_spec <- trimws(as.character(OPT$pre))
if (!nzchar(pre_spec)) {
  # 各 GEO 数据集已由 01 标准化 -> none；extra 是外部矩阵 -> auto
  pre_spec <- paste(c(rep("none", length(gses)), rep("auto", length(extra_in))), collapse = ";")
}

# 临床表: 默认取各 GSE 的 per_gse/clinical_<GSE>.csv（注意 01 是"前缀在前"的命名），
# 按 --inputs 的顺序拼成逗号分隔列表交给 03。03 会据此产出合并临床表 + 分组候选报告。
# 不存在的条目留空（03 会把该数据集的临床列全部记为 NA），extra 数据集没有临床表。
n_all <- length(all_csv)
clin_list <- if (nzchar(OPT$clinical_files)) {
  split_list(OPT$clinical_files)
} else {
  p <- file.path(per_dir, paste0("clinical_", gses, ".csv"))
  ifelse(file.exists(p), p, "")
}
if (length(clin_list) < n_all) clin_list <- c(clin_list, rep("", n_all - length(clin_list)))
if (length(clin_list) > n_all) clin_list <- clin_list[seq_len(n_all)]
n_clin <- sum(nzchar(clin_list))
if (n_clin) {
  msg("  临床表: ", n_clin, "/", n_all, " 个数据集有临床信息 -> 03 会产出合并临床表与分组候选报告")
} else {
  msg("  ! 没有任何数据集找到 clinical_<GSE>.csv，03 无法产出合并临床表。")
  msg("    下游要向用户索要分组时会缺少依据，建议先跑 --stage=normalize。")
}

args <- c(S03,
          paste0("--inputs=", paste(all_csv, collapse = ",")),
          paste0("--names=", paste(all_names, collapse = ",")),
          paste0("--clinical_files=", paste(clin_list, collapse = ",")),
          paste0("--pre=", pre_spec),
          paste0("--dir=", mp_dir),
          paste0("--prefix=", OPT$prefix),
          paste0("--sep=", OPT$file_sep),
          paste0("--gene_case=", OPT$gene_case),
          paste0("--toupper_scan=", OPT$toupper_scan),
          paste0("--match=", OPT$match),
          paste0("--na=", OPT$na),
          paste0("--method=", OPT$method),
          paste0("--primary=", OPT$primary),
          paste0("--combat_prior=", OPT$combat_prior),
          paste0("--combat_mean_only=", OPT$combat_mean_only),
          paste0("--combat_mod=", OPT$combat_mod),
          paste0("--group_na=", OPT$group_na),
          paste0("--pca_top=", OPT$pca_top),
          paste0("--pca_npc=", OPT$pca_npc),
          paste0("--pca_color=", OPT$pca_color),
          paste0("--label_samples=", OPT$label_samples),
          paste0("--ellipse=", OPT$ellipse),
          paste0("--density=", OPT$density),
          paste0("--density_max=", OPT$density_max),
          paste0("--cor_heatmap=", OPT$cor_heatmap),
          paste0("--cor_max=", OPT$cor_max),
          paste0("--fig_style=", OPT$fig_style),
          paste0("--fig_width=", OPT$fig_width),
          paste0("--fig_height=", OPT$fig_height),
          paste0("--hm_width=", OPT$hm_width),
          paste0("--hm_height=", OPT$hm_height),
          paste0("--save_txt=", OPT$save_txt),
          paste0("--save_raw=", OPT$save_raw))
for (kv in list(c("--combat_ref=", OPT$combat_ref), c("--batch=", OPT$batch),
                c("--batch_file=", OPT$batch_file), c("--group_col=", OPT$group_col),
                c("--fig_title=", OPT$fig_title), c("--fig_width_mm=", OPT$fig_width_mm),
                c("--fig_height_mm=", OPT$fig_height_mm), c("--pca_width=", OPT$pca_width),
                c("--pca_height=", OPT$pca_height), c("--box_width=", OPT$box_width),
                c("--box_height=", OPT$box_height))) {
  if (nzchar(kv[2])) args <- c(args, paste0(kv[1], kv[2]))
}
if (nzchar(group_file_use)) args <- c(args, paste0("--group_file=", group_file_use))

merged_csv <- file.path(mp_dir, paste0(OPT$prefix, "_merged.csv"))
logf3 <- file.path(log_dir, "03_merge.log")
r3 <- run_child(S03, args, logf3, tmp_ascii)

if (!file.exists(merged_csv)) {
  msg("")
  msg("!! 03 步失败（退出码 ", r3$status, "，用时 ", round(r3$secs, 1), "s，未生成 ",
      basename(merged_csv), "）")
  msg("---- 03 日志尾部 ----")
  msg(paste0("  ", strsplit(tail_log(logf3, 40L), "\n", fixed = TRUE)[[1]]))
  msg("--------------------")
  stop("03 步未产出合并矩阵。完整日志: ", logf3)
}

## ------------------------------ 汇总 --------------------------------------
tbl <- read.csv(file.path(mp_dir, paste0(OPT$prefix, "_batch_stats.csv")), check.names = FALSE)
n_smp <- length(read_samples_of(merged_csv))
n_gene <- tryCatch({
  con <- file(merged_csv, "r"); on.exit(close(con), add = TRUE)
  length(readLines(con, warn = FALSE)) - 1L
}, error = function(e) NA_integer_)

cat("\n=============== 全流程完成 ===============\n")
cat("预标准化    : ", length(gses), " 个 GEO 数据集", sep = "")
if (length(extra_in)) cat(" + ", length(extra_in), " 个外部矩阵", sep = "")
cat("\n", sep = "")
cat("合并矩阵    : ", n_gene, " 基因 x ", n_smp, " 样本\n", sep = "")
cat("矫正方法    : ", OPT$method, "\n", sep = "")
cat("批次数      : ", nrow(tbl), " -> ", paste(tbl$batch, collapse = ", "), "\n", sep = "")
cat("\n主要产物:\n")
cat("  合并矩阵   : ", merged_csv, "\n", sep = "")
cat("  样本批次表 : ", file.path(mp_dir, paste0(OPT$prefix, "_batch_map.csv")), "\n", sep = "")
cat("  匹配报告   : ", file.path(mp_dir, paste0(OPT$prefix, "_overlap_report.txt")), "\n", sep = "")
cat("  PCA 前后   : ", file.path(mp_dir, paste0(OPT$prefix, "_pca_before")), " / _pca_after\n", sep = "")
if (nzchar(group_file_use)) {
  cat("  分组表     : ", group_file_use, "\n", sep = "")
}
cat("  各 GSE 产物: ", per_dir, "/\n", sep = "")

# 临床信息是"向用户索要分组"的凭据: 有它就一定打出来，别让它埋在目录里
clin_merged <- file.path(mp_dir, paste0(OPT$prefix, "_clinical.csv"))
clin_report <- file.path(mp_dir, paste0(OPT$prefix, "_grouping_candidates.txt"))
has_clin <- file.exists(clin_merged)
if (has_clin) {
  cat("\n临床信息 (做差异分析前必须先给用户看): \n")
  cat("  合并临床表 : ", clin_merged, "\n", sep = "")
  cat("  分组候选报告: ", clin_report, "\n", sep = "")
  cat("  说明: 各数据集的临床列与取值词表都不一样，这两个文件把 N 份临床信息汇总到\n")
  cat("        一起，用户据此才能判断\"用哪一列、哪些取值算 case/control\"。\n")
}
cat("  分组模板   : ", file.path(mp_dir, paste0(OPT$prefix, "_group_template.csv")),
    " (全部样本已列好, 用户可直接填写)\n", sep = "")

cat("\n下一步 (第二步差异分析, 分组须由用户确认):\n")
cat("  bash run_geo.sh 02_deg_plots.R --expr=", merged_csv,
    " --group_file=", if (nzchar(group_file_use)) group_file_use else "<分组文件>",
    " --contrast=<组2>-<组1>\n", sep = "")
if (has_clin) {
  cat("\n!! 不要自行决定分组: 分组是生物学判断。\n")
  cat("   先把 ", basename(clin_report), " 与 ", basename(clin_merged), "\n", sep = "")
  cat("   交给用户，明确问他三件事:\n")
  cat("     ① 每个数据集用哪一列做分组\n")
  cat("     ② 该列的哪些取值算 case、哪些算 control、哪些样本要剔除\n")
  cat("     ③ 对比方向 (contrast 的分子分母)\n")
  cat("   或者更省事: 让他直接给一张自己的分组表 (两列 样本名,分组; 表头可有可无)\n")
  cat("     -> 02 步 --group_file=<他的表>；样本名与矩阵列名一致即可，\n")
  cat("        未覆盖的样本会被 02 步自动剔除。也可以让他填上面那张分组模板。\n")
  if (!nzchar(group_file_use)) {
    cat("   拿到答复后跑 (03_merge_batches.R 在技能脚本目录下):\n")
    cat("     bash run_geo.sh 03_merge_batches.R --inputs=", paste(all_csv, collapse = ","),
        " --names=", paste(all_names, collapse = ","), " --dir=", mp_dir,
        " --group_file=<用户确认后的分组表> [--method=...]\n", sep = "")
    cat("   或直接用合并临床表的某一列:\n")
    cat("     02_deg_plots.R --expr=", merged_csv, " --clinical=", clin_merged,
        " --group_from=<列名> [--group_keep=A,B]\n", sep = "")
  } else {
    cat("   若用户想换一个分组: --stage=group --group_cols=... --group_maps=...\n")
    cat("   (只重建分组表，不重算 ComBat)，或直接换 --group_file=<他的表>\n")
  }
} else if (!nzchar(group_file_use)) {
  cat("\n!! 尚未确定分组，且没有合并临床表可供用户参考。\n")
  cat("   请先跑 --stage=normalize 得到各 GSE 的 clinical_<GSE>.csv，再要分组；\n")
  cat("   也可以直接让用户给一张分组表 -> --group_file=<他的表>。\n")
}
