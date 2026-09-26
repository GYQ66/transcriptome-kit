#!/usr/bin/env Rscript
# ===========================================================================
# 03_merge_batches.R  ——  【第三步】多数据集合并 + 批次矫正
#
# 输入: 多个已完成标准化的表达矩阵 (基因 x 样本)。典型来源:
#         - 01_geo_normalize.R 的产物 <prefix>.csv / <prefix>.txt
#         - 用户已有的 bulk 矩阵 (含 TCGA 表达量，如 TPM/FPKM 文本)
# 输出: 合并矩阵 (可选 ComBat / quantile 批次矫正) + 批次 QC 图表
#
# 本步骤对应生信常规流程:
#   逐数据集预处理 -> 取共同基因 -> cbind 合并 -> ComBat 批次矫正
#   参考实现: <本地参考脚本>/7.R
#   --pre=auto 的判定规则刻意复刻 7.R，本脚本在其上与该实现数值一致。
#
# 用法:
#   # A) 三个数据集，默认 auto 预处理 + ComBat 矫正
#   Rscript 03_merge_batches.R --inputs=A.txt,B.txt,C.txt --dir=./out --prefix=merged
#
#   # B) 逐数据集指定预处理 (分号分隔，顺序对应 --inputs)
#   Rscript 03_merge_batches.R --inputs=A.txt,B.txt,C.txt \
#           --pre='log2p1;quantile;quantile' --names=TCGA,GSE30219,GSE74777
#
#   # C) 提供分组信息，ComBat 保护生物学变异 (mod = ~group)
#   Rscript 03_merge_batches.R --inputs=A.txt,B.txt,C.txt \
#           --group_file=group.csv --combat_mod=group --method=both
# ===========================================================================

options(stringsAsFactors = FALSE, warn = 1)

# 兼容性: 中文用户名 / 部分 R 构建下，输出缓冲会在(偶发的)崩溃时整段丢失，
# 导致脚本“无任何输出、且无产物”。此处把 message 包一层强制 flush.console()，
# 既不影响结果，又能让每步进度即时落盘、崩溃点可被定位。
msg <- function(...) { base::message(...); flush.console() }

## --------------------- Nature 期刊风格绘图组件 (可选) -----------------------
# 与 02_deg_plots.R 保持同一套主题/配色/导出约定，改编自开源项目 nature-skills
# (github.com/Yuan1z0825/nature-skills, Apache-2.0) 的 nature-figure 技能。
NATURE_PALETTE <- c(
  blue_main      = "#0F4D92",   # 深蓝 —— 主方法/关键对象
  blue_secondary = "#3775BA",   # 中蓝
  red_strong     = "#B64342",   # 强调红
  neutral_light  = "#CFCECE",   # 中性浅灰
  neutral_mid    = "#767676",   # 中性中灰
  neutral_dark   = "#4D4D4D"    # 中性深灰
)
# 批次/数据集配色: 期刊常用定性色板
BATCH_COLS <- c("#0F4D92", "#B64342", "#2E8B57", "#B8860B", "#6A3D9A",
                "#3775BA", "#8C564B", "#E377C2", "#17BECF", "#CFCECE")

theme_nature <- function(base_size = 7, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line    = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks   = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title   = ggplot2::element_text(size = base_size),
      axis.text    = ggplot2::element_text(size = base_size - 0.5, colour = "black"),
      legend.title = ggplot2::element_text(size = base_size - 0.5),
      legend.text  = ggplot2::element_text(size = base_size - 0.7),
      legend.key   = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(size = base_size + 0.5, face = "bold", hjust = 0.5),
      panel.grid   = ggplot2::element_blank()
    )
}

theme_classic_std <- function(base_size = 14) ggplot2::theme_bw(base_size = base_size)

# 多格式导出: PDF/PNG 始终输出; nature 风格下额外输出 SVG + 600dpi TIFF
# 注意: 默认 pdf() 设备在未注册 Arial 时会报 "invalid font type"，必须用 cairo_pdf;
#       PNG/TIFF 优先用 ragg (能解析系统字体)。
save_figure <- function(plot, file_base, width, height, style = "classic") {
  pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height, device = pdf_dev)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
                    dpi = 300, device = ragg::agg_png)
  } else {
    ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width, height = height,
                    dpi = 300, type = "cairo")
  }
  if (identical(style, "nature")) {
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggplot2::ggsave(paste0(file_base, ".svg"), plot, width = width, height = height,
                      device = svglite::svglite)
    } else {
      msg("  ! 未安装 svglite，跳过 SVG 输出 (可选: install.packages('svglite'))；",
          "矢量图可由 .pdf 替代")
    }
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(paste0(file_base, ".tiff"), plot, width = width, height = height,
                      dpi = 600, device = ragg::agg_tiff)
    } else {
      msg("  ! 未安装 ragg，跳过 TIFF 输出 (可选: install.packages('ragg'))")
    }
  }
  invisible(NULL)
}

# 把英寸参数统一成 c(宽, 高)，mm 优先 (期刊单栏约 89 mm)
fig_size <- function(w_in, h_in, w_mm, h_mm) {
  c(if (nzchar(w_mm)) as.numeric(w_mm) / 25.4 else as.numeric(w_in),
    if (nzchar(h_mm)) as.numeric(h_mm) / 25.4 else as.numeric(h_in))
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
  inputs = "", names = "", sep = "auto", pre = "auto", gene_case = "asis",
  toupper_scan = "TRUE", match = "intersect", min_expr = "0", na = "error",
  dup = "max", method = "combat", primary = "combat",
  combat_prior = "TRUE", combat_mean_only = "FALSE", combat_ref = "",
  combat_mod = "group",
  batch = "", batch_file = "", group_file = "", group_col = "", group_na = "error",
  clinical_files = "",
  dir = ".", prefix = "merged", save_txt = "TRUE", save_raw = "TRUE",
  pca_top = "0", pca_npc = "10", pca_color = "batch",
  label_samples = "FALSE", ellipse = "TRUE",
  density = "TRUE", density_max = "200", cor_heatmap = "FALSE", cor_max = "150",
  fig_style = "classic", fig_title = "",
  fig_width = "7.5", fig_height = "5.5",
  fig_width_mm = "", fig_height_mm = "",
  pca_width = "", pca_height = "", box_width = "", box_height = "",
  hm_width = "9", hm_height = "8",
  help = "FALSE"
))

usage <- function() {
  cat("
03_merge_batches.R —— 多数据集合并 + 批次矫正

用法: Rscript 03_merge_batches.R --inputs=<矩阵1,矩阵2,...> [选项]

输入:
  --inputs=A.txt,B.txt,C.txt   要合并的矩阵 (基因 x 样本)。支持 .txt(制表符)/.csv，第一列为基因名。
  --names=TCGA,GSE30219,...    数据集显示名 (默认取文件名去扩展名)，用于批次标签与图表。
  --sep=auto                   分隔符: auto|csv|tab (默认按扩展名判断)

预处理 (逐数据集):
  --pre=auto                   auto|none|log2|log2p1|quantile|log2pq|log2p1pq
                               可用分号逐个指定: --pre='log2p1;quantile;quantile'
    auto  = 复刻 7.R 规则: 列名以 TCG 开头 -> log2(x+1);
            以 GSM 开头 / 其他 -> (按需 log2) + normalizeBetweenArrays
    log2     = log2(x)，<=0 置 NA          log2p1 = log2(x+1)
    quantile = normalizeBetweenArrays      log2pq / log2p1pq = 对数化后再分位数归一化
  --min_expr=0                 合并后过滤行均值 <= 该值的基因 (0=不过滤)

基因匹配:
  --gene_case=asis|upper       基因名大小写: asis(默认, 忠于 7.R) | upper(统一大写)
  --toupper_scan=TRUE          额外报告“若统一大写，交集能多出多少基因”(不改数据)
  --match=intersect|union      共同基因取交集(默认) | 并集(缺失补 NA)
  --dup=max|first              同名基因去重: max=保留行均值最大(默认) | first

批次矫正:
  --method=combat|quantile|both|none   默认 combat (= sva::ComBat, 同 7.R)
  --primary=combat|quantile     method=both 时哪一份作为主产物 <prefix>_merged.csv
  --combat_prior=TRUE          ComBat par.prior
  --combat_mean_only=FALSE     ComBat mean.only (只校正均值)
  --combat_ref=                参考批次(批次标签)，该批次不被改动
  --combat_mod=group|none      group=有 --group_file 时用 mod=model.matrix(~group)
                               保护生物学变异(默认); none=不保护

批次与分组:
  --batch=TCGA,GEO1,GEO2       自定义批次标签(长度须等于数据集数)；默认每个数据集=一个批次
  --batch_file=file.csv        样本->批次 两列表(覆盖 --batch)
  --group_file=file.csv        样本->分组 两列表(用于 --combat_mod=group 与分组着色)
  --group_col=                 分组文件列名(多于两列时指定)
  --group_na=error|drop        分组缺失样本: error(默认) | drop(从合并结果剔除)
  --clinical_files=A,B,C       各数据集的临床表(01 步的 clinical_<GSE>.csv)，逗号分隔、
                               与 --inputs 一一对应(某个数据集没有就留空)。
                               产出 <prefix>_clinical.csv (合并临床表, 样本名为行名,
                               可直接给 02 步 --clinical=) 与 <prefix>_grouping_candidates.txt
                               (各数据集的分组候选列与取值, 用来向用户索要分组)。

QC 图表:
  --pca_top=0                  PCA 用变异最大的前 N 个基因 (0=全部)
  --pca_npc=10                 PCA 求前几个主成分 (只需画 PC1/PC2，默认 10 足够)
  --pca_color=batch|group      PCA 着色依据
  --label_samples=FALSE        PCA 是否标注样本名 (仅 classic 风格)
  --ellipse=TRUE               是否画 95% 置信椭圆
  --density=TRUE               是否输出密度图
  --density_max=200            密度图最多画多少条样本曲线(过多会卡且看不清)
  --cor_heatmap=FALSE          是否额外输出样本相关性热图
  --cor_max=150                相关性热图最多用多少样本

输出与风格:
  --dir=.                     输出目录        --prefix=merged   文件名前缀
  --save_txt=TRUE              是否额外输出制表符版 .txt (02 步只吃 .csv 也能跑)
  --save_raw=TRUE              是否输出未矫正矩阵文件 (<prefix>_merged_raw.*)
                               两者都不影响 <prefix>.RData (里面始终保留全部矩阵)
  --fig_style=classic|nature  绘图风格 (nature 见 references/nature-figure-style.md)
  --fig_width=7.5 --fig_height=5.5      图尺寸(英寸); 亦可用 --fig_width_mm/--fig_height_mm
  --pca_width/--pca_height --box_width/--box_height   单独覆盖某类图尺寸
  --hm_width=9 --hm_height=8            相关性热图尺寸(英寸)

产出 (在 --dir 下):
  <prefix>_merged.csv/.txt        矫正后合并矩阵 (第二步 02_deg_plots.R 的输入)
  <prefix>_merged_raw.csv         cbind 合并、未矫正的矩阵 (QC 对照)
  <prefix>_batch_map.csv          样本 -> 数据集/批次/分组 映射 (先核对这张表)
  <prefix>_clinical.csv           合并临床信息 (需 --clinical_files; 样本名为行名,
                                  可直接给 02 步 --clinical=)
  <prefix>_grouping_candidates.txt 各数据集的分组候选列与取值 -> 用来向用户索要分组
  <prefix>_group_template.csv     分组模板 (sample,group; 全部样本已列好, 已推导出
                                  分组时预填) -> 用户填好/或用户自己的分组表都能用
                                  --group_file= 直接喂给 02 步
  <prefix>_overlap_report.txt     各数据集基因数、逐级交集、大小写敏感性
  <prefix>_batch_stats.csv        每个批次矫正前后的均值/中位数/标准差
  <prefix>_pca_before/.after      批次效应 PCA (按批次/分组着色)
  <prefix>_boxplot_before/.after  样本分布箱线图 (按批次着色)
  <prefix>_density_before/.after  样本密度曲线
  <prefix>.RData                  合并矩阵与批次/分组信息
", sep = "")
}

if (isTRUE(OPT$help) || identical(tolower(as.character(OPT$help)), "true")) { usage(); quit(status = 0) }

# 输出目录自动创建 (否则 write.csv 会失败)
if (!dir.exists(OPT$dir)) {
  dir.create(OPT$dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(OPT$dir)) stop("无法创建输出目录: ", OPT$dir)
}

## ------------------------------ 依赖检查 ----------------------------------
need <- c("limma", "ggplot2")
hard_missing <- need[!vapply(need, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
if (length(hard_missing)) {
  stop("缺少必需包: ", paste(hard_missing, collapse = ", "),
       "\n  安装: BiocManager::install(c('limma')); install.packages('ggplot2')")
}
if (tolower(OPT$method) %in% c("combat", "both") && !requireNamespace("sva", quietly = TRUE)) {
  stop("--method=", OPT$method, " 需要 sva 包。\n",
       "  安装: BiocManager::install('sva')   或改用 --method=quantile")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  msg("  提示: 未安装 data.table，大文件读取会明显变慢 (install.packages('data.table'))")
}
if (!requireNamespace("RSpectra", quietly = TRUE)) {
  msg("  提示: 未安装 RSpectra，PCA 将退回 prcomp（样本数多时明显更慢，",
      "15913 基因 x 967 样本实测要 5~10 分钟）。建议 install.packages('RSpectra')")
}

## ------------------------------ 工具函数 ----------------------------------

# 读矩阵: 第一列为基因名, 其余为样本。返回 numeric matrix。
# 不能用 read.csv(row.names=1): 基因名可能是字面量 "NA"，会被 na.strings 转成 NA，
# 随后报 "missing values in 'row.names' are not allowed"。故整体读入再单独取行名。
# na.strings 刻意不含 "NA"，把 "NA" 当行名字面量保留；数值列的空值靠 as.numeric 转 NA。
read_matrix_robust <- function(path, sep = "auto") {
  if (!file.exists(path)) stop("找不到输入矩阵: ", path)
  if (identical(sep, "auto")) sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "csv" else "tab"
  d <- NULL
  if (requireNamespace("data.table", quietly = TRUE)) {
    d <- tryCatch(
      as.data.frame(data.table::fread(
        path, sep = if (identical(sep, "csv")) "," else "\t",
        header = TRUE, check.names = FALSE, data.table = FALSE,
        na.strings = c("", "NaN"), showProgress = FALSE)),
      error = function(e) NULL)
  }
  if (is.null(d)) {
    d <- if (identical(sep, "csv")) {
      read.csv(path, header = TRUE, check.names = FALSE, na.strings = "")
    } else {
      read.table(path, header = TRUE, sep = "\t", check.names = FALSE, na.strings = "")
    }
  }
  if (ncol(d) < 2L) stop("表达矩阵至少需要 2 列（基因名 + 至少 1 个样本）: ", path)
  rn <- as.character(d[[1]])
  m  <- suppressWarnings(matrix(as.numeric(as.matrix(d[, -1, drop = FALSE])), nrow = nrow(d)))
  rownames(m) <- rn
  colnames(m) <- colnames(d)[-1]
  if (anyNA(rownames(m))) stop("基因名中有 ", sum(is.na(rownames(m))), " 个 NA，请检查输入矩阵第一列: ", path)
  m
}

# 与 01 步一致的对数化启发式: 判定“看起来还没 log”
needs_log2 <- function(m) {
  ex <- m[is.finite(m)]
  if (!length(ex)) return(FALSE)
  qx <- as.numeric(stats::quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
  (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
}

qsumm <- function(m) {
  qx <- stats::quantile(m[is.finite(m)], c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  sprintf("min=%.2f Q1=%.2f 中位=%.2f Q3=%.2f Q99=%.2f max=%.2f",
          qx[1], qx[2], qx[3], qx[4], qx[5], qx[6])
}

# 同名基因去重: max=保留行均值最大的那一行 (等价于“保留整体表达量最高的探针”)
dedupe_rows <- function(m, how = "max") {
  if (!anyDuplicated(rownames(m))) return(m)
  n0 <- nrow(m)
  if (identical(how, "max")) {
    m <- m[order(rowMeans(m, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]
  }
  m <- m[!duplicated(rownames(m)), , drop = FALSE]
  msg("      同名基因去重: ", n0, " -> ", nrow(m), " 行 (策略 ", how, ")")
  m
}

## ------------------- 合并临床信息 (供用户确定下游分组) ----------------------
# 为什么需要这一步: 合并矩阵只是数值拼起来，各数据集的临床注释还在各自的
# clinical_<GSE>.csv 里。下游做差异分析必须由用户指定分组，而用户要看的就是
# "每个数据集有哪些临床列、取值怎么分布、同一取值在别的数据集里叫什么"。
# 所以这里把 N 份临床表并成一份，并生成一份候选报告。
#
# 列名策略(重要): 各数据集往往来自不同研究，列名与**取值词表**都不一样，
# 例如 source_name_ch1 在 A 里是 "normal control"、在 B 里是
# "synovial tissue from healthy joint"。因此:
#   * 所有(提供了临床表的)数据集都有的同名列 -> 保留原名
#     (值仍是各数据集自己的词表，需人工核对后再统一映射)
#   * 只出现在部分数据集的列 -> 加 "<数据集>::" 前缀，避免同名互相覆盖
# 值一律按字符处理并 trim: 这张表是给人看的，不在这里做数值推断。
read_clinical_file <- function(path) {
  if (!nzchar(path)) return(NULL)
  if (!file.exists(path)) {
    msg("  ! --clinical_files 里的文件不存在，已跳过该数据集: ", path)
    return(NULL)
  }
  d <- tryCatch(read.csv(path, header = TRUE, check.names = FALSE, na.strings = ""),
                error = function(e) NULL)
  if (is.null(d) || ncol(d) < 2L) {
    msg("  ! 临床表格式不对 (至少需要 样本名列 + 1 个表型列)，已跳过: ", path)
    return(NULL)
  }
  rn <- as.character(d[[1]]); d <- d[, -1L, drop = FALSE]
  if (anyDuplicated(rn)) {
    msg("  ! 临床表里有重复样本名，已加后缀区分: ", path)
    rn <- make.unique(rn)
  }
  rownames(d) <- rn
  d[] <- lapply(d, function(x) trimws(as.character(x)))
  d
}

# 候选分组列: 与 01 步同一套启发式 —— 取值 2..min(10, n/2) 个、每类至少 2 个。
# 返回 table (值 -> 计数) 或 NULL。
clin_candidates <- function(v) {
  v <- trimws(as.character(v))
  v[is.na(v) | v %in% c("NA", "null", "NULL", "---")] <- ""
  tb <- table(v[nzchar(v)])
  if (length(tb) < 2L) return(NULL)
  n <- sum(nzchar(v))
  if (length(tb) > min(10L, max(2L, floor(n / 2)))) return(NULL)
  if (min(tb) < 2L || max(tb) > n - 2L) return(NULL)
  tb
}

# 某个取值当前落到了哪个统一分组 (用来审计 --group_maps 的映射是否符合预期)
clin_grp_ann <- function(val, v, g) {
  if (is.null(g)) return("")
  gg <- unique(g[v == val & !is.na(v)])
  gg <- gg[!is.na(gg) & nzchar(gg)]
  if (!length(gg)) return("   -> (未分组/已剔除)")
  if (length(gg) == 1L) return(paste0("   -> ", gg))
  paste0("   -> 混合(", paste(sort(gg), collapse = "/"), ")")
}

# 逐数据集预处理。返回 list(mat=, action=)
preprocess_one <- function(m, spec) {
  spec <- tolower(trimws(spec))
  if (!spec %in% c("auto", "none", "log2", "log2p1", "quantile", "log2pq", "log2p1pq")) {
    stop("未知的 --pre 取值: ", spec, " (可选 auto|none|log2|log2p1|quantile|log2pq|log2p1pq)")
  }
  if (identical(spec, "auto")) {
    p3 <- substr(colnames(m)[1], 1, 3)
    unlogged <- needs_log2(m)
    if (identical(p3, "TCG")) {
      # TCGA 表达量 (counts/TPM/FPKM): 只取对数，不再分位数归一化 (同 7.R)
      spec <- if (unlogged) "log2p1" else "none"
      note <- sprintf("auto[TCG]: 判定%s -> %s", if (unlogged) "未取对数" else "已取对数",
                      if (unlogged) "log2(x+1)" else "不做处理")
    } else {
      spec <- if (unlogged) "log2pq" else "quantile"
      note <- sprintf("auto[%s]: 判定%s -> %s + normalizeBetweenArrays", p3,
                      if (unlogged) "未取对数" else "已取对数",
                      if (unlogged) "log2(x)" else "仅")
    }
  } else {
    note <- paste0("指定 ", spec)
  }

  if (spec %in% c("log2", "log2p1")) {
    if (identical(spec, "log2p1")) { m[m < 0] <- NA; m <- log2(m + 1) }
    else { m[m <= 0] <- NA; m <- log2(m) }
  } else if (spec %in% c("quantile", "log2pq", "log2p1pq")) {
    if (spec %in% c("log2pq", "log2p1pq")) {
      if (identical(spec, "log2p1pq")) { m[m < 0] <- NA; m <- log2(m + 1) }
      else { m[m <= 0] <- NA; m <- log2(m) }
    }
    m <- limma::normalizeBetweenArrays(m)
  }
  list(mat = m, action = note)
}

## -------------------- [1/8] 读取各数据集 ------------------------------------
cat("\n================ 第三步: 多数据集合并 + 批次矫正 ================\n")
if (!nzchar(OPT$inputs)) {
  msg("缺少 --inputs。")
  usage()
  quit(status = 1)
}
inputs <- trimws(strsplit(OPT$inputs, ",", fixed = TRUE)[[1]])
inputs <- inputs[nzchar(inputs)]
if (length(inputs) < 2L) {
  stop("--inputs 至少需要 2 个矩阵 (收到 ", length(inputs), " 个)。只处理一个矩阵无需本步骤。")
}
n_ds <- length(inputs)

dsets <- if (nzchar(OPT$names)) trimws(strsplit(OPT$names, ",", fixed = TRUE)[[1]]) else character(0)
if (length(dsets) && length(dsets) != n_ds) {
  stop("--names 个数 (", length(dsets), ") 与 --inputs 个数 (", n_ds, ") 不一致")
}
if (!length(dsets)) dsets <- sub("\\.[^.]*$", "", basename(inputs))
dsets <- make.unique(dsets)

pres <- trimws(strsplit(OPT$pre, ";", fixed = TRUE)[[1]]); pres <- pres[nzchar(pres)]
if (length(pres) == 1L) pres <- rep(pres, n_ds)
if (length(pres) != n_ds) {
  stop("--pre 个数 (", length(pres), ") 必须是 1 或与 --inputs 个数 (", n_ds, ") 相等")
}

msg("[1/8] 读取 ", n_ds, " 个表达矩阵 ...")
mats <- vector("list", n_ds)
n_gene_raw <- integer(n_ds)
for (i in seq_len(n_ds)) {
  m <- read_matrix_robust(inputs[i], OPT$sep)
  mats[[i]] <- m
  n_gene_raw[i] <- nrow(m)
  msg(sprintf("  %-16s %6d 基因 x %4d 样本   [%s]", dsets[i], nrow(m), ncol(m), basename(inputs[i])))
}

## -------------------- [2/8] 逐数据集预处理 ---------------------------------
msg("[2/8] 逐数据集预处理 (--pre) ...")
pre_actions <- character(n_ds)
n_gene_used <- integer(n_ds)
for (i in seq_len(n_ds)) {
  before <- mats[[i]]
  res <- preprocess_one(before, pres[i])
  msg("  ", dsets[i], ": ", res$action)
  msg("      对数化前 ", qsumm(before))
  mats[[i]] <- res$mat
  msg("      处理后   ", qsumm(mats[[i]]))
  pre_actions[i] <- res$action
  rm(before); if (i %% 2 == 0) invisible(gc(FALSE))
}

## -------------------- [3/8] 基因名统一 + 交集 ------------------------------
msg("[3/8] 基因名统一并取共同基因 ...")
if (identical(OPT$gene_case, "upper")) {
  msg("  --gene_case=upper: 全部基因名转为大写后匹配")
} else {
  msg("  --gene_case=asis: 基因名原样匹配 (忠于 7.R)；交集偏少时可改用 --gene_case=upper")
}

clean_mat <- function(m, dset) {
  rn <- trimws(as.character(rownames(m)))
  if (identical(OPT$gene_case, "upper")) rn <- toupper(rn)
  bad <- is.na(rn) | !nzchar(rn) | rn %in% c("---", "NA", "na", "NULL")
  if (any(bad)) {
    msg("      ", dset, ": 剔除 ", sum(bad), " 个无名/占位基因名 (空, ---, NA)")
    m <- m[!bad, , drop = FALSE]; rn <- rn[!bad]
  }
  rownames(m) <- rn
  dedupe_rows(m, OPT$dup)
}
for (i in seq_len(n_ds)) {
  mats[[i]] <- clean_mat(mats[[i]], dsets[i])
  n_gene_used[i] <- nrow(mats[[i]])
}
invisible(gc(FALSE))

# 逐级交集，保持第 1 个数据集的基因顺序 (与 7.R 的 Reduce(intersect, ...) 一致)
common_steps <- integer(n_ds)
common <- rownames(mats[[1]]); common_steps[1] <- length(common)
if (n_ds > 1) for (i in 2:n_ds) {
  common <- intersect(common, rownames(mats[[i]]))
  common_steps[i] <- length(common)
}
if (!length(common)) {
  stop("共同基因为 0 个。请检查基因命名是否一致，或改用 --gene_case=upper / --match=union")
}
if (identical(OPT$match, "union")) {
  genes_final <- unique(unlist(lapply(mats, rownames)))
  msg("  --match=union: 使用并集 ", length(genes_final), " 个基因 (交集 ", length(common), " 个)")
} else {
  genes_final <- common
}

# 大小写敏感性探查: 只报告，不改数据
upper_gain <- NA_integer_
if (isTRUE(as.logical(OPT$toupper_scan)) && identical(OPT$gene_case, "asis")) {
  cu <- Reduce(intersect, lapply(mats, function(m) toupper(rownames(m))))
  upper_gain <- length(cu) - length(common)
  rm(cu); invisible(gc(FALSE))
}

## -------------------- [4/8] 合并 -------------------------------------------
msg("[4/8] 合并矩阵 (--match=", OPT$match, ") ...")
subs <- lapply(mats, function(m) {
  out <- matrix(NA_real_, nrow = length(genes_final), ncol = ncol(m),
                dimnames = list(genes_final, colnames(m)))
  hit <- genes_final %in% rownames(m)
  out[hit, ] <- m[match(genes_final[hit], rownames(m)), , drop = FALSE]
  out
})
ds_counts <- vapply(subs, ncol, integer(1))
ds_vec    <- rep(dsets, times = ds_counts)
merged_raw <- do.call(cbind, subs)
rm(subs, mats); invisible(gc(FALSE))

if (anyDuplicated(colnames(merged_raw))) {
  n_dup <- sum(duplicated(colnames(merged_raw)))
  msg("  ! 样本名有 ", n_dup, " 个重复 (不同数据集出现同名样本)，已加后缀区分。")
  msg("    这会影响 --group_file / --batch_file 的匹配，建议先在上游把样本名改唯一。")
  colnames(merged_raw) <- make.unique(colnames(merged_raw))
}

if (as.numeric(OPT$min_expr) > 0) {
  keep <- rowMeans(merged_raw, na.rm = TRUE) > as.numeric(OPT$min_expr)
  msg("  --min_expr=", OPT$min_expr, ": 过滤 ", sum(!keep), " 个低表达基因")
  merged_raw <- merged_raw[keep, , drop = FALSE]
}

na_n <- sum(is.na(merged_raw))
if (na_n > 0) {
  kind <- tolower(as.character(OPT$na))
  if (identical(kind, "error")) {
    worst <- names(sort(rowSums(is.na(merged_raw)), decreasing = TRUE))[seq_len(min(5, nrow(merged_raw)))]
    stop("合并矩阵含 ", na_n, " 个 NA (占 ", round(100 * na_n / length(merged_raw), 2), "%)。\n",
         "  sva::ComBat 不接受 NA。缺失最多的基因: ", paste(worst, collapse = ", "), "\n",
         "  处理方式:\n",
         "    --na=rowmean  用该基因的行均值填补 (最常用)\n",
         "    --na=zero     置为 0\n",
         "    或 --match=intersect 只保留所有数据集都测到的基因")
  } else if (identical(kind, "rowmean")) {
    rm_ <- rowMeans(merged_raw, na.rm = TRUE); rm_[!is.finite(rm_)] <- 0
    idx <- which(is.na(merged_raw), arr.ind = TRUE)
    merged_raw[idx] <- rm_[idx[, 1]]
    msg("  --na=rowmean: 填补 ", na_n, " 个 NA")
  } else if (identical(kind, "zero")) {
    merged_raw[is.na(merged_raw)] <- 0
    msg("  --na=zero: 填补 ", na_n, " 个 NA")
  } else {
    stop("未知的 --na 取值: ", OPT$na, " (可选 error|rowmean|zero)")
  }
} else {
  msg("  合并矩阵无 NA")
}
msg("  合并结果: ", nrow(merged_raw), " 基因 x ", ncol(merged_raw), " 样本")

## -------------------- [5/8] 批次 + 分组 ------------------------------------
msg("[5/8] 构建批次与分组信息 ...")
batch <- factor(ds_vec, levels = unique(ds_vec))

if (nzchar(OPT$batch)) {
  bl <- trimws(strsplit(OPT$batch, ",", fixed = TRUE)[[1]])
  if (length(bl) != n_ds) stop("--batch 个数 (", length(bl), ") 须等于数据集数 (", n_ds, ")")
  batch <- factor(rep(bl, times = ds_counts), levels = unique(rep(bl, times = ds_counts)))
  msg("  使用自定义批次标签: ", paste(levels(batch), collapse = ", "))
}

if (nzchar(OPT$batch_file)) {
  bf <- read.csv(OPT$batch_file, header = TRUE, check.names = FALSE)
  if (ncol(bf) < 2L) stop("--batch_file 需要两列: 样本名, 批次")
  idx <- match(colnames(merged_raw), as.character(bf[[1]]))
  if (all(is.na(idx))) stop("--batch_file 第一列与合并矩阵样本名无一匹配")
  bv <- as.character(bf[[2]])[idx]
  if (anyNA(bv)) {
    stop("--batch_file 未覆盖 ", sum(is.na(bv)), " 个样本，例如: ",
         paste(utils::head(colnames(merged_raw)[is.na(bv)]), collapse = ", "))
  }
  batch <- factor(bv, levels = unique(bv))
  msg("  使用 --batch_file 的批次划分: ", paste(levels(batch), collapse = ", "))
}
msg("  批次数: ", nlevels(batch), " -> ", paste(levels(batch), collapse = ", "))

group <- NULL
if (nzchar(OPT$group_file)) {
  gf <- read.csv(OPT$group_file, header = TRUE, check.names = FALSE)
  if (ncol(gf) < 2L) stop("--group_file 需要两列: 样本名, 分组")
  gcol <- if (nzchar(OPT$group_col)) {
    j <- match(OPT$group_col, colnames(gf))
    if (is.na(j)) stop("--group_col='", OPT$group_col, "' 不在文件列中: ", paste(colnames(gf), collapse = ", "))
    j
  } else 2L
  idx <- match(colnames(merged_raw), as.character(gf[[1]]))
  group <- as.character(gf[[gcol]])[idx]
  n_miss <- sum(is.na(group))
  msg("  分组文件匹配: ", sum(!is.na(group)), "/", ncol(merged_raw), " 个样本",
      if (n_miss) paste0("，", n_miss, " 个未匹配") else "")
  if (n_miss) {
    if (identical(tolower(as.character(OPT$group_na)), "drop")) {
      msg("  --group_na=drop: 从合并矩阵中剔除 ", n_miss, " 个无分组样本")
      keepc <- !is.na(group)
      merged_raw <- merged_raw[, keepc, drop = FALSE]
      ds_vec <- ds_vec[keepc]
      batch  <- factor(as.character(batch[keepc]), levels = unique(as.character(batch[keepc])))
      group  <- group[keepc]
    } else {
      stop("有 ", n_miss, " 个样本在分组文件中缺失 (例如 ",
           paste(utils::head(colnames(merged_raw)[is.na(group)]), collapse = ", "), ")。\n",
           "  用 --group_na=drop 剔除这些样本，或补齐分组文件。")
    }
  }
  msg("  分组水平: ", paste(names(table(group)), collapse = ", "))
}

## -------------------- [6/8] 批次矫正 ---------------------------------------
msg("[6/8] 批次矫正 (--method=", OPT$method, ") ...")
do_combat   <- tolower(as.character(OPT$method)) %in% c("combat", "both")
do_quantile <- tolower(as.character(OPT$method)) %in% c("quantile", "both")

combat_mod <- NULL
if (do_combat && identical(tolower(as.character(OPT$combat_mod)), "group")) {
  if (is.null(group)) {
    msg("  ! --combat_mod=group 但未提供 --group_file，将不保护生物学变异 (mod=NULL)")
  } else {
    grp <- factor(group)
    # ComBat 内部: design <- cbind(batchmod, mod); check <- apply(design,2,all(x==1)); 去掉全 1 列
    # 因此必须传带截距的 model.matrix(~grp)：截距会被自动丢掉，剩下恰好 1 列(两组时)。
    # 若写成 ~0+grp 会得到 k 列虚拟变量，与 batch 列共线 -> 误报 "covariates are confounded"。
    batchmod <- stats::model.matrix(~ -1 + batch)
    dchk <- as.matrix(cbind(batchmod, stats::model.matrix(~ grp)))
    dchk <- as.matrix(dchk[, !apply(dchk, 2, function(x) all(x == 1)), drop = FALSE])
    if (qr(dchk)$rank < ncol(dchk)) {
      ct <- capture.output(print(table(batch = batch, group = group)))
      stop("分组与批次混杂: 用 mod=~group 时设计矩阵秩不足，ComBat 会拒绝运行。\n",
           "  批次 x 分组 交叉表:\n    ", paste(ct, collapse = "\n    "), "\n",
           "  含义: 某些批次内只有一个分组水平，无法把“批次效应”与“生物学差异”分开估计。\n",
           "  可选做法:\n",
           "    1) --combat_mod=none 不做生物变异保护 (仅当批次间分组构成相似时才安全)\n",
           "    2) 只合并分组构成相近的数据集\n",
           "    3) --method=quantile (不做批次参数估计，不受此限制)")
    }
    combat_mod <- stats::model.matrix(~ grp)
    msg("  ComBat 保护生物学变异: mod = model.matrix(~ group)，分组水平 ",
        paste(levels(grp), collapse = "/"))
    msg("  批次 x 分组 交叉表:")
    for (ln in capture.output(print(table(batch = batch, group = group)))) msg("    ", ln)
  }
}

ref_batch <- if (nzchar(OPT$combat_ref)) OPT$combat_ref else NULL
merged <- merged_combat <- merged_quantile <- NULL

if (do_combat) {
  t0 <- Sys.time()
  merged_combat <- sva::ComBat(
    dat = merged_raw, batch = batch, mod = combat_mod,
    par.prior = isTRUE(as.logical(OPT$combat_prior)),
    mean.only = isTRUE(as.logical(OPT$combat_mean_only)),
    ref.batch = ref_batch)
  rownames(merged_combat) <- rownames(merged_raw); colnames(merged_combat) <- colnames(merged_raw)
  msg("  ComBat 完成，用时 ", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " 秒")
}
if (do_quantile) {
  merged_quantile <- limma::normalizeBetweenArrays(merged_raw)
  rownames(merged_quantile) <- rownames(merged_raw); colnames(merged_quantile) <- colnames(merged_raw)
  msg("  normalizeBetweenArrays 完成")
}

if (do_combat && do_quantile) {
  merged <- if (identical(tolower(as.character(OPT$primary)), "quantile")) merged_quantile else merged_combat
  msg("  --method=both: 主产物取 ", tolower(as.character(OPT$primary)),
      " (另一份另存为 ", OPT$prefix, "_merged_", setdiff(c("combat", "quantile"), tolower(as.character(OPT$primary))), ".csv)")
} else if (do_combat) {
  merged <- merged_combat
} else if (do_quantile) {
  merged <- merged_quantile
} else {
  merged <- merged_raw
  msg("  --method=none: 未做批次矫正，主产物 = 合并原始矩阵")
}

## -------------------- [7/8] QC 图表 ----------------------------------------
msg("[7/8] 生成 QC 图表 ...")
suppressMessages(library(ggplot2))
style <- tolower(as.character(OPT$fig_style))
thm <- if (identical(style, "nature")) theme_nature(7) else theme_classic_std(14)

fig_sz <- fig_size(OPT$fig_width, OPT$fig_height, OPT$fig_width_mm, OPT$fig_height_mm)
pick_sz <- function(w, h) c(if (nzchar(w)) as.numeric(w) else fig_sz[1],
                            if (nzchar(h)) as.numeric(h) else fig_sz[2])
pca_sz <- pick_sz(OPT$pca_width, OPT$pca_height)
box_sz <- pick_sz(OPT$box_width, OPT$box_height)

meta <- data.frame(sample = colnames(merged_raw), dataset = ds_vec,
                   batch = as.character(batch), stringsAsFactors = FALSE)
rownames(meta) <- meta$sample
if (!is.null(group)) meta$group <- group

pal_for <- function(lv) stats::setNames(rep(BATCH_COLS, length.out = length(lv)), lv)

# --- PCA ---
# ⚠ 性能关键: 不要直接用 prcomp(t(m), center=TRUE)。
#   prcomp 会对 n(样本) x p(基因) 的矩阵求出**全部** min(n,p) 个奇异向量，
#   LAPACK dgesdd 的代价约 (n+p)*min(n,p)^2。15913 基因 x 967 样本 时约 3e11 flops，
#   参考 BLAS 下单次就要 5~10 分钟，而这里要跑两次(前后各一次)——实测卡住过。
#   改用 RSpectra::svds 只求前 npc 个奇异向量，代价 O(n*p*npc)，同样的数据 <2 秒。
#   前若干主成分的方差占比与 prcomp 一致(见 references/merge-batch-correction.md)。
pca_of <- function(m, top = 0, npc = 10) {
  ok <- apply(m, 1, function(x) all(is.finite(x)) && stats::sd(x) > 0)
  m <- m[ok, , drop = FALSE]
  if (!nrow(m)) stop("所有基因都是常数或含 NA，无法做 PCA")
  if (top > 0 && nrow(m) > top) {
    m <- m[order(apply(m, 1, stats::var), decreasing = TRUE)[seq_len(top)], , drop = FALSE]
  }
  Xc <- scale(t(m), center = TRUE, scale = FALSE)   # 样本 x 基因，已中心化
  tv <- sum(apply(Xc, 2, stats::var))               # 总方差 = sum(sdev^2)
  npc <- max(2L, min(as.integer(npc), nrow(Xc) - 1L, ncol(Xc)))
  res <- NULL
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    sv <- tryCatch({
      set.seed(42)   # ARPACK 起始向量随机，固定种子以保证结果可复现
      RSpectra::svds(Xc, k = npc)
    }, error = function(e) NULL)
    if (!is.null(sv) && length(sv$d) >= 2L) {
      o <- order(sv$d, decreasing = TRUE)
      u <- sv$u[, o, drop = FALSE]
      # ⚠ RSpectra::svds 的 u / v **不保留 dimnames**（实测 rownames/colnames 全为 NULL），
      #   而 stats::prcomp()$x 是保留的。不补回来的话，下游
      #   `meta[rownames(pca$x), ]` 会取到 0 行，报
      #   "arguments imply differing number of rows: 967, 0"。
      rownames(u) <- rownames(Xc)
      colnames(u) <- paste0("PC", seq_len(ncol(u)))
      res <- list(sdev = sv$d[o] / sqrt(nrow(Xc) - 1L),
                  x = u, total_var = tv, engine = "RSpectra::svds")
    }
  }
  if (is.null(res)) {
    p <- stats::prcomp(Xc, center = FALSE, scale. = FALSE)
    k <- min(npc, length(p$sdev))
    res <- list(sdev = p$sdev[seq_len(k)], x = p$x[, seq_len(k), drop = FALSE],
                total_var = sum(p$sdev^2), engine = "prcomp")
  }
  res
}
pca_top <- as.integer(OPT$pca_top)
t0 <- Sys.time()
pca_before <- pca_of(merged_raw, pca_top, OPT$pca_npc)
pca_after  <- pca_of(merged,     pca_top, OPT$pca_npc)
msg("  PCA 完成 (", pca_before$engine, ")，用时 ",
    round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " 秒")
ve_b <- round(100 * pca_before$sdev^2 / pca_before$total_var, 1)
ve_a <- round(100 * pca_after$sdev^2  / pca_after$total_var, 1)

plot_pca <- function(pca, ve, meta, color_by, ttl) {
  # pca$x 的行顺序与 colnames(输入矩阵) 一致，按行对齐更稳（不依赖 dimnames）
  sc <- pca$x
  if (!identical(rownames(sc), rownames(meta))) {
    if (is.null(rownames(sc)) && nrow(sc) == nrow(meta)) rownames(sc) <- rownames(meta)
    else stop("PCA 得分 (", nrow(sc), " 行) 与样本表 (", nrow(meta), " 行) 对不上")
  }
  df <- data.frame(PC1 = sc[, 1], PC2 = sc[, 2],
                   meta[rownames(sc), , drop = FALSE], stringsAsFactors = FALSE)
  cb <- if (identical(color_by, "group")) "group" else "batch"
  df$grpvar <- df[[cb]]
  # ⚠ 图内文字一律用 ASCII: Windows 下 ggplot2 默认字体族不含中文字形，
  #   中文标签会被渲染成乱码（ragg/cairo 都不例外）。见 pitfalls 5.5 / 7.8。
  legend_name <- if (identical(cb, "batch")) "Batch" else "Group"
  p <- ggplot(df, aes(x = PC1, y = PC2, colour = grpvar)) +
    geom_point(size = if (identical(style, "nature")) 0.9 else 2, alpha = 0.85) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
         title = ttl, colour = legend_name) +
    scale_colour_manual(values = pal_for(levels(factor(df$grpvar))), name = legend_name) +
    thm
  if (isTRUE(as.logical(OPT$ellipse)) && nlevels(factor(df$grpvar)) <= 12 &&
      min(table(df$grpvar)) >= 3) {
    p <- p + stat_ellipse(level = 0.95, linewidth = 0.3, show.legend = FALSE)
  }
  if (isTRUE(as.logical(OPT$label_samples)) && identical(style, "classic") &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(aes(label = sample), size = 2.4, max.overlaps = 30)
  }
  p
}

# 批次效应的客观量化: 批次对 PC1/PC2 解释的方差比例 + ANOVA p 值
batch_r2 <- function(pca, b) {
  res <- t(vapply(1:2, function(k) {
    fit <- stats::lm(pca$x[, k] ~ b)
    c(r2 = summary(fit)$r.squared, p = stats::anova(fit)$`Pr(>F)`[1])
  }, numeric(2)))
  rownames(res) <- c("PC1", "PC2")
  res
}
r2_b <- batch_r2(pca_before, batch)
r2_a <- batch_r2(pca_after,  batch)

pca_color <- tolower(as.character(OPT$pca_color))
save_figure(plot_pca(pca_before, ve_b, meta, pca_color, "PCA - before batch correction"),
            file.path(OPT$dir, paste0(OPT$prefix, "_pca_before")), pca_sz[1], pca_sz[2], style)
save_figure(plot_pca(pca_after, ve_a, meta, pca_color, "PCA - after batch correction"),
            file.path(OPT$dir, paste0(OPT$prefix, "_pca_after")), pca_sz[1], pca_sz[2], style)
if (!is.null(group) && identical(pca_color, "batch")) {
  save_figure(plot_pca(pca_after, ve_a, meta, "group", "PCA - after correction (by group)"),
              file.path(OPT$dir, paste0(OPT$prefix, "_pca_after_group")), pca_sz[1], pca_sz[2], style)
}
msg(sprintf("  批次效应: 批次对 PC1 解释方差 %.1f%% (p=%.3g) -> %.1f%% (p=%.3g)",
            100 * r2_b[1, 1], r2_b[1, 2], 100 * r2_a[1, 1], r2_a[1, 2]))

# --- 箱线图 (按样本分位数直接绘制，避免把大矩阵融成长表) ---
box_q <- function(m, meta) {
  q <- apply(m, 2, function(x) stats::quantile(x, c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
  d <- data.frame(sample = colnames(m), t(q), check.names = FALSE)
  colnames(d) <- c("sample", "ymin", "lower", "middle", "upper", "ymax")
  d$batch <- meta[d$sample, "batch"]
  d
}
# 每个样本一个箱体。样本多时箱体很窄，因此**不画边框** (colour = NA)，否则灰色边框
# 会盖住批次填充色、整张图糊成一片（实测 967 样本时必现）。
plot_box <- function(m, meta, ttl) {
  d <- box_q(m, meta)
  d$sample <- factor(d$sample, levels = d$sample)
  ggplot(d, aes(x = sample, ymin = ymin, lower = lower, middle = middle,
                upper = upper, ymax = ymax, fill = batch)) +
    geom_boxplot(stat = "identity", width = 1, linewidth = 0,
                 outlier.shape = NA, colour = NA) +
    labs(x = NULL, y = "Expression (log2 scale)", title = ttl, fill = "Batch") +
    scale_fill_manual(values = pal_for(levels(batch)), name = "Batch") +
    thm + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
# 每个批次一个箱体（箱体来自该批次内各样本的中位数）。样本数很大时这张图始终清晰，
# 用来回答"各批次的整体水平是否被拉齐了"。
plot_box_batch <- function(m, meta, ttl) {
  d <- box_q(m, meta)
  ggplot(d, aes(x = batch, y = middle, fill = batch)) +
    geom_boxplot(width = 0.5, linewidth = 0.3, outlier.shape = NA, colour = "grey25") +
    labs(x = NULL, y = "Per-sample median (log2 scale)", title = ttl, fill = "Batch") +
    scale_fill_manual(values = pal_for(levels(batch)), name = "Batch") +
    thm   # x 轴标签保持水平: 旋转后容易被图片底边裁掉
}
save_figure(plot_box(merged_raw, meta, "Sample distribution - before batch correction"),
            file.path(OPT$dir, paste0(OPT$prefix, "_boxplot_before")), box_sz[1], box_sz[2], style)
save_figure(plot_box(merged, meta, "Sample distribution - after batch correction"),
            file.path(OPT$dir, paste0(OPT$prefix, "_boxplot_after")), box_sz[1], box_sz[2], style)
save_figure(plot_box_batch(merged_raw, meta, "Per-sample medians by batch - before correction"),
            file.path(OPT$dir, paste0(OPT$prefix, "_boxplot_batch_before")), box_sz[1], box_sz[2], style)
save_figure(plot_box_batch(merged, meta, "Per-sample medians by batch - after correction"),
            file.path(OPT$dir, paste0(OPT$prefix, "_boxplot_batch_after")), box_sz[1], box_sz[2], style)

# --- 密度图 (手动算 density，避免融化大矩阵；样本过多时抽样) ---
if (isTRUE(as.logical(OPT$density))) {
  dmax <- as.integer(OPT$density_max)
  sel <- colnames(merged_raw)
  if (dmax > 0 && length(sel) > dmax) {
    set.seed(1); sel <- sort(sample(sel, dmax))
    msg("  密度图: 样本数 ", ncol(merged_raw), " > --density_max=", dmax,
        "，随机抽 ", dmax, " 条曲线 (set.seed(1) 可复现)")
  }
  dens_df <- function(m) {
    parts <- lapply(sel, function(s) {
      x <- m[, s]; x <- x[is.finite(x)]
      if (length(x) < 3 || stats::sd(x) == 0) return(NULL)
      dn <- stats::density(x, n = 256)
      data.frame(x = dn$x, y = dn$y, sample = s,
                 batch = as.character(meta[s, "batch"]), stringsAsFactors = FALSE)
    })
    do.call(rbind, parts)
  }
  plot_dens <- function(m, ttl) {
    d <- dens_df(m)
    ggplot(d, aes(x = x, y = y, group = sample, colour = batch)) +
      geom_line(linewidth = 0.25, alpha = 0.5) +
      labs(x = "Expression (log2 scale)", y = "Density", title = ttl, colour = "Batch") +
      scale_colour_manual(values = pal_for(levels(batch)), name = "Batch") + thm
  }
  save_figure(plot_dens(merged_raw, "Sample density - before batch correction"),
              file.path(OPT$dir, paste0(OPT$prefix, "_density_before")), box_sz[1], box_sz[2], style)
  save_figure(plot_dens(merged, "Sample density - after batch correction"),
              file.path(OPT$dir, paste0(OPT$prefix, "_density_after")), box_sz[1], box_sz[2], style)
}

# --- 样本相关性热图 (可选) ---
if (isTRUE(as.logical(OPT$cor_heatmap))) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    msg("  ! 未安装 pheatmap，跳过相关性热图 (install.packages('pheatmap'))")
  } else {
    cmax <- as.integer(OPT$cor_max)
    selc <- colnames(merged)
    if (cmax > 0 && length(selc) > cmax) {
      set.seed(1)
      per <- ceiling(cmax / nlevels(batch))
      selc <- unlist(lapply(split(selc, as.character(batch)), function(s)
        if (length(s) <= per) s else sample(s, per)), use.names = FALSE)
      selc <- utils::head(selc, cmax)
      msg("  相关热图: 按批次抽样 ", length(selc), " / ", ncol(merged), " 个样本")
    }
    cc  <- stats::cor(merged[, selc, drop = FALSE], method = "pearson")
    ord <- order(as.character(batch)[match(selc, colnames(merged))])
    ann <- data.frame(batch = as.character(batch)[match(selc, colnames(merged))])
    rownames(ann) <- selc
    draw_it <- function() print(pheatmap::pheatmap(
      cc[ord, ord], annotation_col = ann, annotation_row = ann,
      annotation_colors = list(batch = pal_for(levels(batch))),
      cluster_rows = FALSE, cluster_cols = FALSE,
      show_rownames = FALSE, show_colnames = FALSE, border_color = NA, silent = TRUE,
      color = grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
      main = "Sample correlation (after correction)"))
    # pheatmap 返回 gtable，ggsave 可能产出空白文件而不报错 -> 必须用设备包裹
    hw <- as.numeric(OPT$hm_width); hh <- as.numeric(OPT$hm_height)
    grDevices::pdf(file.path(OPT$dir, paste0(OPT$prefix, "_cor_heatmap.pdf")), width = hw, height = hh)
    draw_it(); grDevices::dev.off()
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(file.path(OPT$dir, paste0(OPT$prefix, "_cor_heatmap.png")),
                    width = hw, height = hh, units = "in", res = 300)
    } else {
      grDevices::png(file.path(OPT$dir, paste0(OPT$prefix, "_cor_heatmap.png")),
                     width = hw, height = hh, units = "in", res = 300, type = "cairo")
    }
    draw_it(); grDevices::dev.off()
  }
}

## -------------------- [8/8] 导出 -------------------------------------------
msg("[8/8] 写出结果 ...")
save_txt <- isTRUE(as.logical(OPT$save_txt))
save_raw <- isTRUE(as.logical(OPT$save_raw))
n_cell <- as.numeric(nrow(merged)) * ncol(merged)
if (save_txt && save_raw && n_cell > 5e6) {
  msg("  提示: 矩阵较大 (", nrow(merged), " x ", ncol(merged), ")，写 4 个文本文件较慢；",
      "不需要未矫正矩阵/制表符版可加 --save_raw=FALSE / --save_txt=FALSE (RData 里始终保留)")
}
write_matrix <- function(m, base) {
  # .csv 是 02_deg_plots.R 的输入格式 (read_matrix 同时支持 .csv 与 .txt)
  # ⚠ 用 data.table::fwrite 而不是 write.csv: 15913 x 967 的矩阵 write.csv 要
  #   分钟级，fwrite 快约一个数量级（与 fread vs read.table 同理）。
  d <- data.frame(ID = rownames(m), m, check.names = FALSE)
  if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fwrite(d, paste0(base, ".csv"))
    if (save_txt) data.table::fwrite(d, paste0(base, ".txt"), sep = "\t")
  } else {
    write.csv(m, paste0(base, ".csv"))
    if (save_txt) {
      write.table(d, file = paste0(base, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}
write_matrix(merged, file.path(OPT$dir, paste0(OPT$prefix, "_merged")))
if (save_raw) {
  write_matrix(merged_raw, file.path(OPT$dir, paste0(OPT$prefix, "_merged_raw")))
}
if (do_combat && do_quantile) {
  write_matrix(merged_combat,   file.path(OPT$dir, paste0(OPT$prefix, "_merged_combat")))
  write_matrix(merged_quantile, file.path(OPT$dir, paste0(OPT$prefix, "_merged_quantile")))
}

# 样本 -> 数据集/批次/分组 映射 (先核对这张表再信结果)
bm <- data.frame(sample = colnames(merged), dataset = ds_vec,
                 batch = as.character(batch), check.names = FALSE)
if (!is.null(group)) bm$group <- group
write.csv(bm, file.path(OPT$dir, paste0(OPT$prefix, "_batch_map.csv")), row.names = FALSE)

# 单独再导一份 "sample,group" 两列文件，可直接喂给 02_deg_plots.R 的 --group_file。
# (_batch_map.csv 的第 2 列是 dataset 不是分组，不能直接当 --group_file 用。)
if (!is.null(group)) {
  write.csv(data.frame(sample = colnames(merged), group = group, check.names = FALSE),
            file.path(OPT$dir, paste0(OPT$prefix, "_group.csv")), row.names = FALSE)
}

# 分组模板: 给用户"填好就能用"的一张表 (两列 sample,group，样本顺序与合并矩阵一致)。
# 已推导出分组时预填，用户改几格即可；没有分组时全空。格式与 01 步的
# <GSE>_group_template.csv 一致。用户也可能根本不改，直接上传他自己手上的分组表。
# NA 写成空串而不是字面量 "NA"，读回来时会按"未分组"处理。
tmpl_group <- if (is.null(group)) rep("", ncol(merged)) else ifelse(is.na(group), "", group)
write.csv(data.frame(sample = colnames(merged), group = tmpl_group, check.names = FALSE),
          file.path(OPT$dir, paste0(OPT$prefix, "_group_template.csv")), row.names = FALSE)
msg("  分组模板    : ", OPT$prefix, "_group_template.csv (", ncol(merged),
    " 个样本已列好; 用户自己给了分组表的话也可以直接用 --group_file= 喂进来)")

## -------------------- 合并临床信息 + 分组候选报告 ---------------------------
# 下游做差异分析必须由用户指定分组，而用户需要看到"每个数据集有哪些临床列、
# 取值怎么分布、同一取值在不同数据集里叫什么"。这里把 N 份临床表并成一份，
# 并生成一份可直接交给用户的候选报告。
clinical <- NULL
cf_out <- file.path(OPT$dir, paste0(OPT$prefix, "_clinical.csv"))
gf_out <- file.path(OPT$dir, paste0(OPT$prefix, "_grouping_candidates.txt"))

clin_paths <- character(0)
if (nzchar(OPT$clinical_files)) {
  clin_paths <- trimws(strsplit(OPT$clinical_files, ",", fixed = TRUE)[[1]])
  if (length(clin_paths) < n_ds) clin_paths <- c(clin_paths, rep("", n_ds - length(clin_paths)))
  if (length(clin_paths) > n_ds) {
    msg("  ! --clinical_files 给了 ", length(clin_paths), " 个，数据集只有 ", n_ds, " 个，多余的忽略")
    clin_paths <- clin_paths[seq_len(n_ds)]
  }
}

if (any(nzchar(clin_paths))) {
  smp <- colnames(merged)
  clins <- lapply(clin_paths, read_clinical_file)
  cols_of <- lapply(clins, function(d) if (is.null(d)) character(0) else colnames(d))
  have <- !vapply(clins, is.null, logical(1))
  uniq <- unique(unlist(cols_of))
  n_have <- if (length(uniq)) {
    vapply(uniq, function(cc) sum(vapply(cols_of[have], function(x) cc %in% x, logical(1))), integer(1))
  } else integer(0)
  # 注意: 这里必须用独立变量名。早期版本复用了 `common`，会覆盖上面"基因交集"的
  # 结果，导致报告里 "最终共同基因" 与 "大小写敏感性" 打印成"共同临床列"的个数
  # （2026-09-17 实测: 18036 个基因被打印成 31）。
  common_cols <- uniq[n_have == sum(have)]

  cdf <- data.frame(dataset = ds_vec, batch = as.character(batch), stringsAsFactors = FALSE)
  if (!is.null(group)) cdf$group <- group
  added <- list()
  for (i in seq_len(n_ds)) {
    d <- clins[[i]]; if (is.null(d)) next
    sel <- which(ds_vec == dsets[i]); if (!length(sel)) next
    j <- match(smp[sel], rownames(d))
    for (cc in colnames(d)) {
      nm <- if (cc %in% common_cols) cc else paste0(dsets[i], "::", cc)
      # 别撞上前面的 dataset/batch/group 三列（临床表里真有可能有叫 group 的列）
      if (nm %in% c("sample", "dataset", "batch", "group")) nm <- paste0(dsets[i], "::", cc)
      v <- rep(NA_character_, length(smp))
      v[sel] <- d[[cc]][j]
      if (!is.null(added[[nm]])) {          # 同名列: 各数据集各自填充
        old <- added[[nm]]; k <- !is.na(v); old[k] <- v[k]; v <- old
      }
      added[[nm]] <- v
    }
  }
  if (length(added)) {
    cdf <- cbind(cdf, as.data.frame(added, stringsAsFactors = FALSE, check.names = FALSE))
  }
  rownames(cdf) <- smp                        # 行名=样本名 -> 02 步 --clinical= 可直接读
  clinical <- cdf
  write.csv(cdf, cf_out, row.names = TRUE)

  ## ---- 候选报告 (给用户看，用来决定分组) ----
  L <- character(0); add <- function(...) L <<- c(L, paste0(...))
  W <- 70
  add("合并矩阵的分组候选报告 —— 供人工确定差异分析的分组")
  add(strrep("=", W))
  add(sprintf("生成时间        : %s", format(Sys.time())))
  add(sprintf("合并矩阵        : %d 基因 x %d 样本", nrow(merged), ncol(merged)))
  add(sprintf("参与合并的数据集: %d 个 -> %s", n_ds, paste(dsets, collapse = ", ")))
  add(sprintf("批次            : %s", paste(levels(batch), collapse = ", ")))
  add(sprintf("合并临床表      : %s_clinical.csv (样本名为行名)", OPT$prefix))
  add(sprintf("统一分组        : %s", if (is.null(group)) "尚未提供 (下游跑 02 之前必须补上)"
                                        else paste(names(table(group)), collapse = ", ")))
  add("")
  add("【这份报告是干什么的】")
  add("  合并后的矩阵只有数值。分组是生物学判断，脚本不会替你选。")
  add("  请把本报告 + <prefix>_clinical.csv 一起交给用户，让他确认三件事:")
  add("    1) 每个数据集用哪一列做分组 (下表 [n] 标出的候选列)")
  add("    2) 该列哪些取值算 case、哪些算 control、哪些样本要剔除")
  add("    3) 对比方向 (contrast 里谁做分子、谁做分母)")
  add("  或者他直接给你一张自己的分组表 (两列 样本名,分组) —— 见下方\"下一步 D\"。")
  add("  拿到答复后再跑 02 步；若只想换分组而不想重算批次矫正，先跑")
  add("  03b_merge_from_geo.R --stage=group 重导出分组表。")
  add("  另外: 用户也可以直接给一张他自己的分组表（两列 样本名,分组）——")
  add("  见下方\"下一步 D\"，那是最省事的回答方式。")
  add("")
  add("---------------- 各数据集的候选分组列 ----------------")
  for (i in seq_len(n_ds)) {
    g <- dsets[i]; d <- clins[[i]]
    sel <- which(ds_vec == g)                  # 该数据集进入最终矩阵的样本
    add("")
    if (is.null(d)) {
      add(sprintf("[%s] 未提供临床表 (本矩阵内 %d 个样本)", g, length(sel)))
      add("      -> 无法参与分组推导，只能靠 --group_file 手工给分组。")
      next
    }
    j    <- match(smp[sel], rownames(d))       # 样本 -> 临床表行号
    n_c  <- nrow(d); n_k <- sum(!is.na(j))
    kept <- rep(FALSE, n_c); kept[stats::na.omit(j)] <- TRUE
    add(sprintf("[%s]  临床表 %d 个样本 | 进入最终矩阵 %d 个%s | 表型列 %d 个",
                g, n_c, n_k,
                if (n_c > n_k) sprintf(" (%d 个未进入矩阵)", n_c - n_k) else "",
                ncol(d)))
    gg_keep <- if (is.null(group)) NULL else group[sel]
    k <- 0L; hi <- character(0)
    for (cc in colnames(d)) {
      # 候选判定用该数据集的**完整词表**(临床表全部样本): 某个取值整类被
      # --group_na=drop 剔除时，只看"幸存样本"会以为这列只有 1 个取值而跳过它，
      # 正是用户最需要看到的那种情况。
      vall <- d[[cc]]
      tb   <- clin_candidates(vall)
      if (is.null(tb)) {
        vv <- trimws(vall); vv <- vv[!is.na(vv) & nzchar(vv)]
        if (length(vv)) {
          nu <- length(unique(vv))
          hi <- c(hi, sprintf("%s(%d 个取值%s)", cc, nu,
                              if (nu < 2L) ", 单一" else ", 过多"))
        }
        next
      }
      k <- k + 1L
      add(sprintf("   [%d] %s   (%d 个唯一值, 合计 %d)", k, cc, length(tb), sum(tb)))
      # 取值名字可能很长（GEO 的 source_name_ch1 常见 40+ 字符），按本列最长值
      # 自适应对齐，否则计数会被挤到右边、整块看起来错位
      wv <- min(52L, max(44L, max(nchar(names(tb)))))
      fmtv <- paste0("        %-", wv, "s %5d%s%s")
      for (nm in names(tb)) {
        n_all <- as.integer(tb[[nm]])
        inmat <- which(vall == nm); n_in <- sum(kept[inmat])
        ann <- if (!n_in) "   -> (未进入矩阵, 已剔除)"
               else clin_grp_ann(nm, d[[cc]][j], gg_keep)
        add(sprintf(fmtv, nm, n_all,
                    if (n_in < n_all) sprintf(" [保留 %d/%d]", n_in, n_all) else "", ann))
      }
    }
    if (k == 0L) {
      add("   (没有取值数适中的列，无法直接作分组；可用 --group_regex 从下面的列提取)")
    }
    if (length(hi)) {
      add(sprintf("   其他列 (不直接作分组): %s", paste(utils::head(hi, 6), collapse = ", ")))
      if (length(hi) > 6L) add(sprintf("   ... 另有 %d 列", length(hi) - 6L))
    }
  }
  add("")
  add("---------------- 跨数据集的同名临床列 ----------------")
  if (length(uniq)) {
    w <- max(nchar(dsets))
    shared <- uniq[n_have >= 2L]
    # 标出"在至少一个数据集里算候选分组列"的列: 这一节大多是平台/联系方式等
    # 元数据，没有标记根本看不出哪几列值得看。
    cand_flag <- vapply(uniq, function(cc) {
      any(vapply(seq_along(clins), function(i) {
        d <- clins[[i]]
        !is.null(d) && cc %in% colnames(d) && !is.null(clin_candidates(d[[cc]]))
      }, logical(1)))
    }, logical(1))
    if (length(shared)) {
      add("  (出现在 >= 2 个数据集里的列; 带 * 的在某个数据集里是候选分组列)")
      add(sprintf("  %-34s %s", "column",
                  paste(sprintf(paste0("%-", w, "s"), dsets), collapse = " ")))
      for (cc in shared) {
        marks <- vapply(cols_of, function(x) if (cc %in% x) "Y" else "-", character(1))
        marks[!have] <- "?"
        nm <- substr(if (isTRUE(cand_flag[cc])) paste0(cc, " *") else cc, 1, 34)
        add(sprintf("  %-34s %s", nm,
                    paste(sprintf(paste0("%-", w, "s"), marks), collapse = " ")))
      }
    } else {
      add("  (没有任何列出现在 2 个以上数据集里)")
    }
    # 数据集专属列只列名字：它们不可能"一个列名走遍所有数据集"
    for (i in seq_len(n_ds)) {
      own <- uniq[n_have == 1L & vapply(uniq, function(cc) cc %in% cols_of[[i]], logical(1))]
      if (!length(own)) next
      own <- ifelse(cand_flag[own], paste0(own, " *"), own)
      add(sprintf("  %s 专属 %d 列: %s", dsets[i], length(own),
                  paste(utils::head(own, 8), collapse = ", ")))
      if (length(own) > 8L) add(sprintf("    ... 另有 %d 列", length(own) - 8L))
    }
    add("  (* = 该列在至少一个数据集里是候选分组列; Y=该数据集有这列; ?=该数据集未提供临床表;")
    add(sprintf("   只出现在部分数据集的列在 %s_clinical.csv 里带 \"<数据集>::\" 前缀)", OPT$prefix))
    add("  提示: 即使同一列名在各数据集都有，取值词表也常常不同，")
    add("        必须逐数据集核对后再决定是否统一映射。")
  } else {
    add("  (没有任何数据集提供临床表)")
  }
  add("")
  if (is.null(group)) {
    add("---------------- 统一分组: 尚未提供 ----------------")
    add("  跑 02 步之前必须先补上。三条路:")
    add("   A) 用 03b 重导出分组表 (复用已合并的矩阵，不重算 ComBat):")
    add("      03b_merge_from_geo.R --stage=group \\")
    add("        --group_cols='<数据集1的列>;<数据集2的列>' \\")
    add("        --group_maps='<取值A>=case,<取值B>=control;<同理>'")
    add("   B) 直接用合并临床表的某一列:")
    add(sprintf("      02_deg_plots.R --expr=%s_merged.csv --clinical=%s_clinical.csv \\",
                OPT$prefix, OPT$prefix))
    add("        --group_from=<列名> [--group_regex=<re>] [--group_keep=A,B]")
    add("   C) 让用户直接给一张他自己的分组表 (两列 样本名,分组):")
    add(sprintf("      可直接填 %s_group_template.csv，或提供他手上的表 -> --group_file=<表>",
                OPT$prefix))
  } else {
    add("---------------- 统一分组 (合并时已提供) ----------------")
    tb <- table(group)
    for (nm in names(tb)) add(sprintf("  %-24s %5d 个样本", nm, as.integer(tb[[nm]])))
    if (any(is.na(group))) add(sprintf("  %-24s %5d 个样本", "(未分组)", sum(is.na(group))))
    add("  对照上文每个候选列后面的 \"-> 分组名\" 即可核对映射是否符合预期。")
    add("  ! 该分组同时被用作 ComBat 的 mod=~group。若换一个分组做差异分析，")
    add("    批次矫正本身不必重算，但要让用户知道 ComBat 保护的仍是原分组。")
  }
  add("")
  add("---------------- 下一步 ----------------")
  add("  A) 用统一分组做差异分析:")
  add(sprintf("     02_deg_plots.R --expr=%s_merged.csv --group_file=%s_group.csv --contrast=<组2>-<组1>",
              OPT$prefix, OPT$prefix))
  add("  B) 换一个分组 (不重算批次矫正): 03b_merge_from_geo.R --stage=group ...")
  add(sprintf("  C) 直接用合并临床表某列: 02_deg_plots.R --clinical=%s_clinical.csv --group_from=<列>",
              OPT$prefix))
  add("  D) 用户自己给的分组表 (两列 样本名,分组; 表头可有可无):")
  add(sprintf("     02_deg_plots.R --expr=%s_merged.csv --group_file=<用户的分组表> --contrast=<组2>-<组1>",
              OPT$prefix))
  add(sprintf("     可直接让他填 %s_group_template.csv (全部样本已列好)，", OPT$prefix))
  add("     或提供他手上已有的表。样本名要与合并矩阵列名一致；")
  add("     没覆盖到的样本会被 02 步自动剔除 (不是错误)。")
  add("     合并阶段也能用: 03_merge_batches.R --group_file=<用户的分组表>，")
  add("     这样 ComBat 的 mod 会跟着保护这个分组。")
  add(strrep("=", W))
  writeLines(L, gf_out)

  msg("  合并临床表  : ", basename(cf_out), " (", nrow(cdf), " 样本 x ", ncol(cdf), " 列)")
  msg("  分组候选报告: ", basename(gf_out))
  msg("  ! 分组是生物学判断: 跑 02 步之前必须把这两份文件交给用户，")
  msg("    让他确认用哪一列、哪些取值算 case/control、对比方向。")
}

# 批次统计
bstat <- do.call(rbind, lapply(levels(batch), function(b) {
  i <- which(as.character(batch) == b)
  data.frame(batch = b, n_sample = length(i),
             mean_before = mean(merged_raw[, i]), median_before = stats::median(merged_raw[, i]),
             sd_before = stats::sd(merged_raw[, i]),
             mean_after = mean(merged[, i]), median_after = stats::median(merged[, i]),
             sd_after = stats::sd(merged[, i]), stringsAsFactors = FALSE)
}))
write.csv(bstat, file.path(OPT$dir, paste0(OPT$prefix, "_batch_stats.csv")), row.names = FALSE)

# 交集 / 大小写敏感性报告
# 用 sprintf("%.3g") 而非 formatC(): formatC(format="g") 会按宽度右对齐填充空格，
# 在表里表现为 "p=   0" 这种多余空白。
fmt_p <- function(p) if (is.null(p) || !is.finite(p)) "NA" else sprintf("%.3g", p)
lbl <- sprintf("  %-16s %8d %8d  %s", dsets, n_gene_used, ds_counts, pre_actions)
stp <- sprintf("  前 %d 个数据集交集: %d 个基因", seq_len(n_ds), common_steps)
# 批次统计表手工排版: 默认 print() 在 80 列控制台上会折成两段
stat_hdr <- sprintf("  %-16s %6s  %-9s %-9s %-9s | %-9s %-9s %s",
                    "batch", "n", "mean", "median", "sd", "mean", "median", "sd")
stat_sub <- sprintf("  %-16s %6s  %-9s %-9s %-9s | %-9s %-9s %s",
                    "", "", "(矫正前)", "", "", "(矫正后)", "", "")
stat_rows <- sprintf("  %-16s %6d  %-9.3f %-9.3f %-9.3f | %-9.3f %-9.3f %.3f",
                     bstat$batch, bstat$n_sample,
                     bstat$mean_before, bstat$median_before, bstat$sd_before,
                     bstat$mean_after,  bstat$median_after,  bstat$sd_after)
rep_lines <- c(
  "多数据集合并 - 基因匹配与批次报告",
  strrep("=", 68),
  sprintf("生成时间              : %s", format(Sys.time())),
  sprintf("合并策略              : %s", OPT$match),
  sprintf("基因名大小写          : %s", OPT$gene_case),
  sprintf("批次矫正方法          : %s", OPT$method),
  sprintf("  ComBat 参数         : par.prior=%s, mean.only=%s, ref=%s",
          OPT$combat_prior, OPT$combat_mean_only, ifelse(is.null(ref_batch), "无", ref_batch)),
  sprintf("  ComBat 保护生物变异 : %s", if (is.null(combat_mod)) "否 (mod=NULL)" else "是 (mod=~group)"),
  "",
  sprintf("  %-16s %8s %8s  %s", "数据集", "基因数", "样本数", "预处理"),
  sprintf("  %-16s %8s %8s  %s", strrep("-", 16), strrep("-", 8), strrep("-", 8), strrep("-", 44)),
  lbl,
  "",
  "逐级交集 (保持第 1 个数据集的基因顺序):",
  stp,
  "",
  sprintf("最终共同基因          : %d", length(common)),
  sprintf("最终使用基因          : %d (--match=%s)", length(genes_final), OPT$match),
  sprintf("合并矩阵规模          : %d 基因 x %d 样本", nrow(merged), ncol(merged)),
  sprintf("批次数                : %d (%s)", nlevels(batch), paste(levels(batch), collapse = ", ")),
  "",
  "基因名大小写敏感性:",
  if (is.na(upper_gain)) "  未扫描 (--toupper_scan=FALSE 或已是 --gene_case=upper)"
  else sprintf("  当前大小写交集 = %d；统一为大写后 = %d (%+d)",
               length(common), length(common) + upper_gain, upper_gain),
  if (!is.na(upper_gain) && upper_gain > 0)
    "  -> 数据集跨平台(RNA-seq 与芯片混用)时建议加 --gene_case=upper 以保留更多基因" else "",
  "",
  "批次效应量化 (批次对主成分解释的方差比例, ANOVA by batch):",
  sprintf("  %-5s  %-18s %s", "", "矫正前", "矫正后"),
  sprintf("  %-5s  %-18s %s", "PC1",
          sprintf("%.1f%% (p=%s)", 100 * r2_b[1, 1], fmt_p(r2_b[1, 2])),
          sprintf("%.1f%% (p=%s)", 100 * r2_a[1, 1], fmt_p(r2_a[1, 2]))),
  sprintf("  %-5s  %-18s %s", "PC2",
          sprintf("%.1f%% (p=%s)", 100 * r2_b[2, 1], fmt_p(r2_b[2, 2])),
          sprintf("%.1f%% (p=%s)", 100 * r2_a[2, 1], fmt_p(r2_a[2, 2]))),
  "",
  "批次统计 (各批次全部基因的均值/中位数/sd; 矫正后应彼此接近):",
  stat_hdr, stat_sub, stat_rows,
  "",
  "下一步:",
  sprintf("  bash run_geo.sh 02_deg_plots.R --expr=%s_merged.csv --group_file=<分组文件> --contrast=<组2>-<组1>",
          OPT$prefix),
  sprintf("  分组文件可用 %s_group.csv / %s_group_template.csv，或用户自己的分组表",
          OPT$prefix, OPT$prefix)
)
if (!is.null(clinical)) {
  rep_lines <- c(rep_lines,
    "",
    "临床信息 (供用户确定差异分析的分组):",
    sprintf("  %s_clinical.csv            合并临床表 (样本名为行名, 02 步 --clinical= 可直接读)",
            OPT$prefix),
    sprintf("  %s_grouping_candidates.txt 各数据集候选分组列与取值, 交给用户选分组",
            OPT$prefix))
}
if (is.null(group)) {
  rep_lines <- c(rep_lines,
    "",
    "!! 尚未提供分组: 做差异分析(02 步)之前必须向用户索要分组,",
    "   把上面的临床表与候选报告交给他, 不要自行决定。")
}
writeLines(rep_lines, file.path(OPT$dir, paste0(OPT$prefix, "_overlap_report.txt")))

save(merged, merged_raw, merged_combat, merged_quantile, batch, group, meta,
     clinical, pca_before, pca_after, dsets, ds_counts, inputs, OPT,
     file = file.path(OPT$dir, paste0(OPT$prefix, ".RData")))

cat("\n================ 第三步完成 ================\n")
cat("合并矩阵    : ", nrow(merged), " 基因 x ", ncol(merged), " 样本\n", sep = "")
cat("批次数      : ", nlevels(batch), " -> ", paste(levels(batch), collapse = ", "), "\n", sep = "")
cat("矫正方法    : ", OPT$method, "\n", sep = "")
cat("主产物      : ", file.path(OPT$dir, paste0(OPT$prefix, "_merged.csv")), "\n", sep = "")
cat("样本批次映射: ", file.path(OPT$dir, paste0(OPT$prefix, "_batch_map.csv")), " (先核对这张表)\n", sep = "")
cat("分组模板    : ", file.path(OPT$dir, paste0(OPT$prefix, "_group_template.csv")),
    " (用户可直接填, 或改用自己的分组表)\n", sep = "")
cat("匹配报告    : ", file.path(OPT$dir, paste0(OPT$prefix, "_overlap_report.txt")), "\n", sep = "")
if (!is.null(clinical)) {
  cat("合并临床表  : ", cf_out, "\n", sep = "")
  cat("分组候选报告: ", gf_out, "\n", sep = "")
}
cat("\n下一步 (第二步，分组须由用户指定):\n")
cat("  bash run_geo.sh 02_deg_plots.R --expr=", OPT$prefix, "_merged.csv ",
    "--group_file=<分组文件> --contrast=<组2>-<组1>\n", sep = "")
cat("\n!! 不要自行决定分组:\n")
if (!is.null(clinical)) {
  cat("   把 ", basename(cf_out), " 与 ", basename(gf_out), " 交给用户,\n", sep = "")
  cat("   让他确认 ①用哪一列 ②哪些取值算 case/control ③对比方向；\n", sep = "")
} else {
  cat("   先向用户索要分组（本次没传 --clinical_files，没有临床表可参考）；\n", sep = "")
}
cat("   或者让他直接给一张自己的分组表（两列 样本名,分组;\n")
cat("   可填 ", OPT$prefix, "_group_template.csv）-> 02 步 --group_file=，\n", sep = "")
cat("   样本名与矩阵列名一致即可，未覆盖的样本会被 02 步自动剔除。\n")
