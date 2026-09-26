#!/usr/bin/env Rscript
# ===========================================================================
# 02_deg_plots.R  ——  limma 差异分析 + 热图 + 火山图
#
# 输入: 01_geo_normalize.R 产出的标准化表达矩阵 (基因 x 样本)
# 输出: 差异基因表、火山图、热图、RData
#
# 用法:
#   # A) 分组文件 (两列: 样本名, 分组)
#   Rscript 02_deg_plots.R --expr=GSE55584.csv --group_file=group.csv \
#           --contrast=CASE-CONTROL
#
#   # B) 从临床信息表某列正则提取 (如 title 形如 "N1"/"T1")
#   Rscript 02_deg_plots.R --expr=GSE55584.csv --clinical=clinical_GSE55584.csv \
#           --group_from=title --group_regex='^(.)' --contrast=T-N
#
#   # C) 直接给分组向量 (顺序对应表达矩阵列)
#   Rscript 02_deg_plots.R --expr=GSE55584.csv \
#           --group_values=N,N,N,T,T,T --contrast=T-N
# ===========================================================================

options(stringsAsFactors = FALSE, warn = 1)

# 兼容性: 中文用户名 / 部分 R 构建下，输出缓冲会在(偶发的)崩溃时整段丢失，
# 导致脚本“无任何输出、且无产物”。此处把 message 包一层强制 flush.console()，
# 既不影响结果，又能让每步进度即时落盘、崩溃点可被定位。
msg <- function(...) { base::message(...); flush.console() }

## --------------------- Nature 期刊风格绘图组件 (可选) -----------------------
# 火山图 / 热图的主题、配色与导出约定，改编自开源项目 nature-skills
# (github.com/Yuan1z0825/nature-skills, Apache-2.0) 的 nature-figure 技能。
# 仅在 --fig_style=nature 时启用；默认 classic 完全保持原有输出不变。
NATURE_PALETTE <- c(
  blue_main      = "#0F4D92",   # 深蓝 —— 主方法/关键对象
  blue_secondary = "#3775BA",   # 中蓝
  red_strong     = "#B64342",   # 强调红
  neutral_light  = "#CFCECE",   # 中性浅灰
  neutral_mid    = "#767676",   # 中性中灰
  neutral_dark   = "#4D4D4D"    # 中性深灰
)
# 火山图三态配色 (取自 nature-figure 已验证 volcano 模板: 蓝=下调, 红=上调, 灰=不显著)
NATURE_VOLCANO_COLS <- c(UP = "#B2182B", DOWN = "#2166AC", NOT = "#B3B3B3")

# 期刊主题: classic 版式 + Arial + 细轴线 + 无网格 + 5~7pt 字号
theme_nature <- function(base_size = 7, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line    = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.ticks   = ggplot2::element_line(linewidth = 0.35, colour = "black"),
      axis.title   = ggplot2::element_text(size = base_size),
      axis.text    = ggplot2::element_text(size = base_size - 0.5, colour = "black"),
      legend.title = ggplot2::element_blank(),
      legend.text  = ggplot2::element_text(size = base_size - 0.7),
      legend.key   = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(size = base_size + 0.5, face = "bold", hjust = 0.5),
      panel.grid   = ggplot2::element_blank()
    )
}

# 多格式导出: PDF/PNG 始终输出; nature 风格下额外输出 SVG (可编辑文字) + 600dpi TIFF
# 注意: 默认 pdf() 设备在未注册 Arial 时会报 "invalid font type"，必须用 cairo_pdf；
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
  expr = "", clinical = "", group_file = "", group_from = "", group_regex = "",
  group_values = "", group_spec = "", group_keep = "",
  contrast = "", logfc = "1", pval = "0.05",
  p_type = "P.Value", sort_by = "B", top = "25", prefix = "DEG",
  label_genes = "", label_top = "0", fig_style = "classic", fig_title = "",
  vol_style = "same", hm_style = "same", top_n = "10",
  list_styles = "FALSE", ask_style = "FALSE",
  hm_width = "10", hm_height = "8", heatmap_png = "TRUE",
  vol_width = "7", vol_height = "6",
  vol_width_mm = "", vol_height_mm = "", hm_width_mm = "", hm_height_mm = "",
  help = "FALSE"
))

usage <- function() {
  cat("
02_deg_plots.R —— limma 差异分析 + 热图 + 火山图

用法: Rscript 02_deg_plots.R --expr=<矩阵.csv> [分组方式] [选项]

必填
  --expr=<file>            标准化表达矩阵 (基因 x 样本)，.csv 或制表符 .txt

分组方式 (四选一)
  --group_file=<file>      两列 CSV: 样本名, 分组 (表头可有可无)
                           只列出部分样本也可以，未列出的样本会被剔除
  --group_spec='A=GSM1,GSM2;B=GSM3,GSM4'   直接把分组写在命令行里，组名=样本列表，
                           多组用分号分隔 (最直观，推荐交互场景使用)
  --clinical=<file>        01 步产出的 clinical_*.csv
      --group_from=<col>   取哪一列，如 title / characteristics_ch1 / source_name_ch1
      --group_regex=<re>   带一个捕获组的正则，如 '^(.)' 或 '.*(Tumor|Normal).*'
  --group_values=case,case,control,...   按表达矩阵列顺序直接给出

  --group_keep=A,B         只保留指定的分组做分析 (分组多于 2 个时用)

分析选项
  --contrast=T-N           对比式，必须是 make.names 清洗后的分组名 (默认 第2组-第1组)
  --logfc=1                |log2FC| 阈值
  --pval=0.05              p 值阈值
  --p_type=P.Value         用哪一列卡阈值: P.Value | adj.P.Val
  --sort_by=B              topTable 排序依据: B | P | logFC | none
  --top=25                 热图展示的 top 基因数
  --label_genes=MMP1,VIT   火山图上标注的基因，逗号分隔
  --label_top=0            火山图自动标注 p 值最小的前 N 个显著基因 (0=关闭)
  --prefix=DEG             输出文件名前缀

绘图风格
  --list_styles=TRUE       只列出全部可用风格 (火山图 + 热图, 含依赖与适用条件)，不画图
  --fig_style=classic      绘图风格: classic (默认, 原版主题) | nature (期刊风格)
                           nature: 期刊配色 + 图例计数 + 细轴线 + 额外导出 SVG/600dpi TIFF
  --vol_style=same         火山图风格覆盖 (不写则跟随 --fig_style):
                           same      跟随 --fig_style (classic/nature)
                           enhanced  EnhancedVolcano 四区分色 (需 EnhancedVolcano 包, 未装自动回退 same)
                           gradient  渐变气泡: 颜色/大小随显著性连续渐变 (仿 ggVolcano::gradual_volcano)
                           rainbow   五段彩虹渐变: 深蓝->青->黄->橙红->深红, 大小也随 -log10P
                           tophits   Top-Hits 距离排序 (仿 VolcaNoseR): 按 |logFC|+|-log10P| 选 Top N 标注
  --top_n=10               tophits 风格标注的 Top N 基因数 (也用于 gradient/rainbow 的标签数)
  --hm_style=same          热图风格覆盖 (不写则跟随 --fig_style):
                           same      跟随 --fig_style (classic/nature, pheatmap)
                           complex   ComplexHeatmap 期刊款: 顶部分组注释条 + 发散色带 + 行列自动聚类
                           numbers   pheatmap 数值标注款: 单元格内标注 z-score
                           (仅当样本数 <= 20 时可选; > 20 自动回退 fig_style 并提示)
                           tile      ggplot2 geom_tile 款: 长数据分面/自由定制, 无聚类,
                           顶部分组色条, 适合大样本矩阵 (> 20 样本推荐)
  --fig_title=<文本>       图标题; nature 风格默认不加标题
  --hm_width=10 --hm_height=8    热图尺寸(英寸); 亦可用 --hm_width_mm/--hm_height_mm
  --vol_width=7 --vol_height=6   火山图尺寸(英寸); 亦可用 --vol_width_mm/--vol_height_mm
  --heatmap_png=TRUE       热图是否额外输出 PNG

输出 (<prefix>*)
  <prefix>_sample_group_map.csv   样本->分组映射，先核对这张表再信结果
  <prefix>_all.csv         全部基因差异分析结果 (含 change 列)
  <prefix>_DEG.csv         仅显著差异基因
  <prefix>_volcano.pdf/png 火山图  (nature 风格额外: .svg/.tiff)
  <prefix>_heatmap.pdf/png top 基因热图  (nature 风格额外: .svg/.tiff)
  <prefix>.RData           expr, group, nrDEG, DEG
  <prefix>_for_enrichment.txt  交接清单 (key=value)，下游 enrichment-analysis 技能
                         用 --geo_dir=<本目录> --geo_prefix=<prefix> 直接读
")
  invisible(NULL)
}

# --list_styles=TRUE: 打印全部风格菜单后直接退出 (供用户挑选, 不画图)
# 注意: 必须在 help 检查之前 (list_styles 不需要 --expr)
if (isTRUE(as.logical(OPT$list_styles))) {
  cat("
================ 出图风格菜单 ================
[火山图 --vol_style=]
  same      跟随 --fig_style (默认)
  classic   原版主题: 三分色 + 虚线阈值 + 可选基因标注
  nature    期刊风格: 红蓝灰配色 + 图例计数 + 细轴线 (改编自 nature-skills)
  enhanced  EnhancedVolcano 四区分色 + 引线标注  [需 EnhancedVolcano 包, 未装自动回退]
  gradient  渐变气泡: 颜色/大小随显著性连续渐变  (仿 ggVolcano::gradual_volcano)
  rainbow   五段彩虹渐变: 深蓝->青->黄->橙红->深红, 大小随 -log10P
  tophits   Top-Hits 距离排序: 按 |logFC|+|-log10P| 选 Top N 标注  (仿 VolcaNoseR)
            配套: --top_n=10 (tophits/gradient/rainbow 的标注数)

[热图 --hm_style=]
  same      跟随 --fig_style (默认)
  classic   原版 pheatmap: RdYlBu 色带 + 分组注释
  nature    期刊 pheatmap: 蓝白红发散色带 + 去边框
  complex   ComplexHeatmap 期刊款: 顶部分组注释条 + 行聚类 + 按组分块  [需 ComplexHeatmap]
  numbers   数值标注款: 单元格内标 z-score  [仅 <= 20 样本可用]
  tile      geom_tile 款: ggplot2 长数据, 顶部分组色条, 无聚类  [大样本推荐]

[通用 --fig_style=]
  classic | nature  (决定默认主题与导出格式; nature 额外出 SVG/600dpi TIFF)
==============================================
")
  quit(save = "no", status = 0L)
}

# --ask_style=TRUE: 交互提示版菜单——打印"请先选择"引导语后退出。
# 实际交互由调用方 (AI/用户) 完成: 调用方把菜单呈给用户，拿到选择后带
# --vol_style/--hm_style 重跑。脚本本身不做 stdin 交互 (批处理场景不阻塞)。
if (isTRUE(as.logical(OPT$ask_style))) {
  cat("
================ 请先选择出图风格 ================
出图前请把风格菜单呈给用户挑选，再带 --vol_style/--hm_style 重跑:
  火山图: same|classic|nature|enhanced|gradient|rainbow|tophits   (--vol_style)
  热图  : same|classic|nature|complex|numbers|tile                (--hm_style)
          numbers 仅 <=20 样本可用; tile 适合大样本; complex 需 ComplexHeatmap
  标注数: --top_n=10 (tophits/gradient/rainbow)
  组合  : 不指定时跟随 --fig_style=classic|nature
完整菜单: --list_styles=TRUE
==================================================
")
  quit(save = "no", status = 0L)
}

if (isTRUE(OPT$help) || OPT$expr == "") {
  usage()
  quit(status = if (OPT$expr == "") 1L else 0L, save = "no")
}

LOGCUT  <- as.numeric(OPT$logfc)
PCUT    <- as.numeric(OPT$pval)
TOPN    <- as.integer(OPT$top)
PREFIX  <- OPT$prefix

## ---------------------------- 环境自检 --------------------------------------
if (!dir.exists(tempdir()) || file.access(tempdir(), 2L) != 0L) {
  stop("R 临时目录不可用: ", tempdir(),
       "\nWindows 中文用户名机器上的常见故障。启动 R 前先设置纯 ASCII 临时目录:",
       "\n  bash:        先把 TMP/TEMP/TMPDIR 指向纯 ASCII 目录（run_geo.sh 已自动处理）",
       "\n  PowerShell:  mkdir D:\\Rtmp; $env:TMP='D:\\Rtmp'; $env:TEMP='D:\\Rtmp'; $env:TMPDIR='D:\\Rtmp'")
}

for (pkg in c("limma", "ggplot2", "pheatmap")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("缺少 R 包: ", pkg, "\n请先安装: install.packages(c('ggplot2','pheatmap')); ",
         "BiocManager::install('limma')")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

## ---------------------------- 1. 读表达矩阵 --------------------------------
msg("[1/5] 读取表达矩阵 ...")
if (!file.exists(OPT$expr)) stop("找不到表达矩阵: ", OPT$expr)
read_matrix <- function(path) {
  # 不能用 row.names = 1: 基因名可能是字面量 "NA"，会被 na.strings 转成 NA
  # 而报 "missing values in 'row.names' are not allowed"。故先整体读入再单独取行名。
  d <- if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read.csv(path, header = TRUE, check.names = FALSE, na.strings = "")
  } else {
    read.table(path, header = TRUE, sep = "\t", check.names = FALSE, na.strings = "")
  }
  if (ncol(d) < 2L) stop("表达矩阵至少需要 2 列（基因名 + 至少 1 个样本）")
  rn <- as.character(d[[1]])
  m  <- suppressWarnings(matrix(as.numeric(as.matrix(d[, -1, drop = FALSE])),
                               nrow = nrow(d)))
  rownames(m) <- rn
  colnames(m) <- colnames(d)[-1]
  if (anyNA(rownames(m))) {
    stop("基因名中有 ", sum(is.na(rownames(m))), " 个 NA，请检查输入矩阵第一列")
  }
  m
}
expr <- read_matrix(OPT$expr)
samples <- colnames(expr)
msg("  ", nrow(expr), " 基因 x ", ncol(expr), " 样本")

## ------------------------------ 2. 分组 ------------------------------------
msg("[2/5] 构建分组 ...")
build_group <- function() {
  # D) 命令行内联分组: 组名=样本1,样本2;组名2=样本3
  if (nzchar(OPT$group_spec)) {
    parts <- trimws(strsplit(OPT$group_spec, ";", fixed = TRUE)[[1]])
    parts <- parts[nzchar(parts)]
    if (!length(parts)) stop("--group_spec 解析为空，正确形式: A=GSM1,GSM2;B=GSM3")
    g <- rep(NA_character_, length(samples))
    for (p in parts) {
      kv <- strsplit(p, "=", fixed = TRUE)[[1]]
      if (length(kv) < 2L) stop("--group_spec 每段必须是 '组名=样本1,样本2'，收到: ", p)
      gname <- trimws(kv[1])
      sms <- trimws(strsplit(paste(kv[-1L], collapse = "="), ",", fixed = TRUE)[[1]])
      sms <- sms[nzchar(sms)]
      idx <- match(sms, samples)
      if (all(is.na(idx))) {
        stop("分组 '", gname, "' 的样本名与表达矩阵列名无一匹配: ",
             paste(sms, collapse = ", "))
      }
      if (any(is.na(idx))) {
        msg("  ! 分组 '", gname, "' 中 ", sum(is.na(idx)), " 个样本未匹配，已忽略: ",
                paste(sms[is.na(idx)], collapse = ", "))
      }
      g[idx[!is.na(idx)]] <- gname
    }
    return(g)
  }
  # C) 直接给分组向量
  if (nzchar(OPT$group_values)) {
    g <- trimws(strsplit(OPT$group_values, ",")[[1]])
    if (length(g) != length(samples)) {
      stop("--group_values 有 ", length(g), " 个，表达矩阵有 ", length(samples), " 列，数量不匹配")
    }
    return(g)
  }
  # A) 分组文件
  if (nzchar(OPT$group_file)) {
    if (!file.exists(OPT$group_file)) stop("找不到分组文件: ", OPT$group_file)
    g1 <- read.csv(OPT$group_file, header = TRUE,  check.names = FALSE)
    g0 <- read.csv(OPT$group_file, header = FALSE, check.names = FALSE)
    g  <- if (sum(samples %in% g1[[1]]) >= sum(samples %in% g0[[1]])) g1 else g0
    idx <- match(samples, as.character(g[[1]]))
    if (all(is.na(idx))) stop("分组文件第一列与表达矩阵样本名无一匹配。")
    return(as.character(g[[2]])[idx])
  }
  # B) 临床表 + 正则
  if (nzchar(OPT$group_from)) {
    if (!nzchar(OPT$clinical)) stop("使用 --group_from 时必须同时给 --clinical")
    if (!file.exists(OPT$clinical)) stop("找不到临床信息文件: ", OPT$clinical)
    cl <- read.csv(OPT$clinical, header = TRUE, check.names = FALSE, row.names = 1)
    if (!OPT$group_from %in% colnames(cl)) {
      stop("临床表中无列 '", OPT$group_from, "'，可选: ", paste(colnames(cl), collapse = ", "))
    }
    v <- as.character(cl[[OPT$group_from]])
    v <- v[match(samples, rownames(cl))]
    if (nzchar(OPT$group_regex)) {
      v <- sub(OPT$group_regex, "\\1", v)
      msg("  正则 '", OPT$group_regex, "' 提取结果: ", paste(unique(v), collapse = ", "))
    }
    return(trimws(v))
  }
  stop("必须指定分组方式: --group_file / (--clinical + --group_from) / --group_values")
}

group <- build_group()

# 剔除未分组样本
drop <- is.na(group) | !nzchar(group) | group %in% c("NA", "NULL")
if (any(drop)) {
  msg("  ! 剔除未分组样本 ", sum(drop), " 个: ", paste(samples[drop], collapse = ", "))
  expr <- expr[, !drop, drop = FALSE]; group <- group[!drop]; samples <- colnames(expr)
}
# 只保留指定分组 (分组多于 2 个时用来挑出要对比的两组)
if (nzchar(OPT$group_keep)) {
  keepg <- trimws(strsplit(OPT$group_keep, ",", fixed = TRUE)[[1]])
  miss  <- setdiff(keepg, unique(group))
  if (length(miss)) stop("--group_keep 中的 '", paste(miss, collapse = ","),
                         "' 不在分组中。可用: ", paste(sort(unique(group)), collapse = ", "))
  dropg <- !(group %in% keepg)
  if (any(dropg)) {
    msg("  ! 按 --group_keep 剔除 ", sum(dropg), " 个样本")
    expr <- expr[, !dropg, drop = FALSE]; group <- group[!dropg]; samples <- colnames(expr)
  }
}
if (length(unique(group)) < 2L) stop("有效分组少于 2 个，无法做差异分析。")

# 把样本->分组映射落盘，便于人工核对 (这一步错了后面全白做)
write.csv(data.frame(sample = samples, group = group, check.names = FALSE),
          paste0(PREFIX, "_sample_group_map.csv"), row.names = FALSE)

# 分组名清洗成合法变量名 (makeContrasts 不接受空格/短横等)
lv  <- sort(unique(group))
lv2 <- make.names(lv)
if (!identical(lv, lv2)) {
  msg("  ! 分组名含非法字符，已转换: ", paste(sprintf("%s -> %s", lv, lv2), collapse = ", "))
  group <- lv2[match(group, lv)]
}
# 按分组排序，热图才好分块
ord <- order(group)
expr <- expr[, ord, drop = FALSE]; group <- group[ord]; samples <- colnames(expr)

msg("  分组概览:")
print(table(group))
for (g0 in sort(unique(group))) {
  ss <- samples[group == g0]
  cat("    ", g0, " (n=", length(ss), "): ", paste(head(ss, 4), collapse = ", "),
      if (length(ss) > 4L) ", ..." else "", "\n", sep = "")
}
cat("    完整映射 -> ", PREFIX, "_sample_group_map.csv\n", sep = "")

## --------------------------- 3. limma 差异分析 ------------------------------
msg("[3/5] limma 差异分析 ...")
design <- model.matrix(~0 + factor(group))
colnames(design) <- sort(unique(group))
rownames(design) <- samples

CONTRAST <- OPT$contrast
if (!nzchar(CONTRAST)) {
  if (ncol(design) != 2L) stop("分组数 > 2，必须用 --contrast= 指定对比，如 A-B")
  CONTRAST <- paste0(colnames(design)[2], "-", colnames(design)[1])
  msg("  --contrast 未指定，默认: ", CONTRAST)
}
bad <- setdiff(strsplit(CONTRAST, "-", fixed = TRUE)[[1]], colnames(design))
if (length(bad)) stop("对比式 '", CONTRAST, "' 中的 '", paste(bad, collapse = ","),
                      "' 不在分组中。可用分组: ", paste(colnames(design), collapse = ", "))

fit  <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, makeContrasts(contrasts = CONTRAST, levels = design))
fit2 <- eBayes(fit2)
nrDEG <- na.omit(topTable(fit2, coef = 1, n = Inf, sort.by = OPT$sort_by))
msg("  得到 ", nrow(nrDEG), " 个基因的结果")

## --------------------------- 4. 阈值与上下调 --------------------------------
PCOL <- if (OPT$p_type %in% colnames(nrDEG)) OPT$p_type else "P.Value"
DEG  <- nrDEG
DEG$neglogP <- -log10(DEG[[PCOL]])
DEG$change  <- factor(
  ifelse(DEG[[PCOL]] < PCUT & abs(DEG$logFC) > LOGCUT,
         ifelse(DEG$logFC > LOGCUT, "UP", "DOWN"), "NOT"),
  levels = c("UP", "DOWN", "NOT"))
n_up <- sum(DEG$change == "UP"); n_down <- sum(DEG$change == "DOWN")
cat(sprintf("  |logFC| > %.3g 且 %s < %.3g : 上调 %d, 下调 %d\n", LOGCUT, PCOL, PCUT, n_up, n_down))

write.csv(DEG, paste0(PREFIX, "_all.csv"))
write.csv(DEG[DEG$change != "NOT", ], paste0(PREFIX, "_DEG.csv"))

## ------------------------------ 5. 火山图 -----------------------------------
msg("[4/5] 绘制火山图 ...")
STYLE <- tolower(trimws(if (nzchar(OPT$fig_style)) OPT$fig_style else "classic"))
if (!STYLE %in% c("classic", "nature")) {
  msg("  ! 未知 --fig_style=", STYLE, "，回退为 classic")
  STYLE <- "classic"
}
# 火山图风格: --vol_style 未指定 (same) 时跟随 --fig_style
VSTYLE <- tolower(trimws(if (nzchar(OPT$vol_style)) OPT$vol_style else "same"))
if (VSTYLE == "same") VSTYLE <- STYLE
if (!VSTYLE %in% c("classic", "nature", "enhanced", "gradient", "rainbow", "tophits")) {
  msg("  ! 未知 --vol_style=", VSTYLE, "，回退为 ", STYLE)
  VSTYLE <- STYLE
}
if (VSTYLE == "enhanced" && !requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  msg("  ! --vol_style=enhanced 需 EnhancedVolcano 包 (BiocManager::install('EnhancedVolcano'))，未安装，回退为 ", STYLE)
  VSTYLE <- STYLE
}
# ENV-1: EnhancedVolcano 1.28.x 与 ggplot2 >= 4.0 不兼容 (上游问题):
# 数据点全部挤在 x≈0、y 轴被截断, 与真实分布完全不符; 不降级 ggplot2, 仅提示绕行
if (VSTYLE == "enhanced" && requireNamespace("EnhancedVolcano", quietly = TRUE) &&
    utils::packageVersion("ggplot2") >= "4.0") {
  msg("  ! EnhancedVolcano 与 ggplot2 >= 4.0 存在兼容性问题 (上游): 出图可能为空/坐标异常;")
  msg("    建议改用 classic / gradient / rainbow / tophits, 详见 BUG记录_20260926.md ENV-1")
}
msg("  绘图风格: fig=", STYLE, " / volcano=", VSTYLE)
TOPN_LAB <- suppressWarnings(as.integer(OPT$top_n)); if (is.na(TOPN_LAB) || TOPN_LAB < 1) TOPN_LAB <- 10L

n_not <- sum(DEG$change == "NOT")
if (VSTYLE == "nature") {
  COLS <- NATURE_VOLCANO_COLS; PT_SIZE <- 0.7; PT_ALPHA <- 0.55
} else {
  COLS <- c(UP = "#C31E1F", DOWN = "#1F6FC3", NOT = "#898989")
  PT_SIZE <- 1.75; PT_ALPHA <- 0.4
}
# 图例带样本/基因计数 (nature 约定)
LEG <- c(UP  = sprintf("Up (%d)", n_up),
         DOWN = sprintf("Down (%d)", n_down),
         NOT = sprintf("Not significant (%d)", n_not))

ttl <- sprintf("Cutoff: |logFC| > %s, %s < %s\nUP = %d, DOWN = %d",
               LOGCUT, PCOL, PCUT, n_up, n_down)
TITLE <- OPT$fig_title
if (!nzchar(TITLE) && VSTYLE != "nature") TITLE <- ttl   # nature 默认不加标题

# ========== 新风格分支: enhanced / gradient / rainbow / tophits ==========
# 这四款各自完整接管绘图，画完直接 save_figure 并跳过下方 classic/nature 通用代码。
# 图内文字一律 ASCII (pitfalls 7.8)。
NEWV <- VSTYLE %in% c("enhanced", "gradient", "rainbow", "tophits")

# 标注基因集合: 手动指定 (--label_genes) ∪ top N (--label_top 或 --top_n)
pick_labels <- function(n_top) {
  idx <- rownames(DEG) %in% trimws(strsplit(OPT$label_genes, ",")[[1]])
  sig <- which(DEG$change != "NOT")
  if (length(sig) && n_top > 0) {
    ord_p <- sig[order(DEG[[PCOL]][sig])]
    idx[ord_p[seq_len(min(n_top, length(ord_p)))]] <- TRUE
  }
  idx
}

if (NEWV) {
  # --vol_width* / --vol_height* 对四款新风格同样生效
  base_theme <- if (STYLE == "nature") theme_nature(8) else
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank())

  if (VSTYLE == "enhanced") {
    ## -- V1: EnhancedVolcano 四区分色款 (kevinblighe/EnhancedVolcano) --
    msg("  火山图风格: enhanced (EnhancedVolcano)")
    lab_idx <- pick_labels(if (as.integer(OPT$label_top) > 0) as.integer(OPT$label_top) else TOPN_LAB)
    p_sub <- data.frame(
      logFC = DEG$logFC, neglogP = DEG$neglogP,
      lab = ifelse(lab_idx, rownames(DEG), ""))
    EV_COLS <- if (STYLE == "nature") c("grey30", "#2166AC", "#B2182B", "#B2182B") else
      c("grey30", "#1F6FC3", "#C31E1F", "#C31E1F")
    g <- EnhancedVolcano::EnhancedVolcano(
      p_sub, lab = p_sub$lab, x = "logFC", y = "neglogP",
      pCutoff = PCUT, FCcutoff = LOGCUT,
      xlab = "log2 fold change", ylab = paste0("-log10 ", PCOL),
      title = TITLE, subtitle = NULL, caption = NULL,
      pointSize = if (STYLE == "nature") 0.9 else 2.0,
      labSize = if (STYLE == "nature") 2.2 else 3.2,
      col = EV_COLS, colAlpha = if (STYLE == "nature") 0.55 else 0.45,
      legendPosition = "none",
      drawConnectors = TRUE, widthConnectors = 0.3,
      boxedLabels = FALSE, axisLabSize = if (STYLE == "nature") 8 else 12,
      titleLabSize = if (STYLE == "nature") 9 else 13)
  } else if (VSTYLE == "gradient") {
    ## -- V2: 渐变气泡款 (仿 BioSenior/ggVolcano::gradual_volcano) --
    msg("  火山图风格: gradient (渐变气泡)")
    D <- DEG; D$sig <- -log10(D[[PCOL]])
    D$sig[D$sig > 20] <- 20  # 渐变映射截断，避免极端 p 值把色阶压扁
    sig_rel <- pmin(D$sig / max(D$sig, na.rm = TRUE), 1)
    # 点色: 上调 浅橙->深红, 下调 浅蓝->深蓝, 随显著性加深; 不显著灰
    up_pal  <- colorRampPalette(c("#F4A582", "#B2182B"))(50)
    dn_pal  <- colorRampPalette(c("#92C5DE", "#2166AC"))(50)
    D$col <- ifelse(D$change == "UP",   up_pal[pmax(1, pmin(50, as.integer(sig_rel * 49) + 1L))],
             ifelse(D$change == "DOWN", dn_pal[pmax(1, pmin(50, as.integer(sig_rel * 49) + 1L))],
                    "#BDBDBD"))
    # 点大小: 随显著性从 0.8 渐变到 3.2
    D$size <- ifelse(D$change == "NOT", 1.0, 0.8 + 2.4 * sig_rel)
    lab_idx <- pick_labels(if (as.integer(OPT$label_top) > 0) as.integer(OPT$label_top) else TOPN_LAB)
    g <- ggplot2::ggplot(D, ggplot2::aes(x = logFC, y = neglogP)) +
      ggplot2::geom_point(colour = D$col, size = D$size, alpha = 0.85) +
      ggplot2::geom_hline(yintercept = -log10(PCUT), lty = 4, col = "#555555", lwd = 0.4) +
      ggplot2::geom_vline(xintercept = c(-LOGCUT, LOGCUT), lty = 4, col = "#555555", lwd = 0.4) +
      ggplot2::labs(x = "log2 fold change", y = paste0("-log10 ", PCOL),
                    title = if (nzchar(OPT$fig_title)) TITLE else NULL) +
      base_theme
    if (any(lab_idx) && requireNamespace("ggrepel", quietly = TRUE)) {
      g <- g + ggrepel::geom_text_repel(
        data = D[lab_idx, , drop = FALSE], ggplot2::aes(label = rownames(D)[lab_idx]),
        size = if (STYLE == "nature") 1.9 else 2.8, max.overlaps = Inf,
        segment.size = 0.25, segment.color = "#777777", min.segment.length = 0)
    }
    msg("  ! gradient 风格: 颜色/大小随显著度渐变 (上/下调分色), 无离散图例")
  } else if (VSTYLE == "rainbow") {
    ## -- V5: 五段彩虹渐变款 (深蓝->青->黄->橙红->深红, 中文社区教程) --
    msg("  火山图风格: rainbow (五段彩虹渐变)")
    rainbow_pal <- c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25")
    g <- ggplot2::ggplot(DEG, ggplot2::aes(x = logFC, y = neglogP)) +
      ggplot2::geom_point(ggplot2::aes(size = neglogP, colour = neglogP), alpha = 0.85) +
      ggplot2::geom_hline(yintercept = -log10(PCUT), lty = 4, color = "#999999") +
      ggplot2::geom_vline(xintercept = c(0, 0), lty = 4, color = "#999999") +
      ggplot2::scale_color_gradientn(
        values = seq(0, 1, 0.2), colours = rainbow_pal,
        breaks = signif(seq(0, max(DEG$neglogP, na.rm = TRUE), length.out = 4), 1)) +
      ggplot2::scale_size_continuous(range = c(0.6, 3.2), guide = "none") +
      ggplot2::labs(x = "log2 fold change", y = paste0("-log10 ", PCOL),
                    colour = paste0("-log10 ", PCOL),
                    title = if (nzchar(OPT$fig_title)) TITLE else NULL) +
      base_theme
    lab_idx <- pick_labels(if (as.integer(OPT$label_top) > 0) as.integer(OPT$label_top) else TOPN_LAB)
    if (any(lab_idx) && requireNamespace("ggrepel", quietly = TRUE)) {
      g <- g + ggrepel::geom_text_repel(
        data = DEG[lab_idx, , drop = FALSE],
        ggplot2::aes(label = rownames(DEG)[lab_idx]),
        size = if (STYLE == "nature") 1.9 else 2.8, max.overlaps = Inf,
        segment.size = 0.25, segment.color = "#777777", min.segment.length = 0)
    }
  } else if (VSTYLE == "tophits") {
    ## -- V6: Top-Hits 距离排序款 (仿 VolcaNoseR, Sci Rep 2020) --
    msg("  火山图风格: tophits (Top-Hits 距离排序, N=", TOPN_LAB, ")")
    D <- DEG
    # Manhattan 距离 = |logFC| + |-log10P|; 只在显著点里挑 Top N
    D$dist <- abs(D$logFC) + D$neglogP
    sig_i <- which(D$change != "NOT")
    top_i <- sig_i[order(D$dist[sig_i], decreasing = TRUE)]
    top_i <- head(top_i, TOPN_LAB)
    # 手动指定的基因强制并入
    man_i <- which(rownames(D) %in% trimws(strsplit(OPT$label_genes, ",")[[1]]))
    top_i <- union(top_i, man_i)
    D$lab <- ifelse(seq_len(nrow(D)) %in% top_i, rownames(D), "")
    COLS_T <- if (STYLE == "nature") NATURE_VOLCANO_COLS else
      c(UP = "#99000d", DOWN = "#053061", NOT = "#BDBDBD")
    g <- ggplot2::ggplot(D, ggplot2::aes(x = logFC, y = neglogP, colour = change)) +
      ggplot2::geom_point(alpha = if (STYLE == "nature") 0.55 else 0.5,
                          size = if (STYLE == "nature") 0.9 else 1.8) +
      ggplot2::geom_hline(yintercept = -log10(PCUT), lty = 2, col = "#555555", lwd = 0.4) +
      ggplot2::geom_vline(xintercept = c(-LOGCUT, LOGCUT), lty = 2, col = "#555555", lwd = 0.4) +
      ggplot2::scale_colour_manual(values = COLS_T,
                                   labels = c(UP = sprintf("Up (%d)", n_up),
                                              DOWN = sprintf("Down (%d)", n_down),
                                              NOT = sprintf("NS (%d)", n_not)), drop = FALSE) +
      ggplot2::labs(x = "log2 fold change", y = paste0("-log10 ", PCOL),
                    title = if (nzchar(OPT$fig_title)) TITLE else NULL) +
      base_theme
    if (length(top_i) && requireNamespace("ggrepel", quietly = TRUE)) {
      g <- g + ggrepel::geom_text_repel(
        data = D[top_i, , drop = FALSE], ggplot2::aes(label = lab),
        size = if (STYLE == "nature") 1.9 else 2.8, max.overlaps = Inf,
        fontface = "bold", segment.size = 0.25, segment.color = "#555555",
        min.segment.length = 0)
    }
  }

  VOL_W <- as.numeric(OPT$vol_width);  VOL_H <- as.numeric(OPT$vol_height)
  if (nzchar(OPT$vol_width_mm))  VOL_W <- as.numeric(OPT$vol_width_mm) / 25.4
  if (nzchar(OPT$vol_height_mm)) VOL_H <- as.numeric(OPT$vol_height_mm) / 25.4
  save_figure(g, paste0(PREFIX, "_volcano"), VOL_W, VOL_H, STYLE)
  msg("[4/5] 火山图完成 (", VSTYLE, ")")
} else {
# ========== 通用分支: classic / nature (原有代码, 仅缩进保持可读) ==========

# 图例带样本/基因计数 (nature 约定)
LEG <- c(UP  = sprintf("Up (%d)", n_up),
         DOWN = sprintf("Down (%d)", n_down),
         NOT = sprintf("Not significant (%d)", n_not))

g <- ggplot(DEG, aes(x = logFC, y = neglogP, colour = change)) +
  geom_point(alpha = PT_ALPHA, size = PT_SIZE) +
  scale_colour_manual(values = COLS, labels = LEG, drop = FALSE) +
  geom_vline(xintercept = c(-LOGCUT, LOGCUT), lty = 4, col = "#555555", lwd = 0.4) +
  geom_hline(yintercept = -log10(PCUT), lty = 4, col = "#555555", lwd = 0.4) +
  labs(x = "log2 fold change", y = paste0("-log10 ", PCOL))
if (nzchar(TITLE)) g <- g + labs(title = TITLE)

# 标注基因 = 手动指定 (--label_genes) ∪ 自动 top N (--label_top, 按 p 值最小)
lab_idx <- rep(FALSE, nrow(DEG))
if (nzchar(OPT$label_genes)) {
  want <- trimws(strsplit(OPT$label_genes, ",")[[1]])
  lab_idx <- lab_idx | (rownames(DEG) %in% want)
}
top_lab <- suppressWarnings(as.integer(OPT$label_top))
if (!is.na(top_lab) && top_lab > 0) {
  sig <- which(DEG$change != "NOT")
  if (length(sig)) {
    ord_p <- sig[order(DEG[[PCOL]][sig])]
    lab_idx[ord_p[seq_len(min(top_lab, length(ord_p)))]] <- TRUE
  }
}

if (any(lab_idx)) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    msg("  ! 未安装 ggrepel，跳过基因标注 (install.packages('ggrepel'))")
  } else {
    g <- g + ggrepel::geom_text_repel(
      data = DEG[lab_idx, , drop = FALSE],
      aes(label = rownames(DEG)[lab_idx]),
      size = if (STYLE == "nature") 1.8 else 3,
      box.padding = unit(0.4, "lines"), point.padding = unit(0.5, "lines"),
      segment.size = 0.25, segment.color = "#777777",
      min.segment.length = 0, max.overlaps = Inf, show.legend = FALSE)
  }
}

if (STYLE == "nature") {
  g <- g + theme_nature() + theme(legend.position = "bottom")
} else {
  g <- g + theme_bw(base_size = 14) +
    theme(plot.title = element_text(size = 11, hjust = 0.5), legend.title = element_blank())
}

VOL_W <- as.numeric(OPT$vol_width);  VOL_H <- as.numeric(OPT$vol_height)
if (nzchar(OPT$vol_width_mm))  VOL_W <- as.numeric(OPT$vol_width_mm) / 25.4
if (nzchar(OPT$vol_height_mm)) VOL_H <- as.numeric(OPT$vol_height_mm) / 25.4
save_figure(g, paste0(PREFIX, "_volcano"), VOL_W, VOL_H, STYLE)
}  # end NEWV else-branch

## ------------------------------- 6. 热图 ------------------------------------
msg("[5/5] 绘制热图 ...")
TOPN <- min(TOPN, nrow(DEG))
top_genes <- rownames(DEG)[seq_len(TOPN)]
mat <- expr[top_genes, , drop = FALSE]
# pheatmap scale="row" 遇到 sd=0 或全 NA 的行会报错，先剔除
keep <- apply(mat, 1, function(v) { s <- sd(v, na.rm = TRUE); !is.na(s) && s > 0 })
if (!all(keep)) msg("  ! 热图剔除 ", sum(!keep), " 个零方差/全 NA 基因")
mat <- mat[keep, , drop = FALSE]

# 热图风格: --hm_style 未指定 (same) 时跟随 --fig_style
HSTYLE <- tolower(trimws(if (nzchar(OPT$hm_style)) OPT$hm_style else "same"))
if (HSTYLE == "same") HSTYLE <- STYLE
if (!HSTYLE %in% c("classic", "nature", "complex", "numbers", "tile")) {
  msg("  ! 未知 --hm_style=", HSTYLE, "，回退为 ", STYLE)
  HSTYLE <- STYLE
}
if (HSTYLE == "complex" && !requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  msg("  ! --hm_style=complex 需 ComplexHeatmap 包，未安装，回退为 ", STYLE)
  HSTYLE <- STYLE
}
# numbers 仅 <= 20 样本可用: 超限自动回退并提示
if (HSTYLE == "numbers" && ncol(mat) > 20) {
  msg("  ! --hm_style=numbers 仅适合 <=20 样本 (当前 ", ncol(mat),
      " 列, 数值会重叠不可读)，回退为 ", STYLE)
  HSTYLE <- STYLE
}
msg("  热图风格: ", HSTYLE)

ann   <- data.frame(group = factor(group, levels = sort(unique(group))))
rownames(ann) <- colnames(mat)
gaps  <- cumsum(rle(as.character(group))$lengths)
gaps  <- gaps[-length(gaps)]

# 风格相关参数 (classic/nature 数值型热图与原版一致)
if (STYLE == "nature") {
  pal      <- c("#0F4D92", "#B64342", "#42949E", "#E28E2C", "#9A4D8E", "#3775BA")  # 冷暖交替, 相邻组易分
  hm_ramp  <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100)  # 发散 蓝-白-红
  hm_border <- NA            # 去边框
  hm_fs     <- 6; hm_fs_row <- 6
} else {
  pal      <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#A6761D")
  hm_ramp  <- colorRampPalette(c("#1F6FC3", "white", "#C31E1F"))(100)
  hm_border <- "grey60"
  hm_fs     <- 10; hm_fs_row <- 9
}
gp_col <- setNames(rep(pal, length.out = nlevels(ann$group)), levels(ann$group))

hm_w <- as.numeric(OPT$hm_width); hm_h <- as.numeric(OPT$hm_height)
if (nzchar(OPT$hm_width_mm))  hm_w <- as.numeric(OPT$hm_width_mm) / 25.4
if (nzchar(OPT$hm_height_mm)) hm_h <- as.numeric(OPT$hm_height_mm) / 25.4

if (HSTYLE %in% c("complex", "numbers")) {
  # ===== 新风格 1: complex —— ComplexHeatmap 期刊款 =====
  # 顶部分组注释条 (与 classic/nature 的 pheatmap annotation_col 等价) +
  # 发散色带 + 行列自动聚类; 组间 gaps 在 ComplexHeatmap 里由聚类自然呈现。
  if (HSTYLE == "complex") {
    msg("  ComplexHeatmap 期刊款: 顶部分组注释条 + 自动聚类")
    suppressPackageStartupMessages(library(ComplexHeatmap))
    suppressPackageStartupMessages(library(circlize))
    mat_sc <- t(scale(t(mat)))                       # 行内 z-score (等价 pheatmap scale="row")
    mat_sc[mat_sc >  2.5] <-  2.5; mat_sc[mat_sc < -2.5] <- -2.5  # 截断到色带范围
    ramp_funs <- circlize::colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#F7F7F7", "#B2182B"))
    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      Group = ann$group,
      col = list(Group = gp_col),
      annotation_name_side = "left",
      annotation_name_gp = grid::gpar(fontsize = if (STYLE == "nature") 7 else 9))
    ht <- ComplexHeatmap::Heatmap(
      mat_sc, name = "Z-score",
      col = ramp_funs,
      top_annotation = ha_top,
      cluster_rows = TRUE, cluster_columns = FALSE,   # 列保持分组顺序 (与 pheatmap 版一致)
      column_split = ann$group,                        # 按组分块, 代替 pheatmap 的 gaps_col
      show_column_dend = FALSE,
      row_names_gp = grid::gpar(fontsize = if (STYLE == "nature") 6 else 9),
      column_names_gp = grid::gpar(fontsize = if (STYLE == "nature") 6 else 8),
      show_column_names = FALSE,
      row_dend_width = grid::unit(if (STYLE == "nature") 0.8 else 1.2, "cm"),
      heatmap_legend_param = list(
        title = "Z-score",
        legend_gp = grid::gpar(fontsize = if (STYLE == "nature") 6 else 8),
        title_gp = grid::gpar(fontsize = if (STYLE == "nature") 7 else 9,
                              fontface = "plain")))
    draw_hm <- function() ComplexHeatmap::draw(
      ht, heatmap_legend_side = "right",
      padding = unit(c(2, 12, 2, 2), "mm"))
  } else {
    # ===== 新风格 2: numbers —— pheatmap 数值标注款 =====
    # 与 classic/nature 同一套 pheatmap, 只是打开 display_numbers。
    msg("  pheatmap 数值标注款: display_numbers=TRUE")
    draw_hm <- function() pheatmap::pheatmap(mat,
           scale = "row",
           color = hm_ramp,
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_col = ann,
           annotation_colors = list(group = gp_col),
           gaps_col = gaps,
           border_color = hm_border,
           fontsize = hm_fs, fontsize_row = hm_fs_row,
           display_numbers = TRUE,          # 单元格内标注 z-score
           number_format = "%.1f",
           number_color = "grey30",
           fontsize_number = max(6, hm_fs - 2),
           show_colnames = FALSE, show_rownames = TRUE,
           annotation_legend = TRUE)
  }
} else if (HSTYLE == "tile") {
  # ===== 新风格 6: tile —— ggplot2 geom_tile 款 =====
  # 长数据 + geom_tile, 顶部分组色条 (geom_tile 实现), 无聚类; 大样本矩阵推荐。
  msg("  ggplot2 geom_tile 款: 分组色条 + 无聚类")
  hm_fs_tile <- if (STYLE == "nature") 6 else 8
  # 行内 z-score (与 pheatmap scale="row" 同口径)
  mat_sc <- t(scale(t(mat)))
  lim_q <- max(abs(quantile(mat_sc, probs = 0.02, na.rm = TRUE)), 1)  # 2% 分位截断防极端值压色阶
  mat_sc[mat_sc >  lim_q] <-  lim_q; mat_sc[mat_sc < -lim_q] <- -lim_q
  # 长格式
  long <- data.frame(
    gene = factor(rep(rownames(mat_sc), times = ncol(mat_sc)), levels = rev(rownames(mat_sc))),
    sample = factor(rep(colnames(mat_sc), each = nrow(mat_sc)), levels = colnames(mat_sc)),
    z = as.vector(mat_sc),
    group = rep(ann$group, each = nrow(mat_sc)),
    check.names = FALSE)
  long$yhm  <- as.integer(long$gene) + 1
  brk_y <- seq(1, length(rownames(mat_sc)), by = 1)
  g_hm <- ggplot2::ggplot(long, ggplot2::aes(x = sample, y = yhm, fill = z)) +
    ggplot2::geom_tile(colour = NA) +
    ggplot2::scale_fill_gradient2(
      low = if (STYLE == "nature") "#2166AC" else "#1F6FC3",
      mid = if (STYLE == "nature") "#F7F7F7" else "white",
      high = if (STYLE == "nature") "#B2182B" else "#C31E1F",
      limits = c(-lim_q, lim_q), name = "Z-score") +
    ggplot2::scale_y_continuous(
      breaks = brk_y, labels = levels(long$gene), expand = c(0.01, 0)) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL,
                  title = if (nzchar(OPT$fig_title)) TITLE else NULL) +
    ggplot2::theme_minimal(hm_fs_tile) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = hm_fs_tile, colour = "black"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = hm_fs_tile),
      legend.text = ggplot2::element_text(size = hm_fs_tile))
  # 顶部分组色条: 用 ggh4x 式双图拼接太重, 这里用 patchwork (已装) 上下拼
  grp_bar <- ggplot2::ggplot(
    data.frame(sample = unique(long$sample),
               grp = long$group[match(unique(long$sample), long$sample)]),
    ggplot2::aes(x = sample, y = 1, fill = grp)) +
    ggplot2::geom_tile(colour = NA) +
    ggplot2::scale_fill_manual(values = gp_col, name = "Group") +
    ggplot2::scale_x_discrete(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_void(hm_fs_tile) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = hm_fs_tile),
                   legend.title = ggplot2::element_text(size = hm_fs_tile),
                   legend.key.height = grid::unit(2, "mm"))
  g_hm <- patchwork::wrap_plots(grp_bar, g_hm, ncol = 1, heights = c(0.04, 1))
  draw_hm <- function() print(g_hm)
} else {
  # classic / nature: 原版 pheatmap
  draw_hm <- function() {
  pheatmap(mat,
           scale = "row",
           color = hm_ramp,
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_col = ann,
           annotation_colors = list(group = gp_col),
           gaps_col = gaps,
           border_color = hm_border,
           fontsize = hm_fs, fontsize_row = hm_fs_row,
           show_colnames = FALSE, show_rownames = TRUE,
           annotation_legend = TRUE)
  }
}

hm_pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf
hm_pdf_dev(paste0(PREFIX, "_heatmap.pdf"), width = hm_w, height = hm_h)
draw_hm()
invisible(dev.off())

if (isTRUE(as.logical(OPT$heatmap_png))) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(paste0(PREFIX, "_heatmap.png"), width = hm_w, height = hm_h,
                  units = "in", res = 150)
    draw_hm(); invisible(dev.off())
  } else {
    png(paste0(PREFIX, "_heatmap.png"), width = hm_w * 150, height = hm_h * 150,
        res = 150, type = "cairo")
    draw_hm(); invisible(dev.off())
  }
}

# nature 风格额外导出 SVG (可编辑文字) + 600dpi TIFF
if (STYLE == "nature") {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(paste0(PREFIX, "_heatmap.svg"), width = hm_w, height = hm_h)
    draw_hm(); invisible(dev.off())
  } else {
    msg("  ! 未安装 svglite，跳过热图 SVG 输出 (矢量图可用 .pdf 替代)")
  }
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(paste0(PREFIX, "_heatmap.tiff"), width = hm_w, height = hm_h,
                   units = "in", res = 600)
    draw_hm(); invisible(dev.off())
  } else {
    msg("  ! 未安装 ragg，跳过热图 TIFF 输出")
  }
}

save(expr, group, nrDEG, DEG, file = paste0(PREFIX, ".RData"))

# ---------------------------------------------------------------------------
# 交接清单：给下游 enrichment-analysis 技能用
#   <prefix>_for_enrichment.txt   （key=value 文本，'#' 开头是注释行）
# 富集侧只要给 --geo_dir=<本目录> --geo_prefix=<prefix>，就能自动对上
# ORA / GSEA 的输入文件与列名（gene/logFC/p/change），不必手写一堆参数。
# ---------------------------------------------------------------------------
{
  gtab  <- table(group)
  glv   <- names(gtab)
  kv    <- function(k, v) paste0(k, "=", paste(v, collapse = ""))
  absp  <- function(f) normalizePath(f, winslash = "/", mustWork = FALSE)
  cparts <- strsplit(CONTRAST, "-", fixed = TRUE)[[1]]
  man <- c(
    "# geo-microarray-analysis / 02_deg_plots.R 生成的下游交接清单",
    "# 供 enrichment-analysis 技能读取（--geo_dir 指向本文件所在目录即可）；也能人工核对",
    "# 注：本流程不解析芯片物种，species 一律按 human 交接；",
    "#     小鼠芯片请在富集侧显式加 --species=mouse",
    kv("skill", "geo-microarray-analysis"),
    kv("prefix", PREFIX),
    kv("dir", absp(".")),
    kv("created", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    kv("contrast", CONTRAST),
    if (length(cparts) == 2L) kv("contrast_num", cparts[1]) else "# contrast_num 略（对比式不是两段）",
    if (length(cparts) == 2L) kv("contrast_den", cparts[2]) else "# contrast_den 略",
    kv("groups", paste(glv, collapse = ";")),
    kv("n_samples", length(group)),
    kv("n_per_group", paste(as.integer(gtab), collapse = ";")),
    kv("n_up", n_up),
    kv("n_down", n_down),
    kv("deg_file", absp(paste0(PREFIX, "_DEG.csv"))),
    kv("all_file", absp(paste0(PREFIX, "_all.csv"))),
    kv("group_map_file", absp(paste0(PREFIX, "_sample_group_map.csv"))),
    kv("rdata_file", absp(paste0(PREFIX, ".RData"))),
    kv("gene_col", "rownames"),
    kv("logfc_col", "logFC"),
    kv("p_col", PCOL),
    kv("change_col", "change"),
    kv("deg_p_type", PCOL),
    kv("deg_p", format(PCUT, scientific = FALSE)),
    kv("deg_logfc", format(LOGCUT, scientific = FALSE)),
    kv("species", "human")
  )
  writeLines(man, paste0(PREFIX, "_for_enrichment.txt"), useBytes = TRUE)
  cat("交接清单    : ", PREFIX, "_for_enrichment.txt",
      "  (下游 enrichment-analysis 用 --geo_dir/--geo_prefix 直接读)\n", sep = "")
}

cat("\n================ 完成 ================\n")
cat("对比        : ", CONTRAST, "\n")
cat("差异基因    : 上调 ", n_up, " / 下调 ", n_down, "\n", sep = "")
cat("绘图风格    : ", STYLE, "\n", sep = "")
FIG_EXT <- if (STYLE == "nature") "_volcano.pdf/png/svg/tiff" else "_volcano.pdf/png"
HM_EXT  <- if (STYLE == "nature") "_heatmap.pdf/png/svg/tiff" else "_heatmap.pdf/png"
cat("输出文件    : ", paste0(PREFIX, "_sample_group_map.csv"), ", ",
    paste0(PREFIX, "_all.csv"), ", ", paste0(PREFIX, "_DEG.csv"), ", ",
    paste0(PREFIX, FIG_EXT), ", ", paste0(PREFIX, HM_EXT), ", ",
    paste0(PREFIX, ".RData"), "\n", sep = "")
