#!/usr/bin/env Rscript
## ===========================================================================
## 03_enrich_plot.R —— 富集结果「美化主图」
##
## 版式来自参考代码 kegg美化.R（gground + ggprism）：
##   左侧按 GO 本体 / KEGG 分组的圆角色块标签，色块右边是通路名与基因名（斜体），
##   再往左是基因数量圆点，右侧是 -log10(p.adjust) 的圆角柱，底部一条 x 轴线。
##
## 与参考代码的差异（都是为了稳，不是改设计）：
##   1. 分类顺序与配色**按名字绑定**：KEGG/MF/CC/BP 各自固定颜色，不因缺某一类而错位。
##      参考代码用裸向量 pal 配合 levels=rev(c('BP','CC','MF','KEGG'))，
##      一旦某个本体没有显著条目，颜色就会整体平移一格。
##   2. factor(Description, levels = Description) 改成 unique(Description)：
##      同一 Description 在 BP/CC/MF 里重复出现时，前者会在 R 4.x 直接报
##      "duplicated levels" 而崩掉。
##   3. 过滤后某类为空时不会整体崩掉，而是打印警告并退化为"该类的全部条目"。
##
## 只读 01_enrich_ora.R 产出的 <prefix>_{GO,KEGG,GMT}_<direction>.csv。
## ===========================================================================

t0 <- Sys.time()

.get_script_dir <- function() {
  a <- commandArgs(FALSE)
  i <- grep("^--file=", a)
  if (length(i)) return(dirname(sub("^--file=", "", a[i[1]])))
  j <- grep("\\.[Rr]$", a)
  if (length(j)) return(dirname(a[j[1]]))
  getwd()
}
.find_lib <- function() {
  cand <- character(0)
  d <- Sys.getenv("ENRICH_SKILL_SCRIPTS", "")
  if (nzchar(d)) cand <- c(cand, file.path(d, "lib_enrich_common.R"))
  cand <- c(cand, file.path(.get_script_dir(), "lib_enrich_common.R"),
            file.path(getwd(), "lib_enrich_common.R"))
  for (p in unique(cand)) if (file.exists(p)) return(p)
  stop("找不到 lib_enrich_common.R。已试：\n  ", paste(unique(cand), collapse = "\n  "),
       "\n请用 run_enrich.sh 启动本脚本。", call. = FALSE)
}
source(.find_lib())

USAGE <- '
03_enrich_plot.R —— 富集结果美化主图（gground + ggprism 版式）

输入：
  --in_dir=./enrich       01_enrich_ora.R 的 --outdir
  --prefix=               ORA 输出文件的前缀
  --direction=all         画哪一套：all / up / down / list / deg
  --db=GO,KEGG            要合并进同一张图的表：GO / KEGG / GMT
  --top=5                 挑几条（0 = 不限制，配合 --category 用）

挑法（三种模式）：
  --mode=even             默认（2026-09-25 起）。均衡配额：整张图固定 --even_total 条
                          （默认 20），在实际出现的分类（BP/CC/MF/KEGG，叠 GMT 时
                          还有各集合）间平均分配——四类齐全时各 5 条；某分类不足
                          配额时，富余名额按 BP→CC→MF→KEGG 顺序轮转匀给还有存货
                          的分类，总条数仍尽量凑满。此模式下 --top / --max_terms
                          都不生效。
  --mode=per_ontology     每个分类（BP/CC/MF/KEGG）各取 p.adjust 最小的前 --top 条
  --mode=global           全库按 p.adjust 排名取前 --top 条（"只看排名前 10"就用它）

条数上限：
  --even_total=20         even 模式的总条数（默认 20）。even 模式下总条数就由它决定。
  --max_terms=15          仅 per_ontology / global 模式生效的硬上限，默认 15。挑条 +
                          去重之后如果超过上限，按 p.adjust 从小到大截到前 15 条
                          （打印告警）。0 = 不限制（不推荐用于 SCI 发表）
  GSEA 结果不进本图：本脚本只画 ORA（01 步）的富集表；读表时一旦发现 GSEA 列
  （NES/setSize）直接报错退出。GSEA 结果单独用 02 步的 NES 条形图展示
  （<prefix>_GSEA_top_NES.pdf/png），不与 ORA 美化主图合并。

按功能方向出图（配合 04_enrich_summary.R 的归纳结果）：
  --func_csv=<04 产物>    04_enrich_summary.R 出的 <prefix>_func_pathways_<方向>.csv
  --category=<功能大类>   只画这个功能大类（不指定则保留全部大类）
  --use_representative=TRUE  只画去冗余后的代表通路（默认；04 已挑好）

挑条目与配色：
  --fig_filter=padj       过滤口径：padj | p | qvalue | sig | none（默认 padj）
  --fig_thr=0.05          过滤阈值。某类过滤后为空时会自动退化为"取该类全部条目"并告警
  --dedup_padj=TRUE       同一 p.adjust 只留 Count 最多的一条（参考代码的行为；
                          用了 --func_csv 时默认关闭，因为 04 已按基因重叠挑过代表）
  --pal_idx=3             内置配色 1/2/3（3 = 参考代码最终使用的那套）
  --pal=                 自定义 4 色，逗号分隔，顺序对应 KEGG,MF,CC,BP

版式：
  --width=12 --height=8   英寸
  --format=pdf,png        pdf / png / tiff
  --dpi=600               png / tiff 的分辨率
  --x_break=2             x 轴刻度间隔
  --bar_width=0.5         左侧"基因数圆点/分类色块"的参考宽度（数据单位）
  --bar_height=0.6        彩色柱在 y 方向的高度（参考代码写死 0.6）
  --gene_pos=below        基因名放彩色柱**下方**（默认）还是柱**内部**（inside = 参考代码原样）
  --gene_dy=             基因名的 y 偏移（行单位）；不给就按 --bar_height 自动算，
                         所以调柱高时基因名会跟着柱子一起动
  --cat_size=3.88         左侧分类标签（BP/CC/MF/KEGG）的字号
  --cat_label_chars=0     分类标签每行最多几个字符（0 = 自动按色块可用宽度算；
                         自动时 "KEGG" 通常会折成两行，色块会往左加宽到刚好包住）
  --label_genes=TRUE      是否在通路名下画基因名（斜体）
  --gene_size=3.5         基因名字号
  --auto_gene_size=TRUE   按最长基因名字符数自动缩字号，避免冲出画布
  --auto_width=FALSE      改成自动加宽画布、保持字号（与上一项二选一）
  --max_genes=0           每条通路最多画几个基因名（0 = 全部，忠于参考代码；
                          基因名很长时设 15~25 可避免文字冲出画布）
  --legend=right|bottom|none   图例位置（默认 right，与参考代码一致）
  --title=                图标题，默认无（参考代码无标题）
  --out=                  输出文件路径（不含扩展名），默认 <in_dir>/<prefix>_enrich_<direction>

示例：
  bash run_enrich.sh 03_enrich_plot.R --in_dir=./enrich --prefix=GSE62452_T --direction=all
'

opt <- parse_args()
if (any(c("help", "h") %in% names(opt)) || !length(opt)) {
  cat(USAGE); quit(status = 0)
}

in_dir <- oa(opt, "in_dir", NULL)
prefix <- oa(opt, "prefix", NULL)
if (is.null(in_dir) || is.null(prefix)) {
  stop("必须给 --in_dir=<ORA 输出目录> 与 --prefix=<前缀>。", call. = FALSE)
}
direction <- oa(opt, "direction", "all")
dbs       <- toupper(oa_list(opt, "db", c("GO", "KEGG")))
top_n     <- oa_num(opt, "top", 5)
max_terms <- oa_num(opt, "max_terms", 15)
fig_filter<- tolower(oa(opt, "fig_filter", "padj"))
fig_thr   <- oa_num(opt, "fig_thr", 0.05)
dedup     <- oa_lgl(opt, "dedup_padj", TRUE)
pal_idx   <- oa_num(opt, "pal_idx", 3)
pal_user  <- oa_list(opt, "pal", NULL)
width     <- oa_num(opt, "width", 12)
height    <- oa_num(opt, "height", 8)
fmts      <- tolower(oa_list(opt, "format", c("pdf", "png")))
dpi       <- oa_num(opt, "dpi", 600)
x_break   <- oa_num(opt, "x_break", 2)
bar_width <- oa_num(opt, "bar_width", 0.5)
bar_height<- oa_num(opt, "bar_height", 0.6)
cat_size  <- oa_num(opt, "cat_size", 3.88)
cat_chars <- oa_num(opt, "cat_label_chars", 0)
gene_pos  <- tolower(oa(opt, "gene_pos", "below"))
gene_dy0  <- oa_num(opt, "gene_dy", NA_real_)
lab_genes <- oa_lgl(opt, "label_genes", TRUE)
max_genes <- oa_num(opt, "max_genes", 0)
legend_pos <- tolower(oa(opt, "legend", "right"))
mode       <- tolower(oa(opt, "mode", "even"))
even_total <- oa_num(opt, "even_total", 20)
func_csv   <- oa(opt, "func_csv", NULL)
category   <- oa(opt, "category", NULL)
use_rep    <- oa_lgl(opt, "use_representative", TRUE)
title      <- oa(opt, "title", NULL)
out_base  <- oa(opt, "out", NULL)
if (is.null(out_base) || !nzchar(out_base)) {
  if (!is.null(category) && nzchar(category)) {
    ## 只替换文件名非法字符，保留中文（去掉 LC_ALL 后 R 读写中文路径正常，
    ## 这样产物名自带功能方向，比一串下划线可读）
    out_base <- file.path(in_dir, sprintf("%s_enrich_%s_%s", prefix, direction,
                                          gsub("[/\\\\:*?\"<>|[:space:]]+", "_", category)))
  } else if (mode == "global") {
    out_base <- file.path(in_dir, sprintf("%s_enrich_%s_top%s", prefix, direction,
                                          if (top_n > 0) as.integer(top_n) else "all"))
  } else {
    out_base <- file.path(in_dir, sprintf("%s_enrich_%s", prefix, direction))
  }
}
if (!mode %in% c("even", "per_ontology", "global")) {
  stop("--mode 只能是 even（均衡配额，默认）/ per_ontology（每个分类各取前 N）或 global（全库取前 N），收到：",
       mode, call. = FALSE)
}

## ---- GSEA 防混入检查（★ SCI 主图只画 ORA）---------------------------------
## GSEA 结果有 NES / setSize 列，ORA 表没有。混进来会把"排序富集"和"超代表"
## 两种口径画在同一张图里，审稿人一眼就能挑出来。GSEA 单独用 02 步的 NES 条形图。
gsea_cols <- c("NES", "setSize")
for (db in dbs) {
  f_try <- file.path(in_dir, sprintf("%s_%s_%s.csv", prefix, db, direction))
  if (!file.exists(f_try)) next
  hdr <- names(utils::read.csv(f_try, nrow = 1, check.names = FALSE,
                               stringsAsFactors = FALSE, quote = "\"", comment.char = ""))
  hit <- intersect(gsea_cols, hdr)
  if (length(hit)) {
    stop("检测到 GSEA 结果混入：", basename(f_try), " 含 GSEA 专有列（",
         paste(hit, collapse = ", "), "）。\n",
         "  美化主图只画 ORA 结果（01_enrich_ora.R 产物）。\n",
         "  GSEA 结果单独展示：用 02_enrich_gsea.R 的 NES 条形图 ",
         "(<prefix>_GSEA_top_NES.pdf/png)，或 --db= 去掉 GSEA 来源的表。",
         call. = FALSE)
  }
}
## 用了 04 步的功能归类表时，默认不要再用 p.adjust 撞车去重 ——
## 04 已经按基因重叠挑过代表通路了，再按 p 值去重会把不同的通路误删。
dedup_explicit <- "dedup_padj" %in% names(opt)
if (!dedup_explicit && !is.null(func_csv)) dedup <- FALSE

dir_create(dirname(out_base))
start_log(paste0(out_base, "_plot.log"))
on.exit(close_log(), add = TRUE)

print_banner("03_enrich_plot.R", "富集结果美化主图")
print_args(opt)
note("")
if (mode == "even") {
  if ("top" %in% names(opt)) {
    note("⚠ even 模式下 --top 不生效：总条数由 --even_total=", as.integer(even_total),
         " 均分到各分类。")
  }
  if ("max_terms" %in% names(opt)) {
    note("⚠ even 模式下 --max_terms 不生效：总条数由 --even_total=",
         as.integer(even_total), " 决定。")
  }
}

for (p in c("ggplot2", "gground", "ggprism", "dplyr")) {
  if (!requireNamespace(p, quietly = TRUE)) {
    stop("缺少 R 包：", p, "\n  gground 的安装方式：devtools::install_github(\"dxsbiocc/gground\")",
         call. = FALSE)
  }
}
suppressPackageStartupMessages({
  library(ggplot2); library(gground); library(ggprism); library(dplyr)
})

## ---- 读入各张表 -----------------------------------------------------------
read_one <- function(db) {
  f <- file.path(in_dir, sprintf("%s_%s_%s.csv", prefix, db, direction))
  if (!file.exists(f)) {
    alt <- file.path(in_dir, sprintf("%s_%s_%s_sig.csv", prefix, db, direction))
    if (file.exists(alt)) {
      note("  ", db, "：找不到 ", basename(f), "，改用 ", basename(alt))
      f <- alt
    } else {
      note("  ", db, "：没有这个文件，跳过（", basename(f), "）")
      return(NULL)
    }
  }
  d <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "\"", comment.char = "")
  if (!nrow(d)) return(NULL)
  if (!"ONTOLOGY" %in% names(d)) d$ONTOLOGY <- db
  ## KEGG/GMT 表里的 ONTOLOGY 可能为空；GMT 用 collection。
  ## ★ 一定要把 NA 兜住：后面 `dat[dat$ONTOLOGY == k, ]` 里 k 若是 NA，
  ##   R 会返回一行全 NA 的"假数据"，figures 里多一条空通路。
  if (db == "KEGG") d$ONTOLOGY <- "KEGG"
  if (db == "GMT" && "collection" %in% names(d)) {
    d$ONTOLOGY <- ifelse(is.na(d$collection) | !nzchar(d$collection), "GMT", d$collection)
  }
  bad <- is.na(d$ONTOLOGY) | !nzchar(as.character(d$ONTOLOGY))
  if (any(bad)) d$ONTOLOGY[bad] <- db
  if (!"geneName" %in% names(d)) d$geneName <- d$geneID
  d$geneName[is.na(d$geneName)] <- ""
  note(sprintf("  %-5s <- %s  (%d 条被检验条目)", db, basename(f), nrow(d)))
  d
}

note("读取富集表（方向 = ", direction, "）：")
parts <- lapply(dbs, read_one)
parts <- parts[!vapply(parts, is.null, TRUE)]
if (!length(parts)) {
  stop("没有读到任何富集表。请先跑 01_enrich_ora.R，或检查 --prefix / --direction / --db。",
       call. = FALSE)
}
dat <- do.call(rbind, lapply(parts, function(d) {
  d[, intersect(c("ONTOLOGY", "ID", "Description", "pvalue", "p.adjust", "qvalue",
                  "geneID", "geneName", "Count"), names(d)), drop = FALSE]
}))

## 补齐可能缺的列
for (cc in c("pvalue", "p.adjust", "qvalue", "Count", "geneName", "geneID", "ID")) {
  if (!cc %in% names(dat)) dat[[cc]] <- NA
}
dat$p.adjust <- suppressWarnings(as.numeric(dat$p.adjust))
dat$pvalue   <- suppressWarnings(as.numeric(dat$pvalue))
dat$qvalue   <- suppressWarnings(as.numeric(dat$qvalue))
dat$Count    <- suppressWarnings(as.numeric(dat$Count))
dat$Description <- as.character(dat$Description)
dat$ONTOLOGY <- as.character(dat$ONTOLOGY)
dat <- dat[!is.na(dat$p.adjust) & !is.na(dat$Description), , drop = FALSE]
if (!nrow(dat)) stop("读到的表里没有可用的 p.adjust / Description。", call. = FALSE)

## ---- 过滤 -----------------------------------------------------------------
## --fig_filter 用的是"人话"（padj / p / qvalue），要映射回表里的列名。
filter_col <- switch(fig_filter,
                     padj = "p.adjust", padjval = "p.adjust", p.adjust = "p.adjust",
                     p = "pvalue", pval = "pvalue", pvalue = "pvalue",
                     q = "qvalue", qvalue = "qvalue",
                     fdr = "p.adjust",
                     NA_character_)
if (fig_filter != "none" && fig_filter != "sig" && is.na(filter_col)) {
  stop("--fig_filter 只能是 padj / p / qvalue / sig / none，收到：", fig_filter,
       call. = FALSE)
}
if (!is.na(filter_col) && !filter_col %in% names(dat)) {
  stop("表里没有列 ", filter_col, "，无法按 --fig_filter=", fig_filter, " 过滤。",
       call. = FALSE)
}
note("")
note(sprintf("过滤口径：%s%s <= %s",
             fig_filter,
             if (!is.na(filter_col)) paste0("（列 ", filter_col, "）") else "",
             fig_thr))
if (fig_filter == "none") {
  note("  不筛选，全部条目参与挑选")
} else if (fig_filter == "sig") {
  if (!"sig" %in% names(dat)) {
    note("  表里没有 sig 列，退化为 p.adjust <= ", fig_thr)
    dat <- dat[dat$p.adjust <= fig_thr, , drop = FALSE]
  } else {
    dat <- dat[as.logical(dat$sig), , drop = FALSE]
  }
} else {
  v <- dat[[filter_col]]
  dat <- dat[!is.na(v) & v <= fig_thr, , drop = FALSE]
}

## ---- 可选的"按功能方向出图"（配合 04_enrich_summary.R）---------------------
if (!is.null(func_csv) && nzchar(func_csv)) {
  if (!file.exists(func_csv)) stop("找不到 --func_csv 指定的文件：", func_csv, call. = FALSE)
  fc <- utils::read.csv(func_csv, check.names = FALSE, stringsAsFactors = FALSE,
                        quote = "\"", comment.char = "")
  if (!"ID" %in% names(fc)) stop("--func_csv 的表里没有 ID 列：", func_csv, call. = FALSE)
  if (!"func_cat" %in% names(fc)) stop("--func_csv 的表里没有 func_cat 列（应由 04_enrich_summary.R 产出）：",
                                       func_csv, call. = FALSE)
  before <- nrow(dat)
  dat$func_cat <- fc$func_cat[match(dat$ID, fc$ID)]
  dat$func_sub <- fc$func_sub[match(dat$ID, fc$ID)]
  dat$representative <- if ("representative" %in% names(fc)) {
    as.logical(fc$representative[match(dat$ID, fc$ID)])
  } else TRUE
  note("")
  note("已按功能归类表筛选：", basename(func_csv))
  if (!is.null(category) && nzchar(category)) {
    dat <- dat[!is.na(dat$func_cat) & dat$func_cat == category, , drop = FALSE]
    note("  只保留功能大类「", category, "」：", nrow(dat), " / ", before, " 条")
  } else {
    note("  （未指定 --category，全部功能大类都保留）")
  }
  if (use_rep) {
    n0 <- nrow(dat)
    dat <- dat[!is.na(dat$representative) & dat$representative, , drop = FALSE]
    note("  只保留代表通路（去冗余后）：", nrow(dat), " / ", n0, " 条")
  }
  if (!nrow(dat)) {
    avail <- sort(unique(fc$func_cat))
    stop("按功能大类筛选后一条通路都不剩。\n",
         if (!is.null(category) && nzchar(category))
           paste0("该方向下没有「", category, "」的通路。可用的功能大类：\n  ",
                  paste(avail, collapse = "\n  ")) else "",
         call. = FALSE)
  }
}

## 逐分类检查：某类被筛空时退化为"该类全部条目"（并明确告警，绝不静默）
cats_present <- unique(dat$ONTOLOGY)
cats_present <- cats_present[!is.na(cats_present) &
                               nzchar(as.character(cats_present))]
if (!length(cats_present) && nrow(dat)) {
  note("⚠ 分类列（ONTOLOGY）全是空值，退回按 --db 给的库名分类。")
  dat$ONTOLOGY <- rep(dbs, length.out = nrow(dat))
  cats_present <- unique(dat$ONTOLOGY)
}
if (!length(cats_present)) {
  note("⚠ ", fig_filter, " <= ", fig_thr, " 之后一条都不剩。")
  note("  -> 自动退化为不筛选，用全部条目出图。如果这不是你要的，",
       "请放宽 --fig_thr 或改 --fig_filter=none。")
  dat <- do.call(rbind, lapply(parts, function(d) {
    d[, intersect(c("ONTOLOGY", "ID", "Description", "pvalue", "p.adjust", "qvalue",
                    "geneID", "geneName", "Count"), names(d)), drop = FALSE]
  }))
  for (cc in c("pvalue", "p.adjust", "qvalue", "Count", "geneName", "geneID", "ID")) {
    if (!cc %in% names(dat)) dat[[cc]] <- NA
  }
  dat$p.adjust <- suppressWarnings(as.numeric(dat$p.adjust))
  dat$Count <- suppressWarnings(as.numeric(dat$Count))
  dat$ONTOLOGY <- as.character(dat$ONTOLOGY)
  dat <- dat[!is.na(dat$p.adjust) & !is.na(dat$Description), , drop = FALSE]
  cats_present <- unique(dat$ONTOLOGY)
  cats_present <- cats_present[!is.na(cats_present) &
                                 nzchar(as.character(cats_present))]
  note("  现在共有分类：", paste(cats_present, collapse = ", "))
} else {
  all_cats <- unique(unlist(lapply(parts, function(d) as.character(d$ONTOLOGY))))
  dropped <- setdiff(all_cats, cats_present)
  if (length(dropped)) {
    note("  提示：", paste(dropped, collapse = ", "),
         " 在全表里有条目，但没有一条通过过滤 → 这张图里不会出现这些分类。")
  }
}
if (!nrow(dat)) stop("过滤后没有任何条目，无法出图。", call. = FALSE)

## ---- 挑条目 -----------------------------------------------------------------
## even（默认，2026-09-25 起）：均衡配额。总条数固定 --even_total（默认 20），
##   在实际出现的分类间平均分（BP/CC/MF/KEGG 四类齐全时各 5 条）；某分类库存
##   不足配额时，富余名额按 BP→CC→MF→KEGG 顺序轮转匀给还有存货的分类。
## per_ontology：每个分类各取 p.adjust 最小的前 N 条（参考代码行为）；
## global：全库按 p.adjust 排名取前 N 条。
## --top=0 表示不限制条数（配合 --category 用，仅 per_ontology/global 生效）。
cat_canon <- c(intersect(c("BP", "CC", "MF", "KEGG"), cats_present),
               sort(setdiff(cats_present, c("BP", "CC", "MF", "KEGG"))))
take <- function(d) if (top_n > 0) head(d, top_n) else d
if (mode == "even") {
  ## 去重先做：配额要在去重后的池子上算，否则后面去重会把条数打下去
  pool <- dat[order(dat$p.adjust, -dat$Count, dat$Description), , drop = FALSE]
  if (dedup) {
    b0 <- nrow(pool)
    pool <- pool[!duplicated(pool$p.adjust), , drop = FALSE]
    if (nrow(pool) < b0) {
      note(sprintf("  去重：p.adjust 相同的条目从 %d 条并到 %d 条（每档只留 Count 最多的）",
                   b0, nrow(pool)))
    }
  }
  n_cat      <- length(cat_canon)
  quota_base <- floor(even_total / n_cat)
  rem        <- even_total %% n_cat
  quota      <- setNames(rep(quota_base, n_cat), cat_canon)
  if (rem > 0) quota[cat_canon[seq_len(rem)]] <- quota[cat_canon[seq_len(rem)]] + 1
  stock  <- vapply(cat_canon, function(k) sum(pool$ONTOLOGY == k), 0L)
  take_i <- pmin(quota, stock)
  leftover <- even_total - sum(take_i)
  short <- cat_canon[stock < quota]
  note(sprintf("挑法：均衡配额（--mode=even）：总条数固定 %d，%d 个分类各配 %d 条%s",
               as.integer(even_total), n_cat, as.integer(quota_base),
               if (rem > 0) sprintf("（前 %d 个分类多 1 条）", as.integer(rem)) else ""))
  note("  配额：", paste(sprintf("%s %d", cat_canon, quota), collapse = " / "))
  note("  库存：", paste(sprintf("%s %d", cat_canon, stock), collapse = " / "))
  if (length(short)) {
    note("  ⚠ ", paste(sprintf("%s 只有 %d 条（配额 %d）", short, stock[short], quota[short]),
                      collapse = "；"),
         " → 富余名额按 BP→CC→MF→KEGG 顺序轮转匀给还有存货的分类")
  }
  guard <- 0L
  while (leftover > 0 && guard < 10000L) {
    guard <- guard + 1L
    progressed <- FALSE
    for (k in cat_canon) {
      if (leftover <= 0) break
      if (take_i[[k]] < stock[[k]]) {
        take_i[[k]] <- take_i[[k]] + 1
        leftover <- leftover - 1
        progressed <- TRUE
      }
    }
    if (!progressed) break     # 所有分类库存都见底了
  }
  sel <- do.call(rbind, lapply(cat_canon, function(k) {
    d <- pool[pool$ONTOLOGY == k, , drop = FALSE]
    head(d, take_i[[k]])
  }))
  if (sum(take_i) < even_total) {
    note("⚠ 去重后全部库存只有 ", sum(take_i), " 条，不足 --even_total=",
         as.integer(even_total), "，有多少画多少。")
  }
  note("  实际：", paste(sprintf("%s %d", cat_canon, take_i), collapse = " / "),
       "（合计 ", sum(take_i), " 条）")
} else if (mode == "global") {
  sel <- dat[order(dat$p.adjust, -dat$Count), , drop = FALSE]
  sel <- take(sel)
  note(sprintf("挑法：全库按 p.adjust 排名取前 %s 条",
               if (top_n > 0) as.integer(top_n) else "全部"))
} else {
  sel <- do.call(rbind, lapply(cats_present, function(k) {
    d <- dat[dat$ONTOLOGY == k, , drop = FALSE]
    d <- d[order(d$p.adjust, -d$Count), , drop = FALSE]
    take(d)
  }))
  note(sprintf("挑法：每个分类各取 p.adjust 最小的前 %s 条",
               if (top_n > 0) as.integer(top_n) else "全部"))
}
sel <- sel[!is.na(sel$p.adjust) & !is.na(sel$Description), , drop = FALSE]
if (!nrow(sel)) stop("挑完一条都不剩（检查 --top / --db / 过滤口径）。", call. = FALSE)
if (mode != "even" && dedup) {
  before <- nrow(sel)
  sel <- sel[order(sel$p.adjust, -sel$Count, sel$Description), , drop = FALSE]
  sel <- sel[!duplicated(sel$p.adjust), , drop = FALSE]
  if (nrow(sel) < before) {
    note(sprintf("  去重：p.adjust 相同的条目从 %d 条并到 %d 条（每档只留 Count 最多的）",
                 before, nrow(sel)))
  }
}

## ---- SCI 发表硬上限：整张图最多 --max_terms 条 ------------------------------
## 仅 per_ontology / global 模式生效。even 模式的总条数由 --even_total（默认 20）
## 决定，不走这个截断。
## 挑条（per_ontology/global）+ 去重之后再统一截一次，保证无论哪种挑法、
## 哪个功能方向，最终图上的通路数都不会超过上限（默认 15）。
if (mode != "even" && max_terms > 0 && nrow(sel) > max_terms) {
  note(sprintf("⚠ 条目数 %d 超过 SCI 发表上限 --max_terms=%d，按 p.adjust 从小到大截断",
               nrow(sel), as.integer(max_terms)))
  sel <- sel[order(sel$p.adjust, -sel$Count), , drop = FALSE]
  cut_ids <- sel$Description[(max_terms + 1L):nrow(sel)]
  sel <- head(sel, max_terms)
  note(sprintf("  截掉 %d 条（p.adjust 最大的一部分），例如：%s",
               length(cut_ids),
               paste(head(as.character(cut_ids), 3), collapse = " | ")))
}

## 分类显示顺序：KEGG 在下、BP 在上（与参考代码 levels = rev(c('BP','CC','MF','KEGG')) 一致）
canon <- c("BP", "CC", "MF", "KEGG")
extra <- setdiff(unique(sel$ONTOLOGY), canon)
ord_lv <- rev(c(intersect(canon, unique(sel$ONTOLOGY)), sort(extra)))
sel$ONTOLOGY <- factor(sel$ONTOLOGY, levels = ord_lv)
sel <- sel[order(sel$ONTOLOGY, sel$p.adjust, -sel$Count), , drop = FALSE]
sel$Description <- factor(sel$Description, levels = unique(sel$Description))
sel <- tibble::rowid_to_column(sel, "index")

note("")
note(sprintf("入选 %d 条（%s）：", nrow(sel),
             if (mode == "even") paste0("均衡配额，总 ", as.integer(even_total), " 条")
             else if (top_n > 0) paste0("每类最多 ", top_n, " 条") else "不限条数"))
tab <- as.data.frame(table(sel$ONTOLOGY))
print(tab, row.names = FALSE)
if (dedup) {
  note(sprintf("  （已按 p.adjust 撞车去重；想保留全部请加 --dedup_padj=FALSE，入选中现有 %d 条）",
               nrow(sel)))
}

## 落一份入选清单，便于核对
sel_out <- as.data.frame(sel)
cols_out <- c("index", "ONTOLOGY", "func_cat", "func_sub", "ID", "Description",
              "p.adjust", "Count", "geneName")
write_csv_out(sel_out[, intersect(cols_out, names(sel_out)), drop = FALSE],
              paste0(out_base, "_selected_terms.csv"))

## 基因名太长时会冲出画布：可选截断（默认 0 = 不截，与参考代码一致）
if (max_genes > 0) {
  sel$geneName <- vapply(strsplit(as.character(sel$geneName), "/", fixed = TRUE),
                         function(g) {
                           g <- g[nzchar(g)]
                           if (!length(g)) return("")
                           if (length(g) <= max_genes) return(paste(g, collapse = "/"))
                           paste0(paste(g[seq_len(max_genes)], collapse = "/"),
                                  " ... (+", length(g) - max_genes, ")")
                         }, "", USE.NAMES = FALSE)
  note("基因名已截断到每条通路最多 ", max_genes, " 个")
}

## ---- 版式参数（对齐参考代码）----------------------------------------------
w <- bar_width                       # 参考代码里的 width <- 0.5
xaxis_max <- max(-log10(sel$p.adjust)) + 1

## 高度自适应：条目多时加高、条目少时收矮，否则文字会挤在一起、柱子又会胖得难看。
## ★ 必须在算文字位置之前定下来 —— 基因名要按"一行有多高"来定位。
## 显式给了 --height 就完全按用户给的来（不做自适应）。
h <- if ("height" %in% names(opt)) height else max(4.6, 0.34 * nrow(sel) + 2.2)
if (!("height" %in% names(opt)) && abs(h - height) > 0.05) {
  note(sprintf("画布高度自适应：%d 条 -> %.2f 英寸（默认 %s；想固定就显式写 --height=）",
               nrow(sel), h, format(height)))
}

## ---- 基因名字号：长了就缩（或加宽画布）-------------------------------------
## 经验换算：geom_text 的 size 单位是 mm，实际 pt = size * 2.845；
## 基因名以大写字母、数字、斜杠为主，平均字宽约 0.4 em → 0.0158 * size 英寸/字符。
## （按 0.5 em 估会高估 25%，把字号压得过小。）
gene_size  <- oa_num(opt, "gene_size", 3.5)
auto_gene  <- oa_lgl(opt, "auto_gene_size", TRUE)
auto_width <- oa_lgl(opt, "auto_width", FALSE)
size_use   <- gene_size
CHARS_IN   <- 0.0158          # 英寸 / 字符（size = 1 时）
PANEL_FRAC <- 0.78            # 面板可用横向空间占画布宽度的比例
                              # （右侧要留图例 Category/Count 的位置）
if (auto_gene && lab_genes) {
  chars_max <- suppressWarnings(max(nchar(as.character(sel$geneName)), 1L, na.rm = TRUE))
  if (!is.finite(chars_max) || is.na(chars_max)) chars_max <- 1L
  needed_in <- chars_max * CHARS_IN * gene_size
  avail_in  <- width * PANEL_FRAC
  if (needed_in > avail_in) {
    if (auto_width) {
      width <- needed_in / PANEL_FRAC
      note(sprintf("基因名最长 %d 字符，画布自动加宽到 %.2f 英寸（字号保持 %.1f）",
                   chars_max, width, gene_size))
    } else {
      size_use <- max(2.0, gene_size * avail_in / needed_in)
      note(sprintf("基因名最长 %d 字符，字体自动缩到 %.2f（原 %.1f）以装进 %.0f 英寸画布",
                   chars_max, size_use, gene_size, width))
      if (size_use <= 2.01) {
        note("  已到下限 2.0，仍可能压线。建议改用 --max_genes=20（截断并标注 '+n'）")
        note("  或 --auto_width=TRUE（加宽画布、保持字号）。")
      }
    }
  }
}

## ---- 分类标签：折行，并让彩色框宽到刚好包住它 --------------------------------
## 参考代码把色块写死成 0.5 个数据单位宽（x 从 -1.5 到 -1），x 轴一长就装不下
## "KEGG"，会被面板边界裁成 "EGG"。这里按"这块宽能放几个字符"把标签折行
## （KEGG -> 两行），并把色块**往左**加宽到刚好包住文字（右边界不动，
## 所以左侧基因数圆点和 x=0 的位置都不受影响）。
CHAR_IN_CAT <- 0.0245         # 全大写分类标签比基因名宽，英寸/字符（size = 1 时）
wrap_label <- function(s, n) {
  s <- as.character(s)
  if (is.na(s) || !nzchar(s) || n < 1L) return(s)
  if (nchar(s) <= n) return(s)
  ch <- strsplit(s, "")[[1]]
  paste(vapply(split(ch, ceiling(seq_along(ch) / n)), paste, "", collapse = ""),
        collapse = "\n")
}
cats <- levels(sel$ONTOLOGY)
char_in_cat <- CHAR_IN_CAT * cat_size
box_xmax <- -2 * w
box_units <- w
lab_lines <- stats::setNames(cats, cats)
fit_chars <- 2L
for (it in 1:3) {              # 迭代收敛：块宽 <- 需要的宽度 <- 块宽决定的数据单位换算
  unit_in   <- (width * PANEL_FRAC) / (xaxis_max - (box_xmax - box_units))
  fit_chars <- max(2L, floor((w * unit_in) / char_in_cat))
  if (cat_chars > 0) fit_chars <- max(1L, as.integer(cat_chars))
  lab_lines <- stats::setNames(vapply(cats, function(k) wrap_label(k, fit_chars), ""),
                               cats)
  widest <- max(vapply(strsplit(unname(lab_lines), "\n", fixed = TRUE),
                       function(v) max(nchar(v)), 1L))
  box_units <- max(w, widest * char_in_cat * 1.15 / unit_in)
}
box_xmin <- box_xmax - box_units
unit_in  <- (width * PANEL_FRAC) / (xaxis_max - box_xmin)
row_in   <- h / nrow(sel)
note("")
note(sprintf("分类标签：每行最多 %d 字符（%s）→ 色块宽 %.2f 个数据单位（参考代码是 %.2f）",
             fit_chars,
             paste(sprintf("%s => %s", cats, gsub("\n", "/", unname(lab_lines))),
                   collapse = "  "),
             box_units, w))

rect_data <- sel %>%
  group_by(.data$ONTOLOGY) %>%
  summarise(n = dplyr::n(), .groups = "drop") %>%
  mutate(xmin = box_xmin, xmax = box_xmax,
         ymax = cumsum(n),
         ymin = dplyr::lag(ymax, default = 0) + 0.6,
         ymax = ymax + 0.4)
rect_data$lab <- unname(lab_lines[as.character(rect_data$ONTOLOGY)])

## ---- 基因名的位置：贴在彩色柱"下面"，跟着柱高走 -----------------------------
## 参考代码用 vjust = 2.6 把基因名压在柱子里（位置只跟字号有关，和柱高无关）。
## 这里改成显式给 y：柱下沿(bar_h/2) + 半个字高 + 一点空隙。
## 于是改 --bar_height 时文字会跟着柱子一起动。
gene_dy <- gene_dy0
if (is.na(gene_dy)) {
  line_in      <- size_use * 2.845 / 72 * 1.28     # 一行文字的实际高度（含行距），英寸
  text_h_units <- line_in / row_in                 # 换算成 y 方向的单位（1 行 = 1 单位）
  gene_dy      <- bar_height / 2 + text_h_units / 2 + 0.05
}
note(sprintf("基因名位置：%s（柱高 %.2f，y 偏移 %.3f 行 = %.3f 英寸）",
             if (gene_pos == "below") "彩色柱下方" else "彩色柱内部（参考代码原样）",
             bar_height, gene_dy, gene_dy * row_in))

pal <- build_pal(levels(sel$ONTOLOGY), pal_user, pal_idx)
note("")
note("配色（按分类名绑定）：")
for (k in names(pal)) note(sprintf("  %-6s %s", k, pal[[k]]))

## ---- 作图 -----------------------------------------------------------------
p <- ggplot(sel, aes(x = -log10(.data$p.adjust), y = .data$index, fill = .data$ONTOLOGY)) +
  geom_round_col(aes(y = .data$Description), width = bar_height, alpha = 0.8) +
  geom_text(aes(x = 0.05, label = .data$Description), hjust = 0, size = 5)
if (lab_genes) {
  if (gene_pos == "inside") {
    ## 参考代码原样：vjust 固定倍数（位置只跟字号走）
    p <- p + geom_text(aes(x = 0.1, label = .data$geneName, colour = .data$ONTOLOGY),
                       hjust = 0, vjust = 2.6, size = size_use, fontface = "italic",
                       show.legend = FALSE)
  } else {
    ## 默认：贴在彩色柱下沿（跟随 --bar_height）
    p <- p + geom_text(aes(x = 0.1, y = .data$index - gene_dy, label = .data$geneName,
                           colour = .data$ONTOLOGY),
                       hjust = 0, size = size_use, fontface = "italic",
                       show.legend = FALSE)
  }
}
p <- p +
  geom_point(aes(x = -w, size = .data$Count), shape = 21) +
  geom_text(aes(x = -w, label = .data$Count)) +
  ## Count 图例最多 4 档 —— 档位太多时图例竖向排不下（条目少、画布矮的图会被裁掉）
  scale_size_continuous(
    name = "Count", range = c(5, 12),
    breaks = { cc <- suppressWarnings(range(sel$Count, na.rm = TRUE))
               if (all(is.finite(cc)) && diff(cc) > 0)
                 unique(round(seq(cc[1], cc[2], length.out = 4))) else NULL }) +
  geom_round_rect(aes(xmin = .data$xmin, xmax = .data$xmax,
                      ymin = .data$ymin, ymax = .data$ymax, fill = .data$ONTOLOGY),
                  data = rect_data, radius = unit(2, "mm"), inherit.aes = FALSE) +
  ## 分类标签用我们折好行的 lab 列（r$lab 里已经带了 "\n"）
  geom_text(aes(x = (.data$xmin + .data$xmax) / 2, y = (.data$ymin + .data$ymax) / 2,
                label = .data$lab),
            data = rect_data, inherit.aes = FALSE, size = cat_size, lineheight = 0.95) +
  ## 底部 x 轴线。用 annotate 而不是 geom_segment：
  ## geom_segment 会把标量 aes 按"数据行数"广播，抛
  ## "All aesthetics have length 1, but the data has N rows" 警告。
  annotate("segment", x = 0, y = 0, xend = xaxis_max, yend = 0, linewidth = 1.5) +
  labs(y = NULL, title = title) +
  scale_fill_manual(name = "Category", values = pal) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = seq(0, xaxis_max, x_break),
                     expand = expansion(c(0, 0))) +
  theme_prism() +
  theme(axis.text.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(),
        legend.position = legend_pos,
        plot.title = element_text(size = 12)) +
  ## clip = "off"：允许文字越过面板边界绘制（分类标签、基因名都可能贴到边上）。
  coord_cartesian(clip = "off")

for (ext in fmts) {
  f <- paste0(out_base, ".", ext)
  if (ext == "pdf") {
    grDevices::pdf(f, width = width, height = h)
    print(p); grDevices::dev.off()
  } else {
    ggplot2::ggsave(f, p, width = width, height = h, dpi = dpi)
  }
  note("已出图：", f)
}

note("")
note("共 ", nrow(sel), " 条通路；长宽 = ", width, " x ", round(h, 2), " 英寸")
note("挑法 = ", mode,
     if (mode == "even") paste0("（总 ", as.integer(even_total), " 条均衡配额）")
     else if (top_n > 0) paste0("（前 ", top_n, " 条）") else "（不限条数）",
     if (!is.null(category) && nzchar(category)) paste0("；功能大类 = ", category) else "")
note("耗时：", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " 秒")
