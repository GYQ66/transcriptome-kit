#!/usr/bin/env Rscript
## ===========================================================================
## 02_enrich_gsea.R —— GSEA（排序基因列表富集）
##
## 与 ORA 的根本区别（这一点不能混）：
##   ORA  （01_enrich_ora.R）：输入是**筛过的显著基因列表**，问"这些基因里哪些通路被过度代表"
##   GSEA （本脚本）        ：输入是**全部基因**，按 logFC 排序、**不做任何筛选**，
##                            问"某条通路的基因是否系统性地堆在排序列表的两端"
##   所以本脚本**不接受** --ora_p / --ora_logfc 之类的筛选参数。要筛就不叫 GSEA 了。
##
## 支持的基因集：GO(BP/CC/MF 分别跑再合并)、KEGG、MSigDB GMT(H/C1..C9 或具体文件)
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
02_enrich_gsea.R —— GSEA（全部基因按 logFC 排序后做富集，不做任何筛选）

输入：
  --geo_dir=<目录>        ★ 直接从 geo-microarray-analysis 的差异分析产物接手。
  --geo_prefix=<前缀>     geo 02_deg_plots.R 的 --prefix（不给就自动找目录里唯一的一个）。
                          GSEA 吃 <前缀>_all.csv（全部基因），列名自动对上。
  --geo_manifest=<文件>   直接给 <前缀>_for_enrichment.txt
  --deg=<差异分析表>     必须包含全部基因（不是只含显著基因的表）。
                         （--deg 显式给出时优先于 --geo_dir 的自动识别）
  --logfc_col=logFC      排序列（默认自动识别 logFC / log2FC / stat / t）
  --stat_col=            想用 t 值而不是 logFC 排序时指定它（优先级高于 --logfc_col）
  --gene_col=<列名>      基因名所在列（默认自动识别）
  --id_type=SYMBOL       SYMBOL | ENTREZID | ENSEMBL
  --species=human|mouse  --orgdb_sqlite=<路径>  直接给 OrgDb 的 .sqlite（没装包时用）

基因集：
  --go=TRUE              跑 GO（BP/CC/MF 分别跑再合并，带 ONTOLOGY 列）
  --kegg=auto            auto | online | gmt | none
  --gmt_sets=H,C2        要跑的 MSigDB 集合（配合 --gmt_dir）
  --gmt=<文件>           直接指定 .gmt
  --gmt_dir=<目录>
  --min_gs=10 --max_gs=500

可复现性（fgsea 的 P 值来自随机置换，必须固定种子）：
  --seed=1234
  --nperm=10000          简化置换次数（nPermSimple）
  --eps=0                GSEA 的 p 值精度；0 = 算到最精确

输出：
  --outdir=./enrich --prefix=
  --pvalue=0.05 --padj=0.25   sig 列口径（GSEA 的 padj 惯例阈值是 0.25，不是 0.05）
  --plot=TRUE            出一张"显著基因集 NES 条形图"
  --plot_top=20          条形图画多少个
  --curve=0              额外给前 N 个基因集画 GSEA 富集曲线（0 = 不画）
  --curve_which=both|pos|neg   曲线挑哪些：both / 只要激活 / 只要抑制

示例：
  bash run_enrich.sh 02_enrich_gsea.R --deg=GSE62452_T_all.csv --species=human \\
    --outdir=./enrich --gmt_sets=H --plot=TRUE
'

opt <- parse_args()
if (any(c("help", "h") %in% names(opt)) || !length(opt)) {
  cat(USAGE); quit(status = 0)
}

deg_path <- oa(opt, "deg", NULL)
outdir <- oa(opt, "outdir", "./enrich")

## ---- 与 geo-microarray-analysis 对接 --------------------------------------
geo <- NULL
geo_mode <- !is.null(oa(opt, "geo_dir", NULL)) ||
  !is.null(oa(opt, "geo_manifest", NULL)) ||
  (!is.null(oa(opt, "geo_prefix", NULL)) && nzchar(oa(opt, "geo_prefix", "")))
if (geo_mode) {
  geo <- resolve_geo_inputs(oa(opt, "geo_dir", NULL), oa(opt, "geo_prefix", NULL),
                            oa(opt, "geo_manifest", NULL), verbose = FALSE)
  ## GSEA 一定要吃"全部基因"的表，所以优先用 _all.csv
  if (is.null(deg_path)) deg_path <- geo$all %||% geo$deg
  if (is.null(oa(opt, "prefix", NULL)) || !nzchar(oa(opt, "prefix", ""))) opt$prefix <- geo$prefix
  if (!("species" %in% names(opt)) && !is.null(geo$species)) opt$species <- geo$species
  if (!("logfc_col" %in% names(opt)) && !is.null(geo$logfc_col)) opt$logfc_col <- geo$logfc_col
  if (!("gene_col" %in% names(opt)) && !is.null(geo$gene_col)) opt$gene_col <- geo$gene_col
  if (!("outdir" %in% names(opt)) && !is.null(geo$dir)) opt$outdir <- file.path(geo$dir, "enrich")
}
if (is.null(deg_path)) {
  stop("必须给 --geo_dir=（从 GEO 差异分析接手）或 --deg=<含全部基因的差异分析表>。",
       call. = FALSE)
}
outdir <- oa(opt, "outdir", outdir)
prefix <- oa(opt, "prefix", NULL)
if (is.null(prefix) || !nzchar(prefix)) {
  prefix <- sub("\\.(csv|txt|tsv)$", "", basename(deg_path), ignore.case = TRUE)
}
sp        <- resolve_species(oa(opt, "species", "human"), oa(opt, "orgdb", NULL),
                             oa(opt, "kegg_organism", NULL))
id_type   <- toupper(oa(opt, "id_type", "SYMBOL"))
run_go    <- oa_lgl(opt, "go", TRUE)
kegg_mode <- tolower(oa(opt, "kegg", "auto"))
gmt_tokens<- oa_list(opt, "gmt_sets", character(0))
gmt_files <- oa_list(opt, "gmt", character(0))
gmt_dirs  <- oa_list(opt, "gmt_dir", character(0))
min_gs    <- oa_num(opt, "min_gs", 10)
max_gs    <- oa_num(opt, "max_gs", 500)
seed      <- oa_num(opt, "seed", 1234)
nperm     <- oa_num(opt, "nperm", 10000)
eps       <- oa_num(opt, "eps", 0)
sig_p     <- oa_num(opt, "pvalue", 0.05)
sig_padj  <- oa_num(opt, "padj", 0.25)
do_plot   <- oa_lgl(opt, "plot", TRUE)
plot_top  <- oa_num(opt, "plot_top", 20)
curve_n   <- oa_num(opt, "curve", 0)
curve_which <- tolower(oa(opt, "curve_which", "both"))
to_pdf    <- oa_lgl(opt, "pdf", TRUE)

dir_create(outdir)
start_log(file.path(outdir, paste0(prefix, "_02_GSEA.log")))
on.exit(close_log(), add = TRUE)

print_banner("02_enrich_gsea.R", "GSEA（全部基因按排序值富集）")
print_args(opt)
if (!is.null(geo)) print_geo_handoff(geo)
check_env()
note("")

for (p in c("clusterProfiler", "fgsea")) {
  if (!requireNamespace(p, quietly = TRUE)) stop("缺少 R 包：", p, call. = FALSE)
}
suppressPackageStartupMessages({ library(clusterProfiler) })
set.seed(seed)
note("随机种子已固定：", seed, "（fgsea 的 P 值依赖置换，不固定则不可复现）")

## ---- 读表 & 建排序向量 ----------------------------------------------------
deg <- read_table_robust(deg_path, "差异分析表")
gene_col <- pick_gene_col(deg, oa(opt, "gene_col", NULL))
stat_col <- oa(opt, "stat_col", NULL)
if (is.null(stat_col) || !nzchar(stat_col)) {
  stat_col <- pick_col(deg, oa(opt, "logfc_col", NULL),
                       c("logfc", "log2fc", "log2foldchange", "stat", "t", "zscore", "wald"))
}
if (is.na(stat_col) || is.na(gene_col)) {
  stop("无法识别基因列或排序列，请用 --gene_col / --logfc_col 指定。\n可用列名：",
       paste(names(deg), collapse = ", "), call. = FALSE)
}
note("[输入] ", basename(deg_path), "  ", nrow(deg), " 行")
note("[基因列] ", gene_col, "   [排序列] ", stat_col, "（**不做任何筛选**）")

g <- as.character(deg[[gene_col]])
s <- suppressWarnings(as.numeric(deg[[stat_col]]))
ok <- !is.na(g) & nzchar(g) & !is.na(s) & is.finite(s)
if (sum(!ok)) note(sprintf("剔除 %d 行（基因名为空或排序值为 NA/Inf）", sum(!ok)))
g <- g[ok]; s <- s[ok]

## 重复基因：默认保留 |stat| 最大的一条
dupg <- unique(g[duplicated(g)])
if (length(dupg)) {
  keep <- order(-abs(s))
  g2 <- g[keep]; s2 <- s[keep]
  first <- !duplicated(g2)
  note(sprintf("重复基因 %d 个，保留 |%s| 最大的一条", length(dupg), stat_col))
  g <- g2[first]; s <- s2[first]
}

gl <- setNames(s, g)
gl <- sort(gl, decreasing = TRUE)
if (anyDuplicated(names(gl))) stop("排序向量仍有重复基因名，请检查 --gene_col。", call. = FALSE)
note(sprintf("排序向量：%d 个基因，范围 %.3f ~ %.3f，中位数 %.3f",
             length(gl), min(gl), max(gl), stats::median(gl)))
writeLines(sprintf("%s\t%.6f", names(gl), as.numeric(gl)),
           file.path(outdir, paste0(prefix, "_GSEA_ranked_list.txt")), useBytes = TRUE)

## 是否需要 entrez 版
need_entrez <- (kegg_mode != "none") || run_go
glE <- NULL; sym_map <- NULL
if (need_entrez && id_type != "ENTREZID") {
  orgdb <- load_orgdb(sp$orgdb, oa(opt, "orgdb_sqlite", NULL))
  res <- suppressWarnings(suppressMessages(
    try(clusterProfiler::bitr(names(gl), fromType = id_type, toType = "ENTREZID",
                              OrgDb = orgdb, drop = TRUE), silent = TRUE)))
  if (inherits(res, "try-error") || !nrow(res)) {
    stop("ID 转换失败：一个基因都没转成 ENTREZID。检查 --id_type / --species。", call. = FALSE)
  }
  res <- res[!duplicated(res[[id_type]]), , drop = FALSE]
  sym_map <- setNames(as.character(res[[id_type]]), as.character(res$ENTREZID))
  v <- gl[as.character(res[[id_type]])]
  names(v) <- as.character(res$ENTREZID)
  v <- v[!is.na(v)]
  glE <- sort(v, decreasing = TRUE)
  note(sprintf("ENTREZID 版排序向量：%d / %d（丢失 %d）",
               length(glE), length(gl), length(gl) - length(glE)))
} else if (id_type == "ENTREZID") {
  glE <- gl
}

## ---- 三个富集分支 ---------------------------------------------------------

to_df <- function(obj, ontology = NULL, coll = NULL) {
  if (is.null(obj)) return(NULL)
  if (inherits(obj, "try-error")) return(NULL)
  df <- try(as.data.frame(obj), silent = TRUE)
  if (inherits(df, "try-error") || is.null(df) || !nrow(df)) return(NULL)
  if (!"ID" %in% names(df)) df$ID <- rownames(df)
  if (!is.null(ontology)) df$ONTOLOGY <- ontology
  if (is.null(df$ONTOLOGY)) df$ONTOLOGY <- "GSEA"
  if (!is.null(coll)) { df$collection <- coll; df$ONTOLOGY <- coll }
  ## GSEA 的 leading edge 列里是 entrez，换成 symbol 更好读
  if (!is.null(sym_map) && "core_enrichment" %in% names(df)) {
    conv <- function(x) {
      if (is.na(x) || !nzchar(x)) return(NA_character_)
      ids <- strsplit(x, "/", fixed = TRUE)[[1]]
      hit <- sym_map[ids]; hit[is.na(hit)] <- ids[is.na(hit)]
      paste(hit, collapse = "/")
    }
    df$core_symbol <- vapply(as.character(df$core_enrichment), conv, "", USE.NAMES = FALSE)
  }
  gsea_sig_mark(df, sig_p, sig_padj)
}

gsea_sig_mark <- function(df, p, padj) {
  pv <- suppressWarnings(as.numeric(df$pvalue))
  pa <- suppressWarnings(as.numeric(df$p.adjust))
  df$sig <- !is.na(pv) & pv <= p & !is.na(pa) & pa <= padj
  df$direction <- ifelse(is.na(df$NES), NA_character_,
                         ifelse(df$NES > 0, "activated", "suppressed"))
  want <- c("ONTOLOGY", "collection", "ID", "Description", "setSize", "enrichmentScore",
            "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge",
            "core_enrichment", "core_symbol", "direction", "sig")
  df[, c(intersect(want, names(df)), setdiff(names(df), want)), drop = FALSE]
}

res_all <- list()

## --- GO ---
if (run_go) {
  note("")
  note("---- GSEA: GO ----")
  orgdb <- load_orgdb(sp$orgdb, oa(opt, "orgdb_sqlite", NULL))
  gs <- setNames(list(), character(0))
  for (o in c("BP", "CC", "MF")) {
    r <- suppressMessages(try(
      clusterProfiler::gseGO(geneList = glE, OrgDb = orgdb, keyType = "ENTREZID",
                             ont = o, minGSSize = min_gs, maxGSSize = max_gs,
                             pvalueCutoff = 1, verbose = FALSE,
                             eps = eps, nPermSimple = nperm, seed = TRUE),
      silent = TRUE))
    d <- to_df(r, o)
    if (!is.null(d)) {
      note(sprintf("  GO:%-3s 被检验 %5d 条，显著 %4d 条", o, nrow(d), sum(d$sig)))
      gs[[o]] <- d
    } else {
      note(sprintf("  GO:%-3s 无结果", o))
    }
  }
  if (length(gs)) {
    dgo <- do.call(rbind, gs)
    res_all$GO <- dgo
    write_csv_out(dgo, file.path(outdir, paste0(prefix, "_GSEA_GO.csv")))
    write_csv_out(dgo[dgo$sig, , drop = FALSE],
                  file.path(outdir, paste0(prefix, "_GSEA_GO_sig.csv")))
  }
}

## --- KEGG ---
kegg_src <- kegg_mode
if (kegg_mode != "none" && !is.null(glE)) {
  note("")
  note("---- GSEA: KEGG ----")
  ek <- NULL
  if (kegg_mode %in% c("auto", "online")) {
    ek <- suppressMessages(try(
      clusterProfiler::gseKEGG(geneList = glE, organism = sp$kegg, keyType = "kegg",
                               minGSSize = min_gs, maxGSSize = max_gs,
                               pvalueCutoff = 1, verbose = FALSE,
                               eps = eps, nPermSimple = nperm, seed = TRUE),
      silent = TRUE))
    if (inherits(ek, "try-error")) {
      note("  KEGG 在线 GSEA 失败：", attr(ek, "condition")$message)
      ek <- NULL
    }
  }
  dk <- NULL
  if (!is.null(ek)) {
    dk <- to_df(ek, "KEGG"); kegg_src <- "online"
  } else if (kegg_mode %in% c("auto", "gmt")) {
    kg <- try(resolve_gmt("KEGG", gmt_dirs, character(0)), silent = TRUE)
    if (!inherits(kg, "try-error") && !is.null(kg)) {
      r <- suppressMessages(try(
        clusterProfiler::GSEA(geneList = gl, TERM2GENE = kg[, c("term", "gene"), drop = FALSE],
                              minGSSize = min_gs, maxGSSize = max_gs,
                              pvalueCutoff = 1, verbose = FALSE,
                              eps = eps, nPermSimple = nperm, seed = TRUE),
        silent = TRUE))
      dk <- to_df(r, "KEGG", "KEGG"); kegg_src <- "gmt(MSigDB KEGG_LEGACY)"
      note("  ⚠ 回退到 MSigDB 的 KEGG 子集（KEGG_LEGACY，186 条），",
           "不含 PI3K-Akt / TNF 通路。")
    }
  }
  if (!is.null(dk)) {
    note(sprintf("  KEGG[%s] 被检验 %5d 条，显著 %4d 条", kegg_src, nrow(dk), sum(dk$sig)))
    res_all$KEGG <- dk
    write_csv_out(dk, file.path(outdir, paste0(prefix, "_GSEA_KEGG.csv")))
    write_csv_out(dk[dk$sig, , drop = FALSE],
                  file.path(outdir, paste0(prefix, "_GSEA_KEGG_sig.csv")))
  }
}

## --- 自定义 GMT ---
if (length(gmt_tokens) || length(gmt_files)) {
  note("")
  note("---- GSEA: MSigDB ----")
  t2g <- resolve_gmt(gmt_tokens, gmt_dirs, gmt_files)
  if (!is.null(t2g)) {
    kind <- gmt_id_kind(t2g, names(gl), names(glE %||% gl))
    glg <- if (kind == "symbol") gl else glE
    r <- suppressMessages(try(
      clusterProfiler::GSEA(geneList = glg, TERM2GENE = t2g[, c("term", "gene"), drop = FALSE],
                            minGSSize = min_gs, maxGSSize = max_gs,
                            pvalueCutoff = 1, verbose = FALSE,
                            eps = eps, nPermSimple = nperm, seed = TRUE),
      silent = TRUE))
    d <- to_df(r, NULL, NULL)
    if (!is.null(d)) {
      cm <- t2g[!duplicated(t2g$term), c("term", "collection")]
      d$collection <- cm$collection[match(d$ID, cm$term)]
      d$ONTOLOGY <- d$collection
      d <- gsea_sig_mark(d, sig_p, sig_padj)
      note(sprintf("  MSigDB 被检验 %5d 条，显著 %4d 条", nrow(d), sum(d$sig)))
      res_all$GMT <- d
      write_csv_out(d, file.path(outdir, paste0(prefix, "_GSEA_GMT.csv")))
      write_csv_out(d[d$sig, , drop = FALSE],
                    file.path(outdir, paste0(prefix, "_GSEA_GMT_sig.csv")))
    }
  }
}

if (!length(res_all)) stop("三个分支都没有产出任何 GSEA 结果，检查输入。", call. = FALSE)
saveRDS(res_all, file.path(outdir, paste0(prefix, "_GSEA_results.rds")))

## ---- 汇总 -----------------------------------------------------------------
note("")
note(strrep("=", 72))
note("GSEA 汇总（sig: pvalue<=", sig_p, " & p.adjust<=", sig_padj, "）")
note(strrep("=", 72))
sm <- do.call(rbind, lapply(names(res_all), function(k) {
  d <- res_all[[k]]
  data.frame(db = k, tested = nrow(d), sig = sum(d$sig),
             activated = sum(d$sig & d$NES > 0), suppressed = sum(d$sig & d$NES < 0),
             stringsAsFactors = FALSE)
}))
print(sm, row.names = FALSE)
write_csv_out(sm, file.path(outdir, paste0(prefix, "_GSEA_summary.csv")))
writeLines(capture.output(print(sm, row.names = FALSE)),
           file.path(outdir, paste0(prefix, "_GSEA_summary.txt")), useBytes = TRUE)

## ---- 出图：显著基因集 NES 条形图 ------------------------------------------
if (do_plot) {
  ## 三个库的结果列集不完全一样（GMT 多 collection、GO 多 ONTOLOGY、
  ## core_symbol 只在有 ID 映射时才有），直接 rbind 会报
  ## "numbers of columns of arguments do not match" —— 先对齐列。
  keep_cols <- c("ID", "Description", "setSize", "enrichmentScore", "NES",
                 "pvalue", "p.adjust", "qvalue", "direction")
  parts2 <- lapply(names(res_all), function(k) {
    d <- res_all[[k]]
    d <- d[as.logical(d$sig), , drop = FALSE]
    if (!nrow(d)) return(NULL)
    for (m in setdiff(keep_cols, names(d))) d[[m]] <- NA
    d <- d[, keep_cols, drop = FALSE]
    d$db <- k
    d
  })
  parts2 <- parts2[!vapply(parts2, is.null, TRUE)]
  if (!length(parts2)) {
    note("所有分支都没有显著基因集，跳过条形图（可放宽 --pvalue / --padj）。")
  } else {
    all <- do.call(rbind, parts2)
    all$NES <- suppressWarnings(as.numeric(all$NES))
    all <- all[!is.na(all$NES), , drop = FALSE]
    all <- all[order(-abs(all$NES)), , drop = FALSE]
    top <- head(all, plot_top)
    top$lab <- paste0(top$db, " | ", top$Description,
                      ifelse(nchar(top$Description) > 60,
                             paste0(" (", top$setSize, " genes)"), ""))
    top$lab <- factor(top$lab, levels = rev(unique(top$lab)))
    top$sign <- ifelse(top$NES > 0, "activated (NES > 0)", "suppressed (NES < 0)")
    pal2 <- c("activated (NES > 0)" = "#B2182B", "suppressed (NES < 0)" = "#2166AC")
    p <- ggplot2::ggplot(top, ggplot2::aes(x = .data$NES, y = .data$lab, fill = .data$sign)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey35") +
      ggplot2::scale_fill_manual(values = pal2, name = NULL) +
      ggplot2::labs(x = "Normalized enrichment score (NES)", y = NULL,
                    title = paste0("GSEA top ", nrow(top), " gene sets by |NES|")) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(face = "bold"),
                     axis.text.y = ggplot2::element_text(size = 8))
    for (ext in c("pdf", "png")) {
      f <- file.path(outdir, paste0(prefix, "_GSEA_top_NES.", ext))
      if (ext == "pdf") ggplot2::ggsave(f, p, width = 9, height = max(3, 0.26 * nrow(top) + 1.6))
      else ggplot2::ggsave(f, p, width = 9, height = max(3, 0.26 * nrow(top) + 1.6), dpi = 300)
      note("已出图：", f)
    }
  }
}

## ---- 可选：GSEA 富集曲线 ---------------------------------------------------
if (curve_n > 0) {
  need <- c("enrichplot")
  if (!requireNamespace(need, quietly = TRUE)) {
    note("缺 enrichplot，跳过 --curve。")
  } else {
    for (k in names(res_all)) {
      d <- res_all[[k]][res_all[[k]]$sig, , drop = FALSE]
      if (!nrow(d)) next
      d <- d[order(-abs(d$NES)), , drop = FALSE]
      if (curve_which == "pos") d <- d[d$NES > 0, , drop = FALSE]
      if (curve_which == "neg") d <- d[d$NES < 0, , drop = FALSE]
      d <- head(d, curve_n)
      if (!nrow(d)) next
      ids <- d$ID
      pal <- ifelse(d$NES > 0, "#B2182B", "#2166AC")
      f <- file.path(outdir, paste0(prefix, "_GSEA_curve_", k, ".pdf"))
      grDevices::pdf(f, width = 8, height = 2.4 * nrow(d) + 2)
      for (i in seq_along(ids)) {
        pp <- try(enrichplot::gseaplot2(res_all[[k]], geneSetID = ids[i],
                                        color = pal[i], title = ids[i]), silent = TRUE)
        if (!inherits(pp, "try-error")) print(pp)
      }
      grDevices::dev.off()
      note("已出图：", f)
    }
  }
}

note("")
note("产出目录：", normalizePath(outdir, mustWork = FALSE))
note("耗时：", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " 秒")
