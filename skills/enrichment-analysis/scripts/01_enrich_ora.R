#!/usr/bin/env Rscript
## ===========================================================================
## 01_enrich_ora.R —— 过表达富集分析（ORA）：GO / KEGG / 自定义 MSigDB 基因集
##
## 输入：一份差异基因列表（直接给基因名），或一份差异分析表（按显著性筛出基因）。
## 输出：每个方向（all / up / down）一张完整富集表 + 一张只含显著条目的表。
##
## 关键设计（和参考代码 kegg美化.R 的差异，都是故意的）：
##   1. 富集一律用 cutoff = 1 计算，把"所有被检验的条目"都留下来，
##      再用 `sig` 列标出是否通过 pvalue/p.adjust/qvalue 三道阈值。
##      好处：出图时想放宽/收紧阈值不用重跑（KEGG 每次跑都要联网下载）。
##   2. 默认穷尽 GO 的 BP/CC/MF（ont="ALL"），并额外提供 MSigDB 自定义集合。
##   3. 默认背景集 = clusterProfiler 默认（与参考代码一致）；
##      --universe=deg 可改用"所有被检验的基因"作背景（更严格，推荐发文章用）。
## ===========================================================================

t0 <- Sys.time()

## ---- 定位并加载公共库 ------------------------------------------------------
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
01_enrich_ora.R —— 过表达富集分析（GO / KEGG / MSigDB GMT）

输入（三选一）：
  --geo_dir=<目录>        ★ 直接从 geo-microarray-analysis 的差异分析产物接手。
  --geo_prefix=<前缀>     geo 02_deg_plots.R 的 --prefix（不给就自动找目录里唯一的一个）。
                          会自动认出 ORA 吃 <前缀>_DEG.csv、GSEA 吃 <前缀>_all.csv，
                          并自动对上 gene/logFC/p/change 列名；默认**不再二次筛选**，
                          all/up/down 三套直接按 geo 的 change 列拆（无符号问题）。
  --geo_manifest=<文件>   直接给 <前缀>_for_enrichment.txt（geo 02 会生成）
  --gene_list=<文件>     直接给基因列表（每行一个，或用逗号/空白分隔）。给了它就跑一套名为
                         --list_name（默认 list）的结果。
  --deg=<差异分析表>     给差异分析表（csv/txt），按 --ora_p / --ora_p_type / --ora_logfc
                         筛出显著基因，再按 --direction 拆成 all / up / down 三套。
                         （--deg 显式给出时优先于 --geo_dir 的自动识别）

物种与 ID：
  --species=human|mouse  默认 human（也接受 hsa/mmu/人/小鼠）
  --orgdb_sqlite=<路径>  直接给 OrgDb 的 .sqlite（没装包时用，见 references/pitfalls.md 5.2）
  --id_type=SYMBOL       SYMBOL | ENTREZID | ENSEMBL | UNIPROT 等（默认 SYMBOL）
  --gene_col=<列名>      差异表里基因名所在列（默认自动识别）
  --logfc_col=logFC      差异表里 logFC 列（默认自动识别）
  --p_col=               差异表里 p 值列（默认自动识别 adj.P.Val -> P.Value）

筛显著基因（只影响 ORA 的输入基因，不影响 GSEA）：
  --ora_p=0.05           显著性阈值
  --ora_p_type=auto      auto | adj.P.Val | P.Value
  --ora_logfc=0          |logFC| 阈值（0 = 不筛）
  --direction=all,up,down  要跑哪几套
  --split_by=logfc       logfc（按 --ora_logfc 拆 up/down）| change（按 geo 的 change 列拆，
                         UP/DOWN/NOT，geo 模式下自动启用）

富集方法与阈值：
  --ont=ALL              GO 本体：ALL | BP | CC | MF
  --kegg=auto            auto | online | gmt | none
  --gmt_sets=            例如 H,C2 或 KEGG,GOBP（需 --gmt_dir 或默认目录里有 *.gmt）
  --gmt=                直接指定 .gmt 文件（逗号分隔）
  --gmt_dir=            GMT 所在目录
  --universe=none        none（clusterProfiler 默认，与参考代码一致）| deg（用所有被检验基因）
  --min_gs=10 --max_gs=500   基因集大小范围
  --pvalue=0.05 --padj=0.05 --qvalue=0.05   sig 列的判定口径（富集本身用 cutoff=1 算）

输出：
  --outdir=./enrich --prefix=<自动取输入文件名>
  --list_name=list       给 --gene_list 时这套结果的命名
  --export_sig=TRUE      是否额外导出只含显著条目的 *_sig.csv
  --save_rdata=TRUE      是否保存 <prefix>_ORA.RData

示例：
  bash run_enrich.sh 01_enrich_ora.R --deg=GSE62452_T_all.csv --species=human \\
    --outdir=./enrich --gmt_sets=H
'

## ---- 参数 -----------------------------------------------------------------
opt <- parse_args()
if (any(c("help", "h") %in% names(opt)) || !length(opt)) {
  cat(USAGE); quit(status = 0)
}

outdir   <- oa(opt, "outdir", "./enrich")
deg_path <- oa(opt, "deg", NULL)
list_path<- oa(opt, "gene_list", NULL)

## ---- 与 geo-microarray-analysis 对接（--geo_dir / --geo_prefix / --geo_manifest）----
geo <- NULL
geo_mode <- !is.null(oa(opt, "geo_dir", NULL)) ||
  !is.null(oa(opt, "geo_manifest", NULL)) ||
  (!is.null(oa(opt, "geo_prefix", NULL)) && nzchar(oa(opt, "geo_prefix", "")))
if (geo_mode) {
  geo <- resolve_geo_inputs(oa(opt, "geo_dir", NULL), oa(opt, "geo_prefix", NULL),
                            oa(opt, "geo_manifest", NULL), verbose = FALSE)
  ## 显式 --deg / --gene_list 优先，geo 只负责填空
  if (is.null(deg_path) && is.null(list_path)) {
    src <- tolower(oa(opt, "geo_ora_source", "deg"))
    deg_path <- if (src == "all") geo$all else geo$deg
  }
  if (is.null(oa(opt, "prefix", NULL)) || !nzchar(oa(opt, "prefix", ""))) {
    opt$prefix <- geo$prefix
  }
  if (!("species" %in% names(opt)) && !is.null(geo$species)) opt$species <- geo$species
  if (!("p_col" %in% names(opt)) && !is.null(geo$p_col)) opt$p_col <- geo$p_col
  if (!("logfc_col" %in% names(opt)) && !is.null(geo$logfc_col)) opt$logfc_col <- geo$logfc_col
  if (!("gene_col" %in% names(opt)) && !is.null(geo$gene_col)) opt$gene_col <- geo$gene_col
  ## geo 那次筛基因用的 p 口径（P.Value / adj.P.Val）跟着带过来，日志才好对照
  if (!("ora_p_type" %in% names(opt)) && !is.null(geo$deg_p_type)) opt$ora_p_type <- geo$deg_p_type
  ## geo 的 DEG 表已经筛过了：默认不再二次筛选，改按 change 列拆 up/down
  if (!("ora_p" %in% names(opt))) opt$ora_p <- "1"
  if (!("ora_logfc" %in% names(opt))) opt$ora_logfc <- "0"
  if (!("split_by" %in% names(opt)) && !is.null(geo$change_col)) opt$split_by <- "change"
  ## outdir 默认落在 geo 产物目录下的 enrich/
  if (!("outdir" %in% names(opt)) && !is.null(geo$dir)) opt$outdir <- file.path(geo$dir, "enrich")
}

if (is.null(deg_path) && is.null(list_path)) {
  stop("必须给 --geo_dir=（从 GEO 差异分析接手）、--gene_list=<基因列表>、"
       , "--deg=<差异分析表> 三者之一。", call. = FALSE)
}
outdir <- oa(opt, "outdir", outdir)
prefix <- oa(opt, "prefix", NULL)
if (is.null(prefix) || !nzchar(prefix)) {
  src <- if (!is.null(deg_path)) deg_path else list_path
  prefix <- sub("\\.(csv|txt|tsv)$", "", basename(src), ignore.case = TRUE)
}

sp       <- resolve_species(oa(opt, "species", "human"),
                            oa(opt, "orgdb", NULL), oa(opt, "kegg_organism", NULL))
id_type  <- toupper(oa(opt, "id_type", "SYMBOL"))
directions <- oa_list(opt, "direction", c("all", "up", "down"))
ont      <- toupper(oa(opt, "ont", "ALL"))
kegg_mode<- tolower(oa(opt, "kegg", "auto"))
ora_p    <- oa_num(opt, "ora_p", 0.05)
ora_logfc<- oa_num(opt, "ora_logfc", 0)
universe_mode <- tolower(oa(opt, "universe", "none"))
min_gs   <- oa_num(opt, "min_gs", 10)
max_gs   <- oa_num(opt, "max_gs", 500)
sig_p    <- oa_num(opt, "pvalue", 0.05)
sig_padj <- oa_num(opt, "padj", 0.05)
sig_q    <- oa_num(opt, "qvalue", 0.05)
export_sig <- oa_lgl(opt, "export_sig", TRUE)
save_rdata <- oa_lgl(opt, "save_rdata", TRUE)
list_name  <- oa(opt, "list_name", "list")
gmt_tokens <- oa_list(opt, "gmt_sets", character(0))
gmt_files  <- oa_list(opt, "gmt", character(0))
gmt_dirs   <- oa_list(opt, "gmt_dir", character(0))

dir_create(outdir)
start_log(file.path(outdir, paste0(prefix, "_01_ORA.log")))
on.exit(close_log(), add = TRUE)

print_banner("01_enrich_ora.R", "ORA 富集分析（GO / KEGG / MSigDB）")
print_args(opt)
if (!is.null(geo)) print_geo_handoff(geo)
check_env()
note("")

## ---- 需要哪些 R 包 ---------------------------------------------------------
need <- c("clusterProfiler", "AnnotationDbi")
for (p in need) if (!requireNamespace(p, quietly = TRUE)) {
  stop("缺少 R 包：", p, call. = FALSE)
}
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(AnnotationDbi)
})
note("物种：", sp$label, "  OrgDb=", sp$orgdb, "  KEGG=", sp$kegg)
orgdb <- load_orgdb(sp$orgdb, oa(opt, "orgdb_sqlite", NULL))

## ---- 组装基因集 -----------------------------------------------------------
## sets: named list，名字就是输出的方向标签
sets <- list()
meta <- data.frame()

if (!is.null(list_path)) {
  gl <- readLines(list_path, warn = FALSE)
  gl <- unlist(strsplit(gl, "[,\t;[:space:]]+"))
  gl <- trimws(gl); gl <- gl[nzchar(gl)]
  gl <- unique(gl)
  if (!length(gl)) stop("--gene_list 文件里没有解析到任何基因：", list_path, call. = FALSE)
  sets[[list_name]] <- gl
  note(sprintf("[基因集] %-10s <- %s（直接给的列表，%d 个基因）",
               list_name, basename(list_path), length(gl)))
}

if (!is.null(deg_path)) {
  deg <- read_table_robust(deg_path, "差异分析表")
  gene_col  <- pick_gene_col(deg, oa(opt, "gene_col", NULL))
  logfc_col <- pick_col(deg, oa(opt, "logfc_col", NULL),
                        c("logfc", "log2fc", "log2foldchange", "logfc", "fc"))
  pt <- tolower(oa(opt, "ora_p_type", "auto"))
  ## 列解析优先级：--p_col（显式列名） > --ora_p_type 给的列名 > auto 候选
  ## （注意 --p_col 以前只是写在帮助里、实际没接上，2026-09-18 补上）
  p_col <- pick_col(deg, oa(opt, "p_col", NULL), character(0))
  if (is.na(p_col)) {
    p_col <- if (pt == "auto") {
      pick_col(deg, NULL, c("adj.p.val", "padj", "p.adjust", "fdr", "adj.pvalue",
                            "p.value", "pvalue", "p_val"))
    } else {
      pick_col(deg, oa(opt, "ora_p_type", NULL), character(0))
    }
  }
  if (is.na(gene_col)) stop("无法识别基因名列，请用 --gene_col 指定。", call. = FALSE)
  if (is.na(p_col))    stop("无法识别 p 值列，请用 --p_col 指定。可用列名：",
                            paste(names(deg), collapse = ", "), call. = FALSE)
  note(sprintf("[列识别] 基因名 = %s   显著性 = %s   logFC = %s",
               gene_col, p_col, if (is.na(logfc_col)) "（未找到）" else logfc_col))
  genes_all <- unique(as.character(deg[[gene_col]]))
  genes_all <- genes_all[!is.na(genes_all) & nzchar(genes_all)]

  pv <- suppressWarnings(as.numeric(deg[[p_col]]))
  keep <- !is.na(pv) & pv <= ora_p
  note(sprintf("[输入表] %s  共 %d 行 / %d 个唯一基因", basename(deg_path),
               nrow(deg), length(genes_all)))
  if (ora_p >= 1) {
    note(sprintf("[筛选] 不二次筛选（--ora_p=1）：表里 %d 行全部进入下一步", nrow(deg)))
  } else {
    note(sprintf("[显著性] %s <= %s  ->  %d 个基因", p_col, format(ora_p), sum(keep)))
  }
  g_sig <- unique(as.character(deg[[gene_col]][keep]))

  ## ---- 拆 all / up / down -------------------------------------------------
  ## split_by=logfc（默认）：用 |logFC| >= --ora_logfc 拆
  ## split_by=change：用 geo 的 change 列（UP / DOWN / NOT）拆 —— geo 模式下自动启用，
  ##                  这样 up/down 与 geo 那次差异分析的口径**逐字一致**
  split_by <- tolower(oa(opt, "split_by", "logfc"))
  if (!split_by %in% c("logfc", "change")) {
    stop("--split_by 只能是 logfc 或 change，收到：", split_by, call. = FALSE)
  }
  chg_col <- pick_col(deg, oa(opt, "change_col", NULL),
                      c("change", "regulation", "updown", "direction"))
  if (split_by == "change" && is.na(chg_col)) {
    note("⚠ 指定了 --split_by=change，但表里没有 change 列，退回按 logFC 拆。")
    split_by <- "logfc"
  }

  lfc <- NULL
  if (!is.na(logfc_col)) lfc <- suppressWarnings(as.numeric(deg[[logfc_col]]))

  if (split_by == "change") {
    chg <- toupper(trimws(as.character(deg[[chg_col]])))
    keep_chg <- keep & chg %in% c("UP", "DOWN")
    note(sprintf("[拆分口径] change 列「%s」：UP %d / DOWN %d / 其他 %d",
                 chg_col, sum(keep_chg & chg == "UP"), sum(keep_chg & chg == "DOWN"),
                 sum(!chg %in% c("UP", "DOWN", "NOT"))))
    if (sum(keep_chg) == 0) {
      note("⚠ change 列里没有 UP/DOWN（当前取值：",
           paste(unique(chg), collapse = "/"), "），退回按 logFC 拆。")
      split_by <- "logfc"
    }
  }

  if (split_by == "change") {
    chg <- toupper(trimws(as.character(deg[[chg_col]])))
    for (d in directions) {
      d <- tolower(d)
      if (d == "all") {
        sets[["all"]] <- unique(as.character(deg[[gene_col]][keep & chg %in% c("UP", "DOWN")]))
      } else if (d == "up") {
        sets[["up"]] <- unique(as.character(deg[[gene_col]][keep & chg == "UP"]))
      } else if (d == "down") {
        sets[["down"]] <- unique(as.character(deg[[gene_col]][keep & chg == "DOWN"]))
      } else if (d == "deg") {
        sets[["deg"]] <- genes_all
      } else {
        stop("不认识的 --direction 取值：", d, "（可用 all / up / down / deg）", call. = FALSE)
      }
    }
    note(sprintf("  all(%d) = up(%d) + down(%d)  %s",
                 length(sets[["all"]]),
                 if (is.null(sets[["up"]])) 0L else length(sets[["up"]]),
                 if (is.null(sets[["down"]])) 0L else length(sets[["down"]]),
                 if (length(sets[["all"]]) ==
                     (if (is.null(sets[["up"]])) 0L else length(sets[["up"]])) +
                     (if (is.null(sets[["down"]])) 0L else length(sets[["down"]]))) "✓" else "⚠ 对不上"))
  } else {
  ## 注意：--ora_logfc 对 all 也生效，这样 all 恒等于 up ∪ down，三套结果的基因数能对上账。
  ## （若想让 all = 全部显著基因、不管 logFC，就写 --ora_logfc=0。）
  keep_lfc <- keep
  if (!is.null(lfc)) keep_lfc <- keep & !is.na(lfc) & abs(lfc) >= ora_logfc

  for (d in directions) {
    d <- tolower(d)
    if (d == "all") {
      sets[["all"]] <- if (!is.null(lfc)) unique(as.character(deg[[gene_col]][keep_lfc])) else g_sig
    } else if (d %in% c("up", "down")) {
      if (is.null(lfc)) {
        note("⚠ 表里没找到 logFC 列，跳过 ", d, " 这一套（请用 --logfc_col 指定）")
        next
      }
      sel <- keep & !is.na(lfc) &
        (if (d == "up") lfc >= ora_logfc else lfc <= -ora_logfc)
      sets[[d]] <- unique(as.character(deg[[gene_col]][sel]))
    } else if (d == "deg") {
      sets[["deg"]] <- genes_all
    } else {
      stop("不认识的 --direction 取值：", d, "（可用 all / up / down / deg）", call. = FALSE)
    }
  }
  if (!is.null(lfc) && ora_logfc > 0) {
    note(sprintf("[logFC 阈值] |logFC| >= %s（对 all / up / down 都生效）", ora_logfc))
    note(sprintf("  因而 all(%d) = up(%d) + down(%d)  %s",
                 length(sets[["all"]]),
                 if (is.null(sets[["up"]])) 0L else length(sets[["up"]]),
                 if (is.null(sets[["down"]])) 0L else length(sets[["down"]]),
                 if (!is.null(sets[["all"]]) && !is.null(sets[["up"]]) && !is.null(sets[["down"]])) {
                   if (length(sets[["all"]]) == length(sets[["up"]]) + length(sets[["down"]])) "✓" else "⚠ 对不上"
                 } else ""))
  }
  }
  ## 背景集
  universe <- NULL
  if (universe_mode == "deg") {
    universe <- genes_all
    note("[背景集] --universe=deg：使用输入表里全部 ", length(universe), " 个基因作背景")
  } else if (universe_mode != "none") {
    if (file.exists(universe_mode)) {
      universe <- unique(trimws(readLines(universe_mode, warn = FALSE)))
      universe <- universe[nzchar(universe)]
      note("[背景集] 来自文件 ", universe_mode, "：", length(universe), " 个基因")
    } else {
      stop("--universe 只能是 none / deg / 一个基因列表文件路径", call. = FALSE)
    }
  }
} else {
  universe <- NULL
}

sets <- sets[vapply(sets, length, 1L) > 0]
if (!length(sets)) stop("没有任何可用的基因集（筛完是空的）。请放宽 --ora_p / --ora_logfc。",
                        call. = FALSE)

## 基因集过大的告警：ORA 的意义是"少数基因里某条通路被过度代表"，
## 上千个基因做 ORA 基本等于把背景又抄一遍，几乎什么都显著。
big <- names(sets)[vapply(sets, length, 1L) > 2000]
if (length(big)) {
  note("")
  note("⚠ 基因集偏大：",
       paste(sprintf("%s(%d 个)", big, vapply(sets[big], length, 1L)), collapse = "、"))
  note("  ORA 的假设是「一撮基因里某条通路被过度代表」。几千个基因做 ORA 往往只是把")
  note("  全基因组背景重抄一遍，结果几乎什么都显著、也就没有信息量。建议加上")
  note("    --ora_logfc=1          （补一个 |logFC| 阈值，与 filtered_DEGs 那套口径一致）")
  note("    --ora_p_type=P.Value   （或换一个更宽的 p 口径）")
  note("  或者干脆把已经筛好的 DEG 表 / 基因列表喂进来（--gene_list）。")
  note("")
}

## 导出实际使用的基因清单（可复现性的关键）
for (nm in names(sets)) {
  writeLines(sets[[nm]], file.path(outdir, sprintf("%s_ORA_%s_input_genes.txt", prefix, nm)),
             useBytes = TRUE)
}

## ---- ID 转换 --------------------------------------------------------------
to_entrez <- function(genes, tag) {
  if (id_type == "ENTREZID") {
    return(list(entrez = genes, map = NULL, n_in = length(genes), n_out = length(genes)))
  }
  if (id_type == "SYMBOL" && !requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
      !requireNamespace(sp$orgdb, quietly = TRUE)) {
    stop("需要物种注释包做 ID 转换。", call. = FALSE)
  }
  res <- suppressWarnings(suppressMessages(
    try(clusterProfiler::bitr(genes, fromType = id_type, toType = "ENTREZID",
                              OrgDb = orgdb, drop = TRUE), silent = TRUE)
  ))
  if (inherits(res, "try-error") || is.null(res) || !nrow(res)) {
    stop("ID 转换失败：", tag, " 的 ", length(genes), " 个 ", id_type,
         " 一个都没转成 ENTREZID。\n请检查 --id_type 与 --species 是否配对。", call. = FALSE)
  }
  res <- res[!duplicated(res[[id_type]]), , drop = FALSE]
  list(entrez = unique(res$ENTREZID), map = res,
       n_in = length(genes), n_out = length(unique(res$ENTREZID)))
}

## ---- 三种富集 -------------------------------------------------------------

run_go <- function(entrez) {
  if (is.null(entrez) || !length(entrez)) return(NULL)
  suppressMessages(try(
    clusterProfiler::enrichGO(gene = entrez, OrgDb = orgdb, keyType = "ENTREZID",
                              ont = ont, pAdjustMethod = "BH",
                              pvalueCutoff = 1, qvalueCutoff = 1,
                              minGSSize = min_gs, maxGSSize = max_gs,
                              universe = if (is.null(universe)) NULL else universe,
                              readable = FALSE),
    silent = TRUE))
}

run_kegg_online <- function(entrez) {
  suppressMessages(try(
    clusterProfiler::enrichKEGG(gene = entrez, organism = sp$kegg, keyType = "kegg",
                                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                                minGSSize = min_gs, maxGSSize = max_gs,
                                use_internal_data = FALSE),
    silent = TRUE))
}

run_enricher <- function(genes_use, t2g, uni = NULL) {
  suppressMessages(try(
    clusterProfiler::enricher(gene = genes_use, TERM2GENE = t2g[, c("term", "gene"), drop = FALSE],
                              pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                              minGSSize = min_gs, maxGSSize = max_gs,
                              universe = uni),
    silent = TRUE))
}

## GMT 准备（只读一次）
t2g <- NULL
if (length(gmt_tokens) || length(gmt_files)) {
  note("")
  note("[MSigDB 集合]")
  t2g <- resolve_gmt(gmt_tokens, gmt_dirs, gmt_files)
  if (!is.null(t2g)) {
    note(sprintf("  合计 %d 个基因集 / %d 个唯一基因",
                 length(unique(t2g$term)), length(unique(t2g$gene))))
  }
}
## KEGG 离线兜底所需
kegg_gmt <- NULL
ensure_kegg_gmt <- function() {
  if (!is.null(kegg_gmt)) return(kegg_gmt)
  t <- try(resolve_gmt("KEGG", gmt_dirs, character(0)), silent = TRUE)
  if (inherits(t, "try-error") || is.null(t)) return(NULL)
  kegg_gmt <<- t
  t
}

## ---- 主循环 ---------------------------------------------------------------
all_res <- list()
summary_rows <- list()

for (nm in names(sets)) {
  genes <- sets[[nm]]
  note("")
  note(strrep("-", 72))
  note("基因集：", nm, "  （", length(genes), " 个 ", id_type, "）")
  note(strrep("-", 72))

  conv <- to_entrez(genes, nm)
  note(sprintf("ID 转换：%s -> ENTREZID  %d / %d 成功（丢失 %d，%.1f%%）",
               id_type, conv$n_out, conv$n_in, conv$n_in - conv$n_out,
               100 * (conv$n_in - conv$n_out) / conv$n_in))
  entrez <- conv$entrez
  if (!is.null(conv$map)) {
    write_csv_out(conv$map, file.path(outdir, sprintf("%s_id_map_%s.csv", prefix, nm)))
  }
  sym_map <- if (!is.null(conv$map)) {
    setNames(as.character(conv$map[[id_type]]), as.character(conv$map$ENTREZID))
  } else {
    setNames(entrez, entrez)
  }
  ## GMT 用哪一版的基因列表
  genes_symbol <- if (id_type == "SYMBOL") genes else {
    m <- sym_map; names(m) <- as.character(conv$map$ENTREZID)
    inv <- setNames(names(m), as.character(m))
    unname(inv[entrez])
  }
  genes_symbol <- genes_symbol[!is.na(genes_symbol)]

  writeLines(entrez, file.path(outdir, sprintf("%s_ORA_%s_input_entrez.txt", prefix, nm)),
             useBytes = TRUE)

  ## --- GO ---
  ego <- run_go(entrez)
  df_go <- for_go <- NULL
  if (!inherits(ego, "try-error") && !is.null(ego)) {
    ## ★ ont="ALL" 时 enrichGO 自己带 ONTOLOGY 列；
    ##   ont="BP"/"CC"/"MF" 时**没有这一列**，必须自己补上 ——
    ##   否则 ONTOLOGY 全空，03 出图时 `dat[dat$ONTOLOGY == NA, ]` 会造出一行全 NA。
    for_go <- as_result_df(ego, if (ont == "ALL") NULL else ont)
    if (!is.null(for_go)) {
      df_go <- add_gene_name(mark_sig(order_cols(for_go), sig_p, sig_padj, sig_q), sym_map)
    }
  }
  n_go <- if (is.null(df_go)) 0L else sum(df_go$sig)
  note(sprintf("GO(%s)   ：被检验条目 %5d，其中显著 %4d 条",
               ont, if (is.null(df_go)) 0L else nrow(df_go), n_go))
  if (!is.null(df_go)) {
    write_csv_out(df_go, file.path(outdir, sprintf("%s_GO_%s.csv", prefix, nm)))
    if (export_sig) {
      s <- df_go[df_go$sig, , drop = FALSE]
      write_csv_out(s, file.path(outdir, sprintf("%s_GO_%s_sig.csv", prefix, nm)))
    }
  }

  ## --- KEGG ---
  df_kegg <- NULL
  kegg_src <- kegg_mode
  if (kegg_mode != "none") {
    ek <- NULL
    if (kegg_mode %in% c("auto", "online")) {
      ek <- run_kegg_online(entrez)
      if (inherits(ek, "try-error") || is.null(ek) || !nrow(as.data.frame(ek))) {
        if (kegg_mode == "online") {
          note("⚠ KEGG 在线富集失败：", if (inherits(ek, "try-error")) attr(ek, "condition")$message else "无结果")
        }
        ek <- NULL
      }
    }
    if (!is.null(ek)) {
      df_kegg <- add_gene_name(mark_sig(order_cols(as_result_df(ek, "KEGG")),
                                        sig_p, sig_padj, sig_q), sym_map)
      kegg_src <- "online"
    } else if (kegg_mode %in% c("auto", "gmt")) {
      kg <- ensure_kegg_gmt()
      if (!is.null(kg)) {
        kind <- gmt_id_kind(kg, genes_symbol, entrez)
        g_u <- if (kind == "symbol") genes_symbol else entrez
        cand <- if (kind == "entrez") entrez else genes_symbol
        ## 只保留 GMT 里出现过的，避免 enricher 报 "None of the genes are in the gene sets"
        g_u <- intersect(g_u, unique(kg$gene))
        if (length(g_u)) {
          e2 <- run_enricher(g_u, kg, NULL)
          if (!inherits(e2, "try-error") && !is.null(e2)) {
            d <- mark_sig(order_cols(as_result_df(e2, "KEGG")), sig_p, sig_padj, sig_q)
            d$collection <- "KEGG"
            df_kegg <- d
            kegg_src <- "gmt(MSigDB KEGG_LEGACY)"
            note("⚠ 在线 KEGG 不可用，已回退到本地 MSigDB 的 KEGG 子集。")
            note("  注意：MSigDB 的 KEGG 子集是 KEGG_LEGACY（186 条），",
                 "**不含 KEGG_PI3K_AKT_SIGNALING_PATHWAY、KEGG_TNF_SIGNALING_PATHWAY**；",
                 "通路 ID 也不是 hsaXXXXX 形式。")
          }
        }
      } else {
        note("⚠ KEGG 在线失败且本地没有可用的 KEGG GMT，本次跳过 KEGG。",
             "（可用 --gmt_dir 指向含 c2.all.v*.Hs.symbols.gmt 的目录）")
      }
    }
  }
  if (!is.null(df_kegg)) {
    n_k <- sum(df_kegg$sig)
    note(sprintf("KEGG[%s]：被检验条目 %5d，其中显著 %4d 条",
                 kegg_src, nrow(df_kegg), n_k))
    write_csv_out(df_kegg, file.path(outdir, sprintf("%s_KEGG_%s.csv", prefix, nm)))
    if (export_sig) {
      write_csv_out(df_kegg[df_kegg$sig, , drop = FALSE],
                    file.path(outdir, sprintf("%s_KEGG_%s_sig.csv", prefix, nm)))
    }
  } else {
    n_k <- 0L
    note("KEGG         ：跳过（--kegg=none 或不可用）")
  }

  ## --- 自定义 GMT ---
  df_gmt <- NULL
  if (!is.null(t2g)) {
    kind <- gmt_id_kind(t2g, genes_symbol, entrez)
    g_u <- if (kind == "symbol") genes_symbol else entrez
    g_u <- intersect(g_u, unique(t2g$gene))
    note(sprintf("GMT：输入基因按 %s 解释，其中 %d / %d 出现在集合里",
                 kind, length(g_u), length(genes_symbol)))
    if (length(g_u)) {
      e3 <- run_enricher(g_u, t2g, NULL)
      if (!inherits(e3, "try-error") && !is.null(e3)) {
        d <- mark_sig(order_cols(as_result_df(e3, NULL)), sig_p, sig_padj, sig_q)
        coll_map <- t2g[!duplicated(t2g$term), c("term", "collection")]
        d$collection <- coll_map$collection[match(d$ID, coll_map$term)]
        d$ONTOLOGY <- d$collection
        df_gmt <- d
        note(sprintf("MSigDB       ：被检验条目 %5d，其中显著 %4d 条", nrow(d), sum(d$sig)))
        write_csv_out(df_gmt, file.path(outdir, sprintf("%s_GMT_%s.csv", prefix, nm)))
        if (export_sig) {
          write_csv_out(df_gmt[df_gmt$sig, , drop = FALSE],
                        file.path(outdir, sprintf("%s_GMT_%s_sig.csv", prefix, nm)))
        }
      }
    }
  }

  all_res[[nm]] <- list(GO = df_go, KEGG = df_kegg, GMT = df_gmt,
                        entrez = entrez, symbols = genes_symbol, map = conv$map,
                        ego = if (!inherits(ego, "try-error")) ego else NULL)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    set = nm, n_input = conv$n_in, n_entrez = conv$n_out,
    GO_tested = if (is.null(df_go)) 0L else nrow(df_go),
    GO_sig = n_go,
    KEGG_source = kegg_src,
    KEGG_tested = if (is.null(df_kegg)) 0L else nrow(df_kegg),
    KEGG_sig = n_k,
    GMT_tested = if (is.null(df_gmt)) 0L else nrow(df_gmt),
    GMT_sig = if (is.null(df_gmt)) 0L else sum(df_gmt$sig),
    stringsAsFactors = FALSE)
}

## ---- 汇总 -----------------------------------------------------------------
sm <- do.call(rbind, summary_rows)
note("")
note(strrep("=", 72))
note("汇总（显著口径：pvalue<=", sig_p, " & p.adjust<=", sig_padj, " & qvalue<=", sig_q, "）")
note(strrep("=", 72))
print(sm, row.names = FALSE)
write_csv_out(sm, file.path(outdir, paste0(prefix, "_ORA_summary.csv")))
writeLines(capture.output(print(sm, row.names = FALSE)),
           file.path(outdir, paste0(prefix, "_ORA_summary.txt")), useBytes = TRUE)

if (save_rdata) {
  saveRDS(all_res, file.path(outdir, paste0(prefix, "_ORA_results.rds")))
  note("已保存：", file.path(outdir, paste0(prefix, "_ORA_results.rds")))
}

note("")
note("产出目录：", normalizePath(outdir, mustWork = FALSE))
note("下一步出美化图：")
note(sprintf("  bash run_enrich.sh 03_enrich_plot.R --in_dir=%s --prefix=%s --direction=%s",
             outdir, prefix, names(sets)[1]))
note("耗时：", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " 秒")
