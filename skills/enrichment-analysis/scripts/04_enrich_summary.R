#!/usr/bin/env Rscript
## ===========================================================================
## 04_enrich_summary.R —— 把全部显著通路「统计 + 归纳」成功能大类
##
## 解决的问题：01 步会把上百条 p<0.05 的通路一股脑倒出来，03 步默认只画每类
## 前 5 条。中间缺一层"这些通路到底在讲哪几件事"的归纳。
##
## 本脚本做四件事：
##   1. 取全部通过阈值的通路（默认 p.adjust <= 0.05，可用 --p_type=p 换成原始 p）
##   2. 按 references/functional-rules.tsv 的规则表，把每条通路归纳到一个功能大类
##      （规则表是可编辑的制表符文件，按顺序先匹配到的优先）
##   3. 在同一个功能大类内部，用基因重叠（Jaccard）去冗余，挑出"代表通路"
##   4. 按通路数量排序，报告前 N 个功能大类（默认 20），并把排名前 20 的大类
##      单独落成 CSV 表格（<prefix>_func_top20_<方向>.csv）供用户查阅
##
## ★ 本脚本**不出图**。归纳完把报告交给用户，等用户选：
##      A) 均衡配额出图（03 默认：固定 20 条，各分类均分）-> 03_enrich_plot.R --mode=even
##      B) 只画某个功能方向（该方向的代表通路）-> 03_enrich_plot.R --category=<名> --func_csv=<本步产物>
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
04_enrich_summary.R —— 全部显著通路的功能归纳（统计 + 归类 + 去冗余）

输入：
  --in_dir=./enrich           01_enrich_ora.R 的 --outdir
  --prefix=                   ORA 输出文件的前缀
  --direction=all,up,down     要归纳哪几个方向（可多值）
  --db=GO,KEGG,GMT            要纳入归纳的库

阈值：
  --p=0.05                    显著阈值
  --p_type=padj               padj(p.adjust) | p(pvalue) | qvalue | sig(sig 列)

归纳方法：
  --rules=<tsv>               功能分类规则表（默认用技能自带 references/functional-rules.tsv）
  --jaccard=0.5               同一大类内的基因重叠去冗余阈值（0 = 不去冗余）
  --top_cats=20               报告/落表前几个功能大类（默认 20，正好覆盖全部定义大类）
  --outdir=                   产物目录（默认 = --in_dir）

产出：
  <prefix>_func_pathways_<方向>.csv   逐条通路的归类结果（含 代表/冗余 标记）
  <prefix>_func_counts_<方向>.csv     功能大类汇总（通路数、代表数、最小 p、库构成）
  <prefix>_func_top20_<方向>.csv      排名前 20 的大类清单（供用户查阅；「未归类」不参与排名）
  <prefix>_func_counts_by_direction.csv  所有方向合并的汇总表
  <prefix>_func_report.txt            给人看的报告（含下一步两个选项）
  <prefix>_func_rules_used.tsv        本次实际生效的规则表（便于审计/复现）

示例：
  bash run_enrich.sh 04_enrich_summary.R --in_dir=./enrich --prefix=GSE62452_T \\
    --direction=all,up,down --p_type=p --p=0.05
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
directions <- oa_list(opt, "direction", c("all"))
dbs        <- toupper(oa_list(opt, "db", c("GO", "KEGG", "GMT")))
p_thr      <- oa_num(opt, "p", 0.05)
p_type     <- tolower(oa(opt, "p_type", "padj"))
jaccard    <- oa_num(opt, "jaccard", 0.5)
top_cats   <- oa_num(opt, "top_cats", 20)
rules_path <- oa(opt, "rules", NULL)
outdir     <- oa(opt, "outdir", in_dir)

dir_create(outdir)
start_log(file.path(outdir, paste0(prefix, "_04_summary.log")))
on.exit(close_log(), add = TRUE)

print_banner("04_enrich_summary.R", "显著通路的功能归纳")
print_args(opt)
note("")

## ---- 阈值列名 -------------------------------------------------------------
thr_col <- switch(p_type,
                  padj = "p.adjust", padjval = "p.adjust", p.adjust = "p.adjust", fdr = "p.adjust",
                  p = "pvalue", pval = "pvalue", pvalue = "pvalue",
                  q = "qvalue", qvalue = "qvalue",
                  sig = "sig",
                  NA_character_)
if (is.na(thr_col)) {
  stop("--p_type 只能是 padj / p / qvalue / sig，收到：", p_type, call. = FALSE)
}
note(sprintf("显著口径：%s <= %s%s", thr_col, p_thr,
             if (p_type == "padj") "（想按原始 p 就用 --p_type=p）" else ""))

## ---- 规则表 ---------------------------------------------------------------
rules <- load_func_rules(rules_path)
note(sprintf("功能分类规则表：%s（%d 条规则，%d 个功能大类，按顺序先匹配到的优先）",
             attr(rules, "path"), nrow(rules), length(unique(rules$category))))
write_csv_out(rules, file.path(outdir, paste0(prefix, "_func_rules_used.tsv")), row_names = FALSE)

## 中文是双宽字符，sprintf 的 %-22s 按"字符数"补空格会参差不齐，
## 所以自己按显示宽度对齐（nchar(type="width") 在 UTF-8 locale 下 CJK 记 2）。
dwidth <- function(s) {
  s <- as.character(s)
  w <- suppressWarnings(nchar(s, type = "width"))
  w[is.na(w) | w < 0] <- nchar(s)[is.na(w) | w < 0]
  w
}
pad <- function(s, w, right = FALSE) {
  s <- as.character(s)
  cur <- dwidth(s)
  n <- pmax(0L, w - cur)
  if (right) paste0(strrep(" ", n), s) else paste0(s, strrep(" ", n))
}
trunc_cn <- function(s, w) {
  s <- as.character(s)
  if (dwidth(s) <= w) return(s)
  out <- ""
  for (ch in strsplit(s, "")[[1]]) {
    if (dwidth(out) + dwidth(ch) > w - 1) break
    out <- paste0(out, ch)
  }
  paste0(out, "…")
}
tbl_int <- function(tb, k) {
  v <- tb[k]
  if (!length(v) || is.na(v)) 0L else as.integer(v)
}

## ---- 读表 -----------------------------------------------------------------
read_one <- function(db, d) {
  f <- file.path(in_dir, sprintf("%s_%s_%s.csv", prefix, db, d))
  if (!file.exists(f)) return(NULL)
  x <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "\"", comment.char = "")
  if (!nrow(x)) return(NULL)
  if (!"ONTOLOGY" %in% names(x)) x$ONTOLOGY <- db
  x$ONTOLOGY <- as.character(x$ONTOLOGY)
  if (db == "KEGG") x$ONTOLOGY <- "KEGG"
  if (db == "GMT" && "collection" %in% names(x)) {
    x$ONTOLOGY <- ifelse(is.na(x$collection) | !nzchar(x$collection), "GMT", x$collection)
  }
  bad <- is.na(x$ONTOLOGY) | !nzchar(x$ONTOLOGY)
  if (any(bad)) x$ONTOLOGY[bad] <- db
  x$db <- db
  x
}

keep_cols <- c("db", "ONTOLOGY", "ID", "Description", "pvalue", "p.adjust", "qvalue",
               "Count", "geneID", "geneName", "geneRatio", "category", "subcategory",
               "collection")
norm <- function(x) {
  for (cc in setdiff(keep_cols, names(x))) x[[cc]] <- NA
  x <- x[, keep_cols, drop = FALSE]
  x$pvalue   <- suppressWarnings(as.numeric(x$pvalue))
  x$p.adjust <- suppressWarnings(as.numeric(x$p.adjust))
  x$qvalue   <- suppressWarnings(as.numeric(x$qvalue))
  x$Count    <- suppressWarnings(as.numeric(x$Count))
  x$Description <- as.character(x$Description)
  x$ID <- as.character(x$ID)
  x$geneID <- as.character(x$geneID)
  x$geneName <- as.character(x$geneName)
  x
}

## ---- 逐方向归纳 -----------------------------------------------------------
all_counts <- list()
report_blocks <- list()

for (d in directions) {
  parts <- lapply(dbs, function(db) {
    x <- read_one(db, d)
    if (is.null(x)) { note("  [", d, "] ", db, "：没有文件，跳过"); return(NULL) }
    norm(x)
  })
  parts <- parts[!vapply(parts, is.null, TRUE)]
  if (!length(parts)) {
    note("⚠ 方向 ", d, " 一个库都没读到，跳过。")
    next
  }
  dat <- do.call(rbind, parts)
  n_all <- nrow(dat)

  ## 过滤
  if (thr_col == "sig") {
    if (!"sig" %in% names(dat)) {
      stop("表里没有 sig 列，无法用 --p_type=sig。", call. = FALSE)
    }
    v <- as.logical(dat$sig)
  } else {
    v <- dat[[thr_col]]
  }
  ok <- !is.na(v) & v <= p_thr
  sig <- dat[ok, , drop = FALSE]

  note("")
  note(strrep("-", 72))
  note(sprintf("方向 [%s]：被检验 %d 条 -> 通过 %s<=%s 的 %d 条", d, n_all, thr_col, p_thr, nrow(sig)))
  ## 顺手把"原始 p"口径的数量也报出来，方便对照
  if (thr_col == "p.adjust") {
    n_p <- sum(!is.na(dat$pvalue) & dat$pvalue <= p_thr)
    note(sprintf("        （同阈值下按原始 pvalue 算是 %d 条；--p_type=p 可切到该口径）", n_p))
  }
  note(strrep("-", 72))
  if (!nrow(sig)) {
    note("没有通过阈值的通路，跳过该方向。")
    next
  }

  ## 归类：Description 优先，ID 兜底（有些集合的 Description 就是 ID）
  txt <- paste(sig$Description, sig$ID, ifelse(is.na(sig$subcategory), "", sig$subcategory),
               ifelse(is.na(sig$collection), "", sig$collection), sep = " | ")
  cl <- classify_func(txt, rules)
  sig <- cbind(sig, cl)

  ## 大类内部按基因重叠去冗余
  sig$representative <- FALSE
  for (k in unique(sig$func_cat)) {
    idx <- which(sig$func_cat == k)
    if (length(idx) == 1L) { sig$representative[idx] <- TRUE; next }
    if (jaccard <= 0) { sig$representative[idx] <- TRUE; next }
    sig$representative[idx] <- dedup_by_jaccard(sig$geneID[idx], sig$p.adjust[idx], jaccard)
  }

  ## 逐条结果
  out1 <- sig[order(sig$p.adjust, -sig$Count), , drop = FALSE]
  write_csv_out(out1, file.path(outdir, sprintf("%s_func_pathways_%s.csv", prefix, d)))

  ## 大类汇总
  cnt <- do.call(rbind, lapply(unique(sig$func_cat), function(k) {
    s <- sig[sig$func_cat == k, , drop = FALSE]
    rep_s <- s[s$representative, , drop = FALSE]
    rep_s <- rep_s[order(rep_s$p.adjust), , drop = FALSE]
    dbc <- table(s$db)
    o <- order(s$p.adjust)
    data.frame(
      func_cat = k,
      n_pathways = nrow(s),
      n_representatives = sum(s$representative),
      min_padj = suppressWarnings(min(s$p.adjust, na.rm = TRUE)),
      min_padj_rep = if (nrow(rep_s)) suppressWarnings(min(rep_s$p.adjust, na.rm = TRUE)) else NA_real_,
      n_GO = tbl_int(dbc, "GO"),
      n_KEGG = tbl_int(dbc, "KEGG"),
      n_GMT = tbl_int(dbc, "GMT"),
      top_pathways = paste(head(sprintf("%s(%s)", s$Description[o], s$db[o]), 3), collapse = "; "),
      stringsAsFactors = FALSE)
  }))
  cnt <- cnt[order(-cnt$n_pathways, cnt$min_padj), , drop = FALSE]
  cnt$direction <- d
  cnt <- cnt[, c("direction", setdiff(names(cnt), "direction"))]
  write_csv_out(cnt, file.path(outdir, sprintf("%s_func_counts_%s.csv", prefix, d)))
  all_counts[[d]] <- cnt

  ## ★ 排名前 N 的大类清单（默认 20，单独落成 CSV 供用户查阅）。
  ## 「未归类」不是功能方向，不参与排名，也不进这张表。
  cnt_rank0 <- cnt[cnt$func_cat != "未归类", , drop = FALSE]
  top_tab <- head(cnt_rank0, as.integer(top_cats))
  top_tab$rank <- seq_len(nrow(top_tab))
  top_tab <- top_tab[, c("rank", setdiff(names(top_tab), "rank")), drop = FALSE]
  top_f <- file.path(outdir, sprintf("%s_func_top%d_%s.csv", prefix, as.integer(top_cats), d))
  write_csv_out(top_tab, top_f)
  note(sprintf("已写出前 %d 大类清单（供查阅）：", as.integer(top_cats)), top_f)

  ## 报告块
  ## ★「未归类」不是一个功能方向，不能占排名位次 —— 从排名里摘出去单列。
  cnt_rank <- cnt[cnt$func_cat != "未归类", , drop = FALSE]
  unassigned <- cnt[cnt$func_cat == "未归类", , drop = FALSE]
  show <- head(cnt_rank, top_cats)
  lines <- character(0)
  add <- function(...) lines <<- c(lines, paste0(...))
  add("")
  add("================================================================================")
  add(sprintf("【方向 %s】显著通路 %d 条（%s）-> 去冗余后代表通路 %d 条 -> 功能大类 %d 个",
              d, nrow(sig),
              paste(sprintf("%s %d", names(table(sig$db)), as.integer(table(sig$db))), collapse = " / "),
              sum(sig$representative), nrow(cnt_rank)))
  add("================================================================================")
  add(sprintf("按通路数量排名前 %d 个功能方向（「未归类」不参与排名，单列在下方）：", nrow(show)))
  add("")
  add(paste0(pad("排名", 6), pad("功能大类", 26), pad("通路数", 8, TRUE),
             pad("代表数", 8, TRUE), pad("最小p.adjust", 15, TRUE), "  库构成"))
  add(strrep("-", 92))
  for (i in seq_len(nrow(show))) {
    r <- show[i, ]
    comp <- paste(c(if (r$n_GO > 0) sprintf("GO:%d", r$n_GO) else NULL,
                    if (r$n_KEGG > 0) sprintf("KEGG:%d", r$n_KEGG) else NULL,
                    if (r$n_GMT > 0) sprintf("GMT:%d", r$n_GMT) else NULL),
                  collapse = " ")
    add(paste0(pad(i, 6), pad(r$func_cat, 26), pad(r$n_pathways, 8, TRUE),
               pad(r$n_representatives, 8, TRUE),
               pad(formatC(r$min_padj, format = "e", digits = 1), 15, TRUE),
               "  ", comp))
  }
  if (nrow(cnt_rank) > nrow(show)) {
    add(sprintf("      （其后还有 %d 个功能方向，完整清单见 *_func_counts_*.csv）",
                nrow(cnt_rank) - nrow(show)))
  }
  add("")
  add("每个功能方向的代表通路（按显著性排序前 5 条，完整清单见 *_func_pathways_*.csv）：")
  if (!nrow(show)) add("  （没有任何通路落进功能方向，全部未归类）")
  for (k in show$func_cat) {
    s <- sig[sig$func_cat == k & sig$representative, , drop = FALSE]
    s <- s[order(s$p.adjust), , drop = FALSE]
    add(sprintf("  ● %s（%d 条通路 / %d 条代表）", k,
                show$n_pathways[show$func_cat == k],
                show$n_representatives[show$func_cat == k]))
    for (j in seq_len(min(5, nrow(s)))) {
      add(sprintf("      [%s|%s] %-14s %s   p.adj=%s  n=%d  小类:%s",
                  s$db[j], s$ONTOLOGY[j], trunc_cn(s$ID[j], 16), s$Description[j],
                  formatC(s$p.adjust[j], format = "e", digits = 2), as.integer(s$Count[j]),
                  s$func_sub[j]))
    }
    extra <- sig[sig$func_cat == k & !sig$representative, , drop = FALSE]
    if (nrow(extra)) {
      add(sprintf("      （另 %d 条与上面高度重叠，已折叠为冗余）", nrow(extra)))
    }
  }
  ## 未归类：单列，并列出它的前几条，方便判断要不要补规则
  if (nrow(unassigned)) {
    u <- sig[sig$func_cat == "未归类", , drop = FALSE]
    u <- u[order(u$p.adjust), , drop = FALSE]
    add("")
    add(sprintf("【未归类】%d 条通路没被任何规则命中（不参与上面的排名）：", nrow(u)))
    for (j in seq_len(min(8, nrow(u)))) {
      add(sprintf("      [%s|%s] %-14s %s   p.adj=%s",
                  u$db[j], u$ONTOLOGY[j], trunc_cn(u$ID[j], 16), u$Description[j],
                  formatC(u$p.adjust[j], format = "e", digits = 2)))
    }
    if (nrow(u) > 8) add(sprintf("      …还有 %d 条，见 *_func_pathways_*.csv", nrow(u) - 8))
    add("      想让它们归位就编辑 references/functional-rules.tsv 加规则，或用 --rules= 换一份。")
  }
  report_blocks[[d]] <- lines
  note(paste(lines, collapse = "\n"))
}

if (!length(all_counts)) stop("所有方向都没有可归纳的通路。检查 --prefix / --direction / 阈值。", call. = FALSE)

## ---- 汇总 -----------------------------------------------------------------
## ★ 以前这里写成 <prefix>_func_counts_all.csv，与「方向=all」的逐方向文件同名，
##   后者会被这张合并表覆盖掉。改成 _by_direction，两边都留得住。
cnt_all <- do.call(rbind, all_counts)
write_csv_out(cnt_all, file.path(outdir, paste0(prefix, "_func_counts_by_direction.csv")))

## ---- 总体报告 -------------------------------------------------------------
foot <- c(
  "",
  "================================================================================",
  "下一步要你定（本脚本不出图）",
  "================================================================================",
  "",
  "  A) 均衡配额出图（03 的默认模式）：整张图固定 20 条，BP/CC/MF/KEGG 各 5 条；",
  "     某个分类不足配额时，富余名额自动匀给还有存货的分类",
  sprintf("     Rscript 03_enrich_plot.R --in_dir=%s --prefix=%s --direction=%s --db=%s --mode=even",
          outdir, prefix, directions[1], paste(dbs, collapse = ",")),
  "     （想改回按排名挑法：--mode=global --top=10 或 --mode=per_ontology --top=5）",
  "",
  "  B) 只看某一个功能方向：只画上面某一个大类的通路（默认只画代表通路，不重复）",
  sprintf("     Rscript 03_enrich_plot.R --in_dir=%s --prefix=%s --direction=%s --db=%s \\",
          outdir, prefix, directions[1], paste(dbs, collapse = ",")),
  sprintf("       --func_csv=%s --category=<功能大类名> --top=0",
          file.path(outdir, sprintf("%s_func_pathways_%s.csv", prefix, directions[1]))),
  "",
  "  说明：--top=0 表示不限制条数；配合默认的 --use_representative=TRUE，",
  "        同一个功能方向里高度重叠的通路不会重复画。",
  "        想看该方向全部通路（含冗余）就加 --use_representative=FALSE。",
  "        A/B 两条命令里的 --direction 可换成 up / down；",
  "        B 的 --func_csv 也要换成对应方向的 _func_pathways_<方向>.csv。",
  "",
  "★ 归纳完成后请把本报告交给用户，等他选 A 还是 B、选哪个功能大类，再出图。",
  "================================================================================",
  "")

out_txt <- file.path(outdir, paste0(prefix, "_func_report.txt"))
writeLines(c(paste(report_blocks[[1]]), do.call(c, unname(report_blocks[-1])), foot),
           out_txt, useBytes = TRUE)
note("")
note("已写出报告：", out_txt)
note("耗时：", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " 秒")
