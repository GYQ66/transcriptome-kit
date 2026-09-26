#!/usr/bin/env Rscript
# ===========================================================================
# 02_km_batch.R —— 多基因批量 KM 生存分析（每基因独立走两级策略）
#
# 逐基因执行与 01_km_survival.R 完全相同的两级策略：
#   常规分位值分组 -> 若 log-rank p >= --p_threshold 则自动改用最佳截断值
# 输出：每基因一套 PDF/PNG + 分组表，外加一张汇总表与一份多页 PDF。
#
# 用法:
#   bash run_surv.sh 02_km_batch.R --expr=TCGA_CRC_TPM.txt --clin=time.csv \
#        --genes=TIMP1,GSTP1,CDC25C --outdir=./km --prefix=CRC
#   bash run_surv.sh 02_km_batch.R --table=re.csv --time_col=time --event_col=state \
#        --genes_file=merge2.csv --outdir=./km --prefix=CRCscore
# ===========================================================================

# 定位 lib_surv_common.R（脚本可能被 run_surv.sh 复制到 ASCII 临时目录）
find_surv_lib <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  cand <- character(0)
  if (length(f)) cand <- c(cand, dirname(normalizePath(f[1], mustWork = FALSE)))
  if (nzchar(Sys.getenv("SURV_SKILL_DIR"))) cand <- c(cand, Sys.getenv("SURV_SKILL_DIR"))
  cand <- c(cand, getwd(), file.path(getwd(), "scripts"))
  for (d in cand) {
    if (nzchar(d) && file.exists(file.path(d, "lib_surv_common.R"))) return(d)
  }
  stop("找不到 lib_surv_common.R。请用 run_surv.sh 启动脚本，或在 scripts/ 目录下运行。",
       call. = FALSE)
}
SELF_DIR <- find_surv_lib()
source(file.path(SELF_DIR, "lib_surv_common.R"))

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (!is.null(args[["help"]]) || !is.null(args[["h"]])) {
  cat("
02_km_batch.R —— 多基因批量 KM（每基因独立走 常规分组 -> 最佳阈值回退）

输入（二选一）
  A  --expr=<表达矩阵> --clin=<临床表> --genes=A,B,C  或  --genes_file=<每行一个基因>
  B  --table=<单表>    --genes=A,B,C（=列名）        或  --genes_file=

级联传入 01 的全部参数: --time_col --event_col --time_in --time_out --time_div
  --mode --pct --p_threshold --minprop --maxstat_pmethod --style --palette
  --conf_int --risk_table --median_line --break_time --pval_size --legend_pos
  --tumor_only --tumor_codes --id_norm --min_time --event_levels --width --height --dpi

  --style=unified|classic|threshold  默认 unified。整批共用一个风格，
     所以【走了最佳阈值回退】和【没走回退】的基因出的图版式完全一致。

批量专属
  --outdir=.        产物目录（每基因图放 <outdir>/per_gene/）
  --prefix=batch    产物前缀
  --multi_pdf=TRUE  是否额外输出把所有基因拼在一起的多页 PDF
  --png=TRUE        是否输出 PNG
  --skip_failed=TRUE  某基因失败时跳过继续（FALSE 则整批中止）
")
  quit(save = "no", status = 0)
}

# ------------------------------- 参数 -------------------------------------
style      <- arg_chr(args, "style", "unified")
if (!style %in% STYLE_CHOICES) {
  stop("--style 只能是 ", paste(STYLE_CHOICES, collapse = " / "), "，收到: ", style,
       call. = FALSE)
}
mode_req   <- arg_chr(args, "mode", "auto")
pct        <- arg_num(args, "pct", 0.5)
p_th       <- arg_num(args, "p_threshold", 0.05)
minprop    <- arg_num(args, "minprop", 0.3)
pmethod    <- arg_chr(args, "maxstat_pmethod", "Lau92")
min_time   <- arg_num(args, "min_time", 0)
tumor_only <- arg_bool(args, "tumor_only", TRUE)
id_norm    <- arg_chr(args, "id_norm", "auto")
outdir     <- arg_chr(args, "outdir", ".")
prefix     <- arg_chr(args, "prefix", "batch")
ensure_dir(outdir)
per_dir    <- file.path(outdir, "per_gene")
ensure_dir(per_dir)

def_conf   <- style_default_conf_int(style)
conf_int   <- arg_bool(args, "conf_int", def_conf)
risk_table <- arg_bool(args, "risk_table", TRUE)
hr_on_plot <- arg_bool(args, "hr_on_plot", TRUE)
median_line<- arg_chr(args, "median_line", "hv")
pval_size  <- arg_num(args, "pval_size", 6)
width      <- arg_num(args, "width", 8)
height     <- arg_num(args, "height", 6.25)
dpi        <- arg_num(args, "dpi", 300)
risk_h     <- arg_num(args, "risk_table_height", 0.25)
multi_pdf  <- arg_bool(args, "multi_pdf", TRUE)
do_png     <- arg_bool(args, "png", TRUE)
skip_failed<- arg_bool(args, "skip_failed", TRUE)
time_out   <- arg_chr(args, "time_out", "year")
time_in    <- arg_chr(args, "time_in", "auto")
time_div   <- arg_num(args, "time_div", NA_real_)
if (!time_out %in% names(UNIT_DAYS))
  stop("--time_out 只能是 day / month / year，收到: ", time_out, call. = FALSE)
if (!identical(time_in, "auto") && !time_in %in% names(UNIT_DAYS))
  stop("--time_in 只能是 auto / day / month / year，收到: ", time_in, call. = FALSE)

pal <- arg_vec(args, "palette"); if (length(pal) < 2) pal <- NULL
lp  <- arg_vec(args, "legend_pos")
legend_pos <- if (length(lp) == 2) as.numeric(lp) else NULL
break_by <- arg_num(args, "break_time", NA_real_); if (is.na(break_by)) break_by <- NULL

cat("==================== 02_km_batch ====================\n")

# ------------------------- 1. 基因列表 -------------------------------
genes <- arg_vec(args, "genes")
gf <- arg_chr(args, "genes_file", NULL)
if (!is.null(gf)) genes <- unique(c(genes, read_gene_list(gf)))
if (!length(genes)) stop("请用 --genes=A,B,C 或 --genes_file= 指定基因列表", call. = FALSE)
cat(sprintf("[genes] 共 %d 个: %s%s\n", length(genes),
            paste(head(genes, 10), collapse = ", "),
            if (length(genes) > 10) " ..." else ""))

# ------------------------- 2. 数据准备（只做一次）-------------------------
expr_path <- arg_chr(args, "expr"); clin_path <- arg_chr(args, "clin")
tab_path  <- arg_chr(args, "table")
time_col_in <- arg_chr(args, "time_col"); event_col_in <- arg_chr(args, "event_col")
ev_levels <- arg_vec(args, "event_levels")

value_list <- list()   # 每基因的取值向量；统一按 samples 对齐

if (!is.null(tab_path)) {
  tb <- read_table_auto(tab_path, row_names = 1)
  cat(sprintf("[input] 单表模式：%d 行 × %d 列\n", nrow(tb), ncol(tb)))
  tcol <- if (!is.null(time_col_in)) time_col_in else
    intersect(c("OS.time", "time", "OS_time", "survival_time", "time_days"), colnames(tb))[1]
  ecol <- if (!is.null(event_col_in)) event_col_in else
    intersect(c("OS", "state", "status", "event", "vital_status"), colnames(tb))[1]
  if (is.null(tcol) || is.na(tcol)) stop("未探测到时间列，请用 --time_col= 指定", call. = FALSE)
  if (is.null(ecol) || is.na(ecol)) stop("未探测到状态列，请用 --event_col= 指定", call. = FALSE)
  raw_time <- pick_numeric(tb, tcol, "时间列"); raw_event <- pick_numeric(tb, ecol, "状态列")
  samples <- rownames(tb)
  miss <- setdiff(genes, colnames(tb))
  if (length(miss)) {
    msg <- sprintf("单表中找不到这些列: %s", paste(miss, collapse = ", "))
    if (skip_failed) { warning(msg, call. = FALSE); genes <- setdiff(genes, miss) }
    else stop(msg, call. = FALSE)
  }
  for (g in genes) value_list[[g]] <- suppressWarnings(as.numeric(tb[[g]]))
} else {
  if (is.null(expr_path) || is.null(clin_path))
    stop("路径 A 需要 --expr= 与 --clin=", call. = FALSE)
  expr <- read_table_auto(expr_path, row_names = 1)
  cat(sprintf("[input] 表达矩阵：%d 行 × %d 列\n", nrow(expr), ncol(expr)))
  ft <- filter_tumor(colnames(expr), keep_codes = arg_vec(args, "tumor_codes", "0"),
                     enabled = tumor_only)
  if (ft$applied) {
    cat(sprintf("[filter] 保留 sample type ∈ {%s}：%d -> %d 列\n",
                paste(arg_vec(args, "tumor_codes", "0"), collapse = ","),
                ncol(expr), sum(ft$keep)))
    expr <- expr[, ft$keep, drop = FALSE]
  }
  e_ids <- normalize_ids(colnames(expr), id_norm)
  cli <- read_clin_auto(clin_path, id_col = arg_chr(args, "clin_id_col", NULL))
  c_ids <- normalize_ids(rownames(cli), id_norm)
  tcol <- if (!is.null(time_col_in)) time_col_in else
    intersect(c("OS.time", "time", "OS_time", "survival_time", "time_days"), colnames(cli))[1]
  ecol <- if (!is.null(event_col_in)) event_col_in else
    intersect(c("OS", "state", "status", "event"), colnames(cli))[1]
  if (is.null(tcol) || is.na(tcol)) stop("临床表未找到时间列，请用 --time_col= 指定", call. = FALSE)
  if (is.null(ecol) || is.na(ecol)) stop("临床表未找到状态列，请用 --event_col= 指定", call. = FALSE)
  common <- intersect(e_ids, c_ids)
  cat(sprintf("[merge] 交集样本 %d\n", length(common)))
  if (length(common) < 10) stop("交集样本不足（<10）", call. = FALSE)
  ei <- match(common, e_ids); ci <- match(common, c_ids)
  cli_use <- cli[ci, , drop = FALSE]
  raw_time <- pick_numeric(cli_use, tcol, "时间列"); raw_event <- pick_numeric(cli_use, ecol, "状态列")
  samples <- common
  miss <- setdiff(genes, rownames(expr))
  if (length(miss)) {
    msg <- sprintf("表达矩阵中找不到这些基因: %s", paste(head(miss, 20), collapse = ", "))
    if (skip_failed) { warning(msg, call. = FALSE); genes <- setdiff(genes, miss) }
    else stop(msg, call. = FALSE)
  }
  for (g in genes) value_list[[g]] <- suppressWarnings(as.numeric(expr[g, ei]))
}
if (!length(genes)) stop("没有可分析的基因", call. = FALSE)

# ------------------------- 3. 事件 / 时间 -------------------------------
cat(sprintf("[event] 状态列 '%s' 取值分布:\n", ecol)); print(table(raw_event, useNA = "ifany"))
if (length(ev_levels)) {
  raw_event <- as.integer(raw_event %in% ev_levels)
  cat(sprintf("[event] 按 --event_levels={%s} 重编码\n", paste(ev_levels, collapse = ",")))
} else if (!all(raw_event[!is.na(raw_event)] %in% c(0, 1))) {
  warning("状态列取值不是 0/1，已原样使用，如为 1/2 编码请加 --event_levels=2", call. = FALSE)
}

det <- detect_time_unit(raw_time)
if (identical(time_in, "auto")) {
  time_in <- det$unit
  cat(sprintf("[time] 自动推断输入单位 = %s（中位数 %.3g，最大 %.3g）-- 请核对，必要时用 --time_in= 覆盖\n",
              time_in, det$median, det$max))
}
ana_time <- if (!is.na(time_div)) raw_time / time_div else convert_time(raw_time, time_in, time_out)
cat(sprintf("[time] 分析用单位 = %s（中位数 %.3g）\n", unit_label_cn[[time_out]],
            median(ana_time, na.rm = TRUE)))

keep_base <- is.finite(ana_time) & ana_time > min_time & !is.na(raw_event)
cat(sprintf("[clean] 时间/状态有效样本 %d / %d\n", sum(keep_base), length(keep_base)))

# ------------------------- 4. 逐基因分析 -------------------------------
rows <- list(); plots <- list(); failed <- character(0)
for (g in genes) {
  cat(sprintf("\n---------------- [%s] ----------------\n", g))
  v <- value_list[[g]]
  keep <- keep_base & is.finite(v)
  if (sum(keep) < 10) {
    cat(sprintf("[skip] 有效样本仅 %d 例\n", sum(keep)))
    failed <- c(failed, g); next
  }
  t_i <- ana_time[keep]; e_i <- raw_event[keep]; v_i <- v[keep]; s_i <- samples[keep]

  res <- tryCatch(
    two_tier_km(t_i, e_i, v_i, pct = pct, mode = mode_req, p_th = p_th,
                minprop = minprop, pmethod = pmethod, verbose = TRUE),
    error = function(e) e)
  if (inherits(res, "error")) {
    cat(sprintf("[error] %s\n", conditionMessage(res)))
    failed <- c(failed, g)
    if (!skip_failed) stop("基因 ", g, " 分析失败: ", conditionMessage(res), call. = FALSE)
    next
  }
  st <- res$stats
  cat(sprintf("[%s] 采用=%s 阈值=%.6g %s %s %s\n", g, res$mode_used, res$cutoff,
              fmt_p(st$pvalue), fmt_hr(st$hr, st$hr_low, st$hr_hi),
              ifelse(res$significant, "显著", "未达显著")))

  d <- data.frame(time = t_i, event = e_i,
                  group = factor(res$group, levels = c("Low", "High")), value = v_i)
  fit <- survfit(Surv(time, event) ~ group, data = d)
  p_line <- sprintf("%s\n%s", fmt_p(st$pvalue), fmt_hr(st$hr, st$hr_low, st$hr_hi))
  if (identical(style, "threshold"))
    p_line <- sprintf("%s\n%s", fmt_hr(st$hr, st$hr_low, st$hr_hi), fmt_p(st$pvalue))
  ascii_guard(c(g, p_line), "图内文字")

  pl <- make_km_plot(fit, d, g, style = style, unit_out = time_out,
                     break_by = break_by, palette = pal, p_label = p_line,
                     hr_on_plot = hr_on_plot, risk_table = risk_table,
                     conf_int = conf_int, median_line = median_line,
                     legend_pos = legend_pos, pval_size = pval_size, risk_height = risk_h)
  gsafe <- gsub("[^A-Za-z0-9_.-]", "_", g)
  pdf_i <- file.path(per_dir, paste0(gsafe, "_KM.pdf"))
  png_i <- if (do_png) file.path(per_dir, paste0(gsafe, "_KM.png")) else NULL
  tryCatch(save_km(pl, pdf_i, png_i, width = width, height = height, dpi = dpi),
           error = function(e) { cat(sprintf("[plot error] %s\n", conditionMessage(e)));
                                 failed <<- c(failed, g) })
  plots[[g]] <- pl
  write.csv(data.frame(sample = s_i, value = v_i, group = as.character(res$group),
                       time = t_i, event = e_i),
            file.path(per_dir, paste0(gsafe, "_group_map.csv")),
            row.names = FALSE, fileEncoding = "UTF-8")

  ms <- st$median_surv
  rows[[g]] <- data.frame(
    gene = g, mode_used = res$mode_used, mode_requested = mode_req,
    cutoff = res$cutoff, pct = res$pct, minprop = res$minprop,
    p_used = st$pvalue, p_logrank_quantile = res$stats_quantile$pvalue,
    p_logrank_cutpoint = if (is.null(res$stats_cutpoint)) NA else res$stats_cutpoint$pvalue,
    p_maxstat_adjusted = res$maxstat_p_adjusted,
    significant = res$significant, fallback_used = res$fallback_used,
    hr_high_vs_low = st$hr, hr_ci_low = st$hr_low, hr_ci_high = st$hr_hi,
    p_cox = if (is.null(st$cox)) NA else st$cox_p,
    p_text = fmt_p(st$pvalue), hr_text = fmt_hr(st$hr, st$hr_low, st$hr_hi),
    n_total = st$n, n_event = sum(e_i == 1),
    n_low = sum(res$group == "Low"), n_high = sum(res$group == "High"),
    median_surv_low = if (is.null(ms)) NA else unname(ms["Low"]),
    median_surv_high = if (is.null(ms)) NA else unname(ms["High"]),
    time_unit_in = time_in, time_unit_out = time_out,
    cutpoint_error = res$cutpoint_msg,
    figure_pdf = pdf_i, figure_png = ifelse(is.null(png_i), "", png_i),
    stringsAsFactors = FALSE)
}

if (!length(rows)) stop("没有任何基因分析成功", call. = FALSE)
summary_df <- do.call(rbind, rows)
rownames(summary_df) <- NULL
# 显著者优先、再按 p 升序，便于直接看结果
summary_df <- summary_df[order(!summary_df$significant, summary_df$p_used), ]
sum_path <- file.path(outdir, paste0(prefix, "_km_summary.csv"))
write.csv(summary_df, sum_path, row.names = FALSE, fileEncoding = "UTF-8")

# 多页 PDF
multi_path <- NULL
if (multi_pdf && length(plots)) {
  multi_path <- file.path(outdir, paste0(prefix, "_km_multi.pdf"))
  grDevices::pdf(multi_path, width = width, height = height, onefile = TRUE)
  for (g in names(plots)) { cat(sprintf("[multi] page: %s\n", g)); print_km(plots[[g]]) }
  grDevices::dev.off()
}

cat("\n==================== 汇总 ====================\n")
print(summary_df[, c("gene", "mode_used", "cutoff", "p_used", "hr_high_vs_low",
                     "significant", "fallback_used", "n_low", "n_high")], row.names = FALSE)
cat(sprintf("\n共 %d 个基因：显著 %d，未达显著 %d，失败 %d\n",
            nrow(summary_df), sum(summary_df$significant),
            sum(!summary_df$significant), length(failed)))
if (sum(summary_df$fallback_used))
  cat(sprintf("其中 %d 个基因因常规分组不显著而改用最佳截断值（探索性，建议独立队列验证）\n",
              sum(summary_df$fallback_used)))
if (length(failed)) cat(sprintf("失败/跳过: %s\n", paste(failed, collapse = ", ")))
cat(sprintf("\n[out] 汇总: %s\n", sum_path))
if (!is.null(multi_path)) cat(sprintf("[out] 多页图: %s\n", multi_path))
cat(sprintf("[out] 单基因产物目录: %s\n", per_dir))
cat("=============================================\n")
