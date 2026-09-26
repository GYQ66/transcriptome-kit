#!/usr/bin/env Rscript
# ===========================================================================
# 01_km_survival.R —— 单基因 / 风险评分 Kaplan-Meier 生存分析
#
# 两级策略（--mode=auto，默认）：
#   第一级  常规分位值分组（默认中位值，等价于 生存分析-km.R）
#   第二级  当第一级 log-rank p >= --p_threshold(默认 0.05) 时，
#           自动改用 survminer::surv_cutpoint 的**最佳截断值**重新分组
#           （等价于 最佳阈值-km.R）
#   取两者中 p 更小的分组方案作为最终输出，并把两个 p 都写进决策文件。
#
# 两种输入形态（二选一）：
#   A) 表达矩阵 + 临床表：  --expr=TPM.txt --clin=time.csv --gene=TIMP1
#   B) 单表（含时间/状态/评分）：--table=H19.csv --score_col=risk_score \
#                                --time_col=OS.time --event_col=OS
#
# 用法: bash run_surv.sh 01_km_survival.R --expr=... --clin=... --gene=...
# 完整参数: bash run_surv.sh 01_km_survival.R --help
# ===========================================================================

# 定位 lib_surv_common.R：脚本自身目录 -> run_surv.sh 导出的 SURV_SKILL_DIR
# -> 当前目录/scripts。之所以不能只看脚本自身目录：run_surv.sh 在脚本路径含
# 非 ASCII 时会把它复制到 ASCII 临时目录再运行，那时同目录下没有 lib。
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
01_km_survival.R —— 单基因/风险评分 KM 生存分析（常规分组 → 最佳阈值回退）

输入（二选一）
  A  --expr=<表达矩阵>  --clin=<临床表>  --gene=<基因/行名>
  B  --table=<单表>     --score_col=<评分列名>

必需/常用
  --time_col=<列名>        临床表时间列，默认自动探测(time/OS.time/OS_time/...)
  --event_col=<列名>       临床表状态列，默认自动探测(state/OS/status/event/...)
  --time_in=auto|day|month|year   输入时间单位，默认 auto（**推断结果务必核对**）
  --time_out=year|month|day       分析与画图用的时间单位，默认 year
  --time_div=<数值>        直接 t/time_div（复刻参考代码的除法），给了就覆盖 time_in/out
  --mode=auto|median|quantile|cutpoint   分组策略，默认 auto
  --pct=0.5                分位切分点（median 等价于 pct=0.5；三分位写 0.3333）
  --p_threshold=0.05       低于此值即认为常规分组已够，不再回退
  --minprop=0.3            surv_cutpoint 每组最少样本比例
  --maxstat_pmethod=Lau92  算 maxstat 校正 p 用的近似法（Lau94 可能越界，不推荐）
  --style=unified|classic|threshold  绘图风格，默认 unified
                            unified   统一版式，两种分组路径出的图完全一致（推荐）
                            classic   严格复刻 生存分析-km.R
                            threshold 严格复刻 最佳阈值-km.R
                            **风格只影响外观，不影响统计；也绝不随是否回退而改变**
  --palette=色1,色2        覆盖配色（Low,High）
  --hr_on_plot=TRUE        统一版式是否在图内加 HR 行
  --conf_int=FALSE         是否画置信带（threshold 风格默认 TRUE）
  --risk_table=TRUE        是否画风险表
  --median_line=hv|h|v|none  中位生存线
  --break_time=<数值>      横轴刻度间隔，默认年->2 / 月->24 / 天->365
  --pval_size=6  --legend_pos=0.8,0.8
  --width=8 --height=6.25 --dpi=300
  --outdir=.  --prefix=<前缀，默认取 gene/score_col>
  --tumor_only=TRUE --tumor_codes=0 --id_norm=auto|tcga|none
  --min_time=0             剔除时间 <= min_time 的样本
  --event_levels=           事件取值（默认按 0/1 原样使用；如 1,2 编码时写 2 或 Dead）
  --pdf=TRUE --png=TRUE
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
outdir     <- arg_chr(args, "outdir", ".")
ensure_dir(outdir)
do_pdf     <- arg_bool(args, "pdf", TRUE)
do_png     <- arg_bool(args, "png", TRUE)

# 风格默认值（可被显式参数覆盖）
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

pal <- arg_vec(args, "palette")
pal <- if (length(pal) >= 2) pal else NULL
lp  <- arg_vec(args, "legend_pos")
legend_pos <- if (length(lp) == 2) as.numeric(lp) else NULL

cat("==================== 01_km_survival ====================\n")

# ------------------------- 1. 取数：表达矩阵 + 临床表 ----------------------
expr_path <- arg_chr(args, "expr")
clin_path <- arg_chr(args, "clin")
tab_path  <- arg_chr(args, "table")
score_col <- arg_chr(args, "score_col")
gene      <- arg_chr(args, "gene")

time_col_in  <- arg_chr(args, "time_col")
event_col_in <- arg_chr(args, "event_col")

if (!is.null(tab_path)) {
  # ---------------- 路径 B：单表 ----------------
  if (is.null(score_col)) stop("路径 B 需要 --score_col=<列名>", call. = FALSE)
  tb <- read_table_auto(tab_path, row_names = 1)
  cat(sprintf("[input] 单表模式：%d 行 × %d 列 (列名: %s)\n",
              nrow(tb), ncol(tb), paste(head(colnames(tb), 12), collapse = ", ")))
  tcol <- if (!is.null(time_col_in)) time_col_in else
    intersect(c("OS.time", "time", "OS_time", "survival_time", "time_days"), colnames(tb))[1]
  ecol <- if (!is.null(event_col_in)) event_col_in else
    intersect(c("OS", "state", "status", "event", "vital_status"), colnames(tb))[1]
  if (is.null(tcol) || is.na(tcol)) stop("未能探测到时间列，请用 --time_col= 指定", call. = FALSE)
  if (is.null(ecol) || is.na(ecol)) stop("未能探测到状态列，请用 --event_col= 指定", call. = FALSE)
  raw_time  <- pick_numeric(tb, tcol, "时间列")
  raw_event <- pick_numeric(tb, ecol, "状态列")
  value     <- pick_numeric(tb, score_col, "评分列")
  samples   <- rownames(tb)
  gene_label <- if (!is.null(gene)) gene else score_col
} else {
  # ---------------- 路径 A：表达矩阵 + 临床表 ----------------
  if (is.null(expr_path)) stop("需要 --expr= 或 --table=", call. = FALSE)
  if (is.null(clin_path)) stop("路径 A 需要 --clin=<临床表>", call. = FALSE)
  if (is.null(gene)) stop("路径 A 需要 --gene=<基因名>", call. = FALSE)

  expr <- read_table_auto(expr_path, row_names = 1)
  cat(sprintf("[input] 表达矩阵：%d 行 × %d 列\n", nrow(expr), ncol(expr)))

  # 肿瘤样本筛选（必须在 12 位截断之前：截断后第 4 段就没了）
  ft <- filter_tumor(colnames(expr), keep_codes = arg_vec(args, "tumor_codes", "0"),
                     enabled = tumor_only)
  if (ft$applied) {
    cat("[filter] 样本类型分布:\n"); print(table(ft$codes, useNA = "ifany"))
    cat(sprintf("[filter] 仅保留 sample type 首字符 ∈ {%s}：%d -> %d 列\n",
                paste(arg_vec(args, "tumor_codes", "0"), collapse = ","),
                ncol(expr), sum(ft$keep)))
    expr <- expr[, ft$keep, drop = FALSE]
  } else if (tumor_only) {
    cat("[filter] 样本名不含 TCGA 第 4 段样本身份码，跳过肿瘤筛选（--tumor_only 对非 TCGA 命名无效）\n")
  }

  id_norm <- arg_chr(args, "id_norm", "auto")
  e_ids <- normalize_ids(colnames(expr), id_norm)
  cli <- read_clin_auto(clin_path, id_col = arg_chr(args, "clin_id_col", NULL))
  c_ids <- normalize_ids(rownames(cli), id_norm)

  tcol <- if (!is.null(time_col_in)) time_col_in else
    intersect(c("OS.time", "time", "OS_time", "survival_time", "time_days"), colnames(cli))[1]
  ecol <- if (!is.null(event_col_in)) event_col_in else
    intersect(c("OS", "state", "status", "event"), colnames(cli))[1]
  if (is.null(tcol) || is.na(tcol)) stop("临床表未找到时间列，请用 --time_col= 指定。可用: ",
                                        paste(colnames(cli), collapse = ", "), call. = FALSE)
  if (is.null(ecol) || is.na(ecol)) stop("临床表未找到状态列，请用 --event_col= 指定。可用: ",
                                        paste(colnames(cli), collapse = ", "), call. = FALSE)

  if (!gene %in% rownames(expr)) {
    hint <- grep(gene, rownames(expr), ignore.case = TRUE, value = TRUE)[1]
    stop(sprintf("表达矩阵中找不到 '%s'%s", gene,
                 if (!is.na(hint)) paste0("（最接近的是 '", hint, "'）") else ""), call. = FALSE)
  }
  common <- intersect(e_ids, c_ids)
  cat(sprintf("[merge] 表达矩阵样本 %d，临床表样本 %d，交集 %d\n",
              length(e_ids), length(c_ids), length(common)))
  if (length(common) == 0) {
    stop(sprintf("样本名无交集。\n  表达矩阵前 5: %s\n  临床表前 5:   %s\n提示: 用 --id_norm=tcga 或 --id_norm=none 调整命名规范化方式",
                 paste(head(e_ids, 5), collapse = ", "),
                 paste(head(c_ids, 5), collapse = ", ")), call. = FALSE)
  }
  if (length(common) < 10) warning(sprintf("交集样本仅 %d 例，结果仅供参考", length(common)), call. = FALSE)
  dup_e <- sum(duplicated(e_ids)); dup_c <- sum(duplicated(c_ids))
  if (dup_e || dup_c) cat(sprintf("[merge] 注意: 规范化后重复 ID — 表达矩阵 %d 个、临床表 %d 个，重复者只取首次出现\n",
                                   dup_e, dup_c))

  ei <- match(common, e_ids)
  ci <- match(common, c_ids)
  value <- suppressWarnings(as.numeric(expr[gene, ei]))
  cli_use <- cli[ci, , drop = FALSE]
  raw_time  <- pick_numeric(cli_use, tcol, "时间列")
  raw_event <- pick_numeric(cli_use, ecol, "状态列")
  samples   <- common
  gene_label <- gene
}

# ------------------------- 2. 事件编码 -------------------------------
ev_levels <- arg_vec(args, "event_levels")
cat(sprintf("[event] 状态列 '%s' 取值分布:\n", ecol))
print(table(raw_event, useNA = "ifany"))
if (length(ev_levels)) {
  raw_event <- as.integer(raw_event %in% ev_levels)
  cat(sprintf("[event] 按 --event_levels={%s} 重编码: 事件=%d, 删失=%d\n",
              paste(ev_levels, collapse = ","), sum(raw_event == 1), sum(raw_event == 0)))
} else if (!all(raw_event[!is.na(raw_event)] %in% c(0, 1))) {
  warning("状态列取值不是 0/1，已原样使用（1=事件 的假定可能不成立）。如为 1/2 编码请加 --event_levels=2",
          call. = FALSE)
}

# ------------------------- 3. 时间单位 -------------------------------
time_div <- arg_num(args, "time_div", NA_real_)
time_in  <- arg_chr(args, "time_in", "auto")
time_out <- arg_chr(args, "time_out", "year")
if (!time_out %in% names(UNIT_DAYS))
  stop("--time_out 只能是 day / month / year，收到: ", time_out, call. = FALSE)
if (!identical(time_in, "auto") && !time_in %in% names(UNIT_DAYS))
  stop("--time_in 只能是 auto / day / month / year，收到: ", time_in, call. = FALSE)

det <- detect_time_unit(raw_time)
if (identical(time_in, "auto")) {
  time_in <- det$unit
  cat(sprintf("[time] 自动推断输入单位 = %s （中位数 %.3g，最大 %.3g）\n", time_in, det$median, det$max))
  cat("[time] !! 单位推断只是启发式，请核对：肿瘤 OS 通常 1~5 年 = 12~60 月 = 365~1825 天。\n")
  cat(sprintf("[time] !! 若不符请显式指定 --time_in=day|month|year\n"))
} else {
  cat(sprintf("[time] 输入单位由参数指定 = %s（数据中位数 %.3g，最大 %.3g）\n", time_in, det$median, det$max))
}
if (!is.na(time_div)) {
  ana_time <- raw_time / time_div
  cat(sprintf("[time] 使用 --time_div=%g 直接相除（复刻参考代码），横轴标签仍为 %s\n",
              time_div, unit_label_cn[[time_out]]))
} else {
  ana_time <- convert_time(raw_time, time_in, time_out)
}
cat(sprintf("[time] 分析用时间单位 = %s（中位数 %.3g，最大 %.3g）\n",
            unit_label_cn[[time_out]], median(ana_time, na.rm = TRUE), max(ana_time, na.rm = TRUE)))

# ------------------------- 4. 清洗 -------------------------------
keep <- is.finite(ana_time) & ana_time > min_time & !is.na(raw_event) & is.finite(value)
if (sum(!keep)) cat(sprintf("[clean] 剔除 %d 例（时间缺失/<=%g/状态缺失/表达量缺失），保留 %d 例\n",
                            sum(!keep), min_time, sum(keep)))
ana_time <- ana_time[keep]; ev <- raw_event[keep]
value <- value[keep]; samples <- samples[keep]
if (length(ev) < 10) stop("有效样本太少（<10），无法分析", call. = FALSE)

# ------------------------- 5-6. 两级分组并选定方案 -------------------------------
cat("\n")
two <- two_tier_km(ana_time, ev, value, pct = pct, mode = mode_req,
                   p_th = p_th, minprop = minprop, pmethod = pmethod, verbose = TRUE)
mode_used <- two$mode_used
final_group <- two$group
final_cut   <- two$cutoff
st          <- two$stats
st_q        <- two$stats_quantile
st_c        <- two$stats_cutpoint
qres        <- two$qres
cp          <- two$cutpoint
cp_msg      <- two$cutpoint_msg
sig         <- two$significant

cat(sprintf("\n[final] 采用方案 = %s（请求 %s）| 阈值 = %.6g | %s | %s | %s\n",
            mode_used, mode_req, final_cut, fmt_p(st$pvalue),
            fmt_hr(st$hr, st$hr_low, st$hr_hi), ifelse(sig, "显著 (p<阈值)", "未达显著")))

# ------------------------- 7. 绘图 -------------------------------
d <- data.frame(time = ana_time, event = ev,
                group = factor(final_group, levels = c("Low", "High")),
                value = value)
fit <- survfit(Surv(time, event) ~ group, data = d)

p_line <- sprintf("%s\n%s", fmt_hr(st$hr, st$hr_low, st$hr_hi), fmt_p(st$pvalue))
if (!identical(style, "threshold")) {
  p_line <- sprintf("%s\n%s", fmt_p(st$pvalue), fmt_hr(st$hr, st$hr_low, st$hr_hi))
}
ascii_guard(c(gene_label, p_line), "图内文字")

break_by <- arg_num(args, "break_time", NA_real_)
if (is.na(break_by)) break_by <- NULL
prefix <- arg_chr(args, "prefix", NULL)
if (is.null(prefix)) prefix <- paste0(gsub("[^A-Za-z0-9_.-]", "_", gene_label), "_KM")

out_pdf <- if (do_pdf) file.path(outdir, paste0(prefix, "_KM.pdf")) else NULL
out_png <- if (do_png) file.path(outdir, paste0(prefix, "_KM.png")) else NULL
pl <- make_km_plot(fit, d, gene_label, style = style, unit_out = time_out,
                   break_by = break_by, palette = pal, p_label = p_line,
                   hr_on_plot = hr_on_plot, risk_table = risk_table,
                   conf_int = conf_int, median_line = median_line,
                   legend_pos = legend_pos, pval_size = pval_size,
                   risk_height = risk_h)
save_km(pl, out_pdf, out_png, width = width, height = height, dpi = dpi)
cat(sprintf("[out] 图: %s\n", paste(na.omit(c(out_pdf, out_png)), collapse = "  ")))

# ------------------------- 8. 明细与报告 -------------------------------
gmap <- data.frame(sample = samples, value = value, group = as.character(final_group),
                   time = ana_time, event = ev, stringsAsFactors = FALSE)
gmap_path <- file.path(outdir, paste0(prefix, "_group_map.csv"))
write.csv(gmap, gmap_path, row.names = FALSE, fileEncoding = "UTF-8")

ms <- st$median_surv
med_txt <- if (is.null(ms)) "" else paste(sprintf("%s:%.3g", names(ms), ms), collapse = " | ")

kv <- list(
  gene = gene_label, mode_requested = mode_req, mode_used = mode_used,
  cutoff = final_cut, pct = if (identical(mode_used, "cutpoint")) NA else pct,
  method = if (identical(mode_used, "cutpoint")) "maxstat(surv_cutpoint)" else "quantile split",
  minprop = if (identical(mode_used, "cutpoint")) minprop else NA,
  n_total = length(ev), n_event = sum(ev == 1),
  n_low = sum(final_group == "Low"), n_high = sum(final_group == "High"),
  p_logrank_quantile = st_q$pvalue,
  p_logrank_cutpoint = if (is.null(st_c)) NA else st_c$pvalue,
  p_maxstat_adjusted = two$maxstat_p_adjusted,
  maxstat_statistic = two$maxstat_statistic,
  maxstat_pmethod = two$maxstat_pmethod,
  p_used = st$pvalue,
  p_threshold = p_th,
  significant = sig,
  fallback_used = identical(mode_used, "cutpoint") && !identical(mode_req, "cutpoint"),
  hr_high_vs_low = st$hr, hr_low95 = st$hr_low, hr_high95 = st$hr_hi,
  p_cox = if (is.null(st$cox)) NA else st$cox_p,
  median_surv_by_group = med_txt,
  time_unit_in = time_in, time_unit_out = time_out,
  time_div = time_div, min_time = min_time,
  time_col = tcol, event_col = ecol, event_levels = paste(ev_levels, collapse = ","),
  style = style, conf_int = conf_int, risk_table = risk_table,
  figure = paste(na.omit(c(out_pdf, out_png)), collapse = " | "),
  group_map = gmap_path,
  cutpoint_error = cp_msg
)
stats_path <- file.path(outdir, paste0(prefix, "_stats.txt"))
write_kv(kv, stats_path, title = paste0("KM 生存分析统计  ", gene_label))

# 决策文件：把两级 p 值并列，便于复核"是否真的需要回退"
dec <- c(
  sprintf("基因/变量        : %s", gene_label),
  sprintf("样本数           : %d（事件 %d）", length(ev), sum(ev == 1)),
  sprintf("时间单位         : 输入 %s -> 分析 %s（%s）", time_in, time_out,
          if (is.na(time_div)) "按 30.44/365 换算" else sprintf("除以 %g", time_div)),
  "",
  sprintf("第一级 分位值分组: 阈值 %.6g  Low=%d High=%d  log-rank %s  %s",
          qres$cutoff, sum(qres$group == "Low"), sum(qres$group == "High"),
          fmt_p(st_q$pvalue), fmt_hr(st_q$hr, st_q$hr_low, st_q$hr_hi)),
  if (is.null(st_c)) sprintf("第二级 最佳截断值: %s",
                             if (nzchar(cp_msg)) paste0("未执行/失败 — ", cp_msg) else "未执行（第一级已达标）")
  else sprintf("第二级 最佳截断值: 阈值 %.6g (minprop=%.2g)  Low=%d High=%d  log-rank %s  %s",
               cp$cutoff, minprop, cp$n_low, cp$n_high,
               fmt_p(st_c$pvalue), fmt_hr(st_c$hr, st_c$hr_low, st_c$hr_hi)),
  if (!is.na(two$maxstat_p_adjusted))
    sprintf("         maxstat 校正 p: %s (%s) —— 阈值是选出来的，这个 p 才是应报告的",
            formatC(two$maxstat_p_adjusted, format = "e", digits = 2), two$maxstat_pmethod) else
    "         maxstat 校正 p: 未取得（报告时请注明阈值由数据驱动）",
  "",
  sprintf("最终采用         : %s（请求 %s）", mode_used, mode_req),
  sprintf("显著性           : p = %s  %s",
          ifelse(is.na(st$pvalue), "NA", format(st$pvalue, digits = 4)),
          ifelse(sig, "< 阈值，显著", ">= 阈值，未达显著")),
  "",
  "== 必须核对 ==",
  sprintf("1) 时间单位是否为 %s？原始数据中位数 %.3g、最大 %.3g。不符请加 --time_in=",
          time_in, det$median, det$max),
  sprintf("2) 分组样本数 Low=%d / High=%d 是否可接受？（%s）",
          sum(final_group == "Low"), sum(final_group == "High"),
          ifelse(min(table(final_group)) < 0.1 * length(ev), "严重不平衡，警惕并列值", "正常")),
  "3) 逐样本分组见 *_group_map.csv",
  if (identical(mode_used, "cutpoint"))
    paste0("4) 已使用数据驱动的最佳截断值：阈值是在同一批数据上选出来的，属探索性分析。\n",
           "   报告时不要只写 log-rank p，请一并给出 maxstat 校正 p；正式结论需在独立队列复现。") else
    "4) 未使用数据驱动阈值，p 值可直接报告"
)
dec_path <- file.path(outdir, paste0(prefix, "_km_decision.txt"))
writeLines(dec, dec_path, useBytes = TRUE)
cat(sprintf("[out] 统计: %s\n[out] 决策: %s\n[out] 分组: %s\n", stats_path, dec_path, gmap_path))
cat("\n", paste(dec, collapse = "\n"), "\n", sep = "")
cat("=======================================================\n")
