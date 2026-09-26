# ===========================================================================
# lib_surv_common.R —— survival-km-analysis 公共函数库
#
# 内容：参数解析 / 稳健读表 / TCGA 样本名规范化 / 肿瘤样本筛选 /
#       时间单位推断与换算 / 中位值分组 / surv_cutpoint 最佳阈值 /
#       log-rank + Cox HR 统计 / KM 绘图（两种风格）/ 报告落盘
#
# 本文件被 01_km_survival.R 与 02_km_batch.R 共同 source()。
# ===========================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
})

# --------------------------------------------------------------------------
# 1. 参数解析
# --------------------------------------------------------------------------
# 支持 `--key=value` 与裸 `--flag`（视为 TRUE）。
# 只有遇到 `--` 开头的键才认识，避免把值里的 '-' 误判成选项。
parse_args <- function(args) {
  out <- list()
  for (a in args) {
    if (!grepl("^--", a)) next
    a <- sub("^--", "", a)
    if (grepl("=", a, fixed = TRUE)) {
      k <- sub("=.*$", "", a)
      v <- sub("^[^=]*=", "", a)
      out[[k]] <- v
    } else {
      out[[a]] <- "TRUE"
    }
  }
  out
}

arg_chr <- function(a, k, default = NULL) {
  v <- a[[k]]
  if (is.null(v) || !nzchar(v)) return(default)
  v
}
arg_num <- function(a, k, default = NA_real_) {
  v <- a[[k]]
  if (is.null(v) || !nzchar(v)) return(default)
  suppressWarnings(as.numeric(v))
}
arg_bool <- function(a, k, default = FALSE) {
  v <- a[[k]]
  if (is.null(v) || !nzchar(v)) return(default)
  toupper(v) %in% c("TRUE", "T", "1", "YES", "Y")
}
# 逗号分隔的字符串向量
arg_vec <- function(a, k, default = character(0)) {
  v <- a[[k]]
  if (is.null(v) || !nzchar(v)) return(default)
  trimws(strsplit(v, ",", fixed = TRUE)[[1]])
}

# 把 `--a=b` 里形如 a.b 的键还原（R 里用 a.b 表示 a-b），便于传连字符参数
# 用法：args <- parse_args(commandArgs(TRUE)); key_posix(args, "p-threshold")
arg_get2 <- function(a, k1, k2) {
  if (!is.null(a[[k1]])) return(a[[k1]])
  a[[k2]]
}

# --------------------------------------------------------------------------
# 2. 稳健读表
# --------------------------------------------------------------------------
# 自动判别制表符 / 逗号；保留原始列名；不把 # 当注释；引号关闭。
# row_names: 1 = 第一列作行名；0/NULL = 不作行名
# 去掉包裹在字段外的成对引号。
# 必要性：read.table(quote="") 关闭了引号解析（避免基因名/ID 里的引号引发截断），
# 代价是 write.csv 之类产出的表头/行名会带着**字面引号**（列名变成 `"time"` 而不是
# time），随后按列名取值一律失败。实测 2026-09-18: re.csv（write.csv 产物）报
# `时间列 'time' 不在表中。可用列: "time", "state", ...`。
strip_quotes <- function(x) {
  x <- as.character(x)
  x <- sub('^"(.*)"$', "\\1", x)
  x <- sub("^'(.*)'$", "\\1", x)
  x
}

read_table_auto <- function(path, row_names = 1, quiet = FALSE) {
  if (!file.exists(path)) stop("找不到文件: ", path, call. = FALSE)
  l1 <- readLines(path, n = 1, warn = FALSE)
  if (length(l1) != 1L) stop("文件为空: ", path, call. = FALSE)
  n_tab <- lengths(regmatches(l1, gregexpr("\t", l1, fixed = TRUE)))
  n_com <- lengths(regmatches(l1, gregexpr(",", l1, fixed = TRUE)))
  sep <- if (n_tab >= n_com) "\t" else ","
  if (!quiet) cat(sprintf("[read] %s  分隔符=%s\n", basename(path),
                          ifelse(sep == "\t", "TAB", "COMMA")))
  df <- read.table(path, header = TRUE, sep = sep, check.names = FALSE,
                   row.names = if (is.null(row_names) || row_names == 0) NULL else row_names,
                   quote = "", comment.char = "", fill = TRUE,
                   stringsAsFactors = FALSE)
  colnames(df) <- strip_quotes(colnames(df))
  rn <- rownames(df)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(df))))) {
    rownames(df) <- strip_quotes(rn)
  }
  df
}

# 临床表读取：把样本 ID 列取出作行名。ID 列自动探测（名字像 ID 的优先，否则第 1 列）。
# 全程用**位置索引**操作：
#   实测 2026-09-18 —— `time.csv` 的表头是 `,time,state`（首字段为空），
#   read.table 得到列名 c("", "time", "state")。此时 `df[[""]]` 返回 NULL，
#   长度不匹配会报 `invalid 'row.names' length`，整条链断在这里。
read_clin_auto <- function(path, id_col = NULL, quiet = FALSE) {
  raw <- read_table_auto(path, row_names = 0, quiet = quiet)
  if (ncol(raw) < 2) stop("临床表至少需要 2 列（样本 ID + 时间/状态）: ", path, call. = FALSE)
  cn <- colnames(raw)
  idx <- NA_integer_
  if (!is.null(id_col) && nzchar(id_col)) {
    if (id_col %in% cn) idx <- which(cn == id_col)[1]
    else {
      j <- suppressWarnings(as.integer(id_col))
      if (!is.na(j) && j >= 1 && j <= ncol(raw)) idx <- j
    }
    if (is.na(idx)) stop("--clin_id_col=", id_col, " 不在临床表列名中。可用: ",
                         paste(cn, collapse = ", "), call. = FALSE)
  }
  if (is.na(idx)) {
    idlike <- nzchar(cn) &
      grepl("^(id|ID|Id|sample|Sample|sample_id|sampleID|barcode|patient|Patient)$", cn)
    if (any(idlike)) idx <- which(idlike)[1]
    else if (!nzchar(cn[1]) || is.character(raw[[1]]) || is.factor(raw[[1]])) idx <- 1L
    else idx <- 1L
  }
  ids <- strip_quotes(as.character(raw[[idx]]))
  out <- raw[, setdiff(seq_len(ncol(raw)), idx), drop = FALSE]   # 位置索引，避开空列名
  if (ncol(out) < 1) stop("临床表除 ID 列外没有其它列: ", path, call. = FALSE)
  rownames(out) <- ids
  if (!quiet) cat(sprintf("[clin] ID 列 = '%s'，%d 行 × %d 列\n",
                          ifelse(nzchar(cn[idx]), cn[idx], paste0("第", idx, "列")),
                          nrow(out), ncol(out)))
  out
}

# 读基因列表文件。不假定"每行一个基因"——实测本机的 merge2.csv 是
# `"1","CCNB1"` 这种"索引+基因"两列结构（read.csv 的产物），按行读会得到
# `1","CCNB1` 这种垃圾 token。
# 策略：逐行按逗号/制表符切栏、去掉成对引号，取**第一个像基因符号的字段**
# （以字母开头、只含字母数字和 . _ -）。都不像时退回第一个非空字段。
read_gene_list <- function(path) {
  if (!file.exists(path)) stop("找不到基因列表文件: ", path, call. = FALSE)
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  sym_re <- "^[A-Za-z][A-Za-z0-9_.-]*$"
  hdr_re <- "^(gene|Gene|GENE|genes|Genes|id|ID|symbol|Symbol|V1|X|e|x|rowname|row)$"
  out <- character(0)
  for (l in ln) {
    f <- strip_quotes(strsplit(l, "[,\t;]", perl = TRUE)[[1]])
    f <- trimws(f)
    f <- f[nzchar(f) & !grepl("^(NA|\\.|-)$", f)]
    if (!length(f)) next
    hit <- f[grepl(sym_re, f) & nchar(f) >= 2]
    if (length(hit)) out <- c(out, hit[1]) else out <- c(out, f[1])
  }
  out <- unique(out[!grepl(hdr_re, out)])
  out
}

# 取数值列；列名不存在时报出可用列名
pick_numeric <- function(df, col, what = "列") {
  if (is.null(col) || !nzchar(col)) stop("缺少参数：", what, call. = FALSE)
  if (!col %in% colnames(df)) {
    stop(sprintf("%s '%s' 不在表中。可用列: %s", what, col,
                 paste(colnames(df), collapse = ", ")), call. = FALSE)
  }
  suppressWarnings(as.numeric(df[[col]]))
}

# --------------------------------------------------------------------------
# 3. 样本名规范化（TCGA 条码 -> 12 位 participant ID）
# --------------------------------------------------------------------------
looks_tcga <- function(x) {
  x <- as.character(x)
  mean(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}", x)) >= 0.5
}

# 只对 TCGA 条码做 12 位截断 + . -> -；其它命名（GSM*、自定义）原样保留。
normalize_ids <- function(x, mode = "auto") {
  x <- as.character(x)
  if (identical(mode, "none")) return(x)
  if (identical(mode, "tcga") || (identical(mode, "auto") && looks_tcga(x))) {
    x <- gsub("[.]", "-", x)
    x <- substr(x, 1, 12)
  }
  x
}

# 肿瘤样本筛选：TCGA 第 4 段形如 01A/11A，取首字符比对 codes。
# 必须在 12 位截断之前调用（截断后第 4 段就没了）。
sample_type_codes <- function(ids) {
  parts <- strsplit(as.character(ids), "-", fixed = TRUE)
  vapply(parts, function(p) if (length(p) >= 4L) substr(p[4], 1, 1) else NA_character_,
         character(1))
}

filter_tumor <- function(ids, keep_codes = "0", enabled = TRUE) {
  codes <- sample_type_codes(ids)
  if (!enabled) return(list(keep = rep(TRUE, length(ids)), codes = codes, applied = FALSE))
  if (all(is.na(codes))) {
    return(list(keep = rep(TRUE, length(ids)), codes = codes, applied = FALSE))
  }
  # 认不出身份码的样本（非 TCGA 命名 / 已被截断）一律保留，不做判定
  list(keep = is.na(codes) | codes %in% keep_codes, codes = codes, applied = TRUE)
}

# --------------------------------------------------------------------------
# 4. 时间单位
# --------------------------------------------------------------------------
# 单位换算一律经"天"中转：day=1, month=30.44, year=365
# （与参考代码一致：最佳阈值-km.R 用 /30.44；生存分析-km.R 用 /12 即 365/30.44≈11.99）
UNIT_DAYS <- c(day = 1, month = 30.44, year = 365)

unit_label_cn <- c(day = "days", month = "months", year = "years")

# 推断输入时间单位。启发式，**不可当成定论**：
#   max<=25      -> year
#   25<max<=400  -> month
#   max>400      -> day
detect_time_unit <- function(t) {
  t <- t[is.finite(t)]
  mx <- max(t); md <- median(t)
  unit <- if (mx <= 25) "year" else if (mx <= 400) "month" else "day"
  list(unit = unit, median = md, max = mx)
}

# 默认横轴刻度间隔：年->2（参考 1），月->24（参考 2），天->365
default_break_time <- function(unit_out) {
  switch(unit_out, year = 2, month = 24, day = 365, 2)
}

convert_time <- function(t, from, to) {
  if (identical(from, to)) return(t)
  t * UNIT_DAYS[[from]] / UNIT_DAYS[[to]]
}

# --------------------------------------------------------------------------
# 5. 分组
# --------------------------------------------------------------------------
# 按分位数切分。pct=0.5 即中位值；pct=1/3 即三分位下界（与参考注释里的 1/3 一致）。
# 注意：用 `>` 而非 `>=`，与参考代码一致 —— 表达量取值为 0/固定值的基因会让
# 大量样本落在等于阈值的一侧，必须检查返回的分组样本数。
quantile_split <- function(v, pct = 0.5, high_label = "High", low_label = "Low") {
  cutv <- as.numeric(stats::quantile(v, pct, na.rm = TRUE, names = FALSE))
  g <- ifelse(v > cutv, high_label, low_label)
  g[is.na(v)] <- NA_character_
  list(group = g, cutoff = cutv, pct = pct)
}

# 最佳截断值：survminer::surv_cutpoint（底层 maxstat） + surv_categorize
# 返回 ok / msg / cutoff / group / method
#
# 返回结构要点（2026-09-18 实测 survminer 0.5.2）：
#   names(cp) == c(<变量名>, "data", "minprop", "cutpoint")
#   cp$cutpoint 是 data.frame，列为 cutpoint / statistic（标准化 log-rank 统计量）
#   cp[[<变量名>]] 是 maxstat::maxstat.test 的 "maxtest" 对象，但它是用
#     pmethod="none" 调出来的，**其 $p.value 恒为 NA**，而 maxstat 并不导出
#     pvalue()。要拿校正 p 只能自己用 pmethod="Lau92" 重跑（见下方代码内注释）。
#   surv_categorize() 返回的是**字符向量**（不是 factor），取值 low/high 按
#   数值 ≤ / > 切点划分，所以水平顺序要靠调用方自己用 factor() 固定。
best_cutpoint <- function(df, time_col, event_col, var_col, minprop = 0.3,
                          pmethod = "Lau92") {
  fmls <- names(formals(survminer::surv_cutpoint))
  a <- list(data = df, time = time_col, event = event_col,
            variables = var_col, minprop = minprop)
  if ("progressbar" %in% fmls) a$progressbar <- FALSE
  cp <- tryCatch(do.call(survminer::surv_cutpoint, a), error = function(e) e)
  if (inherits(cp, "error")) {
    return(list(ok = FALSE, msg = conditionMessage(cp)))
  }
  cutv <- tryCatch(as.numeric(cp$cutpoint$cutpoint[1]), error = function(e) NA_real_)
  dfc <- tryCatch(survminer::surv_categorize(cp, labels = c("Low", "High")),
                  error = function(e) e)
  if (inherits(dfc, "error")) {
    return(list(ok = FALSE, msg = conditionMessage(dfc)))
  }
  g <- as.character(dfc[[var_col]])
  stat <- tryCatch(as.numeric(cp$cutpoint$statistic[1]), error = function(e) NA_real_)

  # ---- maxstat 校正 p 值 ----------------------------------------------------
  # surv_cutpoint 内部调用 maxstat.test(pmethod="none")，所以 cp[[var]]$p.value
  # 恒为 NA，必须自己再跑一次 maxstat.test 才能拿到"扣除择优偏倚"的 p。
  # 实测（2026-09-18，re.csv 582 例）：
  #   同一组数据用相同 minprop 重跑，切点与统计量 M 逐位复现（SOX9 均为 153.886 / M=2.1028）
  #   pmethod="Lau92"     -> p = 0.2998   （可用）
  #   pmethod="Lau94"     -> p = 1.156    （**越界，不可用**）
  #   pmethod="condMC"/"HL" 极慢，不适合批量
  # 故默认 Lau92；结果越界或非有限时一律置 NA 并说明，不猜。
  adj_p <- NA_real_; adj_pm <- NA_character_; adj_ok <- NA
  for (pm in unique(c(pmethod, "Lau92", "exactGauss"))) {
    r <- tryCatch(maxstat::maxstat.test(
        survival::Surv(df[[time_col]], df[[event_col]]) ~ V,
        data = data.frame(V = df[[var_col]]),
        smethod = "LogRank", pmethod = pm, minprop = minprop),
      error = function(e) e)
    if (inherits(r, "error")) next
    pv <- suppressWarnings(as.numeric(r$p.value))
    if (length(pv) == 1L && is.finite(pv) && pv >= 0 && pv <= 1) {
      adj_p <- pv; adj_pm <- pm
      adj_ok <- isTRUE(all.equal(as.numeric(r$estimate), cutv, tolerance = 1e-6))
      break
    }
  }

  list(ok = TRUE, msg = "", cutoff = cutv, group = g,
       maxstat_statistic = stat, maxstat_p_adjusted = adj_p,
       maxstat_pmethod = adj_pm, maxstat_cutpoint_match = adj_ok,
       n_low = sum(g == "Low", na.rm = TRUE), n_high = sum(g == "High", na.rm = TRUE),
       method = "maxstat(surv_cutpoint)")
}

# --------------------------------------------------------------------------
# 5b. 两级策略（01 与 02 共用，避免两份逻辑漂移）
# --------------------------------------------------------------------------
# 第一级：分位值分组（pct=0.5 即中位值，等价于 生存分析-km.R）
# 第二级：当第一级 log-rank p >= p_th（或强制 mode="cutpoint"）时，
#         用 surv_cutpoint 的最佳截断值重新分组（等价于 最佳阈值-km.R）
# 最终取 p 更小者；两级 p 都返回，写进决策文件供复核。
two_tier_km <- function(time, event, value, pct = 0.5, mode = "auto",
                        p_th = 0.05, minprop = 0.3, pmethod = "Lau92",
                        verbose = TRUE) {
  stopifnot(length(time) == length(event), length(time) == length(value))
  qres <- quantile_split(value, pct)
  st_q <- km_stats(time, event, qres$group)
  if (verbose) {
    n_unbal <- min(sum(qres$group == "Low"), sum(qres$group == "High"), na.rm = TRUE) <
      0.1 * length(time)
    cat(sprintf("[step1] 分位值分组 pct=%.4g  阈值=%.6g  Low=%d High=%d | log-rank %s | %s%s\n",
                pct, qres$cutoff, sum(qres$group == "Low"), sum(qres$group == "High"),
                fmt_p(st_q$pvalue), fmt_hr(st_q$hr, st_q$hr_low, st_q$hr_hi),
                if (n_unbal) "   [警告: 分组严重不平衡]" else ""))
  }

  base_mode <- if (abs(pct - 0.5) < 1e-9) "median" else "quantile"
  use_cut <- identical(mode, "cutpoint")
  if (identical(mode, "auto") && (is.na(st_q$pvalue) || st_q$pvalue >= p_th)) {
    use_cut <- TRUE
    if (verbose) cat(sprintf("[step2] 常规分组 p=%s 未低于 %.4g -> 回退到最佳截断值\n",
                             ifelse(is.na(st_q$pvalue), "NA", format(st_q$pvalue, digits = 4)), p_th))
  } else if (identical(mode, "auto") && verbose) {
    cat(sprintf("[step2] 常规分组 p=%.4g < %.4g -> 无需回退\n", st_q$pvalue, p_th))
  }

  st_c <- NULL; cp <- NULL; cp_msg <- ""
  if (use_cut) {
    dfc <- data.frame(.t = time, .e = event, .v = value)
    cp <- best_cutpoint(dfc, ".t", ".e", ".v", minprop = minprop, pmethod = pmethod)
    if (cp$ok) {
      st_c <- km_stats(time, event, cp$group)
      if (verbose) {
        cat(sprintf("[step2] 最佳截断值=%.6g (minprop=%.2g)  Low=%d High=%d | log-rank %s | %s\n",
                    cp$cutoff, minprop, cp$n_low, cp$n_high,
                    fmt_p(st_c$pvalue), fmt_hr(st_c$hr, st_c$hr_low, st_c$hr_hi)))
        if (!is.na(cp$maxstat_p_adjusted)) {
          cat(sprintf("[step2] maxstat 校正 p = %s (%s)%s —— 已扣除择优偏倚，是这组结果里最该报告的 p\n",
                      formatC(cp$maxstat_p_adjusted, format = "e", digits = 2), cp$maxstat_pmethod,
                      if (isFALSE(cp$maxstat_cutpoint_match)) " [注意: 重跑切点与 surv_cutpoint 不一致]" else ""))
        } else {
          cat("[step2] 未能取得 maxstat 校正 p（maxstat.test 重跑失败或结果越界），报告时请指明阈值是数据驱动的\n")
        }
      }
    } else {
      cp_msg <- cp$msg
      if (verbose) {
        cat(sprintf("[step2] surv_cutpoint 失败: %s\n", cp$msg))
        cat("[step2] 常见原因: 变量取值过于集中（例如大部分为 0），minprop 下找不到合法切点。\n")
        cat("[step2] 可尝试: 调小 --minprop、改 --pct、或先按常规分组出图再人工定阈值。\n")
      }
    }
  }

  mode_used <- base_mode
  if (use_cut && !is.null(st_c)) {
    better <- is.na(st_q$pvalue) || (!is.na(st_c$pvalue) && st_c$pvalue < st_q$pvalue)
    if (better) mode_used <- "cutpoint"
  }
  final_group <- if (identical(mode_used, "cutpoint")) cp$group else qres$group
  final_cut   <- if (identical(mode_used, "cutpoint")) cp$cutoff else qres$cutoff
  st          <- if (identical(mode_used, "cutpoint")) st_c else st_q
  list(mode_used = mode_used, mode_requested = mode, base_mode = base_mode,
       group = final_group, cutoff = final_cut, stats = st, stats_quantile = st_q,
       stats_cutpoint = st_c, qres = qres, cutpoint = cp, cutpoint_msg = cp_msg,
       maxstat_statistic = if (is.null(cp)) NA_real_ else cp$maxstat_statistic,
       maxstat_p_adjusted = if (is.null(cp)) NA_real_ else cp$maxstat_p_adjusted,
       maxstat_pmethod = if (is.null(cp)) NA_character_ else cp$maxstat_pmethod,
       maxstat_cutpoint_match = if (is.null(cp)) NA else cp$maxstat_cutpoint_match,
       pct = pct, minprop = minprop, p_threshold = p_th,
       significant = !is.na(st$pvalue) && st$pvalue < p_th,
       fallback_used = identical(mode_used, "cutpoint") && !identical(mode, "cutpoint"))
}

# --------------------------------------------------------------------------
# 6. 统计
# --------------------------------------------------------------------------
fmt_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) paste0("p = ", formatC(p, format = "e", digits = 2))
  else sprintf("p = %.4f", p)
}

fmt_hr <- function(hr, lo, hi) {
  if (is.na(hr)) return("HR = NA")
  sprintf("HR = %.2f (95%% CI: %.2f-%.2f)", hr, lo, hi)
}

# 返回 log-rank p、Cox HR（**恒为 High vs Low**）、两组中位生存期与样本数
km_stats <- function(time, event, group, levels_ord = c("Low", "High")) {
  g <- factor(group, levels = levels_ord)
  keep <- is.finite(time) & !is.na(event) & !is.na(g)
  time <- time[keep]; event <- event[keep]; g <- droplevels(g[keep])
  n_groups <- nlevels(g)
  res <- list(n = length(time), n_groups = n_groups,
              n_per = as.list(table(g)), ok = TRUE, msg = "")
  if (n_groups < 2L) {
    res$ok <- FALSE; res$msg <- "有效分组少于 2 组"; return(res)
  }
  surv_obj <- Surv(time, event)
  sd <- survdiff(surv_obj ~ g)
  res$chisq  <- sd$chisq
  res$pvalue <- 1 - pchisq(sd$chisq, df = n_groups - 1L)
  res$fit    <- survfit(surv_obj ~ g)
  # 各组中位生存期
  tb <- tryCatch(summary(res$fit)$table, error = function(e) NULL)
  if (!is.null(tb)) {
    med <- as.data.frame(tb)
    res$median_surv <- setNames(med[["median"]],
                                sub("^g=", "", rownames(med)))
  }
  # Cox：levels 为 c("Low","High") 时系数即 High vs Low
  cox <- tryCatch(coxph(surv_obj ~ g), error = function(e) e)
  if (!inherits(cox, "error")) {
    cs <- summary(cox)
    res$cox <- cs
    res$hr     <- cs$conf.int[1, "exp(coef)"]
    res$hr_low <- cs$conf.int[1, "lower .95"]
    res$hr_hi  <- cs$conf.int[1, "upper .95"]
    res$cox_p  <- cs$coefficients[1, "Pr(>|z|)"]
    # 若水平顺序不是 Low,High，则取反，保证 HR 语义恒为 High vs Low
    if (!identical(levels_ord, c("Low", "High"))) {
      res$hr <- 1 / res$hr; tmp <- res$hr_low; res$hr_low <- 1 / res$hr_hi; res$hr_hi <- 1 / tmp
    }
  } else {
    res$hr <- res$hr_low <- res$hr_hi <- NA_real_
    res$msg <- paste0("Cox 未收敛: ", conditionMessage(cox))
  }
  res
}

# --------------------------------------------------------------------------
# 7. 绘图
# --------------------------------------------------------------------------
# 恒用 levels = c("Low","High") 建 fit/统计，HR 语义 = High vs Low。
#
# 三种风格。**风格是纯出图偏好，与分组方式（常规分位值 / 最佳阈值回退）无关** ——
# 同一批里混着"走了回退"和"没走回退"的基因，出的图版式必须完全一样。
#   unified   默认。统一版式。
#   classic   严格复刻 生存分析-km.R（与 unified 的唯一外观差别：风险表标签
#             用默认黑色，不跟曲线配色）
#   threshold 严格复刻 最佳阈值-km.R（teal/brick 配色、置信带开、theme_minimal、
#             面板内左下 HR+p 标注、图例 <基因>_low/_high）
#
# unified 把两份参考代码里互相冲突的外观项**定死**，用户不必每次去挑：
#   配色 / 图例标题 / 主题 / 标注机制 / 置信带 <- 取自主参考 生存分析-km.R
#   风险表标签跟曲线配色 (risk.table.col="strata") <- 取自 最佳阈值-km.R
#   y 轴标题显式写死 <- 免得随 survminer 版本漂移
STYLE_PALETTE <- list(
  unified   = c("MediumSeaGreen", "Firebrick3", "#6E568C", "#223D6C"),
  classic   = c("MediumSeaGreen", "Firebrick3", "#6E568C", "#223D6C"),
  threshold = c("#3090a1", "#bc5148")
)
STYLE_CHOICES <- names(STYLE_PALETTE)

# 各风格的默认置信带：统一版式与主参考一致为关；最佳阈值版式按参考 2 为开
style_default_conf_int <- function(style) identical(style, "threshold")

# 压掉 ggplot2 4.0 下 ggsurvplot 的一条无害提示
#   `Ignoring unknown labels: • colour : "XXX expression"`
# 图例标题实际是正常渲染的，纯粹是 survminer 内部重复调用 labs() 引起。
# 该提示在构建(ggsurvplot)与渲染(print)两个阶段各出现一次，且成套的 cli 输出
# 既可能走 warning 也可能走 message（2026-09-18 实测两种都要拦），故两者都处理。
quiet_gg <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("Ignoring unknown labels", conditionMessage(w))) invokeRestart("muffleWarning")
    },
    message = function(m) {
      if (grepl("Ignoring unknown labels", conditionMessage(m))) invokeRestart("muffleMessage")
    })
}

make_km_plot <- function(fit, d, gene, style = "unified", unit_out = "year",
                         break_by = NULL, palette = NULL, p_label = "",
                         hr_on_plot = TRUE, risk_table = TRUE, conf_int = FALSE,
                         median_line = "hv", legend_pos = NULL,
                         pval_size = 6, risk_height = 0.25) {
  if (!style %in% STYLE_CHOICES)
    stop("未知的 --style: ", style, "（可选 ", paste(STYLE_CHOICES, collapse = " / "), "）",
         call. = FALSE)
  if (is.null(break_by)) break_by <- default_break_time(unit_out)
  if (is.null(palette)) palette <- STYLE_PALETTE[[style]]
  xlab <- paste0("Time (", unit_label_cn[[unit_out]], ")")

  if (identical(style, "threshold")) {
    if (is.null(legend_pos)) legend_pos <- c(0.85, 0.85)
    # 图例位置必须通过 ggsurvplot 的 legend= 传，不能只写在 ggtheme 里：
    # ggsurvplot 默认 legend="top"，会覆盖 ggtheme 中的 legend.position，
    # 结果是图例跑到图外顶部（参考代码里 legend.position=c(0.85,0.85) 的
    # "面板内右上角"效果丢失）。2026-09-18 实测确认。
    pl <- quiet_gg(ggsurvplot(
      fit, data = d,
      surv.median.line = median_line, pval = FALSE,
      conf.int = conf_int, risk.table = risk_table, risk.table.col = "strata",
      xlab = xlab, ylab = "Survival Probability",
      legend.title = "", legend.labs = paste0(gene, c("_low", "_high")),
      legend = legend_pos, font.legend = 10,
      break.x.by = break_by, color = "strata", palette = palette,
      ggtheme = theme_minimal(base_size = 14) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(color = "black")),
      risk.table.height = risk_height))
    pl$plot <- pl$plot +
      annotate("text", x = max(d$time, na.rm = TRUE) * 0.05, y = 0.12,
               label = p_label, size = pval_size, hjust = 0)
  } else {
    # unified / classic：同一段代码，只差风险表标签是否跟曲线配色
    if (is.null(legend_pos)) legend_pos <- c(0.8, 0.8)
    risk_col <- if (identical(style, "unified")) "strata" else "black"
    label <- if (hr_on_plot) p_label else sub("\n.*$", "", p_label)
    pl <- quiet_gg(ggsurvplot(
      fit, data = d,
      conf.int = conf_int, pval = label, pval.size = pval_size,
      legend.title = paste0(gene, " expression"),
      legend.labs = c("Low", "High"),
      legend = legend_pos, font.legend = 10,
      xlab = xlab, ylab = "Survival probability",
      break.time.by = break_by, palette = palette,
      surv.median.line = median_line,
      risk.table = risk_table, risk.table.col = risk_col, cumevents = FALSE,
      risk.table.height = risk_height))
  }
  # ggsurvplot 会把 legend.title 顺带写进风险表的 y 轴标题，挤出一行多余文字
  if (risk_table && !is.null(pl$table)) {
    pl$table <- pl$table + ggplot2::ylab(NULL)
  }
  pl
}

# 把 ggsurvplot 对象写进 PDF / PNG（PNG 优先用 ragg，缺则退回 grDevices）
# 渲染阶段（print.ggsurvplot）同样会吐 "Ignoring unknown labels" 警告，
# 故统一走 print_km() 压掉。
print_km <- function(pl) quiet_gg(print(pl))

save_km <- function(pl, out_pdf = NULL, out_png = NULL, width = 8, height = 6.25,
                    dpi = 300) {
  if (!is.null(out_pdf)) {
    grDevices::pdf(out_pdf, width = width, height = height, onefile = FALSE)
    print_km(pl); grDevices::dev.off()
  }
  if (!is.null(out_png)) {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(out_png, width = width, height = height, units = "in", res = dpi)
    } else {
      grDevices::png(out_png, width = width, height = height, units = "in", res = dpi)
    }
    print_km(pl); grDevices::dev.off()
  }
}

# --------------------------------------------------------------------------
# 8. 杂项
# --------------------------------------------------------------------------
# 图内文字若含非 ASCII，Windows 缺字体会渲染成乱码。用前先检查。
ascii_guard <- function(x, what = "图内文字") {
  bad <- grepl("[^\x01-\x7F]", x)
  if (any(bad)) {
    warning(sprintf("%s 含非 ASCII 字符，Windows 下会渲染成乱码: %s", what,
                    paste(unique(x[bad]), collapse = ", ")), call. = FALSE)
  }
  invisible(any(bad))
}

ensure_dir <- function(d) {
  if (!is.null(d) && nzchar(d) && !dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(d)
}

# 写 key=value 报告
write_kv <- function(pairs, path, title = NULL) {
  lines <- character(0)
  if (!is.null(title)) lines <- c(lines, paste0("# ", title))
  for (k in names(pairs)) {
    v <- pairs[[k]]
    if (length(v) > 1) v <- paste(v, collapse = " | ")
    if (length(v) == 0) v <- ""
    lines <- c(lines, paste0(k, "=", as.character(v)))
  }
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}
