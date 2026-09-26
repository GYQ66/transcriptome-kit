# ===========================================================================
# lib_string_common.R —— string-ppi-network 公共库
#
# 只依赖 base R（+ 可选 igraph / ggraph）。被 01_string_network.R source。
# 提供：参数解析、基因表读取、symbol 别名修正、引擎无关的边表规范化、配色、
#       日志、报告与参数落盘。
#
# 注意：本文件内所有输出到 stdout / 文件的中文均为 UTF-8；
#       图内文字一律用 ASCII（Windows 缺中文字形会渲染成乱码 ASCII）。
# ===========================================================================

# ---------------------------------------------------------------------------
# 参数解析：只支持 --key=value 形式
# ---------------------------------------------------------------------------
parse_args <- function(args, defaults) {
  vals <- defaults
  vals[["help_flag"]] <- FALSE
  for (a in args) {
    if (identical(a, "-h") || identical(a, "--help")) {
      vals[["help_flag"]] <- TRUE
      next
    }
    if (!grepl("^--", a)) next
    a2 <- sub("^--", "", a)
    if (grepl("=", a2, fixed = TRUE)) {
      k <- sub("=.*$", "", a2)
      v <- sub("^[^=]*=", "", a2)
    } else {
      k <- a2
      v <- ""
    }
    vals[[k]] <- v
  }
  vals
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || length(x) == 0) return(default)
  if (is.logical(x)) return(x)
  s <- tolower(trimws(as.character(x[1])))
  if (s %in% c("1", "true", "t", "yes", "y", "on")) return(TRUE)
  if (s %in% c("0", "false", "f", "no", "n", "off", "")) return(FALSE)
  default
}

as_int <- function(x, default) {
  if (is.null(x) || length(x) == 0) return(default)
  s <- trimws(as.character(x[1]))
  if (!nzchar(s)) return(default)
  v <- suppressWarnings(as.integer(s))
  if (is.na(v)) return(default)
  v
}

as_num <- function(x, default) {
  if (is.null(x) || length(x) == 0) return(default)
  s <- trimws(as.character(x[1]))
  if (!nzchar(s)) return(default)
  v <- suppressWarnings(as.numeric(s))
  if (is.na(v)) return(default)
  v
}

# "a,b,c" -> c("a","b,c") 拆成向量；空串 -> character(0)
csv_vec <- function(x) {
  if (is.null(x) || length(x) == 0) return(character(0))
  s <- trimws(as.character(x[1]))
  if (!nzchar(s)) return(character(0))
  parts <- strsplit(s, "[,;]+")[[1]]
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

num_vec <- function(x, default) {
  v <- suppressWarnings(as.numeric(csv_vec(x)))
  v <- v[!is.na(v)]
  if (length(v) == 0) default else v
}

# ---------------------------------------------------------------------------
# 日志：同时写 stdout 和可选的 --log_file
# ---------------------------------------------------------------------------
.LOG_FILE <- NULL
init_log <- function(path = NULL) {
  if (!is.null(path) && nzchar(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    cat("", file = path)
    .LOG_FILE <<- path
  }
}
msg <- function(...) {
  s <- paste0(...)
  cat(s, "\n", sep = "")
  if (!is.null(.LOG_FILE)) cat(s, "\n", sep = "", file = .LOG_FILE, append = TRUE)
}
rule <- function(ch = "-", n = 72) msg(strrep(ch, n))

# ---------------------------------------------------------------------------
# 基因表读取
#   三种输入，优先级：--table  >  --genes_file(+--logfc_file)  >  --genes(+--logfc)
# ---------------------------------------------------------------------------
read_table_any <- function(path, header = TRUE) {
  if (!file.exists(path)) stop("找不到输入文件: ", path)
  l1 <- readLines(path, n = 1, warn = FALSE)
  sep <- if (grepl("\t", l1)) "\t" else if (grepl(",", l1)) "," else "\t"
  df <- tryCatch(
    utils::read.delim(path, sep = sep, header = header, check.names = FALSE,
                      stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM"),
    error = function(e) NULL
  )
  if (is.null(df)) {
    df <- tryCatch(
      utils::read.csv(path, header = header, check.names = FALSE,
                      stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM"),
      error = function(e) stop("读不了输入文件: ", path, " — ", conditionMessage(e))
    )
  }
  df
}

GENE_HEADERS <- c("gene", "genes", "gene_symbol", "genesymbol", "symbol", "symbols",
                  "gene_name", "genename", "geneid", "gene_id", "id", "name", "query",
                  "protein", "protein_id", "entrez", "ensembl")
LOGFC_HEADERS <- c("logFC", "log2FC", "log2FoldChange", "logfc", "log2_fold_change",
                   "avg_log2FC", "log2fc", "fold_change", "FC", "logFC_mean")

# 基因列表文件：兼容三种写法
#   a) 每行一个基因的纯文本（可能带表头，也可能不带）
#   b) 单行逗号/分号分隔的基因串
#   c) 带表头的表格（用 --gene_col 或自动猜列）
read_gene_file <- function(path, gene_col = "") {
  raw <- readLines(path, warn = FALSE, encoding = "UTF-8")
  raw <- raw[nzchar(trimws(raw))]
  if (!length(raw)) stop("基因文件是空的: ", path)

  # b) 单行多基因
  if (length(raw) == 1 && grepl("[,\t;]", raw)) {
    return(csv_vec(gsub("\t", ",", raw)))
  }
  # a) 无逗号/制表符 -> 纯列表
  if (!grepl("[,\t]", raw[1])) {
    if (nzchar(gene_col) && !(tolower(gene_col) %in% GENE_HEADERS))
      stop("--gene_col=", gene_col, " 与纯文本基因列表冲突（该文件只有一列）")
    drop1 <- tolower(trimws(raw[1])) %in% GENE_HEADERS
    g <- trimws(raw[if (drop1) -1 else TRUE])
    return(g[nzchar(g)])
  }
  # c) 表格
  df <- read_table_any(path, header = TRUE)
  gc <- if (nzchar(gene_col)) gene_col else pick_col(df, GENE_HEADERS)
  if (is.null(gc)) {
    if (ncol(df) == 1) {
      df2 <- read_table_any(path, header = FALSE)
      return(trimws(as.character(df2[[1]])))
    }
    stop("基因文件里找不到基因名列，请用 --gene_col= 指定。现有列: ",
         paste(names(df), collapse = ", "))
  }
  g <- trimws(as.character(df[[gc]]))
  g[nzchar(g)]
}

pick_col <- function(df, candidates) {
  lc <- tolower(names(df))
  for (c0 in candidates) {
    i <- which(lc == tolower(c0))
    if (length(i)) return(names(df)[i[1]])
  }
  NULL
}

build_gene_table <- function(p) {
  # 情况 A：--table（含 gene 列，可选 logFC 列）
  if (nzchar(p$table)) {
    df <- read_table_any(p$table)
    gcol <- if (nzchar(p$gene_col)) p$gene_col else pick_col(df, GENE_HEADERS)
    if (is.null(gcol) && ncol(df) == 1) {
      # 只有一列且没有表头名 -> 当成无表头的基因列表
      df <- read_table_any(p$table, header = FALSE)
      gcol <- names(df)[1]
    }
    if (is.null(gcol)) stop("--table 里找不到基因名列，请用 --gene_col= 指定。现有列: ",
                            paste(names(df), collapse = ", "))
    genes <- trimws(as.character(df[[gcol]]))
    lcol <- if (nzchar(p$logfc_col)) p$logfc_col else pick_col(df, LOGFC_HEADERS)
    logfc <- if (is.null(lcol)) rep(NA_real_, length(genes)) else
      suppressWarnings(as.numeric(as.character(df[[lcol]])))
    keep <- nzchar(genes) & !is.na(genes)
    return(data.frame(gene = genes[keep], logFC = logfc[keep], stringsAsFactors = FALSE))
  }
  # 情况 B：--genes_file（+ 可选 --logfc_file）
  if (nzchar(p$genes_file)) {
    genes <- read_gene_file(p$genes_file, p$gene_col)
    if (nzchar(p$logfc_file)) {
      df2 <- read_table_any(p$logfc_file)
      lcol <- if (nzchar(p$logfc_col)) p$logfc_col else pick_col(df2, LOGFC_HEADERS)
      if (is.null(lcol)) stop("--logfc_file 里找不到 logFC 列，请用 --logfc_col= 指定")
      v <- suppressWarnings(as.numeric(as.character(df2[[lcol]])))
      if (length(v) != length(genes)) stop("--logfc_file 行数 (", length(v),
                                           ") 与基因数 (", length(genes), ") 不一致")
      logfc <- v
    } else {
      logfc <- rep(NA_real_, length(genes))
    }
    return(data.frame(gene = genes, logFC = logfc, stringsAsFactors = FALSE))
  }
  # 情况 C：--genes 内联（+ 可选 --logfc 内联）
  genes <- csv_vec(p$genes)
  if (length(genes) == 0) stop("必须用 --genes= / --genes_file= / --table= 之一提供基因列表")
  lf <- num_vec(p$logfc, numeric(0))
  if (length(lf) == 0) {
    logfc <- rep(NA_real_, length(genes))
  } else if (length(lf) == length(genes)) {
    logfc <- lf
  } else {
    stop("--logfc 个数 (", length(lf), ") 与 --genes 个数 (", length(genes),
         ") 不一致。两个向量必须一一对应。")
  }
  data.frame(gene = genes, logFC = logfc, stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# symbol 别名修正：只作用于「映射不上」的基因，不会动本来就合法的 symbol
# 来源：参考代码注意事项第 1 条 + 常见组蛋白/去乙酰化酶旧名
# ---------------------------------------------------------------------------
ALIAS_TABLE <- list(
  # 参考代码里明确踩坑的四个
  "CoREST"     = "RCOR1",
  "COREST"     = "RCOR1",
  "RCOR"       = "RCOR1",
  "LSD1"       = "KDM1A",
  "AOF2"       = "KDM1A",
  "MACROH2A2"  = "H2AFY2",
  "H2AX"       = "H2AFX",
  # 常见旧名 / 俗名
  "MACROH2A1"  = "H2AFY",
  "H2AZ"       = "H2AZ1",
  "H2AFZ"      = "H2AZ1",
  "H2AB1"      = "H2AB1",
  "HDAC1"      = "HDAC1",
  "SIRT1"      = "SIRT1",
  "BMI1"       = "BMI1",
  "EZH2"       = "EZH2",
  "SUZ12"      = "SUZ12",
  "G9A"        = "EHMT2",
  "GLP"        = "EHMT1",
  "SETDB1"     = "SETDB1",
  "DNMT1"      = "DNMT1",
  "TET1"       = "TET1",
  "BRD4"       = "BRD4",
  "MYC"        = "MYC",
  "PARP1"      = "PARP1",
  "ATM"        = "ATM",
  "ATR"        = "ATR",
  "TP53BP1"    = "TP53BP1",
  "53BP1"      = "TP53BP1",
  "MRE11"      = "MRE11",
  "RAD50"      = "RAD50",
  "NBS1"       = "NBN",
  "NBS"        = "NBN",
  "KAP1"       = "TRIM28",
  "TIF1B"      = "TRIM28",
  "HP1A"       = "CBX5",
  "HP1B"       = "CBX1",
  "HP1G"       = "CBX3",
  "SMARCA4"    = "SMARCA4",
  "BRG1"       = "SMARCA4",
  "SNF5"       = "SMARCB1",
  "INI1"       = "SMARCB1",
  "ARID1A"     = "ARID1A",
  "BAF250A"    = "ARID1A",
  "MED12"      = "MED12",
  "CyclinD1"   = "CCND1",
  "CDK4"       = "CDK4",
  "RB"         = "RB1",
  "p53"        = "TP53",
  "P53"        = "TP53",
  "betaactin"  = "ACTB",
  "ACTIN"      = "ACTB",
  "GAPDH"      = "GAPDH",
  "Tubulin"    = "TUBA1B"
)

apply_alias <- function(genes, mode = "auto") {
  # mode: auto（只替换未能映射者，由调用方决定）/ all（全表替换）/ none
  if (identical(mode, "none")) return(list(genes = genes, changed = character(0)))
  changed <- character(0)
  out <- genes
  for (i in seq_along(out)) {
    g <- out[i]
    if (!is.na(g) && g %in% names(ALIAS_TABLE)) {
      newg <- ALIAS_TABLE[[g]]
      if (identical(mode, "all") || identical(mode, "auto")) {
        out[i] <- newg
        changed <- c(changed, paste0(g, " -> ", newg))
      }
    }
  }
  list(genes = out, changed = changed)
}

# ---------------------------------------------------------------------------
# 边表规范化：无论哪个引擎，最终统一成
#   data.frame(from, to, combined_score)  —— from/to 是 STRING_id
# ---------------------------------------------------------------------------
normalize_edges <- function(edges, score_col = c("combined_score", "score", "combinedscore")) {
  if (is.null(edges) || nrow(edges) == 0) {
    return(data.frame(from = character(0), to = character(0),
                      combined_score = numeric(0), stringsAsFactors = FALSE))
  }
  cn <- names(edges)
  fcol <- intersect(c("from", "stringId_A", "protein1", "item_id_a"), cn)[1]
  tcol <- intersect(c("to", "stringId_B", "protein2", "item_id_b"), cn)[1]
  scol <- intersect(score_col, cn)[1]
  if (is.na(fcol) || is.na(tcol)) stop("边表缺少 from/to 列，现有列: ", paste(cn, collapse = ", "))
  if (is.na(scol)) stop("边表缺少置信度列，现有列: ", paste(cn, collapse = ", "))
  d <- data.frame(
    from = as.character(edges[[fcol]]),
    to = as.character(edges[[tcol]]),
    combined_score = suppressWarnings(as.numeric(as.character(edges[[scol]]))),
    stringsAsFactors = FALSE
  )
  d <- d[!is.na(d$from) & !is.na(d$to), , drop = FALSE]
  d <- d[d$from != d$to, , drop = FALSE]
  # 双向去重（无向图）
  if (nrow(d) > 0) {
    key <- apply(cbind(d$from, d$to), 1, function(z) paste(sort(z), collapse = "|"))
    d <- d[!duplicated(key), , drop = FALSE]
  }
  d$combined_score[is.na(d$combined_score)] <- 0
  rownames(d) <- NULL
  d
}

# ---------------------------------------------------------------------------
# 配色
# ---------------------------------------------------------------------------
build_palette <- function(name, custom = "", n = 7) {
  if (nzchar(custom)) {
    cols <- csv_vec(custom)
    if (length(cols) == 0) stop("--colors= 解析不到颜色")
    return(cols)
  }
  name <- tolower(name)
  if (name %in% c("nature", "string", "default", "orig", "original")) {
    # 参考代码原配色（蓝 -> 浅黄）
    return(c("#08306B", "#2171B5", "#1F9BCD", "#41B6C4", "#7FCDBB", "#C7E9B4", "#FFFFCC"))
  }
  hcl_name <- switch(name,
    "viridis"  = "Viridis",
    "plasma"   = "Plasma",
    "magma"    = "Magma",
    "inferno"  = "Inferno",
    "bluyl"    = "BluYl",
    "buyl"     = "BluYl",
    "rdbu"     = "RdBu",
    "bluered2" = "Blue-Red 2",
    "bluered"  = "Blue-Red 2",
    "redblue"  = "Blue-Red 2",
    "spectral" = "Spectral",
    "blues"    = "Blues",
    "reds"     = "Reds",
    "greens"   = "Greens",
    "purples"  = "Purples",
    "dark"     = "Dark 3",
    "sunset"   = "Sunset",
    name
  )
  cols <- tryCatch(grDevices::hcl.colors(n, palette = hcl_name),
                   error = function(e) NULL)
  if (is.null(cols)) {
    msg("[warn] 未知配色 '", name, "'，回退到 nature 配色。可用：nature/viridis/plasma/",
        "rdbu/bluered2/spectral/bluyl/blues/reds/greens/purples/magma/inferno")
    return(c("#08306B", "#2171B5", "#1F9BCD", "#41B6C4", "#7FCDBB", "#C7E9B4", "#FFFFCC"))
  }
  cols
}

# 把超范围的 logFC 压到边界（避免 ggplot 默认 censoring 把点变灰）
squish_oob <- function(x, range, only.finite = TRUE) {
  x[x < range[1]] <- range[1]
  x[x > range[2]] <- range[2]
  x
}

# ---------------------------------------------------------------------------
# 输出格式
# ---------------------------------------------------------------------------
# --format 归一：白名单 pdf/png/svg；非法值告警忽略；无论如何都强制 pdf+png+svg 三格式。
# （历史 bug：--format=jpg / --format=pdf,png 会原样拼进文件名，而 open_dev 的
#   else 分支对未知格式静默走 PDF 设备 -> 扩展名与内容不符的「未知格式」文件。）
normalize_formats <- function(fmt_raw, always = c("pdf", "png", "svg")) {
  req <- tolower(csv_vec(fmt_raw))
  known <- req[req %in% c("pdf", "png", "svg")]
  unknown <- setdiff(req, c("pdf", "png", "svg"))
  if (length(unknown))
    msg("[warn] 未知输出格式 '", paste(unknown, collapse = ","),
        "'，已忽略（只支持 pdf/png/svg；三种格式始终输出）")
  unique(c(always, known))
}

# ---------------------------------------------------------------------------
# 输出设备
# ---------------------------------------------------------------------------
open_dev <- function(path, fmt, width, height, dpi) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fmt <- tolower(fmt)
  if (fmt == "png") {
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(path, width = width, height = height, units = "in", res = dpi,
                    background = "white")
    } else {
      grDevices::png(path, width = width, height = height, units = "in", res = dpi,
                    type = "cairo", bg = "white")
    }
  } else if (fmt == "svg") {
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(path, width = width, height = height)
    } else {
      grDevices::svg(path, width = width, height = height)
    }
  } else if (fmt == "pdf") {
    grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
  } else {
    stop("未知输出格式 '", fmt, "'（只支持 pdf/png/svg；normalize_formats 应先过滤）")
  }
  invisible(path)
}

# 把同一个 ggplot 对象按多个格式逐一落盘，逐个登记字节数。
# 返回成功产出的文件路径向量。设备打开后无论 print 是否出错都保证关闭，
# 否则进程被 139 段错误打断时会留下残缺文件。
save_plot_formats <- function(plot, outdir, prefix, stem, formats, width, height, dpi) {
  outs <- character(0)
  for (f in formats) {
    fig <- file.path(outdir, paste0(prefix, "_", stem, ".", f))
    ok <- TRUE
    tryCatch({
      open_dev(fig, f, width, height, dpi)
      print(plot)
    }, error = function(e) {
      ok <<- FALSE
      msg("[error] 出图失败 (", f, ") : ", conditionMessage(e))
    }, finally = {
      if (grDevices::dev.cur() > 1) grDevices::dev.off()
    })
    if (ok && file.exists(fig) && file.size(fig) > 0) {
      outs <- c(outs, fig)
      msg("已出图            : ", normalizePath(fig, mustWork = FALSE),
          "  (", format(file.size(fig), big.mark = ","), " bytes)")
    } else {
      msg("[error] 产物缺失或为空: ", fig)
    }
  }
  outs
}

# ---------------------------------------------------------------------------
# 落盘：参数 + 报告
# ---------------------------------------------------------------------------
write_run_params <- function(path, p) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- c("key\tvalue")
  for (k in sort(names(p))) {
    if (!nzchar(k) || identical(k, "help_flag")) next
    lines <- c(lines, paste0(k, "\t", paste(as.character(p[[k]]), collapse = ",")))
  }
  writeLines(lines, path)
}
