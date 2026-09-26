## ===========================================================================
## lib_venn_common.R —— Venn 图技能的公共工具（被 01_venn.R source）
##
## 本文件不单独运行。设计原则：绝不吞掉错误、所有输入输出都打印给人看。
##
## 与 enrichment-analysis 的 lib_enrich_common.R 同源，差别只有一处：
## 这里的参数解析**支持同一个 key 出现多次**（--set 要传 2~5 个集合），
## 重复出现时保留全部值（oa_all），取标量则取最后一次出现的值（oa）。
## ===========================================================================

## ---- 脚本自身位置 ---------------------------------------------------------

get_script_dir <- function() {
  a <- commandArgs(FALSE)
  i <- grep("^--file=", a)
  if (length(i)) {
    return(dirname(normalizePath(sub("^--file=", "", a[i[1]]), mustWork = FALSE)))
  }
  j <- grep("\\.[Rr]$", a)
  if (length(j)) {
    return(dirname(normalizePath(a[j[1]], mustWork = FALSE)))
  }
  getwd()
}

locate_lib <- function() {
  cand <- character(0)
  d <- Sys.getenv("VENN_SKILL_SCRIPTS", "")
  if (nzchar(d)) cand <- c(cand, file.path(d, "lib_venn_common.R"))
  cand <- c(cand,
            file.path(get_script_dir(), "lib_venn_common.R"),
            file.path(getwd(), "lib_venn_common.R"))
  cand <- unique(cand)
  for (p in cand) if (file.exists(p)) return(normalizePath(p, mustWork = FALSE))
  stop("找不到 lib_venn_common.R。已试路径：\n  ", paste(cand, collapse = "\n  "),
       "\n请用 run_venn.sh 启动本技能脚本，或把环境变量 VENN_SKILL_SCRIPTS 指向 scripts/ 目录。",
       call. = FALSE)
}

## ---- 命令行参数 -----------------------------------------------------------
## 支持 --key=value / --key value / 裸 --flag（视为 TRUE）。
## 同一个 key 出现多次 -> 值累积成向量（--set 要用）。

parse_args <- function(argv = NULL) {
  if (is.null(argv)) argv <- commandArgs(trailingOnly = TRUE)
  out <- list()
  push <- function(k, v) {
    if (is.null(out[[k]])) out[[k]] <<- v else out[[k]] <<- c(out[[k]], v)
  }
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (grepl("^--[^=]+=", a)) {
      k <- sub("^--", "", sub("=.*$", "", a))
      push(k, sub("^[^=]*=", "", a))
    } else if (grepl("^--?[A-Za-z]", a) && !grepl("^-?[0-9.]", a)) {
      k <- sub("^--?", "", a)
      nxt <- if (i < length(argv)) argv[[i + 1L]] else NA_character_
      if (!is.na(nxt) && !grepl("^--[A-Za-z]", nxt) && !grepl("^-?[0-9.]", nxt)) {
        push(k, nxt); i <- i + 1L
      } else {
        push(k, "TRUE")
      }
    }
    i <- i + 1L
  }
  out
}

## 标量：取最后一次出现的值
oa <- function(o, k, default = NULL) {
  v <- o[[k]]
  if (is.null(v) || !length(v)) return(default)
  v <- v[length(v)]
  if (!nzchar(v)) default else v
}
## 向量：全部出现的值
oa_all <- function(o, k) {
  v <- o[[k]]
  if (is.null(v)) character(0) else v
}
oa_num <- function(o, k, default) {
  v <- oa(o, k, NULL)
  if (is.null(v)) return(default)
  x <- suppressWarnings(as.numeric(v))
  if (is.na(x)) stop("参数 --", k, " 需要是数字，收到：", v, call. = FALSE)
  x
}
oa_int <- function(o, k, default) {
  v <- oa(o, k, NULL)
  if (is.null(v)) return(default)
  x <- suppressWarnings(as.integer(v))
  if (is.na(x)) stop("参数 --", k, " 需要是整数，收到：", v, call. = FALSE)
  x
}
oa_lgl <- function(o, k, default) {
  v <- oa(o, k, NULL)
  if (is.null(v)) return(default)
  toupper(v) %in% c("TRUE", "T", "1", "YES", "Y")
}
## 逗号/分号分隔 -> 向量
oa_vec <- function(o, k, default = NULL) {
  v <- oa(o, k, NULL)
  if (is.null(v)) return(default)
  if (!nzchar(v)) return(default)
  x <- trimws(unlist(strsplit(v, "[,;]")))
  x <- x[nzchar(x)]
  if (!length(x)) default else x
}
oa_numvec <- function(o, k, default = NULL) {
  v <- oa_vec(o, k, NULL)
  if (is.null(v)) return(default)
  x <- suppressWarnings(as.numeric(v))
  if (any(is.na(x))) stop("参数 --", k, " 需要是数字列表（逗号分隔），收到：",
                          oa(o, k), call. = FALSE)
  x
}

## ---- 日志 ----------------------------------------------------------------
## 一律写到 stderr，避免和产物信息混在一起；同时可被 --log 抓进文件里。
.venn_logfile <- NULL
set_logfile <- function(p) {
  .venn_logfile <<- p
  cat("", file = p, append = FALSE)
}
say <- function(...) {
  s <- paste0(...)
  cat(s, "\n", sep = "")
  if (!is.null(.venn_logfile)) cat(s, "\n", sep = "", file = .venn_logfile, append = TRUE)
  invisible(s)
}
warn <- function(...) say("[警告] ", ...)
step <- function(...) say("[venn] ", ...)

## ---- 配色 ----------------------------------------------------------------
## 1 号 = 参考代码 venn.R 用的那套（紫 + 橙）打头，后面补几个互补色给 3~5 组用
VENN_PALS <- list(
  "1" = c("#6e48fb", "#ffa500", "#00A087", "#3C5488", "#F39B7F"),      # 参考代码紫/橙
  "2" = c("#3C5488", "#E64B35", "#00A087", "#4DBBD5", "#F39B7F"),      # npg
  "3" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),      # Set1
  "4" = c("#0099B4", "#925E9F", "#FDAF91", "#AD002A", "#42B540")       # NEJM 风
)
pal_names <- function() {
  paste0(names(VENN_PALS), " (", vapply(VENN_PALS, function(v) paste(v, collapse = " "),
                                        character(1)), ") = ",
         c("参考代码紫/橙", "npg", "Set1", "NEJM"), collapse = "\n  ")
}

## ---- 读基因列表 -----------------------------------------------------------
## 三种扩展名三种读法；返回**字符向量**（已去空、去重前的原始顺序）
read_gene_list <- function(path, col = NULL, verbose = TRUE) {
  if (!file.exists(path)) stop("基因列表文件不存在：", path, call. = FALSE)
  ext <- tolower(sub("^.*\\.", "", basename(path)))
  if (ext %in% c("csv", "tsv", "tab") && !is.null(col)) {
    ## 表格：按列名取
    sep <- if (ext == "csv") "," else "\t"
    tab <- utils::read.csv(path, sep = sep, header = TRUE,
                           check.names = FALSE, stringsAsFactors = FALSE,
                           fileEncoding = "UTF-8-BOM")
    if (!col %in% names(tab)) {
      stop("文件 ", basename(path), " 里没有列 ", col, "；可选列：",
           paste(names(tab), collapse = ", "), call. = FALSE)
    }
    g <- as.character(tab[[col]])
    if (verbose) step("读入 ", basename(path), " 的列 [", col, "]：", length(g), " 个值")
  } else if (ext %in% c("csv", "tsv", "tab")) {
    sep <- if (ext == "csv") "," else "\t"
    tab <- utils::read.csv(path, sep = sep, header = TRUE,
                           check.names = FALSE, stringsAsFactors = FALSE,
                           fileEncoding = "UTF-8-BOM")
    if (ncol(tab) == 0) stop("文件读不出内容：", path, call. = FALSE)
    g <- as.character(tab[[1]])
    if (verbose) step("读入 ", basename(path), " 的第 1 列 [", names(tab)[1], "]：",
                      length(g), " 个值")
  } else {
    ## txt：一行一个，或用逗号/分号/空白分隔都行
    ln <- readLines(path, warn = FALSE)
    ln[1] <- sub("^\ufeff", "", ln[1])                       # 去 BOM
    ln <- gsub("\r$", "", ln)
    if (length(ln) && grepl("^(gene|genes|symbol|gene_symbol|id|geneid)$",
                            tolower(trimws(ln[1])))) {
      if (verbose) step("跳过首行表头：", trimws(ln[1]))
      ln <- ln[-1]
    }
    g <- unlist(strsplit(ln, "[,;\t ]+"))
    if (verbose) step("读入 ", basename(path), "（纯文本）：", length(g), " 个值")
  }
  g <- trimws(g)
  g <- g[nzchar(g) & !is.na(g) & toupper(g) != "NA"]
  g
}

## ---- 集合整理 -------------------------------------------------------------
## sets: named list of character vectors
## 返回整理后的 named list，并打印每个集合的去重前后大小
clean_sets <- function(sets, force_unique = TRUE) {
  out <- list()
  for (nm in names(sets)) {
    g <- sets[[nm]]
    raw <- length(g)
    if (force_unique) g <- unique(g)
    out[[nm]] <- g
    if (force_unique && raw != length(g)) {
      step(sprintf("集合 [%s]：%d 个值 -> 去重后 %d", nm, raw, length(g)))
    } else {
      step(sprintf("集合 [%s]：%d 个基因", nm, length(g)))
    }
  }
  out
}

## ---- 各区交集成员 ---------------------------------------------------------
## 返回 data.frame(region, idx, n, genes)，只含 n>0 的区
## region 形如 "A&B"，idx 形如 "12"（便于排序）
region_table <- function(sets) {
  n <- length(sets)
  ## 展示/导出用的集合名：标签里的换行（\n）换成 "/"，否则会把 CSV 行拆断
  nm <- gsub("\n", "/", names(sets))
  allg <- unique(unlist(sets, use.names = FALSE))
  if (!length(allg)) stop("所有集合都是空的，没有东西可画。", call. = FALSE)
  ## 每个基因属于哪几个集合
  code <- vapply(allg, function(g) {
    hit <- which(vapply(sets, function(s) g %in% s, logical(1)))
    paste(hit, collapse = ",")
  }, character(1), USE.NAMES = FALSE)
  combos <- sort(unique(code))
  rows <- lapply(combos, function(cb) {
    idx <- as.integer(strsplit(cb, ",")[[1]])
    data.frame(region = paste(nm[idx], collapse = "&"),
               idx = paste(idx, collapse = ""),
               n = sum(code == cb),
               genes = paste(allg[code == cb], collapse = ";"),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df <- df[order(nchar(df$idx), df$idx), , drop = FALSE]
  rownames(df) <- NULL
  df
}

## ---- 写 CSV（带 BOM，Excel 打开中文不乱码）--------------------------------
write_csv_bom <- function(df, path) {
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("\xef\xbb\xbf"), con)
  utils::write.csv(df, file = con, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(path)
}
