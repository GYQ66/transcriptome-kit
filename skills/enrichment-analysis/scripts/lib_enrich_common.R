## ===========================================================================
## lib_enrich_common.R —— 富集分析技能的公共工具（被 01 / 02 / 03 三个脚本 source）
##
## 本文件不单独运行。三个脚本通过 locate_lib() 找到并加载它。
## 设计原则：绝不吞掉错误、绝不替用户做生物学判断、所有输入输出都打印给人看。
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
  d <- Sys.getenv("ENRICH_SKILL_SCRIPTS", "")
  if (nzchar(d)) cand <- c(cand, file.path(d, "lib_enrich_common.R"))
  cand <- c(cand,
            file.path(get_script_dir(), "lib_enrich_common.R"),
            file.path(getwd(), "lib_enrich_common.R"))
  cand <- unique(cand)
  for (p in cand) if (file.exists(p)) return(normalizePath(p, mustWork = FALSE))
  stop("找不到 lib_enrich_common.R。已试路径：\n  ", paste(cand, collapse = "\n  "),
       "\n请用 run_enrich.sh 启动本技能脚本，或把环境变量 ENRICH_SKILL_SCRIPTS 指向 scripts/ 目录。",
       call. = FALSE)
}

## ---- 命令行参数 -----------------------------------------------------------
## 支持 --key=value / --key value / 裸 --flag（视为 TRUE）。
## 返回一个 named list，全部是字符串；取值用下面的 oa_* 系列。

parse_args <- function(argv = NULL) {
  if (is.null(argv)) argv <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (grepl("^--[^=]+=", a)) {
      k <- sub("^--", "", sub("=.*$", "", a))
      out[[k]] <- sub("^[^=]*=", "", a)
    } else if (grepl("^--?[A-Za-z]", a) && !grepl("^-?[0-9.]", a)) {
      k <- sub("^--?", "", a)
      nxt <- if (i < length(argv)) argv[[i + 1L]] else NA_character_
      if (!is.na(nxt) && !grepl("^--[A-Za-z]", nxt) && !grepl("^-?[0-9.]", nxt)) {
        out[[k]] <- nxt; i <- i + 1L
      } else {
        out[[k]] <- "TRUE"
      }
    } else {
      out[[".positional"]] <- c(out[[".positional"]], a)
    }
    i <- i + 1L
  }
  out
}

oa <- function(o, k, default = NULL) {
  v <- o[[k]]
  if (is.null(v) || !nzchar(v)) default else v
}
oa_num <- function(o, k, default) {
  v <- oa(o, k, NULL)
  if (is.null(v)) return(default)
  x <- suppressWarnings(as.numeric(v))
  if (is.na(x)) stop("参数 --", k, " 需要是数字，收到：", v, call. = FALSE)
  x
}
oa_lgl <- function(o, k, default) {
  v <- oa(o, k, NULL)
  if (is.null(v)) return(default)
  toupper(v) %in% c("TRUE", "T", "1", "YES", "Y")
}
## 逗号分隔的多值参数；未给则返回 default
oa_list <- function(o, k, default = NULL) {
  v <- oa(o, k, NULL)
  if (is.null(v)) return(default)
  v <- trimws(unlist(strsplit(v, ",")))
  v[nzchar(v)]
}

## ---- 物种 -----------------------------------------------------------------

SPECIES_TABLE <- list(
  human = list(orgdb = "org.Hs.eg.db", kegg = "hsa", label = "人类",
               alias = c("human", "hsa", "hs", "homo_sapiens", "homosapiens",
                         "human(hs)", "人", "人类", "homo sapiens")),
  mouse = list(orgdb = "org.Mm.eg.db", kegg = "mmu", label = "小鼠",
               alias = c("mouse", "mmu", "mm", "mus_musculus", "musmusculus",
                         "小鼠", "鼠", "mus musculus"))
)

resolve_species <- function(x, orgdb_override = NULL, kegg_override = NULL) {
  x0 <- tolower(trimws(x))
  key <- NULL
  for (nm in names(SPECIES_TABLE)) {
    if (x0 == nm || x0 %in% tolower(SPECIES_TABLE[[nm]]$alias)) { key <- nm; break }
  }
  if (is.null(key)) {
    stop("不认识的物种：", x, "\n本技能支持：",
         paste(sprintf("%s(%s)", names(SPECIES_TABLE),
                       vapply(SPECIES_TABLE, function(s) s$orgdb, "")),
               collapse = "、"),
         call. = FALSE)
  }
  s <- SPECIES_TABLE[[key]]
  if (!is.null(orgdb_override) && nzchar(orgdb_override)) s$orgdb <- orgdb_override
  if (!is.null(kegg_override) && nzchar(kegg_override)) s$kegg <- kegg_override
  s$key <- key
  s
}

## OrgDb 的 sqlite 兜底搜索路径。
## 本机没有 Rtools，`R CMD INSTALL` 装不了源码注释包（见 references/pitfalls.md 5.2），
## 所以支持"不装包，直接 loadDb() 一个现成的 sqlite"。
default_orgdb_sqlite_dirs <- function() {
  cand <- c(
    Sys.getenv("ENRICH_ORGDB_DIR", ""),
    file.path(get_script_dir(), "..", "data", "orgdb"),
    file.path(Sys.getenv("HOME", Sys.getenv("USERPROFILE", "")), "orgdb")
  )
  unique(cand[nzchar(cand)])
}

.find_orgdb_sqlite <- function(orgdb, extra = NULL) {
  fn <- paste0(orgdb, ".sqlite")
  cand <- character(0)
  if (!is.null(extra) && nzchar(extra)) {
    cand <- c(cand, if (dir.exists(extra)) file.path(extra, fn) else extra)
  }
  for (d in default_orgdb_sqlite_dirs()) cand <- c(cand, file.path(d, fn))
  for (p in cand) if (file.exists(p)) return(p)
  NULL
}

load_orgdb <- function(orgdb, sqlite = NULL) {
  if (requireNamespace(orgdb, quietly = TRUE)) {
    ok <- suppressPackageStartupMessages(
      require(orgdb, character.only = TRUE, quietly = TRUE)
    )
    if (ok) return(eval(parse(text = orgdb)))
  }
  ## 退路：直接加载 sqlite（AnnotationDbi::loadDb 会按 sqlite 里的元数据还原成 OrgDb）
  db_file <- .find_orgdb_sqlite(orgdb, sqlite)
  if (!is.null(db_file)) {
    if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
      stop("需要 AnnotationDbi 才能加载 ", basename(db_file), call. = FALSE)
    }
    message("[load_orgdb] 包里没有 ", orgdb, "，改为直接加载 sqlite：", db_file)
    db <- AnnotationDbi::loadDb(db_file)
    return(db)
  }
  stop("缺少物种注释包 ", orgdb, "，也没找到可用的 ", orgdb,
       ".sqlite。\n两条路：\n",
       "  A) 装包（需联网；本机无 Rtools，源码安装会报 cmd 语法错误）：\n",
       "     Rscript -e \"BiocManager::install('", orgdb, "', ask=FALSE, update=FALSE)\"\n",
       "  B) 不装包，只下 sqlite：下载\n",
       "     https://bioconductor.org/packages/3.22/data/annotation/src/contrib/", orgdb,
       "_<版本>.tar.gz\n",
       "     解出 inst/extdata/", orgdb, ".sqlite，放到下列任一位置（或用 --orgdb_sqlite= 指定）：\n  ",
       paste(default_orgdb_sqlite_dirs(), collapse = "\n  "),
       call. = FALSE)
}

## ---- 读表 -----------------------------------------------------------------

read_table_robust <- function(path, what = "输入表") {
  if (is.null(path) || !nzchar(path)) stop("没有提供", what, "路径", call. = FALSE)
  if (!file.exists(path)) stop("找不到", what, "：", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("txt", "tsv")) {
    df <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                            quote = "", comment.char = "", na.strings = c("NA", ""))
  } else {
    df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                          quote = "\"", comment.char = "", na.strings = c("NA", ""))
  }
  if (nrow(df) == 0) stop(what, " 是空表：", path, call. = FALSE)
  ## write.csv(row.names = TRUE) 产出的行名列列名是空字符串（表头形如 "","logFC"...）。
  ## R 里 df[[""]] 取不到东西，但 df[["rownames1"]] 可以 —— 统一改名，后续按名访问。
  nm <- names(df)
  bad <- !nzchar(trimws(nm))
  if (any(bad)) names(df)[bad] <- paste0("rownames", which(bad))
  df
}

## 自动找基因名所在列。
## 优先级：显式指定 > 行名列（rownames1）> 常见基因列名 > 第一个"有非空取值"的字符列。
## 返回列名（一定能在 df[[...]] 里取到）。
pick_gene_col <- function(df, gene_col = NULL) {
  nm <- names(df)
  if (!is.null(gene_col) && nzchar(gene_col)) {
    i <- which(tolower(nm) == tolower(gene_col))
    if (!length(i)) {
      stop("--gene_col=", gene_col, " 不在表里。可用列名：",
           paste(nm, collapse = ", "), call. = FALSE)
    }
    return(nm[i[1]])
  }
  cand <- c(grep("^rownames[0-9]*$", nm, value = TRUE),
            nm[tolower(nm) %in% c("gene", "genes", "gene_name", "genename", "symbol",
                                  "gene_symbol", "genesymbol", "geneid", "gene_id",
                                  "entrezid", "entrez_id", "id", "x")])
  cand <- c(cand, nm)
  cand <- unique(cand)
  for (c in cand) {
    v <- df[[c]]
    if (is.null(v)) next
    v <- suppressWarnings(as.character(v))
    if (sum(!is.na(v) & nzchar(v)) > 0) return(c)
  }
  stop("表里找不到任何可用作基因名的列，请用 --gene_col 指定。\n可用列名：",
       paste(nm, collapse = ", "), call. = FALSE)
}

## 按 candidates 的**优先级顺序**（不是表里的列顺序）挑列。
## 例如 --ora_p_type=auto 的候选是 c("adj.P.Val","P.Value",...)，
## 表里同时有 P.Value 和 adj.P.Val 时必须挑 adj.P.Val。
pick_col <- function(df, col = NULL, candidates = character(0), what = "列") {
  nm <- names(df)
  if (!is.null(col) && nzchar(col)) {
    i <- which(tolower(nm) == tolower(col))
    if (!length(i)) {
      stop("找不到", what, " '", col, "'。可用列名：", paste(nm, collapse = ", "),
           call. = FALSE)
    }
    return(nm[i[1]])
  }
  for (cd in candidates) {
    i <- which(tolower(nm) == tolower(cd))
    if (length(i)) return(nm[i[1]])
    ## 容错：列名里含这个关键词（如 "adj.P.Val.Drug"）。
    ## 只对长度 >= 4 的候选做模糊匹配 —— "t"/"id" 这类短词模糊匹配会乱抓列。
    if (nchar(cd) >= 4) {
      i <- which(grepl(tolower(cd), tolower(nm), fixed = TRUE))
      if (length(i)) return(nm[i[1]])
    }
  }
  NA_character_
}

## ---- 写出 -----------------------------------------------------------------

dir_create <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(p)) stop("无法创建目录（常见原因是路径含非 ASCII 字符时的 locale 问题）：",
                           p, call. = FALSE)
  invisible(p)
}

write_csv_out <- function(df, path, row_names = TRUE) {
  dir_create(dirname(path))
  utils::write.csv(df, path, row.names = row_names, fileEncoding = "UTF-8")
  invisible(path)
}

## ---- 输出通道：同时打到屏幕与日志文件 --------------------------------------
## 用法：LOG <- start_log(path); note("...")

.start_log <- new.env(parent = emptyenv())

start_log <- function(path = NULL) {
  if (!is.null(path) && nzchar(path)) {
    dir_create(dirname(path))
    .start_log$con <- file(path, open = "wt", encoding = "UTF-8")
  } else {
    .start_log$con <- NULL
  }
  invisible(path)
}

close_log <- function() {
  if (!is.null(.start_log$con)) { try(close(.start_log$con), silent = TRUE); .start_log$con <- NULL }
  invisible(NULL)
}

note <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  if (!is.null(.start_log$con)) {
    try(writeLines(msg, .start_log$con), silent = TRUE)
  }
  invisible(NULL)
}

## ---- 配色 -----------------------------------------------------------------
## 三套配色都来自参考代码 kegg美化.R。键名顺序与参考代码
## levels = rev(c('BP','CC','MF','KEGG')) 一致，即 pal[1]->KEGG, [2]->MF, [3]->CC, [4]->BP。
## 用「命名向量」而不是裸向量，是为了让每个分类的颜色不随分类数量变化而漂移。
PALETTES <- list(
  pal1 = c("#7bc4e2", "#acd372", "#fbb05b", "#ed6ca4"),
  pal2 = c("#eaa052", "#b74147", "#90ad5b", "#23929c"),
  pal3 = c("#c3e1e6", "#f3dfb7", "#dcc6dc", "#96c38e")
)
CANON_CATS <- c("KEGG", "MF", "CC", "BP")
## 超出 4 个分类（例如叠加了多个 MSigDB 集合）时的备用色
EXTRA_CATS_PAL <- c("#8FA8C8", "#C8A88F", "#9FBF9F", "#BFA3C8", "#C8C88F", "#8FC8C0")

build_pal <- function(cats, pal = NULL, pal_idx = 3) {
  if (is.null(pal) || !length(pal)) {
    nm <- names(PALETTES)[min(max(pal_idx, 1), length(PALETTES))]
    base <- PALETTES[[nm]]
  } else {
    base <- pal
  }
  base_named <- if (!is.null(names(base)) && any(nzchar(names(base)))) {
    base
  } else {
    stats::setNames(base[seq_along(CANON_CATS)], CANON_CATS)
  }
  out <- character(0)
  extra_i <- 0L
  for (c in cats) {
    if (c %in% names(base_named)) {
      out[[c]] <- base_named[[c]]
    } else if (c %in% CANON_CATS) {
      out[[c]] <- PALETTES$pal3[[match(c, CANON_CATS)]]
    } else {
      extra_i <- extra_i + 1L
      out[[c]] <- EXTRA_CATS_PAL[((extra_i - 1L) %% length(EXTRA_CATS_PAL)) + 1L]
    }
  }
  out
}

## ---- MSigDB GMT -----------------------------------------------------------

## 默认 GMT 搜索目录（按顺序找第一个存在的）：
##   1. 环境变量 ENRICH_GMT_DIR
##   2. 本技能的 data/（套件自带 MSigDB GMT）
##   3. 同级套件里的 msigdb-pathway-gene-lookup/data
default_gmt_dirs <- function() {
  cand <- c(
    Sys.getenv("ENRICH_GMT_DIR", ""),
    file.path(get_script_dir(), "..", "data"),
    file.path(get_script_dir(), "..", "..", "msigdb-pathway-gene-lookup", "data")
  )
  cand <- unique(cand[nzchar(cand)])
  cand[vapply(cand, dir.exists, TRUE)]
}

## 语义集合名 -> 匹配哪个 GMT 文件 + 保留哪些 term 前缀
GMT_SEMANTIC <- list(
  H = list(file = "^h\\.", filter = NULL),
  C1 = list(file = "^c1\\.", filter = NULL),
  C2 = list(file = "^c2\\.", filter = NULL),
  C3 = list(file = "^c3\\.", filter = NULL),
  C4 = list(file = "^c4\\.", filter = NULL),
  C5 = list(file = "^c5\\.", filter = NULL),
  C6 = list(file = "^c6\\.", filter = NULL),
  C7 = list(file = "^c7\\.", filter = NULL),
  C8 = list(file = "^c8\\.", filter = NULL),
  C9 = list(file = "^c9\\.", filter = NULL),
  KEGG = list(file = "^c2\\.", filter = "^KEGG_"),
  REACTOME = list(file = "^c2\\.", filter = "^REACTOME_"),
  WP = list(file = "^c2\\.", filter = "^WP_"),
  BIOCARTA = list(file = "^c2\\.", filter = "^BIOCARTA_"),
  PID = list(file = "^c2\\.", filter = "^PID_"),
  GOBP = list(file = "^c5\\.", filter = "^GOBP_"),
  GOCC = list(file = "^c5\\.", filter = "^GOCC_"),
  GOMF = list(file = "^c5\\.", filter = "^GOMF_"),
  IMMUNE = list(file = "^c7\\.", filter = NULL),
  CELLTYPE = list(file = "^c8\\.", filter = NULL),
  PERTURB = list(file = "^c9\\.", filter = NULL),
  ONCOGENIC = list(file = "^c6\\.", filter = NULL)
)
GMT_SEMANTIC$HALLMARK <- GMT_SEMANTIC$H
GMT_SEMANTIC$WIKIPATHWAYS <- GMT_SEMANTIC$WP
GMT_SEMANTIC$IMMUNESIGDB <- GMT_SEMANTIC$IMMUNE
GMT_SEMANTIC$CELLTYPES <- GMT_SEMANTIC$CELLTYPE
GMT_SEMANTIC$MARKER <- GMT_SEMANTIC$CELLTYPE

read_gmt_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines))]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  ## GMT 行格式：term <TAB> url <TAB> gene1 <TAB> gene2 ...
  ## 必须按行长度 rep() 展开，不能直接 unlist —— 有"只有 2 个字段"的行时
  ## data.frame(term=..., gene=...) 会因为行数不一致报
  ## "arguments imply differing number of rows"。
  parts <- parts[vapply(parts, length, 1L) >= 3L]
  if (!length(parts)) stop("GMT 文件里没有可解析的基因集：", path, call. = FALSE)
  term <- vapply(parts, `[`, "", 1L)
  gene <- unlist(lapply(parts, function(p) p[-c(1L, 2L)]), use.names = FALSE)
  ns   <- vapply(parts, function(p) length(p) - 2L, 1L)
  data.frame(term = rep(term, ns), gene = gene, stringsAsFactors = FALSE)
}

## 把用户写的集合名解析成 TERM2GENE 表。
## tokens: c("H","C2") 或某个具体 .gmt 文件路径
resolve_gmt <- function(tokens, gmt_dirs = NULL, gmt_files = NULL) {
  if (is.null(gmt_dirs) || !length(gmt_dirs)) gmt_dirs <- default_gmt_dirs()
  out <- list()
  add <- function(coll, t2g) {
    if (is.null(t2g) || !nrow(t2g)) return(invisible(NULL))
    t2g$collection <- coll
    out[[length(out) + 1L]] <<- t2g
    invisible(NULL)
  }
  ## 1) 显式给出的 .gmt 文件
  for (f in gmt_files) {
    if (!file.exists(f)) stop("找不到 GMT 文件：", f, call. = FALSE)
    coll <- sub("\\.gmt$", "", basename(f), ignore.case = TRUE)
    add(coll, read_gmt_file(f))
  }
  ## 2) 语义集合名
  for (tok in tokens) {
    if (file.exists(tok)) { add(sub("\\.gmt$", "", basename(tok), ignore.case = TRUE), read_gmt_file(tok)); next }
    up <- toupper(tok)
    if (!up %in% names(GMT_SEMANTIC)) {
      stop("不认识的集合名：", tok, "\n可用：",
           paste(sort(names(GMT_SEMANTIC)), collapse = ", "),
           "\n也可以直接给 .gmt 文件路径。", call. = FALSE)
    }
    spec <- GMT_SEMANTIC[[up]]
    if (!length(gmt_dirs)) {
      stop("要按集合名（", tok, "）取 GMT，但没有找到 GMT 目录。\n",
           "请用 --gmt_dir=<含 *.gmt 的目录> 指定，或直接用 --gmt=<文件路径>。",
           call. = FALSE)
    }
    files <- unlist(lapply(gmt_dirs, function(d) list.files(d, pattern = "\\.gmt$",
                                                             full.names = TRUE,
                                                             ignore.case = TRUE)),
                    use.names = FALSE)
    files <- files[grepl(spec$file, tolower(basename(files)))]
    if (!length(files)) {
      stop("在 GMT 目录里找不到集合 ", tok, " 对应的文件（应匹配 '", spec$file, "'）。\n",
           "已搜索目录：", paste(gmt_dirs, collapse = ", "), call. = FALSE)
    }
    files <- sort(files, decreasing = TRUE)   # 版本号大的优先
    f <- files[1]
    t2g <- read_gmt_file(f)
    if (!is.null(spec$filter)) t2g <- t2g[grepl(spec$filter, t2g$term), , drop = FALSE]
    if (!nrow(t2g)) {
      stop("集合 ", tok, " 过滤后为空（filter='", spec$filter, "'），检查 GMT 文件：", f,
           call. = FALSE)
    }
    note(sprintf("  [GMT] %-10s <- %s  (%d 个基因集, %d 条基因条目)",
                 tok, basename(f), length(unique(t2g$term)), nrow(t2g)))
    add(up, t2g)
  }
  if (!length(out)) return(NULL)
  res <- do.call(rbind, out)
  res <- res[!duplicated(paste(res$collection, res$term, res$gene)), , drop = FALSE]
  res
}

## 判定 GMT 里的基因是 symbol 还是 entrez
gmt_id_kind <- function(t2g, symbols, entrez) {
  g <- unique(t2g$gene)
  n <- min(length(g), 20000L)
  g <- g[seq_len(n)]
  sym_hit <- mean(g %in% symbols)
  ent_hit <- mean(g %in% entrez)
  if (sym_hit >= ent_hit) "symbol" else "entrez"
}

## ---- 富集结果的统一整理 ---------------------------------------------------

## 把 enrichResult 变成 data.frame，并补齐 ID / ONTOLOGY 列
as_result_df <- function(obj, ontology = NULL) {
  if (is.null(obj)) return(NULL)
  df <- try(as.data.frame(obj), silent = TRUE)
  if (inherits(df, "try-error") || is.null(df) || !nrow(df)) return(NULL)
  if (!"ID" %in% names(df) && !is.null(rownames(df))) df$ID <- rownames(df)
  if (!"ONTOLOGY" %in% names(df)) df$ONTOLOGY <- ontology %||% NA_character_
  df$ONTOLOGY[is.na(df$ONTOLOGY)] <- ontology %||% NA_character_
  if (!is.null(ontology) && !is.na(ontology)) df$ONTOLOGY <- ontology
  df
}
`%||%` <- function(a, b) if (is.null(a)) b else a

## 打上 sig 标记（口径：pvalue / p.adjust / qvalue 三条同时满足）
mark_sig <- function(df, pvalue, padj, qvalue) {
  if (is.null(df) || !nrow(df)) return(df)
  pv  <- suppressWarnings(as.numeric(df$pvalue))
  pad <- suppressWarnings(as.numeric(df$p.adjust))
  qv  <- suppressWarnings(as.numeric(df$qvalue))
  ok <- rep(TRUE, nrow(df))
  if (!is.na(pvalue)) ok <- ok & !is.na(pv) & pv <= pvalue
  if (!is.na(padj))   ok <- ok & !is.na(pad) & pad <= padj
  if (!is.na(qvalue)) ok <- ok & (is.na(qv) | qv <= qvalue)
  df$sig <- ok
  df
}

## 把 geneID（entrez，以 / 分隔）换成 symbol，另存一列 geneName
add_gene_name <- function(df, sym_map) {
  if (is.null(df) || !nrow(df) || !"geneID" %in% names(df)) return(df)
  conv <- function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_character_)
    ids <- strsplit(x, "/", fixed = TRUE)[[1]]
    hit <- sym_map[ids]
    hit[is.na(hit)] <- ids[is.na(hit)]
    paste(hit, collapse = "/")
  }
  df$geneName <- vapply(as.character(df$geneID), conv, "", USE.NAMES = FALSE)
  df
}

## 统一的列顺序（对齐参考代码输出的 CSV，末尾加 geneName / sig）
order_cols <- function(df) {
  want <- c("ONTOLOGY", "collection", "ID", "Description", "GeneRatio", "BgRatio",
            "RichFactor", "FoldEnrichment", "zScore", "pvalue", "p.adjust",
            "qvalue", "geneID", "geneName", "Count", "sig")
  keep <- c(intersect(want, names(df)), setdiff(names(df), want))
  df[, keep, drop = FALSE]
}

## ---- 与 geo-microarray-analysis 技能的对接 ---------------------------------
## geo 的 02_deg_plots.R 会写一份 <prefix>_for_enrichment.txt（key=value 文本）。
## 有了它，富集侧只要 --geo_dir=<geo 输出目录> [--geo_prefix=<prefix>] 就能自动
## 认出 ORA / GSEA 该吃哪个文件、列名叫什么，不用手写 8 个参数。
## 读不到清单时退化为按文件名约定找 <prefix>_DEG.csv / <prefix>_all.csv。

GEO_MANIFEST_SUFFIX <- "_for_enrichment.txt"

read_kv_manifest <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!grepl("^\\s*#", ln)]
  ln <- ln[nzchar(trimws(ln)) & grepl("=", ln, fixed = TRUE)]
  if (!length(ln)) return(NULL)
  k <- trimws(sub("=.*$", "", ln))
  v <- trimws(sub("^[^=]*=", "", ln))
  keep <- nzchar(k)
  stats::setNames(as.list(v[keep]), k[keep])
}

## 解析 geo 交接信息。
## 返回 list(ok, source, dir, prefix, deg, all, p_col, logfc_col, change_col,
##           contrast, species, groups, n_per_group, ...)
resolve_geo_inputs <- function(geo_dir = NULL, geo_prefix = NULL, manifest = NULL,
                               verbose = TRUE) {
  ## 允许把目录一起塞进 --geo_prefix（如 D:/x/deg/GSE62452_T）
  if ((is.null(geo_dir) || !nzchar(geo_dir)) && !is.null(geo_prefix) && nzchar(geo_prefix)) {
    if (grepl("[/\\\\]", geo_prefix)) {
      geo_dir <- dirname(geo_prefix)
      geo_prefix <- basename(geo_prefix)
    } else {
      geo_dir <- "."     # geo 02 没有 --outdir，产物就落在它的运行目录
    }
  }

  ## 1) 找清单
  cand_man <- character(0)
  if (!is.null(manifest) && nzchar(manifest)) {
    cand_man <- c(cand_man, manifest)
    if (dir.exists(manifest)) {
      cand_man <- c(cand_man, list.files(manifest, pattern = paste0(GEO_MANIFEST_SUFFIX, "$"),
                                         full.names = TRUE))
    }
  }
  if (!is.null(geo_dir) && nzchar(geo_dir) && dir.exists(geo_dir)) {
    if (!is.null(geo_prefix) && nzchar(geo_prefix)) {
      cand_man <- c(cand_man, file.path(geo_dir, paste0(geo_prefix, GEO_MANIFEST_SUFFIX)))
    }
    cand_man <- c(cand_man, list.files(geo_dir, pattern = paste0(GEO_MANIFEST_SUFFIX, "$"),
                                       full.names = TRUE))
  }
  if (!is.null(geo_prefix) && nzchar(geo_prefix) && file.exists(geo_prefix)) {
    cand_man <- c(cand_man, geo_prefix)
  }
  cand_man <- unique(cand_man)
  man <- NULL; man_path <- NULL
  for (p in cand_man) {
    m <- read_kv_manifest(p)
    if (!is.null(m)) { man <- m; man_path <- p; break }
  }

  out <- list(ok = FALSE, source = NA_character_)
  if (!is.null(man)) {
    out$manifest <- man_path
    out$source   <- "manifest"
    out$dir      <- man$dir %||% dirname(man_path)
    out$prefix   <- man$prefix %||% sub(paste0(GEO_MANIFEST_SUFFIX, "$"), "", basename(man_path))
    pick <- function(k) if (!is.null(man[[k]]) && nzchar(man[[k]])) man[[k]] else NULL
    out$deg    <- pick("deg_file")
    out$all    <- pick("all_file")
    out$p_col  <- pick("p_col")
    out$logfc_col <- pick("logfc_col")
    out$change_col <- pick("change_col")
    ## geo 清单里写的是 gene_col=rownames（表示"基因名在行名里"）。
    ## 但 read_table_robust 会把空列名改成 rownames1，直接拿 "rownames" 当列名会报
    ## "--gene_col=rownames 不在表里"。所以这里统一还原成"自动识别"。
    gc <- pick("gene_col")
    if (!is.null(gc) && tolower(trimws(gc)) %in% c("rownames", "row.names", "rowname",
                                                   "row_names", "none", "")) gc <- NULL
    out$gene_col <- gc
    out$contrast <- pick("contrast")
    out$contrast_num <- pick("contrast_num")
    out$contrast_den <- pick("contrast_den")
    out$species  <- pick("species")
    out$groups   <- pick("groups")
    out$n_per_group <- pick("n_per_group")
    out$deg_p    <- pick("deg_p")
    out$deg_logfc <- pick("deg_logfc")
    out$deg_p_type <- pick("deg_p_type")
    out$n_up <- pick("n_up"); out$n_down <- pick("n_down")
  } else {
    ## 2) 没有清单：按文件名约定在目录里找
    d <- if (!is.null(geo_dir) && nzchar(geo_dir)) geo_dir else NULL
    if (is.null(d) || !dir.exists(d)) {
      stop("--geo_dir 指向的目录不存在：", ifelse(is.null(d), "(未提供)", d),
           "\n也可以直接给 --geo_manifest=<..._for_enrichment.txt>。", call. = FALSE)
    }
    out$source <- "filename"
    out$dir <- d
    pat <- if (!is.null(geo_prefix) && nzchar(geo_prefix)) {
      paste0("^", gsub("([.\\\\+*?\\[^\\]$(){}=!<>|:#-])", "\\\\\\1", geo_prefix))
    } else {
      "^(.*?)"
    }
    fdeg <- list.files(d, pattern = paste0(pat, "_DEG\\.csv$"), full.names = TRUE)
    fall <- list.files(d, pattern = paste0(pat, "_all\\.csv$"), full.names = TRUE)
    if (is.null(geo_prefix) || !nzchar(geo_prefix)) {
      if (length(fall) > 1L) {
        stop("目录里有多个 *_all.csv，无法确定用哪一个，请用 --geo_prefix= 指定：\n  ",
             paste(basename(fall), collapse = "\n  "), call. = FALSE)
      }
      ## 用 _all.csv 的前缀去配 _DEG.csv
      if (length(fall)) {
        pfx <- sub("_all\\.csv$", "", basename(fall[1]))
        fdeg <- file.path(d, paste0(pfx, "_DEG.csv"))
      }
    }
    if (!length(fall) && !length(fdeg)) {
      stop("在 ", d, " 里没找到 geo 的 *_all.csv / *_DEG.csv。\n",
           "请确认 --geo_dir 是 geo 02_deg_plots.R 的产物目录（它没有 --outdir，产物落在运行目录）。",
           call. = FALSE)
    }
    out$prefix <- sub("_all\\.csv$", "", basename(fall[1] %||% fdeg[1]))
    out$all <- if (length(fall)) fall[1] else NULL
    out$deg <- if (length(fdeg) && file.exists(fdeg[1])) fdeg[1] else NULL
    out$p_col <- "adj.P.Val"; out$logfc_col <- "logFC"; out$change_col <- "change"
    out$gene_col <- NULL; out$species <- "human"
  }

  ## 3) 补齐 + 校验
  if (is.null(out$deg) || !file.exists(out$deg)) {
    if (!is.null(out$all) && file.exists(out$all)) {
      out$deg <- out$all
      out$deg_is_all <- TRUE
    }
  }
  if (!is.null(out$all) && !file.exists(out$all)) out$all <- NULL
  if (!is.null(out$deg) && file.exists(out$deg) && is.null(out$deg_is_all)) out$deg_is_all <- FALSE
  out$ok <- !is.null(out$deg) && file.exists(out$deg)

  if (verbose) print_geo_handoff(out)
  if (!out$ok) {
    stop("geo 交接失败：没有可用的 ORA 输入文件。", call. = FALSE)
  }
  out
}

## 把 geo 交接情况打给人看（在日志/banner 之后再调用，位置更清楚）
print_geo_handoff <- function(g) {
  if (is.null(g)) return(invisible(NULL))
  note("")
  note("[geo 交接] ", if (isTRUE(g$source == "manifest"))
    paste0("读到清单 ", basename(g$manifest)) else "按文件名约定识别")
  note("  目录   : ", g$dir)
  note("  前缀   : ", g$prefix)
  if (!is.null(g$contrast)) {
    note("  对比   : ", g$contrast,
         if (!is.null(g$groups)) paste0("   （分组 ", g$groups,
                                        if (!is.null(g$n_per_group))
                                          paste0(" 各 ", g$n_per_group, " 例") else "", "）") else "")
  }
  if (!is.null(g$n_up) || !is.null(g$n_down)) {
    note("  geo 判定: 上调 ", g$n_up %||% "?", " / 下调 ", g$n_down %||% "?")
  }
  note("  ORA 输入: ", if (is.null(g$deg)) "(无)" else basename(g$deg),
       if (isTRUE(g$deg_is_all)) "   ← 注意是 _all.csv，会按阈值重新筛" else
         "   （geo 已筛过，默认不再二次筛）")
  note("  GSEA 输入: ", if (is.null(g$all)) "(无 *_all.csv)" else basename(g$all),
       "   （全部基因，按 logFC 排序）")
  if (!is.null(g$p_col)) {
    note("  列名   : gene=", g$gene_col %||% "自动识别",
         "  logFC=", g$logfc_col %||% "自动识别", "  p=", g$p_col,
         if (!is.null(g$change_col)) paste0("  change=", g$change_col) else "")
  }
  if (!is.null(g$contrast_num) || !is.null(g$contrast_den)) {
    note("  方向   : ", g$contrast_num %||% "?", " 减 ", g$contrast_den %||% "?",
         "（logFC 为正 = ", g$contrast_num %||% "?", " 更高）")
  }
  invisible(NULL)
}

## ---- 功能归类规则（把通路名归纳成功能大类）--------------------------------

default_rules_path <- function() {
  cand <- c(Sys.getenv("ENRICH_FUNC_RULES", ""),
            file.path(get_script_dir(), "..", "references", "functional-rules.tsv"),
            file.path(getwd(), "functional-rules.tsv"))
  cand <- unique(cand[nzchar(cand)])
  for (p in cand) if (file.exists(p)) return(normalizePath(p, mustWork = FALSE))
  NULL
}

load_func_rules <- function(path = NULL) {
  if (is.null(path) || !nzchar(path)) path <- default_rules_path()
  if (is.null(path) || !file.exists(path)) {
    stop("找不到功能分类规则表（functional-rules.tsv）。\n",
         "请用 --rules=<tsv> 指定，或确认技能目录下 references/functional-rules.tsv 存在。",
         call. = FALSE)
  }
  r <- utils::read.delim(path, stringsAsFactors = FALSE, quote = "", comment.char = "",
                         check.names = FALSE)
  need <- c("category", "subcategory", "pattern")
  if (!all(need %in% names(r))) {
    stop("规则表必须含 category / subcategory / pattern 三列（制表符分隔）：", path,
         call. = FALSE)
  }
  r <- r[nzchar(trimws(r$pattern)), , drop = FALSE]
  if (!nrow(r)) stop("规则表里没有有效行：", path, call. = FALSE)
  for (i in seq_len(nrow(r))) {
    ok <- tryCatch({grepl(r$pattern[i], "test", perl = TRUE); TRUE},
                   error = function(e) FALSE)
    if (!ok) {
      stop("规则表第 ", i + 1, " 行的正则表达式无效：", r$pattern[i], call. = FALSE)
    }
  }
  attr(r, "path") <- path
  r
}

## 逐行匹配（**先匹配到的优先**，所以规则表的顺序就是优先级）
classify_func <- function(x, rules, fallback = "未归类") {
  x <- as.character(x)
  x[is.na(x)] <- ""
  cat <- rep(fallback, length(x))
  sub <- rep("", length(x))
  hit <- rep("", length(x))
  done <- rep(FALSE, length(x))
  for (i in seq_len(nrow(rules))) {
    idx <- which(!done & grepl(rules$pattern[i], x, ignore.case = TRUE, perl = TRUE))
    if (length(idx)) {
      cat[idx] <- rules$category[i]
      sub[idx] <- rules$subcategory[i]
      hit[idx] <- rules$pattern[i]
      done[idx] <- TRUE
    }
    if (all(done)) break
  }
  data.frame(func_cat = cat, func_sub = sub, func_rule = hit, stringsAsFactors = FALSE)
}

## 基因重叠去冗余（贪心，按显著性从强到弱，保留"代表通路"）
## thr = Jaccard 阈值：新通路与任一已保留通路的重叠 >= thr 就视为冗余丢掉。
dedup_by_jaccard <- function(gene_str, p_adj, thr = 0.5) {
  n <- length(gene_str)
  if (n <= 1) return(rep(TRUE, n))
  gs <- strsplit(as.character(gene_str), "/", fixed = TRUE)
  gs <- lapply(gs, function(g) unique(g[nzchar(g)]))
  ord <- order(p_adj, na.last = TRUE)
  keep <- rep(FALSE, n)
  kept <- integer(0)
  for (i in ord) {
    gi <- gs[[i]]
    if (!length(gi)) { keep[i] <- TRUE; kept <- c(kept, i); next }
    redundant <- FALSE
    for (j in kept) {
      gj <- gs[[j]]
      inter <- length(intersect(gi, gj))
      if (inter == 0L) next
      if (inter / length(union(gi, gj)) >= thr) { redundant <- TRUE; break }
    }
    if (!redundant) { keep[i] <- TRUE; kept <- c(kept, i) }
  }
  keep
}

## ---- 环境自检 -------------------------------------------------------------

check_env <- function() {
  if (is.na(Sys.getenv("TMP", NA)) || !dir.exists(Sys.getenv("TMP"))) {
    note("⚠ 环境变量 TMP 未指向可用目录，R 的 tempdir() 可能出问题。")
    note("  请用 run_enrich.sh 启动脚本，或在 R 启动前把 TMP/TEMP/TMPDIR 指向一个纯 ASCII 目录")
  }
  if (!nzchar(Sys.getenv("COMSPEC", ""))) {
    note("⚠ 环境变量 COMSPEC 为空，R 的 shell() 可能报 \"'/c' not found\"。")
  }
  lc <- Sys.getenv("LC_ALL", "")
  if (nzchar(lc) && grepl("^C(\\.|$)", lc)) {
    note("⚠ 检测到 LC_ALL=", lc, "。Windows 版 R 会退化成 C locale，",
         "此时含非 ASCII 字符的路径无法寻址，请改用 run_enrich.sh 启动。")
  }
  invisible(TRUE)
}

print_banner <- function(script, subtitle) {
  note(strrep("=", 72))
  note(script, " —— ", subtitle)
  note("时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  note("R：", R.version.string, "   libPaths[1]：", .libPaths()[1])
  note(strrep("=", 72))
}

print_args <- function(o) {
  if (!length(o)) { note("（无参数）"); return(invisible(NULL)) }
  for (k in setdiff(names(o), ".positional")) {
    note(sprintf("  %-18s = %s", k, o[[k]]))
  }
  invisible(NULL)
}
