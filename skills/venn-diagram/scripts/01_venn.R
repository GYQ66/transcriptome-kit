## ===========================================================================
## 01_venn.R —— 韦恩图（2~5 个集合），美化版式与 venn.R 一致
##
## 参考代码：<本地参考脚本>/venn.R（VennDiagram::venn.diagram，scaled=FALSE，
##           圈外标签、半透明填充、2 组集合）
##
## 本脚本保留参考代码的全部版式控制量，并修掉参考代码里的几个坑：
##   - cat.col 与 fill 顺序不一致（第二次调用把 fill 换成了橙/紫，cat.col 没换，
##     于是标签颜色和圈的颜色对不上）-> 默认 cat.col = fill
##   - disable.logging=FALSE + filename=NULL 会在当前目录丢一个
##     VennDiagram.<时间戳>.log -> 这里恒 disable.logging=TRUE
##   - 参考代码写了 height/width/resolution/imagetype/margin 一堆画布参数，
##     但 filename=NULL 时 **根本不建绘图设备**，这些全是死参数 ->
##     这里自己开设备，尺寸/分辨率真正生效
##   - 参考代码连着调了两次 venn.diagram，只画了第二个 -> 只画一张
##
## 用法见 --help
## ===========================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    stop("没装 VennDiagram 包。请先 install.packages('VennDiagram')", call. = FALSE)
  }
  library(VennDiagram)
  library(grid)
})

## ---- 定位 lib -------------------------------------------------------------
.get_script_dir <- function() {
  a <- commandArgs(FALSE); i <- grep("^--file=", a)
  if (length(i)) return(dirname(normalizePath(sub("^--file=", "", a[i[1]]), mustWork = FALSE)))
  getwd()
}
.libpath <- local({
  cand <- c(file.path(Sys.getenv("VENN_SKILL_SCRIPTS"), "lib_venn_common.R"),
            file.path(.get_script_dir(), "lib_venn_common.R"),
            file.path(getwd(), "lib_venn_common.R"))
  cand <- cand[nzchar(cand)]
  hit <- cand[file.exists(cand)]
  if (!length(hit)) stop("找不到 lib_venn_common.R，试过：\n  ", paste(unique(cand), collapse = "\n  "),
                         call. = FALSE)
  normalizePath(hit[1], mustWork = FALSE)
})
source(.libpath)

USAGE <- '
01_venn.R —— 韦恩图（2~5 个集合）

== 输入（二选一，或混用）==
  --set=<标签>=<文件>       重复 2~5 次。文件可为 .txt/.csv/.tsv
                            标签 = 用户手动输入的组名，必给；缺了脚本退出码 4
                            并列出待命名的文件（不从文件名自动猜）
                            标签里写 \\n 可换行
  --table=<表> --by=<分组列> --by_levels=UP,DOWN [--gene_col=<基因列>]
                            从一张长表里按分组列拆集合（例：GEO 的 _DEG.csv
                            按 change 列拆 UP/DOWN）；--by_levels 必给（就是
                            手动指定组名和顺序），不给退出码 4
  --gene_col=<列名>         给 csv/tsv 的 --set / --table 指定基因列（默认第 1 列）

== 输出 ==
  --outdir=./venn           产物目录（默认当前目录）
  --prefix=venn             文件名前缀（默认 venn）
  --format=pdf,png          要出几种格式：pdf/png/tiff/svg
  --width=5 --height=5      画布尺寸（英寸）。2~3 组图要正圆就得 W=H
  --dpi=600                 png/tiff 分辨率
  --bg=white                white / transparent
  --export_members=FALSE    TRUE 则额外导出 <前缀>_regions.csv（各区基因清单）

== 版式（默认值 = 参考代码 venn.R）==
  --fill=#6e48fb,#ffa500    手填颜色；不给就用 --pal 调色板
  --pal=1                   内置调色板 1..4（1=参考代码紫/橙）
  --cat_col=auto            auto=跟 --fill 同色；也可手填
  --alpha=0.5               填充透明度
  --col=black               圈线颜色（注意：不是"列名"，列名是 --gene_col）
  --lwd=1.3                 线宽；--lty=1 线型
  --cex=1.5                 圈内数字字号
  --cat_cex=1.3             圈外标签字号
  --cat_pos=auto            标签角度（度，0=12点，顺时针）。
                            auto: 2 组 = -90,90（参考值）；3~5 组 = 包默认值
  --cat_dist=auto           标签离圈的距离。auto: 2~3 组 = 0.15（参考值）；
                            4~5 组 = 包默认值
  --label_count=below       标签里追加集合大小：below / inline / none
  --print_mode=raw          圈内数字：raw（个数）/ percent（%）/ raw,percent
  --sigdigs=3               百分比有效位
  --rot=0                   整体旋转角度
  --main= --sub=            标题 / 副标题
  --scaled=FALSE            按真实比例缩放圈大小（参考代码 = FALSE）
  --euler_d=FALSE           TRUE 时若一个集合是另一个的子集，会自动改画
                            "欧拉图"样式（2~3 组有效）
  --sep_dist=0.05           2 组时间隙（仅 n=2/3）
  --margin=0.12             每边留白占画布的比例（用视口实现，真的生效；
                            参考代码里的 margin=0.25 在 VennDiagram 1.8.2 已失效）
  --autofit=TRUE            TRUE：出图前先低分辨率试渲一次，若圈外标签超出画布
                            就自动加大留白（必要时再放大画布）直到不裁切

== ★ 出图前的组名确认（必须走一遍）==
  组名永远来自用户的手动输入，脚本不做自动识别（文件名、列名都不拿来当组名）：
  - --set= 缺标签、或 --table 模式没给 --by_levels 时：脚本打印"待命名的集合
    清单"并退出（退出码 4，不是报错）——把清单交给用户，逐个问"这个集合在
    图上叫什么"，拿到名字重新拼命令再跑。
  第一次跑（名字齐了）：脚本只打印"每组在图上叫什么"，并存到
            <outdir>/<前缀>_labels.txt，然后**直接退出（退出码 3，不是报错）**，
            不画图。把这些名字拿给用户最后核对一遍（重点看 \\n 折行位置）。
  --labels_confirmed=TRUE
            用户确认后用同一条命令加上这个标志再跑，才会出图。
            如果这次的名字和上次存盘的不一致（改过名），脚本仍然会停下要求重新确认。

== 其它 ==
  --log=<文件>              同时把日志写进文件
  --help                    显示本说明
'

## ---- 参数 ------------------------------------------------------------------
opt <- parse_args()
if (!is.null(opt[["help"]]) || !length(opt)) {
  cat(USAGE); quit(status = 0)
}

outdir <- gsub("\\\\", "/", oa(opt, "outdir", "."))   # 统一成 /，日志好看
prefix <- oa(opt, "prefix", "venn")
formats <- tolower(oa_vec(opt, "format", c("pdf", "png")))
formats <- match.arg(formats, c("pdf", "png", "tiff", "svg"), several.ok = TRUE)
width  <- oa_num(opt, "width", 5)
height <- oa_num(opt, "height", 5)
dpi    <- oa_num(opt, "dpi", 600)
bg     <- oa(opt, "bg", "white")
exports <- oa_lgl(opt, "export_members", FALSE)
force_unique <- oa_lgl(opt, "force_unique", TRUE)

logf <- oa(opt, "log", NULL)
if (!is.null(logf)) set_logfile(logf)

step("=== venn-diagram | R ", R.version.string, " | VennDiagram ",
     as.character(packageVersion("VennDiagram")), " ===")

## ---- 1) 收集集合 -----------------------------------------------------------
sets <- list()
src_names <- character(0)      # 与 sets 一一对应：这个名字是从哪个文件/哪一列来的

## 1a) --set=标签=文件
## 组名（标签）必须由用户手动输入，脚本不从文件名自动猜：
## 先把全部 --set 解析一遍，缺标签的一次性列出 -> 退出码 4，交给用户逐个起名。
set_args <- oa_all(opt, "set")
if (length(set_args)) {
  parsed_sets <- lapply(set_args, function(s) {
    eq <- regexpr("=", s, fixed = TRUE)
    if (eq > 0) {
      list(raw = s, nm = substr(s, 1, eq - 1), fp = substr(s, eq + 1, nchar(s)))
    } else {
      list(raw = s, nm = "", fp = s)
    }
  })
  need_names <- vapply(parsed_sets, function(p) {
    if (nzchar(trimws(p$nm)) || !nzchar(trimws(p$fp))) "" else p$fp
  }, character(1))
  need_names <- unname(need_names[nzchar(need_names)])
  if (length(need_names)) {
    step("")
    step("================ ★ 有集合还没命名（退出码 4，不是报错）================")
    step("组名不从文件名自动识别 —— 请让用户逐个【手动输入】每个集合的显示名称：")
    for (i in seq_along(need_names)) {
      step(sprintf("  第 %d 个集合 <- %s", i, need_names[i]))
    }
    step("拿到名字后拼成 --set=<名称>=<文件>（名称里写 \\n 可折行）再跑。")
    step("========================================================================")
    quit(save = "no", status = 4L)
  }
  for (p in parsed_sets) {
    nm <- gsub("\\\\n", "\n", p$nm)   # 允许用 \n 手动换行
    fp <- p$fp
    if (!nzchar(trimws(fp))) stop("--set 缺少文件路径：", p$raw, call. = FALSE)
    nm <- trimws(nm)
    if (!nzchar(nm)) stop("--set 的标签是空的：", p$raw, call. = FALSE)
    if (nm %in% names(sets)) stop("集合标签重复：", nm, call. = FALSE)
    step("--- 集合 ", length(sets) + 1, "：", gsub("\n", " / ", nm), " <- ", fp)
    sets[[nm]] <- read_gene_list(fp, oa(opt, "gene_col", NULL))
    src_names <- c(src_names, fp)
  }
}

## 1b) --table + --by
tabf <- oa(opt, "table", NULL)
if (!is.null(tabf)) {
  by <- oa(opt, "by", NULL)
  if (is.null(by)) stop("给了 --table 就必须给 --by（按哪一列分组）", call. = FALSE)
  if (!file.exists(tabf)) stop("表不存在：", tabf, call. = FALSE)
  ext <- tolower(sub("^.*\\.", "", basename(tabf)))
  sep <- if (ext == "csv") "," else if (ext %in% c("tsv", "tab")) "\t" else ","
  tab <- utils::read.csv(tabf, sep = sep, header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  if (!by %in% names(tab)) {
    stop("表里没有分组列 ", by, "；可选列：", paste(names(tab), collapse = ", "),
         call. = FALSE)
  }
  gcol <- oa(opt, "gene_col", NULL)
  if (is.null(gcol)) gcol <- names(tab)[1]
  if (!gcol %in% names(tab)) {
    stop("表里没有基因列 ", gcol, "；可选列：", paste(names(tab), collapse = ", "),
         call. = FALSE)
  }
  lv <- oa_vec(opt, "by_levels", NULL)
  grp <- as.character(tab[[by]])
  if (is.null(lv)) {
    lv_all <- sort(unique(grp[nzchar(grp)]))
    step("")
    step("================ ★ 分组还没命名（退出码 4，不是报错）================")
    step("表 ", basename(tabf), " 的分组列 [", by, "] 取值有：",
         paste(lv_all, collapse = " / "))
    step("组名不从表里自动识别 —— 请让用户手动输入：要用哪几个取值、按什么顺序")
    step("（顺序 = 图上从第 1 个圈到第 n 个圈），用 --by_levels=值1,值2 再跑。")
    step("======================================================================")
    quit(save = "no", status = 4L)
  } else {
    miss <- setdiff(lv, unique(grp))
    if (length(miss)) warn("--by_levels 里的这些取值在表中不存在：", paste(miss, collapse = " / "))
  }
  for (v in lv) {
    nm <- v
    if (nm %in% names(sets)) stop("集合标签重复：", nm, call. = FALSE)
    g <- as.character(tab[[gcol]][grp == v])
    step("--- 集合 ", length(sets) + 1, "：", nm, " <- ", basename(tabf), "[", gcol, "]，",
         sum(grp == v), " 行")
    sets[[nm]] <- g
    src_names <- c(src_names, sprintf("%s[%s] (%s=%s)", basename(tabf), gcol, by, v))
  }
}

n <- length(sets)
if (n < 2) stop("至少要 2 个集合（现在 ", n, " 个）。一个集合画饼图，不是 Venn 图。", call. = FALSE)
if (n > 5) {
  stop("VennDiagram 最多支持 5 个集合（现在 ", n,
       " 个）。6 组以上请改用 UpSet 图（UpSetR 包），Venn 图在这个数量上已经不可读。",
       call. = FALSE)
}
sets <- clean_sets(sets, force_unique = force_unique)

## ---- 2) 各区交集（先算清楚，画的和报的是同一份数）-------------------------
reg <- region_table(sets)
step("=== 各区大小（画在圈里的数就是这些）===")
for (i in seq_len(nrow(reg))) {
  step(sprintf("  %-28s %6d", gsub("\n", "/", reg$region[i]), reg$n[i]))
}
step("  合计唯一基因 ", length(unique(unlist(sets, use.names = FALSE))))
step("  集合大小：", paste(sprintf("%s=%d", gsub("\n", "/", names(sets)),
                                  vapply(sets, length, 0L)), collapse = "  "))
step("  交叉核对：各区之和 = ", sum(reg$n), "，合计唯一基因 = ",
     length(unique(unlist(sets, use.names = FALSE))),
     if (sum(reg$n) == length(unique(unlist(sets, use.names = FALSE)))) "  ✓" else "  ⚠ 对不上！")

## ---- 3) 标签 ---------------------------------------------------------------
lab_count <- oa(opt, "label_count", "below")
lab_count <- match.arg(lab_count, c("below", "inline", "none"))
labels <- names(sets)
if (lab_count != "none") {
  labels <- vapply(seq_along(sets), function(i) {
    nm <- names(sets)[i]
    k <- length(sets[[i]])
    if (grepl("\\(\\s*[0-9,]+\\s*\\)\\s*$", nm)) return(nm)  # 用户自己写了 (n)
    if (lab_count == "below") paste0(nm, "\n(", k, ")") else paste0(nm, " (", k, ")")
  }, character(1))
}
step("=== 圈外标签 ===")
for (l in labels) step("  ", gsub("\n", " / ", l))

## ---- 4) 版式参数 -----------------------------------------------------------
## 颜色
fill <- oa_vec(opt, "fill", NULL)
if (is.null(fill)) {
  pal <- oa(opt, "pal", "1")
  if (!pal %in% names(VENN_PALS)) {
    stop("没有 ", pal, " 号调色板。可用：\n  ", pal_names(), call. = FALSE)
  }
  fill <- VENN_PALS[[pal]]
  step("配色：调色板 ", pal, " = ", paste(fill, collapse = ", "))
} else {
  step("配色：--fill 手填 = ", paste(fill, collapse = ", "))
}
if (length(fill) < n) {
  stop("颜色只有 ", length(fill), " 个，需要 ", n, " 个（--fill 或 --pal 都要给够）。",
       "可选调色板：\n  ", pal_names(), call. = FALSE)
}
fill <- fill[seq_len(n)]

## 圈外标签颜色：参考代码第二次调用里 cat.col 与 fill 顺序不一致（标签颜色和圈对不上）
cat_col <- oa(opt, "cat_col", "auto")
if (identical(cat_col, "auto")) {
  cat_col <- fill
  step("标签颜色：cat.col = fill（参考代码第二次调用里这两者顺序不一致，这里已对齐）")
} else {
  cc <- oa_vec(opt, "cat_col", NULL)
  if (length(cc) != n) stop("--cat_col 要给 ", n, " 个颜色，现在 ", length(cc), " 个", call. = FALSE)
  cat_col <- cc
}

alpha <- oa_numvec(opt, "alpha", 0.5)
if (length(alpha) == 1) alpha <- rep(alpha, n)
if (length(alpha) != n) stop("--alpha 要给 1 个或 ", n, " 个值", call. = FALSE)

col <- oa_vec(opt, "col", "black"); if (length(col) == 1) col <- rep(col, n)
lwd <- oa_numvec(opt, "lwd", 1.3); if (length(lwd) == 1) lwd <- rep(lwd, n)
lty <- oa_numvec(opt, "lty", 1);   if (length(lty) == 1) lty <- rep(lty, n)
if (length(col) != n || length(lwd) != n || length(lty) != n) {
  stop("--col/--lwd/--lty 要给 1 个或 ", n, " 个值", call. = FALSE)
}

cex      <- oa_num(opt, "cex", 1.5)        # 参考代码：圈内 1.5
cat_cex  <- oa_num(opt, "cat_cex", 1.3)    # 参考代码：圈外 1.3
rot      <- oa_num(opt, "rot", 0)
sigdigs  <- oa_int(opt, "sigdigs", 3)
sep_dist <- oa_num(opt, "sep_dist", 0.05)
scaled   <- oa_lgl(opt, "scaled", FALSE)   # 参考代码：FALSE
euler_d  <- oa_lgl(opt, "euler_d", FALSE)
autofit  <- oa_lgl(opt, "autofit", TRUE)

## margin：参考代码写了 0.25，但 VennDiagram 1.8.2 里 margin 已经不是 draw.*.venn
## 的参数了（被 ... 吞掉，完全不生效）。本技能自己用视口把它实现出来，
## 含义是「每边留白占画布的比例」——这正好对应参考代码当年想干的事。
margin <- oa_num(opt, "margin", 0.12)
if (margin < 0 || margin > 0.45) stop("--margin 要落在 0~0.45（每边留白占画布的比例）", call. = FALSE)
if (margin >= 0.22) {
  warn("--margin=", margin, " 是「每边留白占画布的比例」，图形只会占画布的 ",
       round((1 - 2 * margin) * 100), "%。参考代码的 margin=0.25 在 VennDiagram 1.8.2 里",
       "已经失效（不再是 draw.*.venn 的参数），本技能用视口把它实现出来了，所以这个值",
       "现在会真的生效。只是想让圈外标签不被裁的话，交给 --autofit 自动算更省事。")
}

print_mode <- oa_vec(opt, "print_mode", "raw")
print_mode <- match.arg(print_mode, c("raw", "percent", "raw,percent", "percent,raw"),
                        several.ok = TRUE)
main <- oa(opt, "main", NULL); if (!is.null(main) && !nzchar(main)) main <- NULL
sub  <- oa(opt, "sub", NULL);  if (!is.null(sub)  && !nzchar(sub))  sub  <- NULL

## 标签角度：包默认值（探过源码）
PKG_CAT_POS <- list("2" = c(-90, 90),            # 参考代码用的就是这对（包默认是 -50,50）
                    "3" = c(-40, 40, 180),
                    "4" = c(-15, 15, 0, 0),
                    "5" = c(0, 287.5, 215, 145, 70))
cat_pos <- oa_numvec(opt, "cat_pos", NULL)
if (is.null(cat_pos)) {
  cat_pos <- PKG_CAT_POS[[as.character(n)]]
  step("标签角度：auto -> ", paste(cat_pos, collapse = ", "),
       if (n == 2) "（参考代码的 -90,90：标签在圈的正左/正右）" else "")
} else if (length(cat_pos) != n) {
  stop("--cat_pos 要给 ", n, " 个角度，现在 ", length(cat_pos), " 个", call. = FALSE)
}

cat_dist <- oa_numvec(opt, "cat_dist", NULL)
if (is.null(cat_dist)) {
  cat_dist <- if (n <= 3) rep(0.15, n) else NULL   # 参考代码 0.15；n>=4 交给包默认
  if (!is.null(cat_dist)) step("标签距离：auto -> ", paste(cat_dist, collapse = ", "), "（参考代码值）")
} else if (length(cat_dist) != n) {
  stop("--cat_dist 要给 ", n, " 个值，现在 ", length(cat_dist), " 个", call. = FALSE)
}

if (n >= 4 && euler_d) warn("euler.d 对 ", n, " 组无效（包里 draw.quad/quintuple 没这个参数），已忽略")
if (n >= 4 && scaled)  warn("scaled 对 ", n, " 组无效，已忽略")

## ---- 4b) ★★ 出图前的硬闸：逐个确认每组的显示名称 ★★ -------------------
## 用户要求：画之前必须让用户看过"每一组在图上叫什么"并确认。
## 做法：① 把清单打出来 + 存成 <前缀>_labels.txt（"上一版清单"）；
##      ② 没给 --labels_confirmed=TRUE 就不画，退出码 3；
##      ③ 给了标志，还要和"上一版清单"逐字一致才画。
## 于是闸门的真正含义是：**一旦组名变了（用户还没看过新名字），闸门自动重新关上**，
## 想画就得再走一遍"打印清单 -> 用户确认 -> 加标志"。
lab_file <- file.path(outdir, paste0(prefix, "_labels.txt"))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

esc <- function(s) gsub("\n", "\\\\n", s)          # 换行显示成字面 \n，一行一个名字
lab_names <- vapply(labels, esc, character(1), USE.NAMES = FALSE)
lab_lines <- vapply(labels, function(l) length(strsplit(l, "\n")[[1]]), integer(1),
                    USE.NAMES = FALSE)
lab_df <- data.frame(
  序号 = seq_along(sets),
  图片中显示的名称 = lab_names,
  行数 = lab_lines,
  集合大小 = vapply(sets, length, integer(1)),
  颜色 = fill[seq_along(sets)],
  名称来源 = if (length(src_names) == n) src_names else rep("", n),
  stringsAsFactors = FALSE, check.names = FALSE
)

step("")
step("==================== ★ 出图前请确认每组的显示名称 ★ ====================")
step("下面每一行 = 图上那一组圈外标签长什么样（\\n 表示换行）：")
for (i in seq_len(nrow(lab_df))) {
  step(sprintf("  %d) %-34s | %5d 个 | %s | 来自 %s",
               lab_df$序号[i], gsub("\n", "\\\\n", lab_df$`图片中显示的名称`[i]),
               lab_df$集合大小[i], lab_df$颜色[i], lab_df$`名称来源`[i]))
}
step("  （标签在图上会渲染成 ", paste(lab_lines, collapse = " / "), " 行；",
     "颜色与集合顺序一一对应）")
step("------------------------------------------------------------------------")
step("这份清单已存到：", lab_file)

## 读回上一次存盘的名称（第 2 列），用来判断"名字有没有变过"
read_prev_labels <- function(f) {
  if (!file.exists(f)) return(NULL)
  ln <- readLines(f, warn = FALSE)
  if (!length(ln)) return(NULL)
  ln[1] <- sub("^\ufeff", "", ln[1])
  ln <- ln[!grepl("^#", ln)]
  ln <- ln[nzchar(ln)]
  if (!length(ln)) return(NULL)
  vapply(ln, function(x) {
    p <- strsplit(x, "\t", fixed = TRUE)[[1]]
    if (length(p) >= 2) p[2] else NA_character_
  }, character(1), USE.NAMES = FALSE)
}
prev_lab <- read_prev_labels(lab_file)
same_lab <- !is.null(prev_lab) && length(prev_lab) == length(lab_names) &&
  all(prev_lab == lab_names)

## 存盘（制表符分隔纯文本，第 2 列是名称本身，脚本自己也要读回来比对）
lab_out <- c(paste("#", "序号", "图片中显示的名称", "行数", "集合大小", "颜色", "名称来源",
                   sep = "\t"),
             sprintf("%d\t%s\t%d\t%d\t%s\t%s",
                     lab_df$序号, lab_df$`图片中显示的名称`, lab_df$行数,
                     lab_df$集合大小, lab_df$颜色, lab_df$`名称来源`))
con <- file(lab_file, open = "w", encoding = "UTF-8")
writeLines(lab_out, con)
close(con)

conf_flag <- oa_lgl(opt, "labels_confirmed", FALSE)
if (!conf_flag || !same_lab) {
  if (!conf_flag) {
    step("★ 【停】还没画图 —— 先把上面这份组名清单拿给用户，逐个确认。")
  } else if (is.null(prev_lab)) {
    step("★ 【停】虽然带了 --labels_confirmed=TRUE，但还没有任何一版清单记录（第一次跑）",
         "——先让用户看这份清单。")
  } else {
    step("★ 【停】本次组名和上一版清单**对不上**（说明名字改过、用户还没看过这一版）",
         "——要重新确认。")
    m <- max(length(prev_lab), length(lab_names))
    for (i in seq_len(m)) {
      a <- if (i <= length(prev_lab)) prev_lab[i] else "<无>"
      b <- if (i <= length(lab_names)) lab_names[i] else "<无>"
      if (!identical(unname(a), unname(b))) {
        step("    第 ", i, " 组：上次 [", a, "] -> 本次 [", b, "]")
      }
    }
  }
  step("★ 确认无误后：同一条命令加 --labels_confirmed=TRUE 再跑一次即可出图。")
  step("★ 要改名：把对应 --set= 的标签换掉，然后再走一遍确认（脚本会重新停下来）。")
  step("★ 退出码 3 = 等待确认，不是报错。")
  step("========================================================================")
  quit(save = "no", status = 3L)
}
step("名称已确认（与 ", basename(lab_file), " 一致），继续出图。")
step("")

## ---- 5) 画 ----------------------------------------------------------------
## 传给 venn.diagram 的 ... 会一路转发给对应的 draw.*.venn，
## 各 draw 函数只认自己的参数，多传的会被 ... 吞掉（不报错但也不生效），
## 所以这里先按目标函数筛一遍，让"参数没生效"在日志里看得见。
draw_fun <- switch(as.character(n), "2" = "draw.pairwise.venn", "3" = "draw.triple.venn",
                   "4" = "draw.quad.venn", "5" = "draw.quintuple.venn")
allowed <- names(formals(getFromNamespace(draw_fun, "VennDiagram")))

extra <- list(fill = fill, alpha = alpha, col = col, lwd = lwd, lty = lty,
              cex = cex, cat.col = cat_col, cat.cex = cat_cex, cat.pos = cat_pos,
              cat.dist = cat_dist, cat.default.pos = "outer",
              rotation.degree = rot, scaled = scaled, euler.d = euler_d,
              sep.dist = sep_dist)
## 值为 NULL 的要整个丢掉：这些 draw.*.venn 里 NULL 不等于"用默认值"，
## 传进去会报 `Unexpected parameter length for "cat.dist"`。
extra <- extra[!vapply(extra, is.null, TRUE)]
dropped <- setdiff(names(extra), allowed)
if (length(dropped)) step("（", n, " 组图的 ", draw_fun, " 不接受这些参数，已丢弃：",
                          paste(dropped, collapse = ", "), "）")
extra <- extra[names(extra) %in% allowed]

step("=== 出图 ===")
step("画布 ", width, "x", height, " 英寸，", paste(formats, collapse = "/"),
     "，dpi=", dpi, "，bg=", bg, "，margin=", margin)
if (n <= 3 && abs(width - height) > 1e-6) {
  warn("2~3 组图在 VennDiagram 里是按画布比例画的：画布不是正方形时圈会被拉成椭圆",
       "（实测 5x5 时长短轴比 1.000，5x10 时会变 2 倍）。要正圆就令 --width = --height。")
}
if (n >= 4) {
  step("说明：4~5 组图的椭圆形状是 VennDiagram 的固有版式（跟画布比例无关），不是被拉伸。")
}

## venn.diagram() 里 disable.logging=TRUE 的真实行为是「日志改打到控制台」，
## 它会把整个参数列表（含全部基因名！）flog.info 出来，几千行噪声。
## 这里直接把那个 logger 的门槛抬到 ERROR，噪声就没了。
invisible(
  if (requireNamespace("futile.logger", quietly = TRUE)) {
    try(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"),
        silent = TRUE)
  }
)

grob <- do.call(venn.diagram, c(
  list(x = sets, filename = NULL, disable.logging = TRUE,
       category.names = labels, na = "remove", force.unique = force_unique,
       print.mode = print_mode, sigdigs = sigdigs, main = main, sub = sub),
  extra))

## ---- 6) 写设备 -------------------------------------------------------------
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(outdir)) stop("建不了产物目录：", outdir, call. = FALSE)

cairo <- isTRUE(capabilities("cairo"))
if (bg == "transparent" && !cairo) {
  warn("本机没有 cairo 设备，png/tiff 不支持透明背景，已改用白色")
  bg_use <- "white"
} else {
  bg_use <- bg
}

## 把 grob 画到当前设备上。margin 用「缩小视口」实现（见上文注释）：
## clip="off"，所以圈外标签画到视口之外也不会被裁，只是多占一点画布。
draw_page <- function(mrg = margin) {
  grid.newpage()
  if (mrg > 0) {
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5,
                                      width = 1 - 2 * mrg, height = 1 - 2 * mrg,
                                      clip = "off"))
    grid.draw(grob)
    grid::popViewport()
  } else {
    grid.draw(grob)
  }
}

## 内容 bbox（npc）。做法：低分辨率渲一张 png -> png::readPNG 读回来 ->
## 找非背景像素的最外圈。返回 c(x0, y0, x1, y1)，不成功返回 NULL。
## 注：R 的 dev.capture() 在 Windows 的 png 设备上不可用（实测报
## "raster capture is not available for this device"），所以走读文件这条路。
content_bbox <- function(w, h, mrg) {
  if (!requireNamespace("png", quietly = TRUE)) return(NULL)
  tf <- tempfile(fileext = ".png")
  on.exit(unlink(tf), add = TRUE)
  ok <- TRUE
  tryCatch({
    grDevices::png(tf, width = w, height = h, units = "in", res = 100,
                   bg = if (bg_use == "transparent") "transparent" else bg_use,
                   type = if (cairo) "cairo" else "windows")
    draw_page(mrg)
    grDevices::dev.off()
  }, error = function(e) {
    ok <<- FALSE
    try(grDevices::dev.off(), silent = TRUE)
  })
  if (!ok || !file.exists(tf)) return(NULL)
  img <- tryCatch(png::readPNG(tf), error = function(e) NULL)
  if (is.null(img) || length(dim(img)) < 2) return(NULL)
  ch <- dim(img)[3]
  if (!is.null(ch) && ch >= 4 && bg_use == "transparent") {
    hit <- img[, , 4] > 0.03
  } else {
    bgc <- grDevices::col2rgb(bg_use)[, 1] / 255
    if (is.null(ch)) {
      d <- abs(img - bgc[1])
    } else {
      d <- abs(img[, , 1] - bgc[1])
      if (ch >= 2) d <- d + abs(img[, , 2] - bgc[2])
      if (ch >= 3) d <- d + abs(img[, , 3] - bgc[3])
    }
    hit <- d > 0.06
  }
  rows <- which(rowSums(hit) > 0); cols <- which(colSums(hit) > 0)
  if (!length(rows) || !length(cols)) return(NULL)
  py <- nrow(hit); px <- ncol(hit)
  c(x0 = (min(cols) - 1) / px, y0 = (min(rows) - 1) / py,
    x1 = max(cols) / px,       y1 = max(rows) / py)
}

## ---- 6b) 自动适配留白：圈外标签能放进画布 -------------------------------
## 为什么加大留白有用、放大画布几乎没用：
##   在 cat.pos = ±90 这类"标签放圈正左/正右"的布局里，VennDiagram 把标签锚点
##   摆在图内 npc ≈ 0.016 / 0.984 的位置（跟 cat.dist 关系很小），而标签是
##   **以锚点为中心**画的 —— 半个标签必然甩到画布外面，这就是参考代码要设
##   margin 的原因。把留白加大 = 图形整体往内缩、锚点跟着往内走。
##   单纯放大画布没用：锚点到画布边缘的距离只占画布宽的 1.6%，涨得太慢。
##
## 做法：直接把 grob 里每个 text 元素捞出来，量它的实际尺寸（英寸），
## 再解出"让所有文字都落在画布内"所需的最小留白。一次算准，不靠试。
collect_texts <- function(g) {
  out <- list()
  walk <- function(x) {
    if (inherits(x, "text")) {
      out[[length(out) + 1L]] <<- x
    } else if (inherits(x, "gList")) {
      for (i in seq_along(x)) walk(x[[i]])
    } else if (inherits(x, "gTree") && !is.null(x$children)) {
      for (i in seq_along(x$children)) walk(x$children[[i]])
    }
  }
  walk(g)
  out
}

## 解出所需留白（每边占画布的比例）。返回 NA 表示算不出来。
## 推导：设留白 m、画布边长为 T（英寸），元素锚点在"图内"占 u（npc），
## 则它在画布上的位置 P(m) = T*u + m*T*(1-2u)（因为图形被缩到 (1-2m) 并居中）。
## 文字以锚点为中心，要求  P(m) ± size/2  落在 [pad, T-pad] 内，解出 m 的下界。
## u 接近 0.5 时（圈里的交集数字）(1-2u)≈0，这类元素跟留白无关，直接跳过。
fit_margin <- function(g, W, H, pad = 0.02) {
  ok <- TRUE
  tryCatch(grDevices::pdf(NULL), error = function(e) ok <<- FALSE)
  if (!ok) return(NA_real_)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  lower_bound <- function(u, size, total) {
    k <- total * (1 - 2 * u)
    if (!is.finite(k) || abs(k) < 1e-9) return(0)
    if (k > 0) (pad + size / 2 - total * u) / k
    else (total * (1 - u) - pad - size / 2) / k
  }
  need <- 0
  for (tx in collect_texts(g)) {
    xs <- suppressWarnings(as.numeric(tx$x)); ys <- suppressWarnings(as.numeric(tx$y))
    if (length(xs) != 1L || is.na(xs) || length(ys) != 1L || is.na(ys)) next
    w <- tryCatch(grid::convertWidth(grid::grobWidth(tx), "in", valueOnly = TRUE),
                  error = function(e) NA_real_)
    h <- tryCatch(grid::convertHeight(grid::grobHeight(tx), "in", valueOnly = TRUE),
                  error = function(e) NA_real_)
    if (is.finite(w)) need <- max(need, lower_bound(xs, w, W))
    if (is.finite(h)) need <- max(need, lower_bound(ys, h, H))
  }
  need <- max(need, 0)
  if (is.finite(need)) need else NA_real_
}

const_max_margin <- 0.36
if (autofit) {
  for (k in 1:4) {
    req <- fit_margin(grob, width, height)
    if (is.na(req)) { warn("算不出圈外标签需要的留白，跳过自动适配"); break }
    if (req <= margin + 0.002) break
    if (req <= const_max_margin) {
      step(sprintf("圈外标签需要 %.3f 的留白（当前 %.3f），自动加大", req, margin))
      margin <- req
      break
    }
    width <- width * 1.25; height <- height * 1.25
    step(sprintf("画布不够大（需要 %.3f 的留白 > 上限 %.2f），画布放大到 %.2f x %.2f 英寸重试",
                 req, const_max_margin, width, height))
  }
  req2 <- fit_margin(grob, width, height)
  if (!is.na(req2) && req2 > margin + 0.002) {
    warn(sprintf("圈外标签需要 %.3f 的留白，超过上限 %.2f —— 图形会明显偏小。",
                 req2, const_max_margin),
         "建议：缩短标签、在标签里用 \\n 折行、调小 --cat_cex，或用 --cat_pos 换标签位置。")
  }
} else {
  step("（--autofit=FALSE：不检查裁切，完全按你给的 --width/--height/--margin 出图）")
}

## 出图前最后核一遍：低分辨率渲一张，看最外圈有没有内容顶到画布边
bb <- content_bbox(width, height, margin)
if (is.null(bb)) {
  step("（量不出内容范围，跳过裁切自检）")
} else {
  ov <- max(0.004 - bb[["x0"]], bb[["x1"]] - 0.996,
            0.004 - bb[["y0"]], bb[["y1"]] - 0.996)
  if (ov > 0) {
    warn(sprintf("内容仍然顶到画布边缘（越界 %.3f）。请加大 --width/--height，或调小 --cat_cex，",
                 ov), "或在标签里用 \\n 折行。")
  } else {
    step(sprintf("裁切自检通过：画布 %.2f x %.2f 英寸，留白 %.3f，内容占画布 %.0f%% x %.0f%%",
                 width, height, margin, (bb[["x1"]] - bb[["x0"]]) * 100,
                 (bb[["y1"]] - bb[["y0"]]) * 100))
  }
}

written <- character(0)
for (fmt in formats) {
  f <- file.path(outdir, paste0(prefix, ".", fmt))
  ok <- TRUE
  tryCatch({
    if (fmt == "pdf") {
      grDevices::pdf(f, width = width, height = height, bg = bg_use)
    } else if (fmt == "png") {
      grDevices::png(f, width = width, height = height, units = "in", res = dpi,
                     bg = bg_use, type = if (cairo) "cairo" else "windows")
    } else if (fmt == "tiff") {
      grDevices::tiff(f, width = width, height = height, units = "in", res = dpi,
                      bg = bg_use, compression = "lzw",
                      type = if (cairo) "cairo" else "windows")
    } else if (fmt == "svg") {
      grDevices::svg(f, width = width, height = height, bg = bg_use)
    }
    draw_page(margin)
    grDevices::dev.off()
  }, error = function(e) {
    ok <<- FALSE
    try(grDevices::dev.off(), silent = TRUE)
    warn("写 ", fmt, " 失败：", conditionMessage(e))
  })
  if (ok) {
    sz <- file.info(f)$size
    step("  -> ", f, "  (", format(sz, big.mark = ","), " bytes)")
    written <- c(written, f)
  }
}
if (!length(written)) stop("一个格式都没写出来，图没生成。", call. = FALSE)

## ---- 7) 可选：导出各区基因清单 --------------------------------------------
if (exports) {
  rf <- file.path(outdir, paste0(prefix, "_regions.csv"))
  write_csv_bom(reg[, c("region", "n", "genes")], rf)
  step("  -> ", rf, "  (各区基因清单，分号分隔)")
}

step("=== 完成：", length(written), " 个文件 ===")
