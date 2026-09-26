# GEO 芯片分析踩坑手册

配套脚本：`scripts/01_geo_normalize.R`、`scripts/02_deg_plots.R`。
本文档记录脚本无法自动规避、需要人工判断的事项，以及常见报错的定位方式。

**本文所有结论均在 R 4.5.2 / Windows 上实机验证过**（主体 2026-09-10；
第 7 节「多数据集合并与批次矫正」为 2026-09-17 新增），不是推测。

**验证范围**（两个真实数据集，两步流程均端到端跑通）：

| GSE | 平台 | 矩阵规模 | 分组 | |logFC|>1 且 P<0.05 |
|---|---|---|---|---|
| GSE55584 | GPL96 | 13237 基因 × 16 样本 | OA 6 / RA 10 | 上调 288 / 下调 257 |
| GSE62452 | GPL6244 | 23307 基因 × 130 样本 | case 69 / control 61 | 上调 177 / 下调 115 |

- GSE55584：RA 侧 CXCL13 / CXCL9 / CCR7 / CD3D / TRBC1 / ITK / LAT 等淋巴细胞浸润
  标志上调，OA 侧 PRELP / PENK / GDF5 / BMP4 等软骨基质基因相对高表达；
  quantile 归一化后各列中位数完全一致（8.1709）。
- GSE62452：上调 top 为 LAMC2 / LAMB3 / TSPAN1 / LAMA3 / ITGB4 / SLC2A1
  （laminin-332 与 GLUT1，胰腺癌浸润前沿经典组合），下调 top 为
  SLC7A2 / AOX1 / IAPP / PDK4（正常胰腺实质与胰岛基因，被肿瘤取代）。

---

## 1. 环境准备

R 路径：随安装而定；`run_geo.sh` 自动探测，或用环境变量 `RSCRIPT` 指定。
依赖包已装齐：GEOquery 2.78.0、limma 3.66.0、Biobase 2.70.0、ggplot2 4.0.1、pheatmap 1.0.13、ggrepel 0.9.6。

换机器时一次性安装：

```r
install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "Biobase"))
install.packages(c("ggplot2", "pheatmap", "ggrepel"))
```

版本要求：R >= 4.0（`stringsAsFactors` 默认 FALSE）。GEOquery 依赖 `Biobase`、`data.table`、`R.utils`。

---

## 1.1 Windows 中文用户名：临时目录必坏（实测，最高频故障）

**现象**：`getGEO` 一进来就炸，报

```
Error: Failed to create directory (tried 5 times), most likely because of lack
of file permissions (directory 'C:/Users' exists but nothing beyond): C:/Users/i+i:
  at #21. mkdirs.default(path, mustWork = TRUE, ...)
  ...
  at #01. GEOquery::getGEO(GSE, destdir = ".", ...)
```

**根因**：用户名是中文（如 `<中文用户名>`），`TMP` / `TEMP` 环境变量带着中文路径，
R 启动时据其建立 `tempdir()`，路径字节被破坏。实测：

```
> tempdir()
"C:\Users\<中文用户名>\AppData\Local\Temp\RtmpGyEQNb"
> dir.exists(tempdir())
[1] FALSE              # 目录根本不存在
> writeLines("x", file.path(tempdir(), "probe.txt"))
Error: cannot open file 'C:\Users\i+i:\AppData\Local\Temp\RtmpIdpNaM/probe.txt': Invalid argument
```

R 的临时目录**完全不可写**，任何依赖 `tempfile()` 的操作（GEOquery 解压、部分绘图）都会挂。

**解法**：启动 R 前把临时目录指到纯 ASCII 路径。

```bash
mkdir -p D:/Rtmp
export TMP=D:/Rtmp TEMP=D:/Rtmp TMPDIR=D:/Rtmp
```

PowerShell 等价写法：

```powershell
mkdir D:\Rtmp
$env:TMP = 'D:\Rtmp'; $env:TEMP = 'D:\Rtmp'; $env:TMPDIR = 'D:\Rtmp'
```

注意 `tempdir()` 在 R 启动时就已经固定，**脚本内部改不掉**，必须在进程外设置。
两个脚本内置了自检，命中时直接打印上面这段指引。

## 1.2 `COMSPEC` 缺失：`'/c' not found`

**现象**：修好临时目录后，`getGEO` 又报

```
Error in system(cmd, intern = intern, ...) : '/c' not found
Calls: Sys.readlink2 -> sapply -> lapply -> FUN -> shell -> system
```

**根因**：`R.utils::Sys.readlink2` 走 `shell()`，而 `shell()` 的默认解释器取
`Sys.getenv("COMSPEC")`。在 Git Bash 这类环境里 `COMSPEC` 常常没被导出，
R 退回到不存在的 `command.com`。

**解法**：

```bash
export COMSPEC='C:\WINDOWS\system32\cmd.exe'
```

脚本 01 已内置兜底：检测到 `COMSPEC` 为空就自动 `Sys.setenv(COMSPEC=...)`。

## 1.3 下载超时：`downloaded length != reported length`

GEO 的 `*_family.soft.gz` 常见 5~15 MB，R 默认 `options(timeout = 60)`，实测 60 秒只下了 800 KB 就被掐断。

```r
options(timeout = 900)   # 脚本已默认 900 秒，可用 --dl_timeout 覆盖
```

脚本内置 3 次重试，并会打印实际下载到的 MB 数。网络实在不行就手动下载丢进工作目录：

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55584/soft/GSE55584_family.soft.gz
```

## 1.4 退出时段错误（可忽略）

加载过任意包的 R 进程，退出时固定 `Segmentation fault`、退出码 139：

```
$ Rscript -e 'library(ggplot2); cat("ok\n")'
ok
$ echo $?
139
```

实测确认：崩溃发生在**所有工作完成之后**，文件、PDF、PNG 均完整正确。
与 locale（`LC_ALL`）无关，与包是否有编译代码无关，换个 shell 也一样。
**判断成败请检查产物文件，不要看退出码**——这一点在写自动化脚本时尤其要注意。

### 1.4.1 脚本自身路径含中文 → R 打开脚本时乱码、偶发段错误 ★（2026-09-11 新发现）

上面 1.1 解决的是**数据/临时文件**路径带中文的问题。还有一个独立坑：
**R 脚本文件自己的路径**含中文时（本机即 `<TOOLKIT>/.../02_deg_plots.R`），
R 在打开脚本文件时会把路径字节破坏成 `C:/Users/i+ i: /...` 之类的无效路径，
表现为：**进程瞬间段错误（退出码 139）、无任何 stdout、无产物、日志为空**——
与"R 正常跑完退出时那段已知段错误"外观相同但后果更糟（产物全无）。

**关键特征**：同样的脚本，复制到一个**纯 ASCII 路径**（如 `D:/code/_rtest/s1/02_run.R`）
再跑就完全正常。所以只要脚本在中文用户名目录下，就会偶发触发。

**解法（已写入 `run_geo.sh`）**：包装器在跑 R 之前，先检测脚本路径是否含非 ASCII
字符（`printf '%s' "$1" | grep -qP '[^\x00-\x7F]'`），若含则 `cp` 到纯 ASCII 的
`$RTMP/_geo_run_tmp.R` 再运行拷贝件。工作目录不变，传给脚本的 `--expr`/`--group_file`
等参数是用户给的路径、不受影响。**用本技能的脚本一律走 `run_geo.sh`，不要直接 `Rscript`**。

## 1.5 实机验证中发现并已修复的代码缺陷

这几处都是"能跑完但不报错"或"报错信息完全指错方向"的坑，值得记下来避免重犯。

### 1.5.1 `strsplit("", sep)[1]` 返回 NA —— 会凭空多出一个 NA 基因

```r
# 错误写法
vapply(strsplit(x, "///", fixed = TRUE), `[`, character(1), 1L)
```

`strsplit("", "///", fixed = TRUE)` 返回的是 **`character(0)`** 而不是 `""`，
取 `[1]` 就得到 `NA_character_`。GPL96 有 **1058 个探针的 Gene Symbol 为空**，
于是这 1058 个全部变成 NA。

更隐蔽的是随后的过滤：

```r
ids <- ids[!ids$symbol %in% c("", "---", "NA"), ]   # NA 会被留下！
```

因为 `NA %in% c("", "NA")` 返回 **FALSE**，取反成 TRUE，NA 反而全部通过。
最终结果是一个名为 `NA` 的假基因混进矩阵，`write.csv` 写出去后再读回来就会报
`missing values in 'row.names' are not allowed`。

**正确写法**（先补分隔符保证至少有一个字段）：

```r
vapply(strsplit(paste0(x, sep), sep, fixed = TRUE), `[`, character(1), 1L)
ids <- ids[!is.na(ids$symbol) & nzchar(ids$symbol), ]
```

顺便一提：`stringr::str_split("", "///", simplify = TRUE)` 对空串返回 `""`，
所以**用户原版用 stringr 的脚本没有这个问题**，是改用 base R 时引入的回归。
脚本现已在赋值 `rownames(dat)` 后加 `anyNA()` 断言兜底。

### 1.5.2 `slot(GPL, "accession")` 在 GEOquery >= 2.70 报错

```
Error in slot(GPL@gpls[[gpl_idx]], "accession") :
  no slot of name "accession" for this object of class "GPL"
```

新版 `GPL` 类只有 `dataTable` 和 `header` 两个 slot，平台号在 `header$geo_accession`：

```r
gpl_accession <- function(g) as.character(slot(g, "header")[["geo_accession"]])
gpl <- slot(slot(g1, "dataTable"), "table")
```

很多网上流传的教程还在用 `GPL@dataTable@table`，在 2.78 上会直接崩。

### 1.5.3 `read.csv(..., row.names = 1)` 撞上名为 `NA` 的基因

基因符号字面量就是 `"NA"` 时，`na.strings` 默认值会把它转成真 NA，
`row.names = 1` 随即报 `missing values in 'row.names' are not allowed`。

规避方式：整体读入 + 单独取第一列做行名，并把 `na.strings` 设为 `""`：

```r
d  <- read.csv(path, header = TRUE, check.names = FALSE, na.strings = "")
rn <- as.character(d[[1]])
m  <- suppressWarnings(matrix(as.numeric(as.matrix(d[, -1])), nrow = nrow(d)))
rownames(m) <- rn
```

### 1.5.4 `pheatmap` 的图只能靠设备包起来

`pheatmap()` 返回的是 gtable，`ggsave()` 存它可能产出**空白文件而不报错**。
务必用设备包裹，并顺带输出 PNG：

```r
pdf("x_heatmap.pdf", width = w, height = h); pheatmap(...); dev.off()
png("x_heatmap.png", width = w*150, height = h*150, res = 150, type = "cairo"); pheatmap(...); dev.off()
```

---

### 1.5.5 `--dup=elementwise` 漏同步 `ids` → 行名长度不匹配 ★（2026-09-17 修复）

**现象**：用 `--dup=elementwise` 时，注释阶段一切正常，紧接着崩在行名赋值上，
报错完全指不出真正原因：

```
[5/6] 探针匹配与基因合并 ...
  剔除无有效基因名的探针 1058 个
  注释到 21225 个探针 / 13237 个基因
Error in dimnames(x) <- dn :
  length of 'dimnames' [1] not equal to array extent
Calls: rownames<-
```

**根因**：`--dup` 的两个分支对 `dat` / `ids` 的一致性处理不对称。

- `max` 分支：`dat` 与 `ids` 是**并行子集**的，两者长度始终相等；
- `elementwise` 分支：用 `aggregate(num, by = list(symbol = ids$symbol))` **重建**了
  `dat`（nrow = 唯一基因数），但 `ids` 仍是**探针级别**（nrow = 探针数）。
  随后那句无条件的 `rownames(dat) <- ids$symbol` 自然长度不匹配。

**解法**：在 elementwise 分支末尾把 `ids` 同步到去重后的顺序（每个基因取首个探针）：

```r
ids <- ids[match(rownames(dat), ids$symbol), , drop = FALSE]
```

**为什么默认路径没暴露**：默认是 `--dup=max`，根本不走这条分支，所以文档里的
GSE55584 / GSE62452 基线一直是好的。只有用 `--dup=elementwise` 去复现历史
`aggregate(. ~ symbol, max)` 的结果时才会撞上。

**顺带一提**：GPL96 上 `01` 得到 **13237** 个基因（与 GSE55584 基线一致）；
某些历史手工脚本在同一平台上会得到 13433 —— 差异来自注释/占位符清理口径
（`01` 额外剔除 `NA` / `null` / `NULL` 这类字面量），不是回归。

## 1.6 解析本地 GEO 文件的四个坑（GSE62452 实机发现，2026-09-10）

第一步已改为直接读用户提供的 `*_series_matrix.txt.gz` 与 `*_family.soft.gz`，
下面四个坑全部在真实数据上撞过并修复。

### 1.6.1 `gzfile(..., encoding = "latin1")` 会让 readLines 提前 EOF ★

**这是最隐蔽的一个**：给连接指定 `encoding` 参数后，`readLines` 请求 20000 行，
实际**只返回 51 行**，且不报错。结果 33364 行的系列矩阵只解析到头部，
脚本报 "系列矩阵里没读到数据表"。

```r
con <- gzfile(p, "rt", encoding = "latin1")
length(readLines(con, n = 20000L))   # 51   ← 文件实际有 33364 行
close(con)

con2 <- gzfile(p, "rt")              # 不指定 encoding
length(readLines(con2, n = 20000L))  # 20000  ✓
```

**解法**：不要传 `encoding=`，非法字节只会产生 warning，用 `warn = FALSE` 屏蔽即可。

```r
read_chunk <- function(con, n = 50000L) {
  suppressWarnings(readLines(con, n = n, warn = FALSE))
}
```

同一个 `encoding="latin1"` 在 soft 文件上却"正常"（因为 33297 行平台表恰好
在一次 50000 行的读取里读完），所以只在系列矩阵上暴露——**极易误判成文件损坏**。

### 1.6.2 系列矩阵的表头与取值都带双引号

```
!Sample_platform_id	"GPL6244"
!Sample_characteristics_ch1	"tissue: Pancreatic tumor"
ID_REF	"GSM1527105"	"GSM1527106"	...
```

因为解析时用了 `quote = ""`（避免把引号当字段包裹符），引号被当**字面量**留下：

- 平台号变成 `"GPL6244"`，与 soft 里的 `GPL6244` 匹配不上 → 报"没读到平台表"
- 列名变成 `"GSM1527105"`，第二步按样本名对分组时**一个都对不上**
- 分组候选列名变成 `"tissue:ch1` 这种畸形

**解法**：统一去引号。

```r
unq <- function(x) sub('^"(.*)"$', '\\1', x)
cl        <- unq(strsplit(body[1], "\t", fixed = TRUE)[[1]])   # 表头
rownames(m) <- make.unique(unq(as.character(d[[1]])))          # 首列
v <- unq(as.character(meta[[tag]][[k]]))                       # 所有 !Sample_ 取值
```

注意数据行的**数值没有引号**，所以 `as.numeric` 一直是对的——问题只出在
字符串字段上，因此现象是"矩阵能算、但名字全错"。

### 1.6.3 `getGEO(filename = <series_matrix>)` 慢到不可用

130 样本的数据集（14 MB 系列矩阵）：`getGEO` 跑了 **7 分钟仍未结束**，
RSS 内存涨到 855 MB，只能强杀。同样数据自写流式解析器
**端到端 63 秒**（其中平台表解析 4.5 秒）。

结论：**本地文件不要走 `getGEO()`**，自己流式读：
只扫描到 `!series_matrix_table_begin` / `!platform_table_end` 就停，
40 MB 的 soft 文件只解析需要的那一段。

```r
# 只取平台表，撞到 !platform_table_end 立即 break
i2 <- grep("^!platform_table_end", ln, useBytes = TRUE)
if (length(i2)) { body <- c(body, ln[seq_len(i2[1] - 1L)]); break }
```

`useBytes = TRUE` 可以绕开编码相关的匹配问题，配合 `^` 锚点很快。

---

## 1.7 R 的 locale 退化成 `C` → 完全无法寻址非 ASCII 路径 ★★（2026-09-17 实测）

**这是本机最隐蔽的一个坑**，而且它是 1.4.1"脚本路径含中文偶发段错误"、以及
`system2` 报 `unable to translate ... to native encoding` 的共同根因。

**现象**：路径里的中文被**悄悄改写**成乱码，所有文件操作失败，且报错信息不指向真因：

```r
> dir.exists("D:/Rtmp/_probe中文")          # 这个目录确实存在（资源管理器、其他工具都能看到）
[1] FALSE
> file("D:/Rtmp/_probe中文/ascii.txt", "r")
Error: cannot open file 'D:/Rtmp/_probed8-f/ascii.txt': Invalid argument
> dir.create("D:/Rtmp/_probeR中文")         # 静默失败，磁盘上什么都没建
[1] FALSE
```

`list.files()` 更是**静默返回空** —— 表现完全像"文件不存在"，极易误判成路径写错。
实测中 `中文` 会被改写成 `d8-f` 之类，`测试数据` 会变成别的字节序列。

**根因**：环境里带着 `LC_ALL=C.UTF-8`，而 Windows 版 R **应用不了 `C.UTF-8`**——
启动时只打印一行警告就继续：

```
During startup - Warning messages:
1: Setting LC_CTYPE=C.UTF-8 failed
```

于是 R 退化到 `C` locale。在 `C` locale 下 R 无法把字符串转成本地编码，
路径在传给文件 API 时就坏了。

**判定方法**：

```r
Sys.getlocale("LC_CTYPE")
# [1] "C"                                   <- 有问题，中文路径全部不可用
# [1] "Chinese (Simplified)_China.utf8"     <- 正常
```

**解法**：去掉 `LC_ALL` / `LANG`，让 R 用 Windows 原生 UTF-8 locale。

| `LC_CTYPE` | `dir.exists(中文路径)` | 读文件内容 | 传中文命令行给子进程 |
|---|---|---|---|
| `C`（默认，被 `LC_ALL=C.UTF-8` 带偏） | **FALSE** | 失败 | `system2` 拒绝 / `processx` 子进程乱码 |
| `Chinese (Simplified)_China.utf8` | **TRUE** | 正常 | 正常 |

- **`run_geo.sh` 已内置**：启动 R 之前 `unset LC_ALL LANG LC_CTYPE LC_COLLATE LC_MONETARY LC_TIME`
- 已跑起来的会话里补救：`Sys.setlocale("LC_CTYPE", "")`
- **派生 R 子进程时**（如 `03b_merge_from_geo.R` 调 01/03）必须先把这几个变量从环境里
  `Sys.unsetenv()` 清掉，否则子进程一启动又回到 `C`

> 经验法则：**本机所有路径尽量保持纯 ASCII**。这既是规避 1.1 的中文用户名 tempdir
> 问题，也是规避本节的问题，一举两得。

---

## 2. 数据获取阶段

### 2.0 优先用用户提供的本地文件（推荐路径）

第一步的标准输入是用户给的三个信息：**GEO 号 + `GSE*_family.soft.gz` +
`GSE*_series_matrix.txt.gz`**。两个文件都齐时脚本**完全不联网**，速度快且可复现。

```bash
bash run_geo.sh 01_geo_normalize.R \
  --gse=GSE62452 \
  --matrix="D:/path/GSE62452_series_matrix.txt.gz" \
  --soft="D:/path/GSE62452_family.soft.gz" \
  --dir=./out --prefix=GSE62452
```

注意 `--matrix` / `--soft` 的路径会**在 `setwd(--dir)` 之前转成绝对路径**，
所以写相对路径也不受工作目录切换影响。缺文件时才回退到联网下载。

GEO 网页上这两个文件的下载位置：`Series Matrix File(s)` 与
`SOFT formatted family file(s)`（后者比 supplemental 文件更全，含平台注释表）。

### 2.1 `family.soft.gz` 不在目录里
`getGEO(GSE, destdir=".")` 只下载**系列矩阵**，不下载 soft 注释文件。
脚本已内置兜底：按 `GSEnnnnnn` 推算 FTP 路径自动下载

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE55nnn/GSE55584/soft/GSE55584_family.soft.gz
```

若公司/校园网屏蔽 NCBI，手动用浏览器下载后放进工作目录即可，脚本会直接读本地文件。

### 2.2 一个 GSE 对应多个平台
`length(gset) > 1` 时脚本会打印所有平台号并选第 `--platform` 个（默认 1）。
**必须人工核对**：用 `annotation(gset[[i]])` 看平台号，并与 GEO 网页上样本所用平台比对。
soft 文件里也常含多个 GPL，脚本优先用 `annotation(eset)` 自动匹配，匹配失败才回退到 `--platform`。

### 2.3 注释列找不到 —— 有些平台根本没有 `Gene Symbol` 列 ★

不同平台列名差异大：`Gene Symbol` / `Symbol` / `GENE_SYMBOL` / `gene_symbol`。
脚本会依次尝试常见名，再用 `grep("symbol", ignore.case=TRUE)` 兜底。

**但 Affymetrix 转录本芯片（ST 系列）连 symbol 列都没有。** GPL6244
（HuGene 1.0 ST）的 12 个列是：

```
ID  GB_LIST  SPOT_ID  seqname  RANGE_GB  RANGE_STRAND  RANGE_START
RANGE_STOP  total_probes  gene_assignment  mrna_assignment  category
```

基因名藏在 `gene_assignment` 里，格式为：

```
NM_001004195 // OR4F4 // olfactory receptor ... // 15q26.3 // 26682
ENST00000328113 // OR4G2P // ... // --- // ---
---
```

即 `登录号 // 基因名 // 描述 // 位置 // EntrezID`，多个赋值用 ` /// ` 分隔，
无基因时整格是 `---`。

**提取方法**（先按 `///` 取第一条，再按 `//` 取第 2 段）：

```r
v1 <- vapply(strsplit(v, "///", fixed = TRUE),
             function(x) if (length(x)) x[1] else "", character(1))
sym <- vapply(strsplit(v1, "//", fixed = TRUE),
              function(x) if (length(x) > 1) trimws(x[2]) else "", character(1))
```

GPL6244 实测：33297 行中 **25293 行有基因名**（其余 8004 行的
`gene_assignment` 是 `---`，多为非编码/未注释转录本）。脚本已自动兜底，
命中时会打印"平台无 Gene Symbol 列，改用 'gene_assignment' 字段提取基因名"。

其他情况：`getGEO(..., AnnotGPL=TRUE)` 得到的表可能用 `Gene ID` 存 Entrez ID，
此时需要 `--annot_col="Gene ID"` 再用 `org.Hs.eg.db` 转符号，本脚本不覆盖这条路；
更严谨的做法是装对应的注释包（如 `hugene10sttranscriptcluster.db`）用 `mapIds()` 映射。

---

## 3. 表达矩阵处理阶段

### 3.1 log2 判断不是万能的
脚本用分位数启发式规则（99% 分位数 > 100、极差 > 50、或 0<Q1<1 且 1<Q3<2）。
**单通道芯片通常可靠，双通道（如 Agilent two-color）会误判**。
跑完务必用输出的 `_boxplot.pdf` 目测：如果中位数落在 0~15 且分布对称，一般已 log；
落在几百到上万，说明没 log。判断错就用 `--force_log=TRUE/FALSE` 强制。

`log2` 前会把 <=0 的值置为 NA，脚本会打印置 NA 的个数；若这个数很大，说明数据可能已做背景校正/已 log，别硬转。

### 3.2 重复基因合并方式
- `--dup=max`（默认）：先按每个探针的行均值降序排，再每个 symbol 取第一个。
  结果等价于"保留整体表达量最高的探针"，复杂度 O(n log n)，两万基因秒级完成。
- `--dup=elementwise`：逐元素取最大值，严格等价于原版 `aggregate(. ~ symbol, max)`。
  用 `aggregate` + 匿名函数时大矩阵会明显变慢，仅在需要严格复现旧结果时用。

### 3.3 原版 `aggregate(. ~ symbol, dat, max)` 的两个隐患
1. `probe_id` 是字符型，也被 `max` 作用，得到字典序最大的探针号（无害，因为随后被丢弃，但纯属浪费）。
2. 默认 `na.action = na.omit`：任一样本含 NA 就整行丢弃，会**静默丢基因**。
   脚本改用 `by=` 形式 + `na.rm=TRUE`，只在全 NA 时才留 NA。

### 3.4 归一化方法怎么选
| 方法 | 适用 |
|---|---|
| `quantile`（默认） | 单通道芯片（Affymetrix、Illumina）首选 |
| `median` | 只需校正中位数偏移，最保守 |
| `cyclicloess` | 样本数少（< 10）且需要精细校正，慢 |
| `scale` | 仅缩放，不校正分布形状 |
| `none` | RNA-seq 的 FPKM/TPM 后续流程，或已归一化的数据 |

`normalizeBetweenArrays` 默认按**列**（样本）归一化，与 limma 的假设一致，不要转置传参。

---

## 4. 差异分析阶段

### 4.1 分组信息是整条流程最容易出错的地方 ★

**这是整条流程唯一必须由人拍板的环节。** 第一步跑完会生成
`<prefix>_grouping_candidates.txt`，中性列出所有"取值 2~10 个、每组至少 2 个样本"
的表型列，附每个取值的样本数（**不包含任何** `--group_from` / `--contrast` 建议）：

```
[1] tissue:ch1   (唯一值 2 个)
      Pancreatic tumor                           69
      adjacent pancreatic non-tumor              61
```

**重要（用户硬性规则）**：第一步只负责输出临床信息表 + 候选列报告（取值与样本数），
**绝不能替用户决定分组、绝不自动写死 `--group_from` / `--contrast`**。把报告交给用户，
**停下来明确询问**"用哪一列、取哪些值、谁减谁"，等用户给出明确分组指令后才跑第二步。

四种分组方式：

| 方式 | 参数 | 适用 |
|---|---|---|
| 内联指定 | `--group_spec='case=GSM1,GSM2;control=GSM3'` | 交互场景首选，无需额外文件 |
| 分组文件 | `--group_file=group.csv`（sample,group 两列） | 用户已有分组表；只列部分样本也行，未列出的会被剔除 |
| 临床表提取 | `--clinical=... --group_from=tissue:ch1` | 从 clinical 表某列取 |
| 直接给向量 | `--group_values=case,case,control` | 样本少且顺序明确 |

分组多于 2 组时用 `--group_keep=A,B` 挑出要对比的两组。

`--group_from` 常配合 `--group_regex='^(Tumor|Normal).*'`（必须**有且仅有一个**捕获组）。

脚本会自动剔除分组为空的样本，打印分组概览表（`table(group)`）与每个分组的样本示例，
并把完整映射写入 `<prefix>_sample_group_map.csv`。
**先核对这张映射表：样本数对不上、分组反了，就别往下看结果。**
`--contrast` 的方向决定了 logFC 的正负含义（`A-B` 表示 A 相对 B）。

### 4.2 `makeContrasts` 的非法变量名
分组名含空格、短横、数字开头时 `makeContrasts` 会报
`Error: unexpected symbol` 或 `object 'xxx' not found`。
脚本已用 `make.names()` 清洗并把映射关系打印出来，此时 `--contrast` 要用**清洗后的名字**
（如 `Tumor.tissue-Normal.tissue`）。脚本会校验对比式中的名字是否存在于分组并报错。

### 4.3 `P.Value` vs `adj.P.Val`
原版脚本用 `P.Value < 0.05`。严格做法是用 FDR：`--p_type=adj.P.Val --pval=0.05`。
发文章推荐后者，探索性分析可用前者。阈值别硬套 `|logFC| > 1`，看数据分布再定。

### 4.4 `topTable` 排序
`sort.by` 默认 `"B"`（log-odds），不是 p 值。想让热图/火山图 top 基因按显著性排就用 `--sort_by=P`。
`--sort_by=none` 保持原始顺序。

---

## 5. 绘图阶段

### 5.1 原版热图的三个问题（脚本已修）
1. `ann_colors` 定义了却**从未传进 `pheatmap()`**，分组配色不生效 → 改为 `annotation_colors = list(group = ...)`。
2. `gaps_col = 8` 是硬编码 → 改为按分组游程 `cumsum(rle(group)$lengths)` 自动算分界。
3. `ggsave("heatmap.pdf", a, ...)` 对 `pheatmap` 返回的 gtable 在部分 ggplot2 版本下不产出内容
   → 改用 `pdf(); pheatmap(); dev.off()`。

### 5.2 `pheatmap(scale="row")` 报 `row scaling: sd = 0`
某基因在所有样本中表达完全相同（或全 NA）时触发。脚本已预先剔除零方差行并打印剔除数量。
若剔除数量异常多，回头检查是不是过滤/归一化把数据搞坏了。

### 5.3 样本数过多时热图看不清
列数 > 50 建议 `show_colnames = FALSE`（脚本默认如此）并加大 `--hm_width`。
`cluster_cols = FALSE` 是硬编码的——因为原始样本按分组排序后才有 `gaps_col` 的意义；
若想看聚类结果，需要手动改脚本并去掉 `gaps_col`。

### 5.4 配色约定
上调 **红 `#C31E1F`**、下调 **蓝 `#1F6FC3`**、不显著灰 `#898989`。
这是生信领域通用约定，也符合中文语境（红=升 / 蓝绿=降）。改配色时注意 `scale_colour_manual`
必须用**命名向量**绑定因子水平，否则换数据后颜色会串位。

### 5.5 中文乱码
Windows 下 ggplot2 默认字体族不含中文。图中如需中文标题，
改用 `showtext` 包，或导出 PDF 后在 Illustrator 里改字；
`pdf()` 设备默认不嵌入中文字体，直接写中文会变成空白或方框。

> 补充（2026-09-17）：实测 `ragg::agg_png` / `cairo_pdf` 缺中文字形时**不是画方框，
> 而是挑到别的字形、渲染成看起来像乱码的 ASCII**，比方框更难察觉。
> **本技能的约定是绘图标签一律只用 ASCII**。详见 7.8。

### 5.6 默认 `pdf()` 设备不认 Arial → `invalid font type` ★（2026-09-11 新发现）
`--fig_style=nature` 用 Arial 系列字体。直接 `ggsave(".pdf", ...)` 走默认 `pdf()`
设备时，其 PostScript 字体库**没有 Arial 条目**，会刷屏
`font family 'Arial' not found in PostScript font database`，最终
`Error in grid.Call.graphics(C_text, ...): invalid font type` 并**中断脚本**
（产物只出到一半）。
- **修复**：PDF 一律用 `grDevices::cairo_pdf`（cairo 能解析系统字体），
  PNG/TIFF 用 `ragg::agg_png` / `ragg::agg_tiff`。脚本已内置。
- SVG：R ≥ 4.5 的 `svg()` 需要 `svglite`，缺失时报
  `The package "svglite" is required to save as SVG.`；脚本会自动跳过，
  **矢量图可用 cairo 出的 PDF 顶替**。

---

## 6. 输出文件清单

**01 步**（在 `--dir` 目录下，`<prefix>` 默认等于 GSE 号）

| 文件 | 内容 |
|---|---|
| `<prefix>.csv` / `<prefix>.txt` | 标准化后表达矩阵（基因 × 样本），02 步的输入 |
| `clinical_<prefix>.csv` | 样本临床/表型信息，用于提取分组 |
| **`<prefix>_grouping_candidates.txt`** | **分组候选列报告，第一步跑完要交给用户看** |
| `<prefix>_group_template.csv` | 分组模板（sample,group），用户填好可直接喂给 02 步 |
| `<prefix>_boxplot.pdf` | 归一化前后箱线图，**用于判断 log2 是否正确** |
| `<prefix>_normalized.RData` | `dat`、`pd` 两个对象 |

**02 步**（`<prefix>` 默认 `DEG`）

| 文件 | 内容 |
|---|---|
| `<prefix>_sample_group_map.csv` | **样本→分组映射，先核对这张表再信结果** |
| `<prefix>_all.csv` | 全部基因的 logFC / P.Value / adj.P.Val / change |
| `<prefix>_DEG.csv` | 仅显著差异基因 |
| `<prefix>_volcano.pdf/.png` | 火山图（`--fig_style=nature` 时另加 `.tiff`，及有 svglite 时的 `.svg`） |
| `<prefix>_heatmap.pdf/.png` | top N 差异基因热图（PNG 可用 `--heatmap_png=FALSE` 关闭；nature 同上另加 `.tiff`/`.svg`） |
| `<prefix>.RData` | `expr`、`group`、`nrDEG`、`DEG` |

**03 步**（`<prefix>` 默认 `merged`）

| 文件 | 内容 |
|---|---|
| **`<prefix>_merged.csv` / `.txt`** | **矫正后合并矩阵，02 步的输入** |
| `<prefix>_merged_raw.csv` / `.txt` | cbind 合并、未矫正（QC 对照） |
| `<prefix>_merged_combat.csv`、`_merged_quantile.csv` | 仅 `--method=both` 时输出 |
| **`<prefix>_batch_map.csv`** | **样本→数据集/批次/分组 映射，先核对这张表** |
| `<prefix>_group.csv` | `sample,group` 两列，给了 `--group_file` 时生成，可直接喂给 02 步 |
| **`<prefix>_clinical.csv`** | **合并临床表**（需 `--clinical_files`；03b 自动填 `per_gse/clinical_<GSE>.csv`）。行名=样本名，列为 `dataset`/`batch`/`group` + 各数据集临床列（跨数据集的列带 `<数据集>::` 前缀）；02 步 `--clinical=` 可直接读 |
| **`<prefix>_grouping_candidates.txt`** | **向用户索要分组的主材料**：逐数据集候选列与取值计数、每个取值 `-> 分组名` 的映射审计、跨数据集同名列矩阵 |
| `<prefix>_overlap_report.txt` | 各数据集基因数、逐级交集、大小写敏感性、批次效应量化、批次统计 |
| `<prefix>_batch_stats.csv` | 每个批次矫正前后的均值/中位数/标准差 |
| `<prefix>_pca_before` / `_after`（+ `_after_group`） | 批次效应 PCA |
| `<prefix>_boxplot_before` / `_after` | 样本分布箱线图（每个样本一个箱体，按批次着色） |
| `<prefix>_boxplot_batch_before` / `_after` | 每个批次一个箱体（批内各样本中位数）——样本数大时仍清晰 |
| `<prefix>_density_before` / `_after` | 样本密度曲线 |
| `<prefix>_cor_heatmap.pdf/.png` | 仅 `--cor_heatmap=TRUE` 时输出 |
| `<prefix>.RData` | `merged`、`merged_raw`、`batch`、`group`、PCA 对象 |

`--fig_style=nature` 时图中另加 `.svg`/`.tiff`。

---

## 7. 多数据集合并与批次矫正阶段（03 步）

配套脚本：`scripts/03_merge_batches.R`。设计说明见 `references/merge-batch-correction.md`。

### 7.1 `prcomp` 会求全部奇异向量 → PCA 卡住数分钟 ★（2026-09-17 实测）

**现象**：03 步跑到 `[7/8] 生成 QC 图表` 后长时间无输出、产物目录为空，
R 进程 CPU 只占 ~40%（看着像死锁，其实在算）。

**根因**：`prcomp(t(m), center = TRUE)` 对 `n(样本) × p(基因)` 的矩阵会求出
**全部 `min(n,p)` 个奇异向量**，LAPACK `dgesdd` 的代价约 `(n+p)·min(n,p)²`。
15913 基因 × 967 样本 时约 `3×10¹¹` flops，属分钟级以上，而脚本前后各要跑一次。
实测：在本数据上 QC 步骤停留数分钟仍未产出全部图表（进程 CPU 仅 ~40%，**不是死锁**）。

**解法**：用 `RSpectra::svds` 只求前几个主成分，代价 `O(n·p·k)`。

| 检验（500 样本 × 3000 基因） | prcomp | svds(k=10) |
|---|---|---|
| 用时 | 6.14 s | **0.31 s（20×）** |
| PC1 方差占比 | 79.6112% | 79.6112% |
| PC2 方差占比 | 3.4421% | 3.4421% |
| 前 10 个 sdev 最大相对差 | — | **3.1e-15** |

全尺寸 967 × 15913：`svds(k=10)` **2.78 s**，PC1 79.659% / PC2 3.508%。
脚本已默认走这条路（缺 `RSpectra` 时退回 `prcomp` 并打印提示）。

#### 7.1.1 `RSpectra::svds` 不保留 `dimnames` ★

换成 `svds` 后立刻踩到：**`svds()$u` 和 `$v` 的 `rownames`/`colnames` 全是 `NULL`**
（实测），而 `prcomp()$x` 是保留 rownames 的。若下游按行名对齐，例如

```r
data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], meta[rownames(pca$x), , drop = FALSE])
```

`meta[NULL, ]` 会返回 **0 行**，报
`arguments imply differing number of rows: 967, 0`。

**解法**（两条都做，互为兜底）：

```r
rownames(u) <- rownames(Xc)                      # 手工补回行名
colnames(u) <- paste0("PC", seq_len(ncol(u)))
# 且下游改成按行位置对齐 + 显式断言，不依赖 dimnames
```

`RSpectra::svds` 的起始向量还是随机的，要固定 `set.seed()` 才能保证结果可复现。

### 7.2 `ComBat` 的 `mod` 必须用带截距的 `model.matrix(~ group)` ★

`sva::ComBat` 内部是：

```r
design <- cbind(batchmod, mod)
check  <- apply(design, 2, function(x) all(x == 1))   # 全 1 列
design <- as.matrix(design[, !check])                 # 丢掉全 1 列
if (qr(design)$rank < ncol(design)) stop("...confounded...")
```

**全 1 的截距列会被自动丢掉**，所以必须传带截距的 `~ group`：截距去掉后正好剩
1 列虚拟变量。若写成 `~ 0 + group` 会得到 2 列虚拟变量，两者之和等于全 1 向量
（= batch 列之和）→ 共线 → 秩亏 → **误报**
`The covariates are confounded!`，而实际并不混杂。

### 7.3 批次与分组混杂

某些批次内只有一个分组水平（例如每个 GSE 只含一种组织）时，"批次效应"与
"生物学差异"不可分。`ComBat` 会抛 `The covariate is confounded with batch!`；
脚本在调用前就主动检测并打出"批次 × 分组"交叉表，给出三条出路
（`--combat_mod=none` / 只合并分组构成相近的数据集 / `--method=quantile`）。

### 7.4 `ComBat` 不接受 NA

```r
if (any(is.na(dat))) stop("Data contains NA. To proceed, please remove or impute missing values.")
```

`--match=union` 时很容易产生 NA。脚本默认 `--na=error`，会打印缺失最多的基因名
并提示 `--na=rowmean` / `--na=zero`，而不是让用户撞 sva 那句无信息量的报错。

### 7.5 读大矩阵：`read.table` 慢一个数量级

7.R 那类写法用 `read.table`。实测同样三个文件：

| 文件 | `read.table` | `fread` |
|---|---|---|
| 19937×553（79 MB） | 26.9 s | **2.4 s** |
| 20824×307（58 MB） | 19.8 s | **1.7 s** |
| 30905×107（28 MB） | 10.5 s | **0.8 s** |

脚本优先用 `data.table::fread`，缺包时退回 `read.table`。
注意 `fread` 的 `na.strings` **不要**包含 `"NA"`，否则名为 `NA` 的基因会被转成真 NA
（同 1.5.3）。

### 7.6 `<prefix>_batch_map.csv` 不能直接当 02 步的 `--group_file`

它的列是 `sample, dataset, batch[, group]`，第 2 列是 `dataset` 不是分组。
02 步的 `--group_file` 取第 2 列当分组，会静默把 `dataset` 当成分组。
脚本另导了 `<prefix>_group.csv`（`sample, group` 两列）供 02 步直接使用。

### 7.7 并发跑多个大 R 任务会因换页显著变慢

15913×967 的合并矩阵，raw / combat / quantile 三份同时驻留时单进程约
**1.6 GB**。本机 16 GB 内存下同时跑两个这类任务会把物理内存吃光，
PCA 一步会从十几秒退化到数分钟（实测踩过，症状与 7.1 混淆）。
**一次只跑一个。**

### 7.8 图内文字必须一律用 ASCII ★（2026-09-17 实测）

**现象**：03 步写中文图标题时，图上出现的是**乱码 ASCII**而不是方框——
标题 `PCA - 批次矫正前` 渲染成 `PCA - f 9f,lg □+f-#e □□`，
图例 `批次` 渲染成 `f 9f,!`。

**根因**：同 5.5（Windows 下 ggplot2 默认字体族不含中文字形），但表现更隐蔽：
`ragg::agg_png` 与 `cairo_pdf` 在缺字形时不是画空框，而是**挑到别的字形**，
看起来像乱码文本，很容易被当成"标题写错了"而不是"字体问题"。

**结论（本技能的约定）**：**绘图函数里一律只用 ASCII 标签**。
- 02 步一直是这么做的：`labs(x = "log2 fold change", y = paste0("-log10 ", PCOL))`，
  图例标题直接 `element_blank()`
- 03 步的 PCA / 箱线图 / 密度图同理，标签用 `"Batch"` / `"Group"` /
  `"Expression (log2 scale)"` / `"Density"` 等英文
- **控制台 message 里可以照常用中文**（实测 `message()` 输出中文正常），
  只有**画到图上的文字**受限

要中文图，得先装 `showtext` 并显式注册中文字体，或导出后在 Illustrator 里改字。

### 7.9 合并矩阵的分组不能想当然 ★（2026-09-17 实测）

合并矩阵的样本来自多个数据集，**分组是跨数据集对齐出来的，比单数据集更容易错**：

1. **列名与取值词表都不同** —— 同一个概念在两份数据里可能连列名都不一样：
   GSE55235 用 `source_name_ch1`（值 `synovial tissue from osteoarthritic joint`），
   GSE55457 用 `clinical status:ch1`（值 `osteoarthritis`）。不存在"一个列名
   走遍所有数据集"这种捷径，必须逐数据集给 `--group_cols` / `--group_maps`。
   **所以必须把 N 份临床信息汇总后交给用户，让他来定。** 03 步的
   `<prefix>_clinical.csv`（合并临床表）与 `<prefix>_grouping_candidates.txt`
   （候选列 + 每个取值 `-> 分组名` 的映射审计）就是干这个的。
2. **别用数据集/批次当分组** —— "A 数据集全是 case、B 全是 control"会让
   `mod=~group` 与批次列共线，ComBat 直接拒绝（7.3）。候选报告里每个取值后面的
   `-> 分组名` 能一眼看出某个分组是否只落在单个数据集里。
3. **`--group_na=drop` 之后，临床表必须跟着裁** —— 否则会给用户看到已经不在矩阵
   里的样本。03 步的合并临床表是以**最终矩阵的样本**为准建的；候选报告里那些被剔
   掉的取值标成 `(未进入矩阵, 已剔除)`、部分被剔的标 `[保留 3/10]`。
4. **候选列判定要用"完整临床表"，不能用幸存样本** —— 某个取值整类被剔除后，这列
   在幸存样本里只剩 1 个取值，会被当成"单一取值"跳过——而这正是用户最需要看到的
   情况（"我把 RA 全剔了，是不是这列就没法用了"）。
5. **`--clinical_files` 的文件名是 `clinical_<GSE>.csv`（前缀在前）** ——
   不是 `<GSE>_clinical.csv`。拼错时 03b 会用 `ifelse(file.exists(..), .., "")`
   静默转成空串，症状是"报告说没有任何数据集找到临床表"，而 `per_gse/` 里明明有。
6. **换分组不必重算 ComBat** —— `03b --stage=group` 只重建分组表，样本范围取自
   合并临床表。注意它**不改变**合并矩阵里 ComBat 用的 `mod`；要让 ComBat 也保护
   新分组，得带 `--group_cols/--group_maps` 重跑 `--stage=merge`。
7. **02 步吃合并临床表时，跨数据集的列名带 `<数据集>::` 前缀** ——
   `--group_from=source_name_ch1` 在合并表里可能不存在（会被加前缀），
   报错信息会列出全部可选列；也可以先用 `group` 这一列（03 已把统一分组写进去）。

### 7.10 图内标签列宽用 `sprintf` 时注意 CJK 与长取值

报告表格用 `sprintf("%-44s")` 对齐时有两个坑：`nchar()` 按**字符数**算宽度，
中日韩字符在等宽终端里占 2 列，中英混排的表格看起来会错位；GEO 的
`source_name_ch1` 取值常超过 44 字符（如
`synovial tissue from rheumatoid arthritis joint` 是 45），超出后计数会被挤到右边。
03 步的报告因此：表格首列名用 ASCII（`column`），取值列宽按本列最长值自适应
（上限 52），其余说明性文字里出现 CJK 不影响排版。

### 7.11 用户自己给的分组表：格式与实测行为（2026-09-17 实测）

向用户索要分组时，**别逼他按你的列名规则回答**——他手上往往已经有一张表
（样本清单、表型汇总、文章补充材料）。只要两列：**第一列样本名、第二列分组**，
表头可有可无，逗号/制表符都行，**不必覆盖全部样本**。
`<prefix>_group_template.csv` 就是给他填的（全部样本已列好）。

实测（14 样本的合并矩阵，只给 8 个样本的分组）：

| 情形 | 行为 |
|---|---|
| Excel 另存的 **UTF-8 BOM** 表 | **可用**。02 的 `--group_file` 按**位置**取列（第 1 列样本名、第 2 列分组），不对列名做匹配，BOM 只影响列名文本 |
| 表头叫什么 | 无所谓，同上（按位置取） |
| 只覆盖部分样本 | 其余样本被**自动剔除**并打印清单：`! 剔除未分组样本 6 个: GSM202, ...`——不是错误 |
| **两列写反** | `Error: 分组文件第一列与表达矩阵样本名无一匹配。` 必须第一列样本名、第二列分组 |
| **中文组名**（UTF-8） | 能正常读入；但 `--contrast` 必须用中文组名，否则 `对比式 'case-control' 中的 'case,control' 不在分组中。可用分组: 病例, 对照` |
| 中文组名但文件是 **GBK**（Excel 的"CSV (逗号分隔)"） | 组名会乱码 → `--contrast` 对不上。**让他另存为"CSV UTF-8"**，或组名一律用 ASCII |
| 组名带空格/短横 | 会被 `make.names` 清洗（`rheumatoid arthritis` → `rheumatoid.arthritis`），`--contrast` 用清洗后的形式 |

**合并阶段（03/03b）用同一张表时行为不同**：默认 `--group_na=error`，样本没覆盖到会
**直接报错**并提示"用 `--group_na=drop` 剔除这些样本，或补齐分组文件"。原因很简单：
合并阶段的 `group` 还要喂给 `ComBat` 的 `mod`，静默丢样本不合适。
所以**用户给的表只覆盖部分样本时，最省事的用法是直接给 02 步**（差异分析），
不要在合并阶段用；要用在合并阶段就补全或显式加 `--group_na=drop`。

拿到表之后先看两处再信结果：02 步打印的剔除清单 + `<prefix>_sample_group_map.csv`
（样本→分组完整映射）。

---

## 8. 快速排错表

| 报错 | 原因 / 处理 |
|---|---|
| `Failed to create directory ... C:/Users/xxx:` | 中文用户名导致 tempdir 不可写，见 1.1 |
| `'/c' not found` | `COMSPEC` 未导出，见 1.2 |
| `downloaded length != reported length` / `Timeout of 60 seconds` | 下载超时，加 `--dl_timeout`，见 1.3 |
| 系列矩阵"没读到数据表"（文件明明完整） | `gzfile(encoding=)` 导致 readLines 提前 EOF，见 1.6.1 |
| 平台号显示成 `"GPL6244"`（带引号）/ 列名带引号 | 系列矩阵取值带引号未剥离，见 1.6.2 |
| `分组文件第一列与表达矩阵列名无一匹配` | 多半是列名带引号，见 1.6.2 |
| `getGEO` 解析本地矩阵卡住不返回 | 别用 getGEO 读本地文件，见 1.6.3 |
| `平台表里找不到基因名列` | Affy ST 类平台无 symbol 列，走 gene_assignment 兜底，见 2.3 |
| `no slot of name "accession" for ... class "GPL"` | GEOquery >= 2.70 结构变化，见 1.5.2 |
| `invalid font type` / `font family 'Arial' not found` | 默认 `pdf()` 不认 Arial，改用 `cairo_pdf`，见 5.6 |
| `The package "svglite" is required to save as SVG` | 缺 svglite，SVG 被跳过，矢量图用 PDF 顶替，见 5.6 |
| `missing values in 'row.names' are not allowed` | 有名为 `NA` 的基因，见 1.5.1 / 1.5.3 |
| 矩阵里出现名为 `NA` 的基因 | 同上，`strsplit` 空串取值陷阱 |
| `找不到对象 'exprs'` | 未加载 Biobase，`library(GEOquery)` 没跑起来 |
| `cannot open file 'GSE..._family.soft.gz'` | 见 2.1，手动下载或检查网络 |
| 注释后基因为 0 行 | 平台选错，加 `--gpl=GPLxxx` 或改 `--platform` |
| `contrasts can be applied only to factors with 2 or more levels` | 分组只有 1 类，见 4.1 |
| `object 'XXX' not found`（makeContrasts 处） | 分组名未 `make.names` 清洗，见 4.2 |
| `row scaling: sd = 0` | 见 5.2 |
| 火山图一片灰 | 阈值太严或数据没归一化，先看 `_boxplot.pdf` |
| 输出 PDF 是空白 | 见 1.5.4，必须用设备包裹而非 `ggsave` |
| 03 步卡在 `[7/8] 生成 QC 图表`、产物为空、CPU 只占 ~40% | `prcomp` 求全部奇异向量太慢，见 7.1；装 `RSpectra` |
| `The covariate is confounded with batch!`（03 步） | 批次与分组混杂，见 7.3 |
| `The covariates are confounded!`（分组明明没混杂） | `mod` 写成了 `~0+group`，见 7.2 |
| `Data contains NA. To proceed, please remove or impute...` | ComBat 不接受 NA，加 `--na=rowmean`，见 7.4 |
| `reference level ref.batch is not one of the levels` | `--combat_ref` 必须是批次标签之一，见 `merge-batch-correction.md` |
| `Note: one batch has only one sample, setting mean.only=TRUE` | 正常提示，非错误 |
| `共同基因为 0 个`（03 步） | 基因命名不一致，试 `--gene_case=upper` 或 `--match=union` |
| 03 步批次效应矫正后 PCA 没改善 | 看 `_overlap_report.txt` 的 `ComBat 保护生物变异` 是否为"是" |
| 03 步 `--method=quantile` 后批次效应量化几乎没降 | 正常现象：分位数归一化去不掉**基因特异**的平台效应（实测 99.3%→98.6%），改用 `--method=combat`，见 `merge-batch-correction.md` 5.1 |
| 03 步 `--group_file` 匹配率低 | 样本名重复被加后缀（脚本会警告），或命名不一致 |
| 把 `<prefix>_batch_map.csv` 当 02 步的 `--group_file` | 第 2 列是 `dataset`，会静默错分组；用 `<prefix>_group.csv`，见 7.6 |
| 合并后没有 `<prefix>_clinical.csv` / `_grouping_candidates.txt` | 03 没收到 `--clinical_files`；03b 应自动填 `per_gse/clinical_<GSE>.csv`（**前缀在前**），见 7.9 |
| 03b 提示"没有任何数据集找到 clinical_<GSE>.csv"但文件其实在 | 文件名拼成了 `<GSE>_clinical.csv`（应在前面），见 7.9 第 5 条 |
| 合并矩阵上做 DEG，不知道该用什么分组 | 必须向用户索要：把 `_clinical.csv` + `_grouping_candidates.txt` 交给他，见 7.9 与 `merge-batch-correction.md` 6.2 |
| `03b --stage=group` 报找不到合并临床表 | 上次 `--stage=merge` 没带临床表；重跑 `--stage=merge` |
| `03b --stage=group` 换了分组，但 RA 样本没回来 | 正常：样本范围取自合并临床表，被 `--group_na=drop` 剔掉的样本不会回来，见 7.9 第 6 条 |
| 02 步报 `临床表中无列 'source_name_ch1'`（合并矩阵场景） | 合并临床表里跨数据集的列名带 `<数据集>::` 前缀；报错会列出全部可选列，也可用 `group` 列 |
| 报告表格里数字列错位 | `sprintf` 按字符数对齐，CJK 占 2 列宽；取值过长时列宽自适应，见 7.10 |
| `分组文件第一列与表达矩阵样本名无一匹配`（用户给的表） | 两列写反了：必须第一列样本名、第二列分组（按位置取，不看表头），见 7.11 |
| 用户的分组表只覆盖部分样本 | 02 步会自动剔除其余样本并打印清单；合并阶段默认 `--group_na=error` 会直接报错，见 7.11 |
| `对比式 'case-control' 中的 ... 不在分组中。可用分组: 病例, 对照` | 分组名是中文；`--contrast` 要用实际组名，见 7.11 |
| 用户表里的中文组名变成乱码 | Excel 存成了 GBK 的 CSV；让他另存为"CSV UTF-8"，或组名用 ASCII，见 7.11 |
