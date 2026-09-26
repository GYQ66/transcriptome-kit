# INSTALL —— 环境与依赖

## 0. 一键环境检查（先跑这个）

```bash
python check_env.py
```

- 检测 **R 是否安装**（含 Windows 多盘符 / macOS / Linux 常见位置）；
- **没装 R 时**：打印官方 + 清华镜像下载地址，**询问是否自动下载**最新安装包
  （默认否，需明确确认；Windows 可再确认后静默安装）；
- 检测内置 GMT 数据完整性；
- `--check-packages`：再用 Rscript 检查全部 R 包依赖，缺失时给出可复制的安装命令；
- `--mirror tuna`：走清华镜像；`--no-download`：只打印地址不询问；
- `--print-download-url`：只解析并打印最新 R 安装包直链。

以下章节是手动安装的完整说明（check_env.py 缺什么会告诉你对应装什么）。

## 1. 基础运行时

### R（geo / enrichment / survival / string / venn 需要）

- 版本 **≥ 4.2**（在 R 4.5.2 上开发与回归验证）。
- 下载：https://cran.r-project.org （Windows 用户选 base binary）。
- `run_*.sh` 包装器按以下顺序探测 Rscript，任一命中即可：
  1. 环境变量 `RSCRIPT`（如 `RSCRIPT="C:/Program Files/R/R-4.4.1/bin/x64/Rscript.exe"`）
  2. PATH 里的 `Rscript`
  3. Windows 常见安装位置（`C:/Program Files/R/R-*/bin/x64/` 等）

### Python（msigdb-pathway-gene-lookup 需要）

- 版本 **≥ 3.8**，**纯标准库，零第三方依赖**。
- 任何系统 Python / conda Python 均可。首次运行自动建索引（约 10 秒，仅一次）。

## 2. R 包依赖（按技能）

安装方式：CRAN 包用 `install.packages()`，Bioconductor 包用
`BiocManager::install()`。国内网络可用镜像：

```r
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
```

### geo-microarray-analysis

| 步骤 | 包 |
|---|---|
| 01 标准化 / 02 差异 | GEOquery, limma, Biobase, ggplot2, pheatmap, ggrepel |
| 03 合并 + 批次矫正 | sva（ComBat）, data.table, RSpectra（推荐，大矩阵 PCA 提速百倍） |
| 可选风格 | ComplexHeatmap + circlize（热图 complex 风格）、patchwork（tile 风格）、EnhancedVolcano（火山 enhanced 风格）、svglite + ragg（SVG/TIFF 导出） |

### enrichment-analysis

| 用途 | 包 |
|---|---|
| 核心 | clusterProfiler, org.Hs.eg.db, enrichplot, DOSE, fgsea |
| 出图 | ggplot2, ggprism, **gground**（GitHub 包：`devtools::install_github("dxsbiocc/gground")`） |
| 小鼠（可选） | org.Mm.eg.db |

### survival-km-analysis

`survival, survminer, ggplot2, ggpubr, ragg`（survminer 自带 maxstat 依赖）。

### string-ppi-network

| 类别 | 包 |
|---|---|
| 必需 | igraph, ggraph, ggplot2, tidygraph, graphlayouts |
| 建议 | ggrepel, ragg, svglite |
| 可选 | STRINGdb（仅 `--engine=stringdb` 引擎用；默认 `api` 引擎不需要） |

### venn-diagram

`VennDiagram`（自带 futile.logger 依赖）、`png`。

### 一键安装脚本（可选）

```r
cran <- c("GEOquery","limma","Biobase","ggplot2","pheatmap","ggrepel",
          "data.table","sva","RSpectra","clusterProfiler","enrichplot","DOSE",
          "fgsea","ggprism","survival","survminer","ggpubr","ragg","igraph",
          "ggraph","tidygraph","graphlayouts","ggrepel","svglite","VennDiagram",
          "png","patchwork","circlize","devtools")
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("org.Hs.eg.db","ComplexHeatmap"))
install.packages(setdiff(cran, rownames(installed.packages())))
if (!requireNamespace("gground", quietly=TRUE))
  devtools::install_github("dxsbiocc/gground")
```

> ⚠️ `enrichment-analysis` 里 `ComplexHeatmap` 属 geo 技能可选包，装不装不影响
> enrichment 本体；`gground` 必须装，否则富集美化主图（03 步）跑不了。

## 3. 内置数据（无需下载）

- `skills/msigdb-pathway-gene-lookup/data/*.gmt`：MSigDB v2026.1.Hs 十大集合
  symbols 版（H + C1…C9），共 35,361 个基因集，约 29 MB。
  - enrichment 的 MSigDB 富集（`--gmt_sets=H` 等）会自动复用这批 GMT；
    也可用环境变量 `ENRICH_GMT_DIR` 指向自己的 GMT 目录。
- KEGG 在线富集（`--kegg=online/auto`）需要联网；离线时自动退到本地 GMT 的
  KEGG 子集（注意那是 KEGG_LEGICY 186 条，不含 2020 年后新增通路）。

## 4. 离线 / 内网环境

| 步骤 | 离线方案 |
|---|---|
| MSigDB 检索 | 完全离线（数据内置） |
| GO/MSigDB 富集 | 完全离线（GMT 内置） |
| KEGG 富集 | 自动退 GMT 子集 |
| GEO 标准化 | 给齐 family.soft.gz + series_matrix.txt.gz 即完全离线 |
| STRING 网络 | `--engine=local` + 官方 `protein.info/links` 文件，或 `--engine=api` 在线 |

## 5. 已知平台注意事项（包装器已自动处理）

- **Windows + 中文用户名**：R 的 `tempdir()` 会因 TMP/TEMP 带中文而损坏 →
  包装器自动把 TMP/TEMP/TMPDIR 指到纯 ASCII 目录（默认 `C:/Rtmp`）。
- **Windows locale 陷阱**：`LC_ALL=C.UTF-8` 会让 Windows R 退化成 C locale，
  含非 ASCII 路径完全无法寻址 → 包装器自动 unset 相关变量。
- **R 退出码 139（Segmentation fault）**：Windows 上加载过包的 R 退出时常见，
  发生在所有工作完成之后，**判断成败看产物文件，不要看退出码**。
- **图内中文**：缺中文字形的环境会把图内文字渲染成乱码，所有技能图内文字一律 ASCII。
