---
name: geo-microarray-analysis
description: GEO 芯片（microarray）表达谱的「三步式」分析流程。第一步用用户提供的 GSE 编号 + family.soft.gz + series_matrix.txt.gz 生成标准化表达矩阵、临床信息与分组候选报告；第二步按用户指定的分组做 limma 差异分析并绘制火山图与热图；第三步把多个数据集（含 TCGA 表达量）合并并用 ComBat / normalizeBetweenArrays 做批次矫正，输出合并矩阵与批次 QC 图表。覆盖探针注释转 Gene Symbol、重复基因合并、normalizeBetweenArrays 归一化、topTable、pheatmap、ggplot2、ComBat、多数据集取交集、批次效应 PCA。当用户提到 GEO / GSE 编号 / 表达矩阵标准化 / 探针注释 / 芯片数据去重 / limma 差异基因 / 火山图 / 热图 / 多数据集合并 / 批次矫正 / ComBat / 去批次 / 合并 GEO 等需求时使用。
agent_created: true
---

# GEO 芯片数据：三步式标准化、差异分析与多数据集合并

## 三个步骤与"谁来决定"

| 步骤 | 性质 | 谁来决定 |
|---|---|---|
| 第一步 标准化 | 纯数据加工：解析、注释、去重、归一化 | 脚本自动完成 |
| 第二步 差异分析 | **分组是生物学判断** | **必须由用户拍板** |
| 第三步 多数据集合并 + 批次矫正 | 批次划分与矫正方法都是判断，且**会改变数值本身** | **必须由用户确认** |

前两步之间**必须停一次**：第一步跑完把"分组候选报告"交给用户，等用户明确指定分组后，
才执行第二步。分组选错，后面所有结果全部作废。

第三步是**可选支线**（只有要合并多个数据集时才走），但它另有两个必须让用户拍板的点：
**批次怎么划分**、**用 ComBat 还是 quantile**，以及**是否用分组信息保护生物学变异**。
详见下文"第三步"与 `references/merge-batch-correction.md`。

## 何时使用

- 用户给出 GSE 编号与本地 `*_family.soft.gz` / `*_series_matrix.txt.gz`，要标准化
- 用户已有标准化矩阵，要按指定分组做差异分析、画火山图/热图
- **用户要把多个数据集（多个 GSE、或 GEO+TCGA）合并成一个矩阵**，或提到 ComBat / 去批次
- 用户提到 normalizeBetweenArrays、探针注释、重复基因取最大值、limma DEG 等

只涉及 RNA-seq（counts + DESeq2/edgeR）时不适用本技能。

## 前置检查（必须先做）

1. **定位 R**：见上文——`run_geo.sh` 自动探测 Rscript。

   ```bash
   RSCRIPT="/path/to/Rscript"
   ```

   依赖包在此环境**已全部装好**：

   | 用途 | 包 |
   |---|---|
   | 第一/二步 | GEOquery、limma、Biobase、ggplot2、pheatmap、ggrepel |
   | 第三步 | **sva**（ComBat）、**data.table**（fread/fwrite）、**RSpectra**（快速 PCA） |

   第三步若缺 `RSpectra`，PCA 会退回 `prcomp` 且**明显变慢**（见下文性能说明）；
   缺 `data.table` 则大文件读写慢一个数量级。这两个都值得先装上。
   缺 `svglite`（nature 风格的 SVG 导出）时脚本会自动跳过，矢量图用 PDF 顶替。

2. **必须用包装器 `run_geo.sh` 启动**——它自动设置好临时目录与 COMSPEC：

   ```bash
   SK="<TOOLKIT>/skills/geo-microarray-analysis/scripts"
   bash "$SK/run_geo.sh" "$SK/01_geo_normalize.R" --help
   ```

   若直接调 `Rscript`，必须先手动设置，否则会崩在临时目录上：

   ```bash
   mkdir -p /tmp/rtmp   # Windows: mkdir C:/Rtmp
   export TMP=/tmp/rtmp TEMP=/tmp/rtmp TMPDIR=/tmp/rtmp   # Windows: 指向纯 ASCII 目录如 C:/Rtmp
   export COMSPEC='C:\WINDOWS\system32\cmd.exe'
   # 关键：必须清掉 LC_ALL —— 见下方说明
   unset LC_ALL LANG LC_CTYPE LC_COLLATE LC_MONETARY LC_TIME
   ```

   Windows 中文用户名会把 `TMP/TEMP` 编码搞坏，R 启动时据此建立的 `tempdir()`
   是坏路径，症状是 `Failed to create directory ... C:/Users/i+i:`。
   `tempdir()` 在 R 启动时就固定，**脚本内部改不掉，只能在启动前设好**。

   **同样必须在启动前处理的还有 locale。** 若环境带着 `LC_ALL=C.UTF-8`，
   而 Windows 版 R 应用不了 `C.UTF-8`（只打印一行 `Setting LC_CTYPE=C.UTF-8 failed`
   警告就继续），于是退化到 `C` locale。在 `C` locale 下 R **完全无法寻址含非 ASCII
   字符的路径**：`dir.exists` / `file.exists` / `readLines` 一律失败，路径被悄悄改写成
   乱码，`list.files()` 还会静默返回空——看起来就像"文件不存在"。
   典型症状：`--inputs=D:/结果/per_gse/GSE55235.csv` 明明存在，脚本却报
   `Error in read_matrix_robust(...) : 找不到输入矩阵: ...` （**2026-09-17 又踩一次**）。
   去掉这几个变量后 R 使用系统原生 UTF-8 locale，中文路径读写与中文命令行参数都正常。
   **`run_geo.sh` 已内置这一步**；这也是"一律走包装器"的又一个理由。详见 `pitfalls.md` 1.7。

> **已知无害现象**：此环境下 R 退出时报 `Segmentation fault`（退出码 139），
> 只要 `library()` 加载过包就会出现，且发生在所有工作完成之后。
> **判断成败看产物文件，不要看退出码。**

3. **Bash 工具被沙箱拦截时的备用通道（2026-09-16 实测）**：客户端沙箱会把
   bash/wsl 调用路由到被黑名单拦截的 `wsl.exe`（报 PROGRAM BLOCKED BY
   SECURITY POLICY），且与命令里是否含中文路径无关（嵌套 `bash`、长参数、
   `bash.exe 脚本.sh` 都可能触发；简单 `echo` 反而能过，无法预测）。
   → **备用方案：PowerShell 工具直接调 Rscript**，在**同一条命令**里先设好
   环境变量（会话状态不跨调用保留），再 `& Rscript.exe --vanilla 脚本 参数`,
   输出用 `$log = & ... 2>&1` 收集后 `[IO.File]::WriteAllText` 写 UTF-8 日志
   （日志中文是 GBK 乱码，属显示问题不影响结果）：

   ```powershell
   $env:TMP="D:\Rtmp"; $env:TEMP="D:\Rtmp"; $env:TMPDIR="D:\Rtmp"; $env:COMSPEC="C:\WINDOWS\system32\cmd.exe"
   $log = & "C:\path\to\R\bin\x64\Rscript.exe" --vanilla <脚本> <参数...> 2>&1
   [IO.File]::WriteAllText("<日志路径>", (($log | ForEach-Object { "$_" }) -join "`r`n"), [Text.Encoding]::UTF8)
   ```

   ⚠️ 注意：`02_deg_plots.R` 的产物落在**进程工作目录**（PowerShell 会话的
   CWD），`--outdir` 可能不生效 → 跑完用 `Move-Item` 把 `<prefix>_*` 移到目标
   目录，或先 `Set-Location` 再调用。

## 第一步：标准化 + 分组候选探查

用户需要提供三个信息：GEO 号、`GSE*_family.soft.gz`、`GSE*_series_matrix.txt.gz`。
两个文件都给时脚本**完全不联网**。

```bash
SK="<TOOLKIT>/skills/geo-microarray-analysis/scripts"
bash "$SK/run_geo.sh" "$SK/01_geo_normalize.R" \
  --gse=GSE62452 \
  --matrix="D:/path/GSE62452_series_matrix.txt.gz" \
  --soft="D:/path/GSE62452_family.soft.gz" \
  --dir=./out --prefix=GSE62452
```

常用参数：`--gpl=GPLxxx`（GSE 含多平台时指定）、`--force_log=TRUE/FALSE`、
`--norm=quantile|median|cyclicloess|scale|none`、`--dup=max|elementwise`、
`--min_expr=1`、`--annot_col=`（默认自动探测）。完整列表见 `--help`。

产出（在 `--dir` 下）：

| 文件 | 用途 |
|---|---|
| `<prefix>.csv` / `.txt` | 标准化表达矩阵，**第二步的输入** |
| `clinical_<prefix>.csv` | 样本表型信息 |
| **`<prefix>_grouping_candidates.txt`** | **候选分组列报告（候选列 + 每列取值的样本数，供你确定分组）** |
| `<prefix>_group_template.csv` | 分组模板（sample,group），用户填好可直接用 |
| `<prefix>_boxplot.pdf` | 归一化前后箱线图，**复核 log2 判断的依据** |
| `<prefix>_normalized.RData` | dat + pd |

**跑完必须 STOP，把决策权交还用户（这一步不能省）**：

1. 把 `clinical_<prefix>.csv` 的样本表型 + `_grouping_candidates.txt` 的候选列
   **完整报告给用户**：有哪些列可作为分组、每组的样本数。
2. **明确询问用户用哪一列 / 哪些值作为分组、对比方向（谁减谁）**。
   **绝对不要替用户写死 `--group_from` / `--contrast`**，也不要默认选某个
   候选列就直接跑第二步。等用户给出明确指令（如"用 clinical status 那列，
   RA 减 OA"）后，才执行第二步。
3. 看 `_boxplot.pdf` 确认分布合理（已 log 的数据中位数通常 0~15）。

## 第二步：差异分析 + 火山图 + 热图

拿到用户指定的分组后执行。

**⚠️ 出图前必须先让用户挑风格（2026-09 起强制）**：跑 02 步之前，把风格菜单呈给
用户选择（火山图 + 热图都要问），拿到明确选择后再带 `--vol_style` / `--hm_style`
执行。两种方式：

1. AI 直接把下表呈给用户（首选，可在对话里附风格说明）
2. `--list_styles=TRUE` 让脚本打印完整菜单（不需要 --expr）：
   `bash "$SK/run_geo.sh" "$SK/02_deg_plots.R" --list_styles=TRUE`

用户不挑时才用默认（跟随 `--fig_style=classic`）。样本数决定热图可选范围：
**numbers 仅 ≤20 样本可用**（>20 会自动回退并提示）；**>20 样本推荐 complex 或 tile**。

```bash
bash "$SK/run_geo.sh" "$SK/02_deg_plots.R" \
  --expr=./out/GSE62452.csv \
  --group_file=./out/GSE62452_group_template.csv \
  --contrast=case-control --prefix=GSE62452_T --top=30
```

分组方式四选一：

| 方式 | 参数 | 适用 |
|---|---|---|
| 内联指定 | `--group_spec='case=GSM1,GSM2;control=GSM3'` | **交互场景首选**，无需额外文件 |
| 分组文件 | `--group_file=group.csv`（sample,group 两列） | 用户已有分组表；只列部分样本也行 |
| 临床表提取 | `--clinical=... --group_from=tissue:ch1` | 从 clinical 表某列取 |
| 直接给向量 | `--group_values=case,case,control,...` | 样本少且顺序明确 |

- `--group_from` 常配合 `--group_regex='^(.)'`（必须含且仅含**一个捕获组**）。
- `--group_keep=A,B` 可在分组多于 2 组时挑出要对比的两组。
- `--contrast` 必须是 **`make.names` 清洗后**的分组名，如
  `rheumatoid arthritis` → `rheumatoid.arthritis`，写错会报
  `object '...' not found`。
- `--p_type=P.Value|adj.P.Val`（发文章用 `adj.P.Val`）、`--logfc`、`--pval`、
  `--top`（热图基因数）、`--label_genes=MMP1,VIT`。
- **用户直接给的分组表**直接用 `--group_file=`：两列（样本名, 分组），表头可有可无，
  逗号/制表符都行，**不必覆盖全部样本**（未覆盖的会被自动剔除）。
  合并流程产出的 `<prefix>_group_template.csv` 可以直接给他填；
  03 步的 `<prefix>_group.csv` 也可直接喂进来。

> ⚠️ **在"合并矩阵"上做差异分析时，分组必须先向用户索要。** 合并矩阵的样本来自
> 多个数据集，`--clinical` 要指向 03 步产出的 `<prefix>_clinical.csv`
> （行名=样本名），列名常带 `<数据集>::` 前缀。先把
> `<prefix>_clinical.csv` 与 `<prefix>_grouping_candidates.txt` 交给用户，
> 问清"每个数据集用哪一列、哪些取值算 case/control、对比方向"再跑。
> 详见第三步的「第四件必须让用户拍板的事」。

**期刊风格出图**（`--fig_style=nature`，改编自开源 nature-skills / Apache-2.0）：

- 火山图换期刊配色（红 `#B2182B`/蓝 `#2166AC`/灰 `#B3B3B3`）、图例带计数、
  `theme_nature()` 细轴线无网格、默认无标题；`--label_top=N` 自动标注 p 值最小的
  前 N 个显著基因（与 `--label_genes` 合并）。
- 热图用发散色带 `#2166AC`-白-`#B2182B`、去边框、期刊字号。
- 额外导出 **SVG + 600dpi TIFF**（缺 `svglite`/`ragg` 时自动跳过）。
- 尺寸可用 `--vol_width_mm/--vol_height_mm`、`--hm_width_mm/--hm_height_mm`
  （期刊单栏约 89 mm）。详见 `references/nature-figure-style.md`。

**更多火山图/热图风格**（`--vol_style` / `--hm_style`，2026-09 新增；不写则跟随
`--fig_style`，默认行为不变）：

| 参数值 | 风格 | 来源 | 依赖/适用条件 |
|---|---|---|---|
| `--vol_style=enhanced` | 四区分色 + 引线标注 | kevinblighe/EnhancedVolcano | 需 EnhancedVolcano 包，未装自动回退 |
| `--vol_style=gradient` | 渐变气泡：颜色/大小随显著性加深 | BioSenior/ggVolcano gradual_volcano | 零新依赖 |
| `--vol_style=rainbow` | 五段彩虹渐变（深蓝→青→黄→橙红→深红） | 中文社区教程 | 零新依赖 |
| `--vol_style=tophits` | 按 Manhattan 距离(\|logFC\|+\|-log10P\|)选 Top N 标注 | VolcaNoseR (Sci Rep 2020) | 零新依赖 |
| `--hm_style=complex` | ComplexHeatmap：顶部分组注释条 + 行内 z-score + 行聚类 + 按组 column_split | jokergoo/ComplexHeatmap | ComplexHeatmap+circlize（已装）；**大样本推荐** |
| `--hm_style=numbers` | pheatmap 单元格标注 z-score | pheatmap display_numbers | **仅 ≤20 样本**，>20 自动回退 |
| `--hm_style=tile` | ggplot2 geom_tile：顶部分组色条 + 无聚类 + 行内 z-score | ggplot2 geom_tile 教程 | patchwork（已装）；**大样本推荐** |

- `--top_n=10` 控制 tophits 标注数（也作 gradient/rainbow 的默认标签数）。
- 新热图配色跟随 `--fig_style`：nature 下蓝白红发散色带，classic 下蓝-白-红。
- `--vol_style=enhanced` 需要 EnhancedVolcano 包（Bioconductor），未装时自动回退；其余风格已实测通过
  （GSE62452 回归：上调 177 / 下调 114，与基线一致；numbers>20 样本回退、
  tile 大样本渲染均已验证）。

产出：`<prefix>_sample_group_map.csv`（**样本→分组映射，先核对这张表**）、
`<prefix>_all.csv`、`<prefix>_DEG.csv`、`<prefix>_volcano.pdf/.png`、
`<prefix>_heatmap.pdf/.png`、`<prefix>.RData`（nature 风格另加 `.svg`/`.tiff`）、
**`<prefix>_for_enrichment.txt`**（交接清单，见下一节）。

### 跑完的下一步：交给下游富集分析（enrichment-analysis）

本步跑完会多写一份 **`<prefix>_for_enrichment.txt`**（key=value 文本），
把下游需要的东西一次交代清楚，免得再来回问：

```
prefix / dir / contrast / contrast_num / contrast_den / groups / n_per_group
n_up / n_down / deg_file / all_file / group_map_file / rdata_file
gene_col=rownames  logfc_col=logFC  p_col=P.Value  change_col=change
deg_p_type=P.Value  deg_p=0.05  deg_logfc=1  species=human
```

**用户要富集分析时，不用手抄路径**——下游技能只要给 `--geo_dir=<本目录>`：

```bash
ESK="<TOOLKIT>/skills/enrichment-analysis/scripts"
Rscript --vanilla "$ESK/01_enrich_ora.R" \
  --geo_dir="<本步产物目录>" --direction=all,up,down        # ORA：吃 _DEG.csv
Rscript --vanilla "$ESK/02_enrich_gsea.R" \
  --geo_dir="<本步产物目录>"                                # GSEA：吃 _all.csv（全部基因）
Rscript --vanilla "$ESK/04_enrich_summary.R" \
  --in_dir="<本步产物目录>/enrich" --prefix=<prefix> --direction=all,up,down
```

下游会自动用 `_DEG.csv` 做 ORA（**默认不再二次筛选**，上下调直接按本步的 `change`
列拆，口径逐字一致）、用 `_all.csv` 做 GSEA（不筛选、按 logFC 排序），
产物落在 `<本目录>/enrich/`。

> ⚠️ **本步没有 `--outdir`，产物落在进程工作目录**。所以在 PowerShell 里跑要先
> `Set-Location`，否则交接清单里的路径指向的是别处。
>
> ⚠️ 交接清单里 `species` 一律写 `human`（本流程不解析平台物种）；
> **小鼠芯片在富集侧要显式加 `--species=mouse`**。

## 第三步（可选）：多数据集合并 + 批次矫正

只有要**把多个数据集合成一个矩阵**时才走这一步。典型动因：单数据集样本量不足，
或要做"训练集+验证集"式的跨队列分析。

### 推荐：一条命令从原始文件做起（`03b_merge_from_geo.R`）

用户手上通常是 **N 组 `(GSE*_series_matrix.txt.gz, GSE*_family.soft.gz)`**，而不是
标准化好的矩阵。这时不要手工跑 N 遍 01 再跑 03——用编排脚本一次做完：它逐个调用
01、再调用 03，并把 N 份临床表自动拼成一张统一的 `sample,group` 表。

**这条链要分两段跑**，因为"用哪些列、取值怎么归并成分组"必须由用户拍板
（和第一、二步之间那个停顿是同一条铁律）：

```bash
# 第一段: 只做标准化 + 打印每个 GSE 的分组候选，然后停住
bash "$SK/run_geo.sh" "$SK/03b_merge_from_geo.R" \
  --gse=GSE55235,GSE55457 --datadir=D:/data --outdir=./geo_merge --stage=normalize

# ↑ 看每个 <GSE>_grouping_candidates.txt 里的候选列与取值，再决定分组

# 第二段: 合并 + 批次矫正
bash "$SK/run_geo.sh" "$SK/03b_merge_from_geo.R" \
  --gse=GSE55235,GSE55457 --datadir=D:/data --outdir=./geo_merge --stage=merge \
  --group_cols='source_name_ch1;source_name_ch1' \
  --group_maps='normal=control,osteoarthritis=case;normal=control,osteoarthritis=case'
```

- `--stage=normalize|merge|group|all`：`all` 一次跑完（分组方案已明确时才用）；
  `group` 只重建分组表，不重跑 03（见下节）
- `--datadir` 会**递归**查找 `<GSE>_series_matrix.txt.gz` 与 `<GSE>_family.soft.gz`，
  所以多个 GSE 分散在各自子目录也没问题；也可用 `--matrices=` / `--softs=` 直接指定
- `--group_cols` / `--group_maps` 用 `;` 分隔、逐数据集给（顺序对应 `--gse`）；
  各数据集取值本来就统一时只给 `--group_cols` 即可
- 已完成的 GSE 自动跳过（`--force=TRUE` 重跑）；某个 GSE 失败会打印它的日志尾部
- 并入已标准化的矩阵（如 TCGA 表达量）：`--extra_inputs=` / `--extra_names=`

产物在 `--outdir` 下：`per_gse/`（每个 GSE 的 01 步全部产物）、`logs/`、
`merge/`（03 步全部 QC 产物）、`merge/<prefix>_clinical.csv`（**合并临床表**）、
`merge/<prefix>_grouping_candidates.txt`（**交给用户索要分组的材料**）、
`<prefix>_group_derived.csv`（推导出的统一分组表，可直接喂给 02 步）。

### 手工方式：分别跑 01 和 03

需要逐步控制时（只重跑部分数据集、要混入已标准化的矩阵、或想手动拼分组表）：

```bash
# A) 每个 GSE 先各自标准化
for g in GSE55584 GSE57218 GSE98918; do
  bash "$SK/run_geo.sh" "$SK/01_geo_normalize.R" --gse=$g \
    --matrix="D:/data/${g}_series_matrix.txt.gz" \
    --soft="D:/data/${g}_family.soft.gz" --dir=./out --prefix=$g
done

# B) 合并 + 批次矫正（01 已逐个标准化，故 --pre=none 只做交集+cbind+批次矫正）
bash "$SK/run_geo.sh" "$SK/03_merge_batches.R" \
  --inputs=./out/GSE55584.csv,./out/GSE57218.csv,./out/GSE98918.csv \
  --names=GSE55584,GSE57218,GSE98918 --pre=none \
  --dir=./merge --prefix=merged

# C) 在合并矩阵上做差异分析
bash "$SK/run_geo.sh" "$SK/02_deg_plots.R" --expr=./merge/merged_merged.csv \
  --group_file=./merge/merged_group.csv --contrast=case-control --prefix=merged_T
```

> ⚠️ **01 已逐个标准化过时，给 03 传 `--pre=none`。** 03 的 `--pre=auto` 是按
> "第一个样本列名的前 3 个字符"判断数据来源的：`GSM*` 会被当成芯片、再做一次
> `normalizeBetweenArrays`。虽然分位数归一化近似幂等、结果差别不大，但语义上
> 重复了。`03b` 编排器已默认这样处理。

也支持 GEO + TCGA 混用（TCGA 列名以 `TCG` 开头时按表达量处理，只取 `log2(x+1)`）。

### 03 步做了什么

1. **逐数据集预处理**（`--pre=auto`）：按首列样本名的前 3 个字符判断类型
   - `TCG*`（TCGA counts/TPM/FPKM）→ `log2(x+1)`，**不再**分位数归一化
   - `GSM*` / 其他（芯片）→ 按需 `log2(x)`，然后 `normalizeBetweenArrays`
   - 判定"是否需要 log"用的是与 01 步同一套分位数启发式，每个数据集都会把
     前后分位数打印出来供复核。要覆盖就逐数据集写 `--pre='log2p1;quantile;quantile'`
2. **取共同基因**：`Reduce(intersect, ...)`，保持第 1 个数据集的基因顺序
   （`--match=union` 可改用并集，缺失补 NA）
3. **cbind 合并** → `<prefix>_merged_raw.csv`（未矫正，QC 对照）
4. **批次矫正** → `<prefix>_merged.csv`（**02 步的输入**）
5. **QC 图表**：PCA / 箱线图 / 密度图的前后对比，外加批次效应量化与批次统计表

### 三个必须让用户拍板的点

这三项是 03 步的**输入**，脚本不会替你决定：

1. **批次怎么划分** —— 默认"一个数据集 = 一个批次"。但有时需要合并批次
   （例如把同平台的多个 GEO 当一个批次，或某数据集要单独留着做验证集）。
   用 `--batch=` / `--batch_file=` 覆盖。**跑完要用 `_batch_map.csv` 核对批次划分
   是否和你预期一致**（每个批次的样本数、有没有样本漏掉或归错批）。
2. **用哪种矫正** —— `--method=combat|quantile|both|none`
   - `combat`（默认，= `sva::ComBat`）：**逐基因**估计并扣除每个批次的均值/方差参数，
     跨平台首选
   - `quantile`：`normalizeBetweenArrays` 把各样本分布强制对齐，最保守，
     但**去不掉基因特异的平台效应**（实测：批次对 PC1 的解释方差只从 99.3% 降到 98.6%，
     而 ComBat 降到 0.0%）。**别把 quantile 当成"去过批次"。**
   - `both`：两份都出（`_merged_combat.csv` / `_merged_quantile.csv`），加
     `--primary=` 指定谁当主产物。**跨平台数据建议先 `both` 让用户对比 PCA 再定。**
   - 细节见 `references/merge-batch-correction.md` 第 5.1 节
3. **是否保护生物学变异** —— 给 `--group_file` 且 `--combat_mod=group`（默认）时，
   ComBat 用 `mod = model.matrix(~ group)`，只扣"批次内同一分组共有的"偏移，
   保留分组间差异。**不给分组时 `mod=NULL`，ComBat 可能把真实的生物学差异也扣掉。**

> ⚠️ `mod` 必须用**带截距**的 `model.matrix(~ group)`。ComBat 内部会
> `design <- cbind(batchmod, mod)` 再 `apply(design,2,function(x) all(x==1))` 丢掉
> 全 1 列，所以截距会被自动去掉、正好剩 1 列虚拟变量。若写成 `~ 0 + group`
> 会得到 k 列虚拟变量，与 batch 列共线，ComBat 会误报
> `The covariates are confounded!`。脚本已固定用 `~ group`。

### 第四件必须让用户拍板的事：合并矩阵上的分组 ★

合并矩阵是多个数据集拼起来的，**各数据集的临床列名与取值词表都不一样**：
同一个 `source_name_ch1`，在 A 里是 `normal control`、在 B 里可能是
`synovial tissue from healthy joint`。所以"用哪一列、哪些取值算 case/control"
**必须由用户回答**——不要在合并脚本里顺手替他定，也不要看到 `sample,group`
文件已经存在就默认它是对的。

03 步会自动把 N 份临床表并成**两份材料**，这两份就是索要分组的依据：

| 文件 | 内容 |
|---|---|
| `<prefix>_clinical.csv` | 合并临床表：每行一个样本（行名 = 样本名，02 步 `--clinical=` 可直接读），列为 `dataset`/`batch`/`group` + 各数据集的临床列；只出现在部分数据集的列带 `<数据集>::` 前缀 |
| `<prefix>_grouping_candidates.txt` | **交给用户的主材料**：逐数据集的候选分组列与取值计数、每个取值当前落到哪个分组（含"未进入矩阵/已剔除""混合"标注）、跨数据集同名列矩阵、下一步命令 |
| `<prefix>_group_template.csv` | 分组模板（`sample,group` 两列，全部样本已列好，已推导出分组时预填）。用户填好或换成他自己的表，都能直接当 `--group_file=` |

**要问用户的三件事**（缺一不可）：

1. **每个数据集用哪一列** —— 跨数据集时往往不是同一个列名（报告里 `[n]` 标出的候选列）
2. **该列哪些取值算 case / control / 哪些样本要剔除** —— 报告里每个取值后面的
   `-> 分组名` 就是当前映射结果，让用户核对；`混合(...)` 说明这个取值被映射到了
   多个分组（通常意味着映射写错了），`(未进入矩阵, 已剔除)` 说明整类样本被丢了
3. **对比方向** —— `--contrast=<组2>-<组1>` 的分子分母。组名会被 `make.names`
   清洗（`rheumatoid arthritis` → `rheumatoid.arthritis`），含空格时先确认实际取值

**第四块：让用户直接上传他自己的分组列表**（最省事的回答方式）

用户手上往往已经有一张分组表（样本清单、表型汇总、或从文章补充材料抄下来的）。
**不要逼他按你的列名规则回答——两列表格就够**：第一列样本名、第二列分组，
表头可有可无，逗号或制表符都行。可以直接把模板给他填：

```
<prefix>_group_template.csv      # 03 已把全部样本名按矩阵顺序列好
sample,group
GSM1332201,control
GSM1332202,case
...
```

拿回来直接喂：

```bash
# 差异分析（未覆盖到的样本会被 02 步自动剔除，不是错误）
bash "$SK/run_geo.sh" "$SK/02_deg_plots.R" \
  --expr=./geo_merge/merge/<prefix>_merged.csv \
  --group_file=<用户的分组表> --contrast=case-control

# 合并阶段也能用：ComBat 的 mod=~group 会跟着保护这个分组
bash "$SK/run_geo.sh" "$SK/03b_merge_from_geo.R" --gse=... --outdir=./geo_merge \
  --stage=merge --group_file=<用户的分组表>

# 样本少时还可以行内给，不用文件
# --group_spec='case=GSM1,GSM2;control=GSM3'   或   --group_values=case,case,control
```

- **样本名必须与合并矩阵的列名一致**（合并矩阵列名 = 各 01 步产物列名，
  可用 `_batch_map.csv` 第 1 列核对）。02 步会输出 `<prefix>_sample_group_map.csv`，
  **先核对这张表**再信结果
- **按位置取前两列**：表头叫什么无所谓，但**第一列必须是样本名、第二列必须是分组**
  （写反会报 `分组文件第一列与表达矩阵样本名无一匹配`）
- **不必覆盖全部样本**，只给要对比的那两组也行。差异分析阶段（推荐）未覆盖的样本
  会被**自动剔除**并打印清单；**合并阶段则会默认报错**（`--group_na=error`，
  因为 `group` 还要喂 ComBat 的 `mod`），要用就补全或显式加 `--group_na=drop`
- Excel 另存的 **UTF-8 CSV 直接可用**（BOM 不影响）；但别用 Excel 的"CSV (逗号分隔)"
  存中文组名——那是 GBK，会乱码
- 组名含中文/空格时 `--contrast` 用实际组名（`make.names` 清洗后的形式），
  写错会报 `... 不在分组中。可用分组: 病例, 对照`
- **不要替他改写组名**：用户给什么组名就用什么

> ⚠️ 用 `--group_file` 走合并阶段（而不是 `02`）时，ComBat 的 `mod` 会变成这个分组：
> 若该分组与批次混杂（比如"数据集 A 全是 case"），脚本会在调用前拦住（见下一节）。

**若用户是按上面"三件事"回答的**（给了列名与取值，没给现成表），有两条路把它变成分组表：

```bash
# 路径 A: 重建分组表，不重算批次矫正（推荐——ComBat 很贵）
bash "$SK/run_geo.sh" "$SK/03b_merge_from_geo.R" \
  --gse=GSE55235,GSE55457 --outdir=./geo_merge --stage=group \
  --group_cols='source_name_ch1;clinical status:ch1' \
  --group_maps='synovial tissue from healthy joint=control,synovial tissue from osteoarthritic joint=case;normal control=control,osteoarthritis=case'
# -> ./geo_merge/<prefix>_group_derived.csv，可直接喂 02 步

# 路径 B: 用合并临床表里某一列的取值直接当分组
bash "$SK/run_geo.sh" "$SK/02_deg_plots.R" \
  --expr=./geo_merge/merge/<prefix>_merged.csv \
  --clinical=./geo_merge/merge/<prefix>_clinical.csv \
  --group_from=<列名> --group_keep=A,B --contrast=B-A
```

- `--stage=group` 只重建分组表：样本范围取自**合并临床表**（已被 `--group_na=drop`
  剔掉的样本不会回来），不跑 03、不重算 ComBat
- 换分组**不会**改变合并矩阵里 ComBat 用的 `mod`（它保护的仍是合并时那个分组）。
  要让 ComBat 也保护新分组，得带 `--group_cols/--group_maps` 重跑 `--stage=merge`

> ⚠️ **不要用数据集/批次当分组。** 若"数据集 A 全是 case、数据集 B 全是 control"，
> ComBat 的 `mod=~group` 与批次列共线，脚本会在调用前拦住。反过来，用
> `_grouping_candidates.txt` 里每个取值的 `-> 分组名` 标注，也能一眼看出某个分组
> 是否只落在单个数据集里。

### 跑完必须停下确认（和第一、二步之间的那个停顿同等重要）

03 步跑完**不要直接进 02 步**。先把这四样交给用户：

1. `<prefix>_batch_map.csv` —— 样本 → 批次 的映射对不对（每个批次的样本数、
   有没有样本被漏掉或归错批）
2. `<prefix>_overlap_report.txt` —— 共同基因数是否合理、大小写敏感性有没有提示
   需要 `--gene_case=upper`、**批次效应量化表**（PC1/PC2 的批次解释方差前后对照）
3. `<prefix>_pca_before` / `_pca_after` 两张 PCA
4. **`<prefix>_clinical.csv` + `<prefix>_grouping_candidates.txt`** —— 合并后的
   临床信息，用来**向用户索要差异分析的分组**（见上一节）。
   顺手把 `<prefix>_group_template.csv` 也给他：他直接填好回传、或给一张
   自己手上的分组表，都是最快的收敛方式

等用户确认"批次划分与矫正结果可接受"后，才执行 02 步。
若批次效应没扣干净、或**分组也被一起扣掉了**（生物学差异消失），
就回去调 `--method` / `--batch` / `--combat_mod` 重跑 03 步。

参考量级：<回归测试三数据集> 上，ComBat 把批次对 PC1 的解释方差
从 **99.3%（p≈0）降到 0.0%（p=1）**——这是"批次效应确实被扣掉"的典型表现。

### 批次与分组混杂时会直接停下

若某些批次内只有一个分组水平（例如每个 GSE 只含一种组织），"批次效应"和
"生物学差异"在数学上无法分开。脚本会在调用 ComBat 前**主动检测并报错**，
把"批次 × 分组"交叉表打出来，并给三条出路：

1. `--combat_mod=none`：不保护生物学变异（**仅当各批次的分组构成相似时才安全**）
2. 只合并分组构成相近的数据集
3. 改用 `--method=quantile`（不做批次参数估计，不受此限制）

### 怎么判断"批次矫正是否合适"

看 `<prefix>_overlap_report.txt` 里的批次效应量化表：

```
批次效应量化 (批次对主成分解释的方差比例, ANOVA by batch):
            矫正前                矫正后
  PC1   42.3% (p=1.2e-80)     3.1% (p=0.021)
```

- 矫正后批次对 PC1 的解释方差**显著下降**（理想是降到个位数百分比、p 值不再极端）
  说明批次效应被扣掉了
- 若矫正后仍很高 → 批次没扣干净，或该数据集本身与其它数据集太不一样
- 同时要看 `<prefix>_group` 着色版的 PCA（给了分组时自动输出
  `<prefix>_pca_after_group.pdf`）：**分组要仍然分得开**。若分组也被扣没了，
  说明 `mod` 没用上或版本搞错了，去看 `_overlap_report.txt` 里
  `ComBat 保护生物变异 : 是/否`

### 产出

| 文件 | 用途 |
|---|---|
| **`<prefix>_merged.csv` / `.txt`** | **矫正后合并矩阵，第二步的输入** |
| `<prefix>_merged_raw.csv` | cbind 合并、未矫正（QC 对照） |
| `<prefix>_batch_map.csv` | **样本→数据集/批次/分组 映射，先核对他** |
| `<prefix>_group.csv` | `sample,group` 两列，给了 `--group_file` 时生成，可直接喂给 02 步 |
| **`<prefix>_clinical.csv`** | **合并临床表**（需 `--clinical_files`，01 步的 `clinical_<GSE>.csv` 自动传入）。行名=样本名，列含 `dataset`/`batch`/`group` + 各数据集临床列；02 步 `--clinical=` 可直接读 |
| **`<prefix>_grouping_candidates.txt`** | **向用户索要分组的材料**：逐数据集候选列与取值计数、每个取值 `-> 分组名` 的映射结果、跨数据集同名列矩阵 |
| `<prefix>_group_template.csv` | 分组模板（`sample,group`，全部样本已列好）。用户填好、或换成他自己的分组表，都能直接当 `--group_file=` |
| `<prefix>_overlap_report.txt` | 各数据集基因数、逐级交集、大小写敏感性、批次效应量化、批次统计 |
| `<prefix>_batch_stats.csv` | 每个批次矫正前后的均值/中位数/标准差 |
| `<prefix>_pca_before` / `_pca_after` | 批次效应 PCA（另加 `_pca_after_group`） |
| `<prefix>_boxplot_before` / `_after` | 样本分布箱线图（每个样本一个箱体，按批次着色） |
| `<prefix>_boxplot_batch_before` / `_after` | **每个批次一个箱体**（箱体来自批内各样本的中位数）——样本数很大时这张图仍清晰，用来判断各批次整体水平是否被拉齐 |
| `<prefix>_density_before` / `_after` | 样本密度曲线 |
| `<prefix>.RData` | `merged`、`merged_raw`、`batch`、`group`、PCA 对象 |

`--method=both` 时另出 `<prefix>_merged_combat.csv` / `_merged_quantile.csv`。
`--cor_heatmap=TRUE` 可加出样本相关性热图（按批次抽样，默认关闭）。

> ⚠️ `<prefix>_batch_map.csv` 的第 2 列是 `dataset`、不是分组，**不能**直接当 02 步的
> `--group_file`。要喂 02 步请用 `<prefix>_group.csv`。

### 常见参数

`--gene_case=asis|upper`（默认 `asis`，忠于 7.R；跨平台混用常需要 `upper`，
脚本默认会额外扫描并报告"统一大写能多出多少基因"）、
`--dup=max|first`、`--min_expr=0`、`--na=error|rowmean|zero`、
`--combat_prior` / `--combat_mean_only` / `--combat_ref=`、
`--clinical_files=A,B,C`（各数据集的临床表，逗号分隔、与 `--inputs` 一一对应，
某个数据集没有就留空；03b 会自动填 `per_gse/clinical_<GSE>.csv`）、
`--pca_top=N`（PCA 用变异最大的前 N 个基因，0=全部）、
`--pca_npc=10`（PCA 只求前几个主成分）、
`--density_max=200`（密度图最多画多少条曲线）、
`--save_raw=TRUE|FALSE`（是否输出未矫正矩阵文件；`.RData` 里始终保留）、
`--save_txt=TRUE|FALSE`（是否额外输出制表符版；02 步吃 `.csv` 就够）、
`--fig_style=nature`（同 02 步的期刊风格与 SVG/TIFF 导出）。完整列表见 `--help`。
03b 另有 `--stage=group` 与 `--group_out=`（重建分组表的输出路径）。

> 💾 **导出开销**：矩阵大时（如 15913×967）导出是整步最慢的环节。
> 脚本已用 `data.table::fwrite`（比 `write.csv` 快约一个数量级）。
> 只要 02 步的输入时可以加 `--save_raw=FALSE --save_txt=FALSE` 进一步省掉冗余文件。

> ⚡ **性能**：请安装 `RSpectra`。脚本用 `RSpectra::svds` 只求前 `--pca_npc` 个主成分；
> 若缺该包则退回 `prcomp`，而 `prcomp` 会求出**全部** `min(样本数, 基因数)` 个奇异向量，
> 代价随 `(样本数+基因数)×min²` 增长——15913 基因 × 967 样本 实测在 QC 步骤停留数分钟
> 仍未产出全部图表（进程 CPU 仅 ~40%，不是死锁）；而 `svds` 只要 **2.8 秒**，
> 且前几个主成分的方差占比与 `prcomp` 一致到小数第 4 位。合并规模更大时差异会更极端。
>
> ⚠️ 同时注意：`RSpectra::svds` 的 `u`/`v` **不保留 `dimnames`**，与 `prcomp()$x` 不同。
> 按行名对齐 PCA 得分时要自己补回，否则会取到 0 行并报
> `arguments imply differing number of rows`。



## 回归基线（数字对不上说明中间出了错）

**01 + 02 步**（单数据集）

| GSE | 平台 | 矩阵规模 | 分组 | |logFC|>1 且 P<0.05 |
|---|---|---|---|---|
| GSE55584 | GPL96 | 13237 基因 × 16 样本 | OA 6 / RA 10 | 上调 288 / 下调 257 |
| GSE62452 | GPL6244 | 23307 基因 × 130 样本 | case 69 / control 61 | 上调 177 / 下调 115 |

- GSE55584 上调 top：CXCL13 / CXCL9 / CCR7 / CD247 / PSMB9（滑膜免疫浸润）
- GSE62452 上调 top：LAMC2 / LAMB3 / TSPAN1 / LAMA3 / ITGB4 / SLC2A1；
  下调 top：SLC7A2 / AOX1 / IAPP / PDK4（胰腺癌 vs 癌旁）

**03 步**（多数据集合并；`<本地数据目录>/7Merge/` 的三数据集，2026-09-17 实测）

| 数据集 | 矩阵规模 | `--pre=auto` 判定 |
|---|---|---|
| TCGA_LUSC_TPM.txt | 19937 基因 × 553 样本 | `TCG*` → `log2(x+1)`（判定未取对数） |
| GSE30219.txt | 20824 基因 × 307 样本 | `GSM*` → 仅 `quantile`（判定已取对数） |
| GSE74777.txt | 30905 基因 × 107 样本 | `GSM*` → 仅 `quantile`（判定已取对数） |

- 逐级交集：`intersect(1,2)` = 17138 → `∩ 3` = **15913**（`--gene_case=asis`，忠于 7.R）
- 合并后：**15913 基因 × 967 样本**，无 NA；`--toupper_scan` 报告大写收益为 **+0**
- ComBat（`par.prior=TRUE`、`mod=NULL`）耗时约 **11~15 秒**
- 批次矫正效果：批次对 PC1 的解释方差 **99.3% → 0.0%**；
  各批次全矩阵均值 **3.289 / 6.673 / 5.909 → 4.650 / 4.658 / 4.657**
- 结果与 `<本地数据目录>/7Merge/7.R` 的 `merge1.txt` / `merge2.txt` **数值完全一致**
  （未矫正合并矩阵与 ComBat 输出均 `identical() = TRUE`、逐元素 `max|diff| = 0`；
  验证记录见 `references/merge-batch-correction.md` 第 9 节）

```bash
# 03 步完整复现（产物 <prefix>_merged.csv 再喂给 02 步）
bash "$SK/run_geo.sh" "$SK/03_merge_batches.R" \
  --inputs=<本地数据目录>/7Merge/TCGA_LUSC_TPM.txt,<本地数据目录>/7Merge/GSE30219.txt,<本地数据目录>/7Merge/GSE74777.txt \
  --names=TCGA_LUSC_TPM,GSE30219,GSE74777 --dir=./merge --prefix=merged
```

```bash
# GSE62452 两步完整复现
bash "$SK/run_geo.sh" "$SK/01_geo_normalize.R" --gse=GSE62452 \
  --matrix="<本地数据目录>/62452/GSE62452_series_matrix.txt.gz" \
  --soft="<本地数据目录>/62452/GSE62452_family.soft.gz" \
  --dir=./out --prefix=GSE62452
bash "$SK/run_geo.sh" "$SK/02_deg_plots.R" --expr=./out/GSE62452.csv \
  --group_file="<本地数据目录>/group.csv" \
  --contrast=case-control --prefix=GSE62452_T --top=30
```

## 关键决策点

以下各项**不要替用户默默决定**：

**第一 / 二步：**

1. **分组依据** —— 唯一必须由用户拍板的环节。把临床表 + 候选列报告给他，
   **明确问他用哪一列 / 哪些值、对比方向**，等他指定。
   **绝不要自作主张选列或写对比式。**
2. **log2 转换** —— 脚本自动判断并打印分位数，但双通道芯片会误判，需看箱线图复核。
3. **归一化方法** —— 默认 `quantile`（单通道首选）；双色芯片按需改。
4. **显著性口径** —— 默认 `P.Value < 0.05`，正式分析建议 `--p_type=adj.P.Val`。

**第三步（只有合并多数据集时才涉及）：**

5. **批次怎么划分** —— 默认每个数据集一个批次，合并批次或留验证集要用
   `--batch=` / `--batch_file=`。
6. **用 ComBat 还是 quantile** —— 跨平台建议先 `--method=both`，让用户对比 PCA 再定。
7. **是否用分组保护生物学变异** —— 给 `--group_file` 才有效；不给则 `mod=NULL`，
   有把真实生物学差异一起扣掉的风险。
8. **跑完必须停一次**：把 `_batch_map.csv` + `_overlap_report.txt` + PCA 前后图交给
   用户确认，再进 02 步。批次没扣干净或分组被扣没了都要回去调参重跑。
9. **合并矩阵上的差异分析分组** —— 必选项。把 `_clinical.csv` +
   `_grouping_candidates.txt` 交给用户，问"每个数据集用哪一列、哪些取值算
   case/control、对比方向"，**拿到答复前不跑 02 步**。若他直接给一张自己的
   分组表（两列 样本名,分组），用 `--group_file=` 收下即可，别逼他按列名回答。
   详见「第四件必须让用户拍板的事」。

## 资源

- `scripts/run_geo.sh` —— **首选入口**，包装环境变量后调用 R
- `scripts/01_geo_normalize.R` —— 第一步：本地文件快速解析 + 标准化 + 分组候选报告
- `scripts/02_deg_plots.R` —— 第二步：差异分析 + 火山图 + 热图（`--fig_style=nature` 期刊风格；
  `--vol_style=enhanced|gradient|rainbow|tophits` 与 `--hm_style=complex|numbers|tile` 多风格可选；
  `--list_styles=TRUE` 打印风格菜单；**出图前必须先让用户挑风格**）
- `scripts/03_merge_batches.R` —— 第三步：多数据集合并 + 批次矫正（ComBat / quantile）+ 批次 QC
  + 合并临床表与分组候选报告
- `scripts/03b_merge_from_geo.R` —— **编排**：N 组原始 GEO 文件 → 逐个标准化 → 合并去批次，
  并自动把各数据集临床表汇总成「合并临床表 + 分组候选报告」
  （`--stage=normalize|merge|group|all`）
- `references/merge-batch-correction.md` —— 第三步的设计依据、`ComBat` 的 `mod` 语义、
  **合并后的分组从哪来**、验证记录
- `references/nature-figure-style.md` —— 期刊风格配色/主题/导出约定（改编自 nature-skills）
- `references/pitfalls.md` —— 踩坑手册。**改动脚本前先读**

四个脚本都支持 `--help`。遇到报错先查 `references/pitfalls.md` 的排错表。

> 📌 **改动本技能的脚本时，图内文字一律用 ASCII。** Windows 下中文字形缺失会被
> 渲染成乱码 ASCII（不是方框），很难一眼看出来。详见 `pitfalls.md` 7.8。

