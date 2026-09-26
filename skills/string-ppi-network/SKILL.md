---
name: string-ppi-network
description: 用 STRING 数据库构建蛋白互作（PPI）网络的参数化流程 —— 基因列表（可带 logFC）→ symbol 映射到 STRING_id → 取互作边（combined_score）→ igraph 建图 → ggraph 出图（节点颜色=logFC，线宽=互作置信度），始终同时输出 PDF+PNG+SVG 三格式。默认生成自包含 HTML 交互式布局编辑器：浏览器里拖动节点摆位 → 导出坐标 CSV → --positions 回读重跑即按摆位出图（manual 布局）。取数支持三个引擎：STRINGdb（Bioconductor 包，与参考代码完全一致）、STRING REST API（不依赖任何 Bioconductor 包，推荐）、本地 STRING 数据文件（完全离线）。参考代码来自「STRINGdb：在R语言中构建STRING蛋白互作网络」（复现 Nature s41586-024-08533-3 Fig.1.d：KBTBD4 突变型 vs 野生型显著下调蛋白的互作网络）。当用户提到 STRING / STRINGdb / 蛋白互作网络 / PPI 网络 / 互作网络 / protein-protein interaction / get_interactions / STRING 映射 / combined_score / 置信度评分 / 互作置信度 / 节点颜色 logFC / 线宽表示置信度 / stress 布局 / ggraph 网络图 / 蛋白网络图 / hub 蛋白 / Cytoscape graphml / 物种 9606 / 用差异蛋白画互作网络 / 拖动节点 / 拖拽摆位 / 交互式布局 / 手动布局 / positions / layout editor / SVG 输出 / 三格式输出 时使用。
agent_created: true
---

# STRING 蛋白互作网络：查询 → 建图 → 可视化

参考代码：`<参考代码路径>`
（来源：微信公众号「华哥生信 / Bioinformation」，复现 Nature Fig.1.d）

参考代码的 5 步全部覆盖，另外把「取数方式」抽成了可切换的引擎：

| 文件 | 干什么 |
|---|---|
| `scripts/01_string_network.R` | **唯一主脚本**：读基因表 → 映射 → 取边 → 建图 → 出图 → 落盘报告（v1.2） |
| `scripts/lib_string_common.R` | 公共库：参数解析、基因表读取、别名表、边表规范化、配色、绘图设备、格式归一 |
| `scripts/lib_layout_editor.R` | 交互式布局编辑器（自包含 HTML）生成器：拖节点 → 导出坐标 CSV |
| `scripts/run_string.sh` | **首选入口**，包装 PATH / TMP / COMSPEC / locale / 路径转换 |

> 参考代码里 `library(tidygraph)` 根本没被用到（且不是常规 CRAN 包），本 skill 不依赖它。
> `tidyverse` 也不需要 —— 只用到 `igraph` + `ggraph` + `ggplot2`。

## 三个取数引擎（关键决策）

| `--engine=` | 取数方式 | 适用场景 | 实测 |
|---|---|---|---|
| **`api`（推荐）** | `https://version11_5.string-db.org/api/tsv/` 的 `get_string_ids` + `network` | 只想尽快出图 | ✅ 13 基因 31 边，秒级 |
| `stringdb` | Bioconductor `STRINGdb` 包，与参考代码逐行对应 | 要复刻参考代码的原路径；小物种很快 | ✅ 已用大肠杆菌 K-12（`--species=511145`）跑通全链路：6/6 映射、12 条边<br>⚠️ 人类数据首次要下 20.7 MB aliases + 72 MB links，慢速网络 ~80 KB/s + 易断 + `load()` 要吞 1200 万行 → **人类建议直接用 `api`** |
| `local` | 读本地 `protein.info.*.txt.gz` + `protein.links.*.txt.gz` | 完全离线 / 内网 | ✅ 逻辑已验证（用合成文件） |
| `auto`（默认） | cache_dir 里已有 STRINGdb 缓存 → 走 `stringdb`；否则走 `api` | 日常 | — |

**`api` 与 `stringdb` 查的是同一个 STRING 11.5 数据库，结果同源**，只是 `stringdb` 多一层本地缓存。
两个引擎都不需要用户自己比对 —— `--engine=auto` 的表现就是「有缓存用包，没缓存用 API」。

### 网络连通性备忘（一台受限网络机器上的实测，供参考）

| 目标 | 结果 |
|---|---|
| `https://string-db.org/api/...`（当前版本 v12） | ❌ `SSL connect error` |
| `https://version11.string-db.org/api/...` | ❌ `SSL connect error` |
| `https://version12_0.string-db.org/api/...` | ❌ `SSL connect error` |
| **`https://version11_5.string-db.org/api/...`** | ✅ **可达**（参考代码用的正是 v11.5，运气好） |
| `https://stringdb-static.org/download/protein.info.v11.5/...`（1.8 MB） | ✅ 可达 |
| `https://stringdb-downloads.org/download/protein.aliases.v11.5/...`（20.7 MB） | ⚠️ 时通时断 |
| `https://bioconductor.org/...` | ❌ `SSL connect error`（BiocManager 装不了东西） |
| `https://mirrors.westlake.edu.cn/bioconductor/...` | ✅ 可达（镜像可用来装 STRINGdb） |

→ 结论：受限网络下**默认用 `--engine=api`**（先自测连通性）；网络通畅时 `auto` 即可。

## 前置检查

1. **定位 R**：任意 R ≥ 4.2（包装器自动探测；也可用环境变量 `RSCRIPT` 指定）。**必须**通过包装器启动：

   ```bash
   SK="<TOOLKIT>/skills/string-ppi-network/scripts"
   bash "$SK/run_string.sh" "$SK/01_string_network.R" --genes=... --outdir=...
   ```

   `run_string.sh` 做四件事，缺一个就跑不起来（与 `run_wgcna.sh` 同源）：
   - 防御性补 PATH（精简 Git 环境可能缺 coreutils）
   - `TMP/TEMP/TMPDIR` 改指纯 ASCII 目录（中文用户名会让 R 的 `tempdir()` 不可写）
   - 导出 `COMSPEC`（否则 R 的 `shell()` 报 `'/c' not found`）
   - `unset LC_ALL LANG ...`（否则 Windows R 退化为 `C` locale，含中文的路径全部寻址失败）
   - 把 POSIX 路径转成 `D:\...`（否则 Windows 版 Rscript 段错误退出，无产物无报错）

2. **依赖**（安装命令见套件根目录 INSTALL.md）：
   - 必需：`igraph`、`ggraph`、`ggplot2`、`tidygraph`、`graphlayouts`（stress 布局）
   - 建议：`ggrepel`（标签防重叠）、`ragg`（PNG 渲染）、`svglite`（SVG 矢量；缺省退回 `grDevices::svg`）
   - 可选：`STRINGdb`（只有 `--engine=stringdb` 用得到）

   > STRINGdb 安装备忘（Bioconductor 主站不可达时）：可从国内镜像下 Windows 二进制本地安装，
   > 依赖包走 CRAN 镜像。

3. **已知无害现象**（不要当故障处理）：
   - `package 'ggraph' was built under R version 4.5.3` 之类的版本警告
   - R 退出码可能是 139（Segmentation fault），发生在所有产物写完之后。
     **判断成败看产物文件，不要看退出码。**

## 标准工作流

### 例 1：复刻参考代码的 Nature 图（13 个基因 + 示例 logFC）

```bash
SK="<TOOLKIT>/skills/string-ppi-network/scripts"
bash "$SK/run_string.sh" "$SK/01_string_network.R" \
  --genes=RCOR3,CoREST,LSD1,RCOR2,MACROH2A2,H2AC7,HDAC2,MDC1,CENPC,UBTF,H2AX,CTCF,MIER3 \
  --logfc=-2.0,-1.5,-1.4,-1.9,-1.8,-0.6,-0.5,-0.7,-0.4,-0.3,-0.6,-0.8,-0.9 \
  --color_limits=-2,0 --engine=api \
  --outdir=./string_out --prefix=kbtbd4
```

回归基线（见 `references/pitfalls.md`）：映射 13/13、原始边 46、过阈值（400）边 31、
单连通分量、边权 417~999、节点 13 个。图内标签用的是**输入时的写法**（CoREST/LSD1/…），
和参考原图一致。

### 例 2：从差异分析结果表出发（最常用）

```bash
bash "$SK/run_string.sh" "$SK/01_string_network.R" \
  --table=D:/data/deg.csv --score_threshold=400 \
  --outdir=./string_out --prefix=deg --format=png
```

`--table` 会自动猜列名：基因列认 `gene/gene_symbol/symbol/id/ID...`，
logFC 列认 `logFC/log2FC/log2FoldChange/avg_log2FC/FC...`。猜不到就用
`--gene_col=` / `--logfc_col=` 显式指定。

### 例 3：只要网络结构，不做 logFC 着色

```bash
bash "$SK/run_string.sh" "$SK/01_string_network.R" \
  --genes_file=my_genes.txt --node_colour="#D73027" --outdir=./string_out
```

不给 logFC 时：`--prune` 自动视为 `none`（不会把节点全剪掉），节点统一用 `--node_colour`。

### 例 4：完全离线（内网机器）

```bash
bash "$SK/run_string.sh" "$SK/01_string_network.R" \
  --genes_file=my_genes.txt --engine=local \
  --info_file=D:/string/9606.protein.info.v11.5.txt.gz \
  --links_file=D:/string/9606.protein.links.v11.5.txt.gz \
  --outdir=./string_out
```

离线文件从 `https://stringdb-static.org/download/` 取。`protein.links.v11.5` 人类约 72 MB。

### 例 5：交互式布局微调（拖节点摆位）

**编辑器默认生成**：不加任何参数的普通运行就产出 `*_layout_editor.html`（连同
pdf/png/svg 最终图）。想按自己摆的布局出图，再走两步：

```bash
SK="<TOOLKIT>/skills/string-ppi-network/scripts"

# 第 1 步（人工）：浏览器打开 *_layout_editor.html，按住节点拖动摆位，
#                  点「导出坐标 CSV」下载 <prefix>_positions.csv
#                  （不想拖也可以直接手改 <prefix>_positions_auto.csv）

# 第 2 步：坐标喂回来，按 manual 布局重新出 pdf/png/svg
bash "$SK/run_string.sh" "$SK/01_string_network.R" \
  --table=D:/data/deg.csv --outdir=./string_out --prefix=deg \
  --positions=./string_out/deg_positions.csv
```

要点：
- 编辑器是**单文件自包含 HTML**（原生 SVG + JS，无 CDN 无依赖），内网/离线都能用；
  节点颜色 = logFC、线宽 = 置信度，与最终图一致。不想要就 `--interactive=0`。
- `--positions` 重跑时编辑器会以你导入的摆位为起点重新生成，可以反复微调迭代。
- 坐标体系与 R 完全一致（SVG 的 y 向下已在编辑器内部翻转处理），「导出 → 回读」无损往返，
  实测往返误差 ~1e-16。
- `--positions` 按 STRING_id 优先、symbol 兜底逐节点匹配；没覆盖到的节点保留自动布局坐标
  （seed 兜底）并告警。报告里会写 `布局来源: manual（覆盖 N/M 个节点）`。

## 参数速查

### 输入（三选一）

| 参数 | 说明 |
|---|---|
| `--genes=A,B,C` | 内联基因列表 |
| `--genes_file=f.csv` | 基因文件（第一列，或 `--gene_col` 指定列） |
| `--table=deg.csv` | 含基因列 + logFC 列的结果表（**推荐**，一次到位） |
| `--logfc=-2,-1.5,...` | 与 `--genes` **一一对应**的内联 logFC，个数不一致会直接报错 |
| `--logfc_file=fc.csv` | 单独文件，行数须与基因数一致 |

### STRING 查询

| 参数 | 默认 | 说明 |
|---|---|---|
| `--species` | `9606` | 也接受 `human` / `mouse` / `rat` / `zebrafish` / `fly` / `worm` / `yeast` |
| `--version` | `11.5` | 参考代码用 11.5；受限网络下可能只有部分版本域名可达 |
| `--score_threshold` | `400` | combined_score 阈值，**0–1000 口径**。400 中高可信、200 低可信 |
| `--network_type` | `functional` | 或 `physical`（只看直接物理结合） |
| `--add_nodes` | `0` | >0 时让 STRING 额外补充 N 个邻居蛋白；此时自动改成 `--prune=none`（除非显式给了 `--prune`） |
| `--engine` | `auto` | `auto` / `stringdb` / `api` / `local` |
| `--timeout` / `--retries` | `120` / `2` | 单次请求超时（秒）/ 重试次数 |

### 基因名处理

| 参数 | 默认 | 说明 |
|---|---|---|
| `--alias_fix` | `auto` | `auto`：**只救映射失败的基因**（两阶段映射，原名先试，失败的才换别名）；`all`：无条件全表替换；`none`：不替换 |

内置别名表覆盖参考代码踩过的 4 个（`CoREST→RCOR1`、`LSD1→KDM1A`、
`MACROH2A2→H2AFY2`、`H2AX→H2AFX`）+ 常见旧名（`53BP1→TP53BP1`、`KAP1→TRIM28`、
`NBS1→NBN`、`BRG1→SMARCA4`、`G9A→EHMT2`、`HP1A→CBX5`、`p53→TP53` 等 50 余条）。

> ⚠️ **实测反直觉的一点**：参考代码注意事项说这 4 个名字「映射会失败」，但在
> **STRING 11.5** 里 `CoREST / LSD1 / MACROH2A2 / H2AX` **全都能直接映射上**
> （13/13 命中，别名表一次都没触发）。所以别照着笔记里的结论去改基因名 ——
> 脚本会实打实告诉你哪几个没映射上，以它为准。

### 建图

| 参数 | 默认 | 说明 |
|---|---|---|
| `--prune` | `logfc` | `logfc`：只保留有 logFC 的节点（**等同参考代码的剪枝**）；`none`：全留 |
| `--keep_isolated` | `1` | 把没有任何边的输入基因也加进图。**参考代码会静默丢掉这些节点**，本 skill 默认补回 |

### 绘图

| 参数 | 默认 | 说明 |
|---|---|---|
| `--layout` | `stress` | `stress`/`fr`/`kk`/`dh`/`circle`/`graphopt`/`lgl`/`mds` |
| `--seed` | `42` | 参考代码就是 `set.seed(42)`，固定种子保证可复现 |
| `--palette` | `nature` | `nature` = 参考图原始 7 色（`#08306B→#FFFFCC`）；另有 `viridis`/`plasma`/`rdbu`/`bluered2`/`spectral`/`bluyl`/`blues`/`reds`/`magma`/`inferno` 等 |
| `--colors` | 空 | 自定义渐变色（逗号分隔 hex），优先于 `--palette` |
| `--color_limits` | `auto` | `auto` 按 logFC 实际范围；**要对齐参考原图请写 `--color_limits=-2,0`**（参考代码就是硬编码 -2,0）。超范围值会被压到边界色，不会变灰 |
| `--node_size` / `--label_size` | `6` / `3` | 同参考代码 |
| `--edge_width_range` | `0.1,1.2` | 线宽映射范围，同参考代码 |
| `--edge_colour` / `--edge_alpha` | `grey70` / `0.9` | 同参考代码 |
| `--legend` | `bottom` | `bottom`/`right`/`left`/`top`/`none` |
| `--color_label` | 空 | 默认图例标题是 `log₂ fold-change`（与参考代码同） |
| `--node_colour` | `#2171B5` | **没有** logFC 时的统一节点色 |

### 输出

| 参数 | 默认 | 说明 |
|---|---|---|
| `--outdir` / `--prefix` | `./string_out` / `ppi` | 产物目录与前缀 |
| `--format` | `pdf,png,svg` | **三格式始终同时输出**。非法值（如 `jpg`）告警并忽略，绝不产出扩展名与内容不符的文件 |
| `--width` / `--height` / `--dpi` | `7` / `7` / `300` | 英寸 / DPI |
| `--interactive` | `1` | 默认生成 HTML 布局编辑器 + 坐标快照；`0` = 关闭（纯一键出图） |
| `--positions` | 空 | 布局编辑器导出的坐标 csv（列：symbol 或 STRING_id + x + y），给了就走 manual 布局 |
| `--graphml` | `0` | `1` 时额外导出 `.graphml`，可直接拖进 Cytoscape |
| `--log_file` | 空 | 把日志同时写进文件 |

## 输出清单

以 `--prefix=ppi` 为例：

| 文件 | 内容 |
|---|---|
| `ppi_ppi_network.pdf` | **主图（PDF）**：节点色=logFC、线宽=STRING 置信度、stress/manual 布局 |
| `ppi_ppi_network.png` | **主图（PNG，300 dpi）**：与 PDF 同图同源 |
| `ppi_ppi_network.svg` | **主图（SVG 矢量）**：与 PDF 同图同源，AI/浏览器可直接编辑 |
| `ppi_layout_editor.html` | 自包含 HTML 布局编辑器（默认生成）：浏览器拖节点 → 导出坐标 |
| `ppi_positions_auto.csv` | 初始自动布局坐标快照（给过 `--positions` 时不再覆盖写） |
| `ppi_nodes.csv` | 节点表：STRING_id、symbol、logFC、degree（按 degree 降序） |
| `ppi_edges.csv` | 边表：from/to、两端 symbol、combined_score、score_norm（原值就是 `score_norm = combined_score / max(combined_score)`） |
| `ppi_unmapped_genes.csv` | **没映射上的基因 + 建议别名**（这些基因不会出现在图里，也不会报错） |
| `ppi_string_report.txt` | 人读报告：映射统计、网络统计、**布局来源**、绘图参数、top10 hub、**交付前必须核对**清单 |
| `ppi_run_params.txt` | 全部实际参数，可复现 |
| `ppi_network.graphml` | 仅 `--graphml=1` 时产出 |
| `_cache/` | API 原始 TSV 响应 / STRINGdb 数据文件缓存 |

## 交付前必须核对

1. **映射成功率**：`映射成功 N / M`，以及 `*_unmapped_genes.csv` 里有没有东西。
   参考代码的坑就是 `removeUnmappedRows=TRUE` 会**静默剔除**，图里少几个节点看不出来。
2. **别名替换**：日志里 `[alias]` 那几行有没有不该有的替换。
3. **阈值**：默认 400。若最终边数很少或网络碎片化，先下调到 200 看是否连通性改善
   （参考代码注意事项第 4 条）。
4. **颜色区间**：`logFC 颜色区间` 那行。`auto` 是实际范围，图表数值不能和参考图的
   `-2~0` 直接对比 —— 要对比就显式 `--color_limits=-2,0`。
5. **孤立节点/连通分量**：日志会报 `连通分量` 与 `孤立节点`。这是真实拓扑，不是脚本问题；
   想补连接要么降阈值、要么 `--add_nodes` 让 STRING 带邻居进来。
6. **logFC 数值来源**：参考代码里的 logFC 是**示例数据**，不是原文献真实值。
   复现真图必须换成自己的差异表达结果 —— 不要在交付里把示例值当结论。
7. **`--add_nodes>0` 时**：图里会出现输入基因之外的蛋白（它们没有 logFC，颜色是灰色），
   正文里要说清楚。
8. **布局来源**：报告里 `布局来源` 一行。为 `manual` 时核对 `--positions` 的节点覆盖数
   是否 = 节点总数（没覆盖到的节点用的是 seed 兜底的自动坐标，不是你拖的位置）。

## 资源

- `scripts/run_string.sh` —— **首选入口**，包装 PATH / TMP / COMSPEC / locale / 路径转换
- `scripts/01_string_network.R` —— 主脚本，支持 `--help`
- `scripts/lib_string_common.R` —— 公共库（参数解析、读表、别名表、边表规范化、配色、
  绘图设备、`normalize_formats` / `save_plot_formats`）
- `scripts/lib_layout_editor.R` —— 交互式布局编辑器（HTML）生成器
- `references/methodology.md` —— 三步取数的口径、score 标度陷阱、配色与布局选择、
  与参考代码的逐条差异
- `references/pitfalls.md` —— 参考代码 8 条注意事项的实测复核 + 已知坑 + **回归基线表**

> 📌 **改动脚本时，图内文字一律用 ASCII。** Windows 下中文字形缺失会被渲染成乱码
> ASCII（不是方框），很难一眼看出来。日志与 CSV 里的中文没问题（UTF-8）。
