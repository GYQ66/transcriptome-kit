---
name: msigdb-pathway-gene-lookup
description: 在本地 MSigDB GMT 文件（Hallmark / KEGG / REACTOME / WikiPathways / GO / 免疫与细胞类型签名等十大集合）中做「通路 → 全部基因」正向检索和「基因 → 参与通路」反向检索，并在检索前先按用户的组织/样本与实验类型推荐该用哪几个集合（H + C1…C9）再由用户拍板。当用户问「某个信号通路包含哪些基因」「XXX 通路里的基因列表」「这个基因参与了哪些通路」「我这个实验该查哪个数据库/哪方面的通路」「把这几条通路导出成 gmt / csv 给 GSEA 用」时使用。关键词：信号通路、通路基因、基因集、GMT、MSigDB、GSEA、Hallmark、KEGG、REACTOME、免疫签名、细胞类型签名、集合选择、该用哪个数据库。
metadata:
  short-description: 本地 MSigDB 十大集合的通路↔基因双向检索 + 集合选择向导
  agent_created: true
---

# MSigDB 通路 ↔ 基因 双向快速检索

`<SKILL_DIR>` 指本 skill 的安装目录（agent 加载时会注入实际路径；手动运行脚本时，把它替换成
`scripts/gmt_lookup.py` 所在目录的上一级即可）。

用这一个脚本回答「某条信号通路里有哪些基因」和「某个基因参与了哪些通路」。

**数据已内置，不依赖任何外部文件**：`<SKILL_DIR>/data/` 里带了 MSigDB v2026.1.Hs 全部 10 个 `*.all` 集合的
symbols 版 GMT（H + C1…C9），共 **35,361 个基因集 / 43,633 个唯一基因 / 4,113,934 条基因条目**，约 29 MB。
底层把这些 GMT 建一次 SQLite 索引（`sets` 表存通路→基因，`gene_index` 表存基因→通路倒排），
之后每次查询都是索引命中，检索本身 <100 ms。

## 零、先定范围（必做）：问清「什么组织 + 什么实验」，再选集合

MSigDB 的十大集合（H + C1…C9）**不是同一类东西**：真信号通路只有 H、C2、C6，
其余是染色体位置（C1）、调控靶标（C3）、计算模块（C4）、本体与表型（C5）、
免疫签名（C7）、细胞类型标记（C8）、扰动签名（C9）。
**集合选错，查出来的东西对用户毫无意义**——所以本 skill 的第一动作不是检索，而是锁定范围。

### 第 1 步：问两件事（用户没说就必须问）

1. **什么组织/样本**：物种 + 组织或细胞（如「人肝癌组织」「小鼠胰岛 β 细胞」「患者外周血 PBMC」）；
2. **什么实验**：技术平台 + 实验设计（如「bulk RNA-seq，对照组 vs 敲低组」「scRNA-seq 做细胞类型注释」
   「ChIP-seq 找 XX 的靶基因」「CRISPR 筛选」。

用 AskUserQuestion **一次把两个问题都问掉**（每问 3~4 个常见选项，允许用户自填），别只问一半就开始查。

> 用户在本轮已经说清了（如「我做人肝癌 bulk RNA-seq，想看该查哪些通路」）→ **不要重复问**，
> 用一句话复述理解（「人肝癌组织、bulk RNA-seq 的差异基因，对吧？」）直接进第 2 步。

### 第 2 步：用向导拿建议（别凭感觉推）

```bash
python "<SKILL_DIR>/scripts/gmt_lookup.py" guide "人肝癌 bulk rna-seq 差异表达"
```

`guide` 按「组织/实验体系」匹配出**主推 + 可叠加**的集合，并为每个集合给出：

- `是信号通路吗` / `内容` / `适合组织` / `适合实验`
- **`适合的场景`：这个库具体适合哪些研究场景（每个集合 4~6 条，这是要讲给用户看的重点）**
- `典型产出`：用了它能得到什么（富集图 / 打分矩阵 / 细胞类型标签 / gmt…）
- `什么时候选它` / `不适合 / 坑`

不带关键词时打印全部 10 个集合；带关键词但想同时看全部 10 个加 `--all`。
完整矩阵见 `references/collection-selection-guide.md`。

### 第 3 步：把建议 + 每个库的场景介绍摆给用户，等他拍板

给用户的推荐里，**每个集合都要带上「它适合什么场景」**（从 `guide` 的 `适合的场景` 清单里挑
2~3 条最贴近用户体系的，别只写一句「推荐这个」），形式建议：

> 你这个体系（人肝癌组织、bulk RNA-seq 差异表达）我建议看这三个：
>
> | 集合 | 里面是什么 | 适合什么场景 | 用了能得到 |
> |---|---|---|---|
> | **H** Hallmark（50 集） | 50 条标志性生物学过程 | 差异基因第一层富集、队列分组比较、跨数据集可比 | 条数少、好解释的富集图 |
> | **C2** 经典通路 | KEGG / REACTOME / WP / PID 具体通路 | 要机制级细节、画通路图、多库交叉验证 | 具体通路的富集结果 |
> | **C6** 致癌签名（可选） | 癌基因激活 / 抑癌失活签名 | 判断哪条致癌通路被激活、肿瘤亚型比较 | 致癌通路激活打分 |
>
> 你想看哪一个或哪几个？定下来我就按这个范围去查通路和基因。

然后**停下来等确认**（AskUserQuestion，multiSelect=true，问「看哪几个集合的通路」）。
在用户点头之前，不要跑 `find` / `gene`。

### 第 4 步：只在用户选定的范围内检索

之后每条 `find` / `show` / `gene` 都带上用户选定的 `--collection`，
并在回复里写明本次范围（例：「本次范围：H + C2 的 REACTOME/KEGG」）。

用户如果完全放权（「你决定」「都看看」「随便」），按这个顺序兜底：

| 兜底选择 | 说明 |
|---|---|
| `--collection PATHWAY` | 经典信号通路（H + KEGG/REACTOME/WP/BIOCARTA/PID，约 3,900 集）——通用默认 |
| `--collection PATHWAY,GOBP` | 再叠 C5 的 GO:BP，需要穷尽式功能注释时 |
| `--collection H,C2,C5,C7,C8` | 真要「都看看」时的五集合折中组合 |

**不要一次把 10 个集合全选**：C5 一个就有 16,283 个集，全选会把检索结果彻底淹没。

## 一、怎么调用

```bash
python "<SKILL_DIR>/scripts/gmt_lookup.py" <子命令> [参数]
```

只需要 **Python 3.8+，纯标准库**，不用装任何第三方包。直接调用系统里的 `python` / `python3` 即可：

```bash
python "<SKILL_DIR>/scripts/gmt_lookup.py" ...
```

若 `python` 不在 PATH，用你机器上任一 Python 3.8+ 解释器的**绝对路径**替换上面的 `python` 即可。

**首次运行**会自动建索引（约 10~16 秒，只需一次），之后每次查询约 2 秒（其中 1.7 秒是解释器启动）。
索引落在 `<SKILL_DIR>/data/.msigdb_index.sqlite3`（约 78 MB，纯生成物，删掉会自动重建）；
若安装目录不可写，会自动改放到用户缓存目录（Windows: `%LOCALAPPDATA%\msigdb-gmt-lookup\cache`；
Linux/macOS: `$XDG_CACHE_HOME/msigdb-gmt-lookup`，默认 `~/.cache/msigdb-gmt-lookup`）。
**GMT 文件有任何增删改，下次查询会自动重建索引**，无需手动 `build`。

想改用自己的 GMT 目录（而不是自带的），优先级从高到低：
`--data-dir <目录>` > 环境变量 `MSIGDB_DIR` > `<SKILL_DIR>/config.json` 里的 `data_dir`。
（config.json 默认留空 = 用自带数据。）

## 二、五个子命令

| 子命令 | 用途 | 典型写法 |
|---|---|---|
| `guide` | **选择向导**：按组织/实验类型推荐该用哪几个集合 | `guide "人肝癌 bulk rna-seq"` |
| `list` | 列出集合、规模、来源库分布 | `list` |
| `find` | **正向**：通路名/关键词 → 基因 | `find "TNF signaling"` |
| `show` | **正向**：完整集合名 → 全部基因 | `show HALLMARK_APOPTOSIS` |
| `gene` | **反向**：基因 → 参与的通路 | `gene TP53` |

```bash
# 第 0 步的选择向导（不带关键词 = 打印全部 10 个集合的适用面）
... guide
... guide "人肝癌 bulk rna-seq 差异表达"
... guide "scRNA-seq 细胞类型注释"

# 正向：一个或多个关键词，一次进程解决（推荐把多个查询合并成一次调用，见「坑」）
... find "TNF signaling" "WNT signaling" apoptosis -n 5

# 限定范围：只在用户选定的集合里查（第 0 步拍板的结果）
... find "PI3K" --collection H,C2
... find apoptosis --collection H            # 只要 Hallmark
... find "cell cycle" --collection C2,KEGG,REACTOME
... find "TNF" --collection PATHWAY          # H + 经典通路（通用默认）
... find "IFN" --collection IMMUNE           # 语义别名，= C7 免疫签名

# 只要名称精确匹配 / 显示全部命中 / 一次打印全部基因
... find MAPK --exact
... find MAPK --all -n 100
... find "TNF signaling" --preview 0          # 0 = 每个基因集打印全部基因

# 精确集合名取全部基因（等价于原来 R 里的 subset(gmt, term=='...')）
... show HALLMARK_ADIPOGENESIS
... show KEGG_MAPK_SIGNALING_PATHWAY --format genes   # 一行一个基因，便于管道

# 反向
... gene TP53                                # 全库
... gene TP53 --collection PATHWAY -n 10     # 只看经典通路
... gene FOXO1 --collection H                # 这个基因在哪些 Hallmark 里
... gene "TP53,EGFR" --mode intersect        # 同时包含两个基因的通路（交集）
... gene "TP53,EGFR" --mode union            # 至少包含一个（并集，默认）
```

### 导出

```bash
... find apoptosis --collection H --csv out.csv     # 长格式：set_id/集合/来源/基因数/匹配/基因
... show HALLMARK_APOPTOSIS --csv genes.csv        # set_id/集合/来源/基因数/基因
... gene TP53 --collection PATHWAY --csv hits.csv  # .../matched_genes/genes
... find "TNF signaling" --gmt tnf_sets.gmt        # 标准 GMT，可直接 clusterProfiler::read.gmt 读入
```

CSV 一律 UTF-8 BOM，Excel 直接双击不乱码。`--gmt` 导出与原始 MSigDB 行**逐字节一致**（已验证）。

### 其他开关

`--collection/-c`（可重复或用逗号）、`--file`（按 GMT 文件名子串过滤）、`--keep-dups`、
`--limit/-n`、`--all`、`--exact`、`--preview N`、`--format table|json|genes`、`--data-dir`、`--db`、`--no-build`。

`--collection` 可用别名：`H`/`HALLMARK`、`C1`…`C9`、`KEGG`、`REACTOME`、`WP`/`WIKIPATHWAYS`、
`BIOCARTA`、`PID`、`GO`/`GOBP`/`GOCC`/`GOMF`、`HP`/`HPO`、`MIR`、`TFT`、`PATHWAY`（= 经典信号通路）。

按「研究方面」选集合时可以用语义别名（与 C1…C9 完全等价，第 0 步推荐里给的就是这些）：

| 语义别名 | = 集合 | 适用场景 |
|---|---|---|
| `IMMUNE` / `IMMUNESIGDB` | C7 | 免疫、炎症、感染、肿瘤免疫微环境 |
| `CELLTYPE` / `CELLTYPES` / `MARKER` | C8 | 单细胞注释、细胞类型标记、组成解卷积 |
| `PERTURB` / `PERTURBATION` / `DEPMAP` | C9 | CRISPR / 敲低 / 药物处理的扰动响应 |
| `ONCOGENIC` | C6 | 致癌通路激活状态 |
| `CANCER` | C6 + C9 | 肿瘤方向的两类签名一起看 |
| `REGULATORY` / `TARGETS` | C3 | TF / miRNA 靶标（含 `TFT`、`MIR`） |
| `POSITIONAL` / `LOCATION` | C1 | 染色体位置 / CNA / 核型 |
| `ONTOLOGY` | C5 | GO + HPO 本体（含 `GO`、`HPO`） |
| `COMPUTATIONAL` / `MODULE` | C4 | 共表达模块 |
| `CP` / `CLASSIC` | C2 的经典通路 | 只要 KEGG/REACTOME/WP/BioCarta/PID |
| `CGP` | C2 的扰动签名 | C2 里 3,574 个非经典通路来源的集 |

## 三、结果怎么给用户

用户要的是**能直接用的基因列表**，所以：

1. **开头先写明本次范围**：用的哪几个集合（「本次范围：H + C2 的 REACTOME/KEGG」），
   以及为什么是这几个；用户后面要复现或换范围时才有的可依。
2. 再给一张概览表（通路名 / 集合 / 来源库 / 基因数），让用户知道命中了哪些；
3. 用户关心的是哪条通路，就把该通路**完整基因列表**贴出来（`--preview 0`）；
4. 同时把 `--csv` 落到用户指定的输出目录，并在回复里写明文件路径。
5. 若用户要做 GSEA/富集，用 `--gmt` 导出子集 GMT 比 CSV 更省事。

默认按「精确 > 前缀 > 词组 > 全词 > 包含 > 近似」分层打分，并给经典通路来源（KEGG/REACTOME/WP/Hallmark 等）加权。
同名多条时默认只返回主版本；用 `--keep-dups` 可看到全部版本，再用 `--file` 定位。

用户在第 0 步选定的集合里查不到东西时（例：在 H/C2 里找肝细胞特异标记），
**不要默默扩到全部 10 个集合**——先说明「这个方向不在你选的集合里」，
再建议该换哪个集合（如 C8 细胞类型签名），由用户决定。

## 四、坑（重要）

- **传给脚本/curl 的路径参数必须是 Windows 形式**：`--csv "D:/out/out.csv"`。
  不要写 Git Bash 的 `/d/out/out.csv` —— 那是给 bash 自带命令用的；
  Python 会把它当成「当前盘符下的 `\d\out\...`」而**静默**写错位置（不报错），
  curl 则会直接报 `curl: (23) client returned ERROR on write`。
- 要新增集合时，官方下载地址规则是
  `https://data.broadinstitute.org/gsea-msigdb/msigdb/release/<版本>.Hs/<文件名>`，
  下好丢进 `data_dir` 即可（详见 `references/msigdb-collections.md`）。
- **部分 Windows 机器 Python 启动可慢至 1~2 s**（与脚本无关）。所以
  **尽量把多个关键词写在一次 `find` 调用里**（`find "A" "B" "C"`），不要为每个词起一次进程。
- **MSigDB 的 KEGG 子集本来就不含 `KEGG_PI3K_AKT_SIGNALING_PATHWAY`、`KEGG_TNF_SIGNALING_PATHWAY`。**
  已跨版本核实 v6.2 / v7.4 / v7.5.1 / v2026.1 的 KEGG 子集都是稳定 186 个集（KEGG_LEGACY），
  从来就没有这两条 —— 不是 2026.1 的变动，也没有老版本可下回来。检索落空是**数据源的覆盖范围**，
  替代集合与完整证据见 `references/msigdb-collections.md`。
- **基因 symbol 必须用 MSigDB 官方写法**（`TP53`、`NFKB1`、`MTOR`），蛋白名/别名查不到。
  反向检索不到时会明确告知哪些 symbol 没找到。
- 连字符/空格会被归一化：`NF-kB`、`NFkB`、`NF KB` 都能命中，`PI3K-Akt` 也能命中
  `..._PI3K_AKT_SIGNALING_PATHWAY`，不需要手工改成下划线。
- 只搜基因集名，**不搜描述文字**（GMT 第 2 列只是 MSigDB 官网 URL，没有信息量）。
- 自带数据只有 `Hs`（人类）symbols 版。要加小鼠 / Entrez 版，建议放到**独立目录**并配独立索引
  （`Entrez` 是数字 ID，和 symbol 混在一个索引里难以辨认）。
- `data/.msigdb_index.sqlite3`（78 MB）是自动生成的，可以安全删除，下次查询会重建。
  要分享/打包这个 skill 时建议把它排除，只带 `data/*.gmt`（29 MB）。

## 五、参考

- `references/collection-selection-guide.md`：**十大集合怎么选**（第 0 步的完整矩阵）——
  每个集合适合什么组织 / 什么实验 / **适合什么场景（逐个展开）**、11 类常见场景的推荐组合、
  问用户的话术模板、常见选错的例子。
- `references/msigdb-collections.md`：10 个集合（H + C1…C9）的含义、各来源库规模、
  「KEGG 子集从不含 PI3K-Akt / TNF」的跨版本核对证据与替代集合、以及如何下载补齐新集合。
