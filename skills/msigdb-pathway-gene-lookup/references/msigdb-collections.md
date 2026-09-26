# MSigDB 集合说明与数据核对

> **该用哪一个集合？** 见 `references/collection-selection-guide.md`（按组织 + 实验类型的选择矩阵）
> 和 `SKILL.md` 第 0 节；本文只管「每个集合里到底有多少东西、数据对不对」。

数据位置：**`<SKILL_DIR>/data/`**，即随 skill 一起分发的目录（本 skill 自带，无需另下）。共 10 个 GMT 文件，
**MSigDB v2026.1.Hs，symbols 版**，约 29 MB，装完即用，不依赖任何外部目录。
所有数字都是在这批文件上实测的，并与 MSigDB 官网 Collections 页交叉核对过。

## 一、当前数据规模（已补齐，与官网总数一致）

| 文件 | 集合 | 基因集数 | 基因条目 | 说明 |
|---|---|---|---|---|
| `h.all.v2026.1.Hs.symbols.gmt` | H | 50 | 7,322 | Hallmark 标志性基因集 |
| `c1.all.v2026.1.Hs.symbols.gmt` | C1 | 302 | 43,707 | 位置型（染色体细胞遗传带） |
| `c2.all.v2026.1.Hs.symbols.gmt` | C2 | 7,670 | 588,281 | 专家审编（CP + CGP） |
| `c3.all.v2026.1.Hs.symbols.gmt` | C3 | 3,714 | 819,042 | 调控靶标（TFT + MIR） |
| `c4.all.v2026.1.Hs.symbols.gmt` | C4 | 1,006 | 98,548 | 计算型（3CA / CGN / CM） |
| `c5.all.v2026.1.Hs.symbols.gmt` | C5 | 16,283 | 1,373,024 | 本体（GO + HPO） |
| `c6.all.v2026.1.Hs.symbols.gmt` | C6 | 189 | 30,586 | 致癌签名 |
| `c7.all.v2026.1.Hs.symbols.gmt` | C7 | 5,219 | 990,387 | 免疫签名（ImmuneSigDB + VAX） |
| `c8.all.v2026.1.Hs.symbols.gmt` | C8 | 866 | 157,462 | 细胞类型签名 |
| `c9.all.v2026.1.Hs.symbols.gmt` | C9 | 62 | 5,575 | 计算扰动签名（DepMap CRISPR / CCLE） |

合计 **35,361 个基因集 / 43,633 个唯一基因 / 4,113,934 条基因条目**。

> **35,361 正好等于官网「Human MSigDB Collections」页给出的全部基因集总数**
> （该页原文：*The 35361 gene sets in the Human Molecular Signatures Database*）。
> 也就是说这 10 个 `*.all` 文件合起来已经等于官网的完整人类 gene set 集合，
> 无需再下载那个「All gene sets」合并包（官网本身也不建议直接用它做分析）。

set_id 全局唯一，无跨文件重复。

### 各集合的子库拆分（`--collection` 用的就是这些）

| 来源库 / 子库 | 集合数 | 别名 |
|---|---|---|
| GO:BP / GO:CC / GO:MF | 7,538 / 1,080 / 1,872 | `GOBP` `GOCC` `GOMF` `GO` |
| HPO | 5,793 | `HPO` `HP` |
| REACTOME | 1,839 | `REACTOME` |
| WikiPathways | 925 | `WP` `WIKIPATHWAYS` |
| KEGG（= LEGACY 186 + MEDICUS 658） | 844 | `KEGG` |
| BioCarta | 292 | `BIOCARTA` |
| PID | 196 | `PID` |
| Hallmark | 50 | `H` `HALLMARK` |
| C2 的 CGP 子库（化学/遗传扰动签名） | 3,555 | — |
| C2 的 CP 子库（经典通路） | 4,115 | — |

> `--collection GO` 横跨 C5 的 BP/CC/MF；`--collection PATHWAY` = 经典信号通路
> （H + C2 的 KEGG/REACTOME/WP/BIOCARTA/PID，约 3,900 个集），做通路分析时最常用。

## 二、关于经典 KEGG 通路名的核对结论

**结论：MSigDB 的 KEGG 子集从来没有包含过 `KEGG_PI3K_AKT_SIGNALING_PATHWAY`
或 `KEGG_TNF_SIGNALING_PATHWAY`。这不是 2026.1 版本的变动，也没有「老版本可以下回来」这回事。**

实测证据（每个文件的集合数 + 对这两条通路名的 grep 计数）：

| 版本 | 文件 | 集合数 | `KEGG_PI3K_AKT_SIGNALING_PATHWAY` | `KEGG_TNF_SIGNALING_PATHWAY` |
|---|---|---|---|---|
| v6.2 (2019) | `c2.cp.kegg.v6.2.symbols.gmt` | 186 | 0 | 0 |
| v7.4 | `c2.cp.kegg.v7.4.symbols.gmt` | 186 | 0 | 0 |
| v7.5.1 | `c2.cp.kegg.v7.5.1.symbols.gmt` | 186 | 0 | 0 |
| v2026.1 | `c2.all.v2026.1.Hs.symbols.gmt`（KEGG 部分） | 186 | 0 | 0 |

MSigDB 的 `KEGG_LEGACY` 子库（官网明示 "considered Legacy gene sets"）长期稳定在 **186 个集**，
是一个固定的、部分覆盖的旧 KEGG 通路快照。2026.1 新增的是 **658 个 `KEGG_MEDICUS_*`** 细粒度集
（186 + 658 = 844，与官网一致），后者按注释流向命名，不含那种「整条 XXX 信号通路」的集。

**所以检索 `PI3K-Akt` / `TNF` 的 KEGG 整条通路会落空——这是数据源本身的覆盖范围，不是脚本问题。**

### 替代方案（都在本库内，推荐按这个顺序取）

PI3K-Akt 方向：

| 集合 | 集合族 | 基因数 |
|---|---|---|
| `HALLMARK_PI3K_AKT_MTOR_SIGNALING` | H | 105 |
| `REACTOME_PI3K_AKT_SIGNALING_IN_CANCER` | C2/REACTOME | 113 |
| `PID_PI3KCI_AKT_PATHWAY` | C2/PID | 按需查看 |
| `REACTOME_PI3K_AKT_ACTIVATION` | C2/REACTOME | 9 |
| `BIOCARTA_AKT_PATHWAY` | C2/BIOCARTA | 按需查看 |

TNF 方向：

| 集合 | 集合族 | 基因数 |
|---|---|---|
| `REACTOME_TNF_SIGNALING` | C2/REACTOME | 57 |
| `HALLMARK_TNFA_SIGNALING_VIA_NFKB` | H | 200 |
| `KEGG_MEDICUS_REFERENCE_TNF_NFKB_SIGNALING_PATHWAY` | C2/KEGG | 16 |
| `KEGG_MEDICUS_REFERENCE_TNF_JNK_SIGNALING_PATHWAY` | C2/KEGG | 16 |

一条命令拿到全部候选：

```bash
... find "PI3K-Akt" "PI3K" "TNF" --collection H,KEGG,REACTOME,PID,BIOCARTA -n 20
... gene FOXO1 --collection PATHWAY        # 反向确认某个基因落在哪些主干通路
```

> 如果你**确实需要** MSigDB 口径的那条整通路（比如要复现别人论文里的
> `KEGG_PI3K_AKT_SIGNALING_PATHWAY`），那只能从别处取基因集（如 KEGG 官网 hsa04151，
> 或另建一份自定义 GMT），放进本目录即可——索引会自动纳入，同名冲突也不会报错
> （默认只返回主版本，`--keep-dups` + `--file` 可查非主版本）。

## 同类坑：`REACTOME_GLYCOLYSIS` 也不存在（v2026.1 已核实）

`show REACTOME_GLYCOLYSIS` 返回「没找到集合」，`find "glycolysis" --collection REACTOME`
**只命中 1 条**（`REACTOME_REGULATION_OF_GLYCOLYSIS_BY_FRUCTOSE_2_6_BISPHOSPHATE_METABOLISM`，12 基因）。
不是脚本或索引问题：Reactome 的糖酵解各步骤现挂在**父通路** `REACTOME_GLUCOSE_METABOLISM`（84 基因）之下，
没有单列一条 "Glycolysis"。查糖酵解走下面这套更直接（都在 H + C2 内）：

| 集合 | 集合族 | 基因数 |
|---|---|---|
| `HALLMARK_GLYCOLYSIS` | H | 200 |
| `KEGG_GLYCOLYSIS_GLUCONEOGENESIS` | C2/KEGG | 62 |
| `WP_GLYCOLYSIS_AND_GLUCONEOGENESIS` | C2/WP | 45 |
| `KEGG_MEDICUS_REFERENCE_GLYCOLYSIS` | C2/KEGG | 25 |
| `KEGG_MEDICUS_REFERENCE_GLUCONEOGENESIS` | C2/KEGG | 5 |
| `REACTOME_GLUCONEOGENESIS` | C2/REACTOME | 26 |
| `MOOTHA_GLYCOLYSIS` / `MOOTHA_GLUCONEOGENESIS` | C2（CGP） | 21 / 31 |

⚠️ **别拿 `REACTOME_GLUCOSE_METABOLISM`（84 基因）当糖酵解集**：它虽然含糖酵解各步骤
（HK1-3/GCK/GCKR/HKDC1 → GPI → PFKL/M/P → ALDOA/B/C → GAPDH/GAPDHS → PGK1/2 → PGAM1/2 →
ENO1-4 → PKLR/PKM → LDHA/B/C → PC/PCK1/2 → FBP1/2 → G6PC1-3 → SLC37A1-4），
但**同时混进了约 30 个核孔蛋白基因**（`NUP35/37/42/43/50/54/58/62/85/88/93/98/107/133/153/155/160/188/205/210/214`、
`NDC1`、`POM121`、`POM121C`、`RAE1`、`RANBP2`、`SEC13`、`SEH1L`、`TPR`、`AAAS`），
占全集的 1/3 以上，做富集会直接污染结果（疑为 MSigDB 侧 Reactome 集构建时的串接问题）。
糖酵解就用上表里的 KEGG / Hallmark / WP 条目。

一条命令拿到全部候选：

```bash
... find "glycolysis" "gluconeogenesis" --collection H,C2 --all -n 40
```

> CSV 里 `source_db` 为空的集不是缺数据：该列由 **set_id 前缀**推断
> （`HALLMARK_`/`KEGG_`/`REACTOME_`/`WP_`/`BIOCARTA_`/`PID_`，见脚本 `SOURCE_PREFIX`），
> CGP 型集（`MOOTHA_*`、`CUI_*` 等）不带这些前缀，所以留空——用 `--collection CGP` 可整体选中它们。

## 三、如何补齐/新增集合

下载地址规则（**已实测可用**）：

```
https://data.broadinstitute.org/gsea-msigdb/msigdb/release/<版本>.Hs/<文件名>
# 例：2026.1.Hs/c9.all.v2026.1.Hs.symbols.gmt
# 老版本无 .Hs 后缀：7.5.1/c2.cp.kegg.v7.5.1.symbols.gmt、6.2/c2.cp.kegg.v6.2.symbols.gmt
```

把下好的 `.gmt` 直接丢进数据目录（默认就是 `<SKILL_DIR>/data/`；若你用
`--data-dir` / `MSIGDB_DIR` / config.json 指向了别处，就放进那里），
下次查询会自动重建索引，不用手动 `build`。

⚠️ **curl 是 Windows 程序，`-o` 参数必须写 Windows 路径**（`C:/data/x.gmt`），
写成 Git Bash 的 `/c/data/x.gmt` 会报 `curl: (23) client returned ERROR on write`。

还想扩充时的常见选项：

| 想加什么 | 文件 | 说明 |
|---|---|---|
| Entrez 版（不是 symbol） | `<集合>.v2026.1.Hs.entrez.gmt` | 基因 ID 会变成数字，与现有 symbol 索引混在一起难以辨认，建议单独目录 + 独立索引 |
| 小鼠 | `<集合>.v2026.1.Mm.symbols.gmt` | 同上，建议单独目录 + 独立索引（用 `--data-dir` 或 `MSIGDB_INDEX`） |
| 官方 SQLite 版 | 官网 "SQLite database" | 含完整元数据（描述、来源、PMID），本脚本不读它 |
| KEGG 官网整通路 | KEGG hsa04151 等 | 自行转成 GMT 后放入 |

## 四、GMT 格式要点（排查问题时用）

每行：`set_id <TAB> url <TAB> gene1 <TAB> gene2 …`

- 第 2 列在本批数据里只是 MSigDB 官网 URL（无信息量），所以检索**只匹配 set_id**。
- `set_id` 统一大写、下划线分隔；基因 symbol 统一大写，字符集干净（`[A-Za-z0-9_.@-]`）。
- 行内无重复基因；set_id 无跨文件重复。
- 脚本对超大集合不设限（最大 `c3` 有 2,000 基因的集）；`--preview 0` 会全量打印。

## 五、相关工具

- `gmt_lookup.py list` —— 任何时候想知道当前索引里到底有什么，跑这个。
- 反向查某个基因（如 `FOXO1`）在经典通路中的位置：`gmt_lookup.py gene FOXO1 --collection PATHWAY`，
  按基因数排序，大集合在前，便于先看到主干通路。
- 导出的 `--gmt` 与原始 MSigDB 行逐字节一致，可直接 `clusterProfiler::read.gmt()` 或
  `GSEA` 的 `.gmt` 输入使用。
