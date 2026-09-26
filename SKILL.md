---
name: transcriptome-kit
description: 转录组 bulk RNA 全流程分析套件的总控入口：把 GEO 芯片标准化 / 差异分析 / 富集分析（ORA+GSEA+美化出图）/ MSigDB 通路检索 / KM 生存分析 / STRING 蛋白互作网络 / 韦恩图 六个子技能串成一条「数据 → 差异 → 功能 → 临床 → 网络 → 汇总图」的流水线。当用户提到 bulk RNA 分析、转录组分析流程、GSE 编号分析、差异基因之后做什么、一整套转录组分析、这套 skill 怎么用、哪个子技能负责哪一步、把 GEO 差异结果继续做富集/生存/网络时使用。具体每一步的执行细节在各子技能的 SKILL.md 里。
---

# Bulk RNA 转录组分析套件（总控）

本目录是 6 个子技能的**总控入口 + 流水线编排说明**。每个子技能的完整用法见
`skills/<名称>/SKILL.md`；本文件只回答三个问题：**整体流程怎么走、每步用谁、
步骤之间怎么交接**。

## 套件地图

| 子技能 | 在流水线里的角色 | 关键输入 → 产出 |
|---|---|---|
| `geo-microarray-analysis` | ① 数据标准化 → ② 差异分析（火山图/热图）→ ③ 多数据集合并去批次（可选） | GSE 原始文件 → 标准化矩阵、`_DEG.csv`、`_all.csv` |
| `enrichment-analysis` | ③ 功能解读：ORA + GSEA + 功能归纳 + 美化主图 | `_DEG.csv` / `_all.csv` → 富集表 + 出版级图 |
| `msigdb-pathway-gene-lookup` | 通路↔基因双向检索（离线，纯 Python） | 通路名/基因名 → 基因列表 / 参与通路 |
| `survival-km-analysis` | ④ 临床关联：KM 生存曲线（常规 → 最佳阈值两级策略） | 表达矩阵 + 临床表 → KM 图 + 决策文件 |
| `string-ppi-network` | ⑤ 蛋白互作网络（三格式输出 + 交互式布局编辑器） | 差异基因（可带 logFC） → PPI 网络图 |
| `venn-diagram` | ⑥ 汇总交集图（2~5 集合，组名必须用户确认） | 基因列表 → 韦恩图 + 各区交集表 |

## 标准流水线（按需取用，不必全跑）

```
GSE 原始文件（family.soft.gz + series_matrix.txt.gz）
   │   ★ 推荐让用户提前从 NCBI 下载好这两个文件再开工：AI 会话内直连
   │     NCBI 很慢且易超时；文件就位后整条流水线离线可跑（直链格式见 README）
   ▼
① geo-microarray-analysis 01_geo_normalize.R     标准化 + 分组候选报告
   │   ★ 停：把分组候选表交给用户，等用户拍板分组与对比方向
   ▼
② geo-microarray-analysis 02_deg_plots.R          limma 差异 + 火山图 + 热图
   │   产出 <prefix>_DEG.csv / <prefix>_all.csv / <prefix>_for_enrichment.txt
   │
   ├─► ③a enrichment-analysis（--geo_dir 一条命令接手）
   │        01 ORA → 02 GSEA → 04 归纳 ★停：用户选出图方式 → 03 美化主图
   │
   ├─► ③b msigdb-pathway-gene-lookup              查通路含哪些基因 / 基因在哪条通路
   │
   ├─► ④ survival-km-analysis                      每个候选基因出 KM 曲线
   │        （表达矩阵 + 临床时间/状态表）
   │
   ├─► ⑤ string-ppi-network --table=<prefix>_DEG.csv   差异蛋白互作网络
   │
   └─► ⑥ venn-diagram                              多来源基因列表取交集（如 2 个 GSE 的 DEG）
```

**两处必须停下来的地方**（子技能里的硬规则，总控同样执行）：

1. ①→② 之间：分组是生物学判断，**必须等用户指定**分组列、取值、对比方向。
2. 富集 04→03 之间：功能归纳报告交给用户后，**必须等用户选**出图方式（A 均衡配额 / B 某功能方向）。

## 步骤间的交接文件（记住这 3 个就够）

| 交接文件 | 谁产出 | 谁消费 | 怎么接 |
|---|---|---|---|
| `<prefix>_for_enrichment.txt` | geo 02 步 | enrichment 全部脚本 | `--geo_dir=<geo 产物目录>`，不用手抄参数 |
| `<prefix>_DEG.csv` | geo 02 步 | string（`--table=`）、venn（`--table= --by=change`）、survival（基因来源） | 各自的 `--table` 参数直接吃 |
| `<prefix>_group_template.csv` | geo 01/03 步 | 用户填好 → geo 02 步 `--group_file=`、survival 的临床表加工 | 两列表格（样本名,分组） |

## 环境要求（详见根目录 INSTALL.md / README.md）

- **首次使用前先跑 `python <套件根>/check_env.py`**：自动检测 R；**没装 R 时打印
  官方/清华镜像下载地址并询问是否自动下载（默认否，必须用户确认）**；装好后加
  `--check-packages` 可补查 R 包依赖并给出安装命令。
- **R ≥ 4.2**：geo / enrichment / survival / string / venn 五个子技能；包装器
  `run_*.sh` 会自动探测 Rscript（环境变量 `RSCRIPT` > PATH > Windows 常见位置）。
- **Python ≥ 3.8（纯标准库）**：msigdb-pathway-gene-lookup 唯一依赖。
- MSigDB GMT 数据已内置在 `skills/msigdb-pathway-gene-lookup/data/`（约 29 MB），
  enrichment 的 MSigDB 富集也会自动找到它（同一套件的 `data/` 目录）。

## 快速上手（最小可用链）

```bash
TK=<本套件根目录>

# ① 标准化（不联网）
bash $TK/skills/geo-microarray-analysis/scripts/run_geo.sh 01_geo_normalize.R \
  --gse=GSE62452 --matrix=GSE62452_series_matrix.txt.gz --soft=GSE62452_family.soft.gz \
  --dir=./out --prefix=GSE62452
# ★ 把 ./out/GSE62452_grouping_candidates.txt 交给用户，等分组指令

# ② 差异分析
bash $TK/skills/geo-microarray-analysis/scripts/run_geo.sh 02_deg_plots.R \
  --expr=./out/GSE62452.csv --group_file=./group.csv \
  --contrast=case-control --prefix=GSE62452_T

# ③ 富集（一条命令接手）
bash $TK/skills/enrichment-analysis/scripts/run_enrich.sh 01_enrich_ora.R \
  --geo_dir=./out --gmt_sets=H --direction=all,up,down
# ★ 04 归纳跑完把报告交给用户，等选 A/B 再出美化主图

# ⑤ PPI 网络
bash $TK/skills/string-ppi-network/scripts/run_string.sh 01_string_network.R \
  --table=./out/GSE62452_T_DEG.csv --engine=api --outdir=./string_out
```

## 设计原则（所有子技能共守）

1. **不替用户做生物学决定**：分组、对比方向、出图方式都要问，脚本有硬闸门的按闸门走。
2. **决策留痕**：每次运行输出决策文件 / 参数快照 / 可核对清单，结果可复现。
3. **口径一致**：上游 `change` 列拆上下调，下游不重新筛——全流程同一套阈值。
4. **图内文字一律 ASCII**：Windows 缺中文字形会把中文渲染成乱码 ASCII。
5. **回归基线**：各子技能 SKILL.md 里带实测基线数字，对不上说明中间出了错。
