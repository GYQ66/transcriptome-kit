# transcriptome-kit

**转录组 bulk RNA 从原始数据到可视化的全流程 Agent 技能套件**（Agent Skill Suite）。

把 6 个经过真实数据回归验证的子技能打包成一套，覆盖：
**GEO 芯片标准化 → 差异分析 → 富集分析（ORA + GSEA + 美化出图）→ MSigDB 通路检索 → KM 生存分析 → 蛋白互作网络 → 韦恩图**。

每个子技能既可以在总控流程里串联使用，也可以单独拿出来用。

## 子技能一览

| # | 子技能 | 一句话 | 语言/依赖 |
|---|---|---|---|
| 1 | [geo-microarray-analysis](skills/geo-microarray-analysis/SKILL.md) | GEO 芯片三步式：标准化 + limma 差异（火山图/热图，多风格）+ 多数据集 ComBat 合并 | R |
| 2 | [enrichment-analysis](skills/enrichment-analysis/SKILL.md) | ORA + GSEA + 功能大类归纳 + gground/ggprism 出版级富集主图 | R |
| 3 | [msigdb-pathway-gene-lookup](skills/msigdb-pathway-gene-lookup/SKILL.md) | MSigDB 十大集合通路↔基因双向检索（数据内置，离线可用） | Python（纯标准库） |
| 4 | [survival-km-analysis](skills/survival-km-analysis/SKILL.md) | KM 生存分析两级策略：中位值分组 → 最佳截断值自动回退 + maxstat 校正 p | R |
| 5 | [string-ppi-network](skills/string-ppi-network/SKILL.md) | STRING 蛋白互作网络：API/Bioconductor/本地三引擎 + 交互式布局编辑器 | R |
| 6 | [venn-diagram](skills/venn-diagram/SKILL.md) | 2~5 集合韦恩图：自动留白防裁切、组名确认闸门、各区交集导出 | R |

整体流程与步骤间交接文件见 [总控 SKILL.md](SKILL.md)。

## 安装

### 方式一：装进支持 Agent Skills 的客户端（WorkBuddy / Claude Code 等）

把本仓库（或解压后的 zip）里的 `skills/` 下各子目录，复制到客户端的用户级技能目录：

- WorkBuddy：`~/.workbuddy/skills/`（Windows: `C:\Users\<你>\.workbuddy\skills\`）
- Claude Code：`~/.claude/skills/`

```bash
git clone https://github.com/<你的用户名>/transcriptome-kit.git
mkdir -p ~/.workbuddy/skills
cp -r transcriptome-kit/skills/* ~/.workbuddy/skills/
```

之后在对话里说「帮我分析 GSE62452」「做富集分析」「画 KM 曲线」等即可自动触发。

### 方式二：作为命令行工具直接用

所有 R 子技能通过各自的 `run_*.sh` 包装器启动，**自动探测 Rscript 路径**并处理
Windows 中文用户名 / locale / 临时目录等环境坑：

```bash
# 依赖安装见 INSTALL.md；R 与包就绪后：
bash skills/venn-diagram/scripts/run_venn.sh \
  --set="处理组=gene_list_A.txt" --set="对照组=gene_list_B.txt" \
  --labels_confirmed=TRUE --outdir=./out

# MSigDB 检索（纯 Python，零依赖）
python skills/msigdb-pathway-gene-lookup/scripts/gmt_lookup.py find "TNF signaling"
```

## 准备 GEO 数据（推荐：先下载好，再丢给 AI）

geo 技能每分析一个 GSE 编号，需要两个原始文件：

| 文件 | 内容 | 用途 |
|---|---|---|
| `GSE{编号}_series_matrix.txt.gz` | 表达矩阵 + 样本信息 | 标准化、差异分析 |
| `GSE{编号}_family.soft.gz` | 样本临床元信息 | 分组、注释 |

**强烈推荐：先自己把这两个文件下载到本地，再把本地路径丢给 AI 分析。**
AI 的运行环境直连 NCBI 的速度往往很慢也不稳定，让它在会话里现下载，
一个数据集可能要等很久甚至超时失败；提前下载好之后，整条流水线
全程离线可跑，又快又稳。

直链格式（以 GSE62452 为例；`GSE62nnn` 的规则：编号去掉末 3 位，
剩余数字原样保留再补 `nnn`）：

```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62452/soft/GSE62452_family.soft.gz
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62452/matrix/GSE62452_series_matrix.txt.gz
```

不想算规律也可以直接打开
[GEO 编号页](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452)
（把 acc= 换成你的编号），点页面下方的 **Download family** 下载
`_family.soft.gz`，在 **Download data files** 里选 **Series Matrix File(s)**
下载 `_series_matrix.txt.gz`。
下载完成后对 AI 说：「用这两个本地文件做 GSE62452 差异分析」即可。

## 环境要求

| 组件 | 版本 | 用于 |
|---|---|---|
| R | ≥ 4.2 | 5 个 R 子技能 |
| Python | ≥ 3.8（纯标准库） | msigdb-pathway-gene-lookup |
| R 包 | CRAN + Bioconductor 若干 | 逐技能清单见 [INSTALL.md](INSTALL.md) |
| 网络 | 仅部分步骤需要 | KEGG 在线富集、STRING API；均可离线回退 |

## 第一步永远是：环境检查

安装（或克隆）后先运行环境检查脚本（纯 Python 标准库，无需任何依赖）：

```bash
python check_env.py
```

它会依次检查 **R → GMT 数据 → （可选）R 包依赖**。**如果没检测到 R，它会打印
官方与清华镜像的下载地址，并询问是否自动下载最新安装包**——不回答就默认不下载，
绝不未经确认就动你的系统：

```
[1/3] 检查 R …
  ✗ 未检测到 R（本套件 6 个技能中有 5 个依赖 R）

R 下载地址（cran）：https://cran.r-project.org/bin/windows/base/
国内镜像（清华）    ：https://mirrors.tuna.tsinghua.edu.cn/CRAN/bin/windows/base/
是否现在自动下载 R 安装包？[y/N]
```

- 回答 `y`：自动解析最新版本直链（91 MB 左右），带进度条下载到「下载」文件夹，
  下载完可再确认是否启动安装（Windows 静默安装 / macOS 打开向导）。
- 回答 `n` 或直接回车：只保留地址，你自己装。
- 装完 R 后重跑 `python check_env.py --check-packages` 还能把缺的 R 包列出来，
  并给出可直接复制的安装命令（如 `install.packages(c('svglite'))`）。
- `--mirror tuna` 强制走清华镜像；`--no-download` 只看地址不询问。

MSigDB GMT 数据（v2026.1.Hs，10 个集合 / 35,361 个基因集）**已内置**，
不需要另下数据。

## 特性

- **流水线交接**：geo 差异分析跑完生成 `_for_enrichment.txt` 交接清单，
  富集侧一条 `--geo_dir` 接手，不用手抄参数；`_DEG.csv` 可直接喂 PPI / 韦恩图。
- **决策留痕**：分组候选报告、KM 决策文件、布局坐标快照、运行参数快照——
  每一步「为什么这么做」都可回查。
- **硬闸门防越权**：分组要用户拍板、韦恩图组名要用户确认、富集出图方式要用户选，
  脚本层面强制停顿，不会悄悄替你做生物学决定。
- **跨平台包装器**：`run_*.sh` 统一处理 Windows 中文用户名、`LC_ALL` locale 陷阱、
  POSIX→Windows 路径转换、Rscript 探测；macOS / Linux 开箱即用。
- **回归基线**：各子技能附实测基线数字（如 GSE62452 上调 177 / 下调 115），
  升级或换机后可快速自检。

## 适用范围

- **适用**：bulk 转录组（GEO 芯片、TCGA 表达量）、基因列表的功能与临床解读。
- **不适用**：scRNA-seq、RNA-seq counts 的 DESeq2/edgeR 流程（geo 技能只吃芯片/表达量矩阵）。

## License

MIT —— 详见 [LICENSE](LICENSE)。
