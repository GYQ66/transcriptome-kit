# 多数据集合并与批次矫正（03 步）

`scripts/03_merge_batches.R` 的配套参考。本文记录设计依据、`ComBat` 的源码级行为、
以及必须由用户拍板的决策点。

**本文所有结论均在 R 4.5.2 / sva 3.58.0 / limma 3.66.0 / Windows 上实机验证过**
（2026-09-17），不是推测。

---

## 1. 定位

```
01_geo_normalize.R  × N   --  每个数据集各自标准化
        |
        v
03_merge_batches.R        --  取共同基因 -> cbind -> 批次矫正 -> 批次 QC
        |
        v
02_deg_plots.R            --  在合并矩阵上做差异分析
```

参考实现：社区常见的"多数据集合并 + ComBat"写法。
`--pre=auto` 的判定规则**刻意复刻**该实现，使本脚本在其上与其**数值一致**
（见第 9 节验证记录）。差异只在于：多了 QC 图表、多了基因命名/NA/混杂的显式检查、
读取改用 `data.table::fread`。

### 1.1 从原始 GEO 文件开始：`03b_merge_from_geo.R`

上面那套流程要求输入是**已标准化**的矩阵。但用户手上通常只有
N 组 `(GSE*_series_matrix.txt.gz, GSE*_family.soft.gz)`。手工做法要跑 N 遍 01
再跑 03，而 N 份临床表还得人工拼成一张分组表（用户历史脚本 `merge处理.R` 就是这么干的，
且极易出错）。`03b_merge_from_geo.R` 把这条链自动化。

**它分两段跑，因为"分组"必须由用户拍板：**

```bash
# 第一段: 逐个 GSE 标准化 + 打印每个 GSE 的分组候选，然后停住
03b_merge_from_geo.R --gse=A,B --datadir=D:/data --outdir=./geo_merge --stage=normalize
# 第二段: 合并 + 批次矫正（分组由用户指定）
03b_merge_from_geo.R --gse=A,B --datadir=D:/data --outdir=./geo_merge --stage=merge \
  --group_cols='列1;列2' --group_maps='取值=case,...;取值=control,...'
# 第三段（可选）: 用户改了差异分析用的分组 -> 只重建分组表，不重算 ComBat
03b_merge_from_geo.R --gse=A,B --outdir=./geo_merge --stage=group \
  --group_cols='列1;列2' --group_maps='...'
```

`--stage=all` 可以一次跑完（分组方案已经明确时）。`--stage=merge` 可反复重跑调参，
不必重做标准化（已完成的 GSE 默认跳过，`--force=TRUE` 才重跑）。

**几个设计要点：**

- **`--datadir` 递归查找** `<GSE>_series_matrix.txt.gz` / `<GSE>_family.soft.gz`，
  也支持 `--matrices=` / `--softs=` 直接指定（与 `--gse` 一一对应，允许留空项）。
  默认**不允许联网**（`--allow_download=FALSE`），找不到文件就把期望的文件名与
  GEO 的 FTP 位置打出来。
- **`--pre` 默认 `none`**：各 GEO 数据集已经被 01 标准化（log2 + 探针注释 + 去重 +
  `normalizeBetweenArrays`），03 侧只该做"取交集 + cbind + 批次矫正"。
  若走 `--extra_inputs=` 并入未标准化的外部矩阵（如 TCGA TPM），那些会默认给 `auto`。
- **分组推导**：从每个 GSE 的 `clinical_<GSE>.csv` 里取指定列，可选先用
  `--group_regexes` 提取，再用 `--group_maps` 把取值映射成统一分组名；
  没被映射到的取值留 NA，交给 03 的 `--group_na=drop` 剔除
  （等价于用户手工脚本里"只保留 OA/Normal"那一步）。
  推导出的统一分组表写到 `<outdir>/<prefix>_group_derived.csv`，可直接喂给 02 步。
- **临床信息自动汇总**：编排器会把各 GSE 的 `per_gse/clinical_<GSE>.csv`
  按 `--inputs` 的顺序传成 03 的 `--clinical_files=`，于是 03 会额外产出
  合并临床表 `<prefix>_clinical.csv` 与分组候选报告
  `<prefix>_grouping_candidates.txt`——**这两份是向用户索要差异分析分组的材料**
  （见第 6.2 节）。也可用 `--clinical_files=` 覆盖自动定位。
- **`--stage=group`**：只重建分组表。样本范围取自**合并临床表**（已被
  `--group_na=drop` 剔掉的样本不会回来），不跑 03、不重算 ComBat。
  注意它**不改变**合并矩阵里 ComBat 用的 `mod`。
- **判成败看产物不看退出码**：本机 R 退出时固定报 Segmentation fault(139)，
  与工作是否完成无关（pitfalls 1.4）。编排器按 `<GSE>.csv` / `<prefix>_merged.csv`
  是否存在判定，失败时打印该子进程日志的尾部。
- **locale 与路径编码**：编排器在最开头自愈 locale（见下），
  并把 `LC_ALL` 等从环境里清掉，保证派生出的 01/03 子进程一启动就是
  `Chinese (Simplified)_China.utf8`。**这一步是中文路径能用的前提**，见 pitfalls 1.7。

---

## 2. 为什么这一步也必须让用户拍板

第一步（标准化）是纯数据加工，脚本可以自动完成；第二步（差异分析）的分组是生物学判断。
第三步介于两者之间，**但它的判断会直接改变矩阵里的数值**，所以同样不能默认：

| 决策 | 为什么不能替用户决定 | 怎么表达 |
|---|---|---|
| 批次怎么划分 | 一个数据集不一定等于一个批次；某数据集也可能要留着当验证集不参与矫正 | `--batch=` / `--batch_file=` |
| 用哪种矫正 | ComBat 与 quantile 的结果、适用场景都不同 | `--method=combat|quantile|both|none` |
| 是否保护生物学变异 | 不给分组时 ComBat 可能把真实的组间差异一起扣掉 | `--group_file=` + `--combat_mod=group` |

**跑完 03 步要把 `<prefix>_batch_map.csv`、`<prefix>_overlap_report.txt` 和
PCA 前后图交给用户确认，再进 02 步。** 另外，若要做差异分析，还要把
`<prefix>_clinical.csv` + `<prefix>_grouping_candidates.txt` 交给用户**索要分组**
（见 6.2 节），拿到答复前不要跑 02 步。

---

## 3. `--pre=auto` 的判定规则

`auto` 看**第一个样本列名的前 3 个字符**，然后结合对数化启发式：

| 首列名前 3 字符 | 判定为 | 处理 |
|---|---|---|
| `TCG` | TCGA 表达量（counts / TPM / FPKM） | 需要时 `log2(x+1)`，**不做** `normalizeBetweenArrays` |
| `GSM` / 其他 | GEO 芯片 | 需要时 `log2(x)`，**然后** `normalizeBetweenArrays` |

对数化启发式（与 01 步同一套，作用在全部有限值上）：

```r
qx <- quantile(x, c(0, .25, .5, .75, .99, 1))
needs_log2 <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) ||
              (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
```

**为什么 TCGA 用 `log2(x+1)` 而不是 `log2(x)`**：TPM/FPKM 可以是 0，
`log2(0) = -Inf` 会污染整个矩阵。`log2(x+1)` 把 0 映射到 0。

**为什么 TCGA 不再做 `normalizeBetweenArrays`**：同一个 TCGA 队列内部各样本
已经是统一流程产出的，样本间可比；而它与 GEO 数据集之间的尺度差异交给后面的
批次矫正处理。这也是 7.R 的选择（`rt1` 只取对数，`rt2`/`rt3` 才做 quantile）。

每个数据集在处理前后都会打印分位数：

```
  TCGA_LUSC_TPM: auto[TCG]: 判定未取对数 -> log2(x+1)
      对数化前 min=0.00 Q1=0.31 中位=6.13 Q3=25.31 Q99=491.32 max=108151.41
      处理后   min=0.00 Q1=0.38 中位=2.83 Q3=4.72 Q99=8.94 max=16.72
```

**复核方法**：看 `处理后` 那一行。已对数化的芯片数据中位数通常在 0~15。
若 `auto` 判错（双色芯片容易误判），用分号逐个覆盖：
`--pre='log2p1;quantile;quantile'`。

---

## 4. 基因匹配

### 4.1 交集顺序

```r
common <- rownames(mats[[1]])
for (i in 2:n) common <- intersect(common, rownames(mats[[i]]))
```

结果**保持第 1 个数据集的基因顺序**，与 7.R 的
`Reduce(intersect, list(...))` 一致——这是数值可比的前提之一。

### 4.2 大小写（跨平台混用的关键）

`--gene_case=asis`（默认）原样匹配；`upper` 统一转大写再匹配。
**默认取 `asis` 是为了忠于 7.R**，但跨平台（RNA-seq 与芯片混用）经常需要 `upper`。

脚本默认会**额外扫描**并把收益写进报告（不改数据）：

```
基因名大小写敏感性:
  当前大小写交集 = 15913；统一为大写后 = 15913 (+0)
```

⚠️ **收益可能是 0，但 `upper` 仍然有实际作用。** 在这三个数据集上，统一大写
既**新增**了匹配（如 `C10orf53` 与 `C10ORF53` 对齐），也**合并/撞车**掉了一些
原本不同的名字，两者相抵，净收益恰好为 0。所以：

- 报告显示收益 **明显为正** → 加 `--gene_case=upper` 能保留更多基因
- 报告显示收益 **为 0** → 两种写法净基因数相同，此时按平台惯例选即可

### 4.3 并集模式

`--match=union` 用全部基因的并集，某数据集没测到的位置补 NA。
**必须配合 `--na=rowmean`**（或 `zero`），否则 ComBat 会因为 NA 拒绝运行。

---

## 5. 批次矫正方法

| `--method` | 实现 | 适用 |
|---|---|---|
| `combat`（默认） | `sva::ComBat` | 跨平台、跨队列合并；能同时校正均值与方差 |
| `quantile` | `limma::normalizeBetweenArrays` | 各数据集分布差异不大；不做参数估计，最保守 |
| `both` | 两份都出 | **跨平台建议先用这个**，让用户对比 PCA 再定 |
| `none` | 不矫正 | 只想要合并矩阵时 |

`--method=both` 产物：主产物 `<prefix>_merged.csv`（由 `--primary=` 决定，
默认 `combat`），另一份为 `<prefix>_merged_combat.csv` / `_merged_quantile.csv`。

**ComBat 在做什么**：对每个基因，先按 `mod` 拟合生物学协变量并取残差 →
标准化 → 用 L/S 模型估计每个批次的加性偏移 `gamma` 与乘性缩放 `delta` →
用**经验贝叶斯**把各批次的估计向共同先验收缩（`par.prior=TRUE`，即参数化先验）
→ 反变换回原尺度。收缩的作用是防止样本量小的批次过拟合。

### 5.1 ComBat vs quantile：同一份数据上的实测差异 ★

在三数据集合并（15913×967，TCGA + 2 个 GEO）上分别矫正后，
用"批次对 PC1 的解释方差"衡量残留批次效应：

| `--method` | 批次对 PC1 解释方差 | 三批次均值（全矩阵） |
|---|---|---|
| `none`（不矫正） | 99.3% | 3.289 / 6.673 / 5.909 |
| `quantile` | **98.6%**（几乎没降） | **4.654 / 4.653 / 4.653**（已拉齐） |
| `combat` | **0.0%** | 4.650 / 4.658 / 4.657（已拉齐） |

**注意 `quantile` 那一行的反常之处**：各批次的均值/中位数/sd 已经被拉得**完全一致**
（因为分位数归一化让每个样本的分布都相同），但 PC1 上的批次信号仍然高达 98.6%。
这正是问题的关键——

**这是"跨平台必须优先 ComBat"的直接证据。** 原因：

- `normalizeBetweenArrays` 分位数归一化把**每个样本的分布**强制变成同一条，
  所以样本级、批次级的均值/分位数会被拉齐；
- 但它**不建模、也不扣除基因特异的批次效应**。同一个基因在 TCGA 样本里排在
  分布的高位、在 GEO 样本里排在中位，两边的"值"仍然系统性不同——
  PC1 恰好就是抓这种基因水平的系统性差异，所以 98.6% 居高不下。
- `ComBat` 是**逐基因**估计每个批次的加性/乘性参数（`gamma`/`delta`）并扣除，
  所以能把 PC1 上的批次信号压到 0。

**结论**：不要看"各批次均值是否拉齐"来判断批次效应是否去掉——分位数归一化
能做到均值拉齐却完全没去批次。要看**批次对主成分的解释方差**（脚本自动算在报告里）。
跨平台（尤其 RNA-seq/TCGA 与芯片混用）合并时**不要只做 quantile 归一化就当成
"去过批次"**，先用 `--method=both` 各出一份，对比 `_pca_after` 与报告里的量化表再决定。

### 5.2 参数

- `--combat_prior=TRUE`（默认，参数化先验）。`FALSE` 用非参数先验，
  批次内基因数少时更稳，但更慢。
- `--combat_mean_only=TRUE` 只校正均值不校正方差。
  **批次数里有单样本批次时 sva 会自动强制 `mean.only=TRUE`** 并打印提示。
- `--combat_ref=<批次标签>` 指定参考批次，该批次不被改动。
  标签必须是 `levels(batch)` 之一，否则 sva 报
  `reference level ref.batch is not one of the levels of the batch variable`。

---

## 6. `ComBat` 的 `mod`：必须用带截距的 `model.matrix(~ group)` ★

这是最容易写错、且报错信息完全指错方向的一处。看 `sva::ComBat` 的实际逻辑：

```r
design <- cbind(batchmod, mod)
check  <- apply(design, 2, function(x) all(x == 1))   # 找出"全 1"的列
design <- as.matrix(design[, !check])                 # 丢掉全 1 列
if (qr(design)$rank < ncol(design)) {
  if (ncol(design) == (n.batch + 1))
    stop("The covariate is confounded with batch! Remove the covariate and rerun ComBat")
  if (ncol(design) > (n.batch + 1))
    ...
    stop("The covariates are confounded! ...")
}
```

关键在这一步：`mod` 里**全 1 的截距列会被自动丢掉**。所以：

| 写法 | design 列数（n.batch=3、两组） | 后果 |
|---|---|---|
| `model.matrix(~ group)` ✅ | 3 个 batch 列 + 1 个 group 虚拟变量（截距被丢） | 正常 |
| `model.matrix(~ 0 + group)` ❌ | 3 个 batch 列 + **2** 个 group 虚拟变量 | 两个虚拟变量之和 = 全 1 向量 = batch 列之和 → 共线 → 秩亏 → **误报 `covariates are confounded!`** |

**结论：脚本固定用 `stats::model.matrix(~ grp)`（带截距）。** 不要"顺手"改成
`~ 0 + grp`，否则会得到一个看起来很吓人、但完全是误报的错误。

### 6.1 批次与分组混杂

若某些批次内**只有一个分组水平**（例如每个 GSE 只含一种组织），则"批次效应"与
"生物学差异"在数学上不可分。上面的 `qr(design)$rank` 检查会失败，ComBat 抛
`The covariate is confounded with batch!`。

脚本在调用 ComBat **之前**就主动复算了同样的检查，报出更可操作的错误：

```
分组与批次混杂: 用 mod=~group 时设计矩阵秩不足，ComBat 会拒绝运行。
  批次 x 分组 交叉表:
       group
  batch  Normal Tumor
  GSE1        0   100
  GSE2      120     0
  GSE3       80     0
  含义: 某些批次内只有一个分组水平，无法把"批次效应"与"生物学差异"分开估计。
  可选做法:
    1) --combat_mod=none 不做生物变异保护 (仅当批次间分组构成相似时才安全)
    2) 只合并分组构成相近的数据集
    3) --method=quantile (不做批次参数估计，不受此限制)
```

**注意第 1 条是有条件的**：`mod=NULL` 时 ComBat 不知道分组，会把"组织间差异"
也当成批次偏移扣掉。只有当各批次的分组构成相似（每个批次都同时有两种组织）时，
这样做才安全。

### 6.2 合并后的分组从哪来：把 N 份临床信息交给用户 ★

合并矩阵的样本来自多个数据集，**每个数据集的临床列名与取值词表都不一样**：

| 数据集 | 列名 | 取值 |
|---|---|---|
| GSE55235 | `source_name_ch1` | `synovial tissue from healthy joint` / `... osteoarthritic joint` / `... rheumatoid arthritis joint` |
| GSE55457 | `clinical status:ch1` | `normal control` / `osteoarthritis` / `rheumatoid arthritis` |

同一个概念在两边叫法完全不同（连列名都不是同一个），所以"哪一列、哪些取值算
case/control"**只能由用户回答**。03 步在 `--clinical_files`（03b 会自动填
`per_gse/clinical_<GSE>.csv`）下产出两份材料：

**① `<prefix>_clinical.csv`** —— 合并临床表，行名 = 样本名（与 01 的
`clinical_<GSE>.csv` 同格式，02 步 `--clinical=` 可直接读）。列包含：

- `dataset` / `batch` / `group`（分组已提供时）
- 各数据集的临床列。**所有提供了临床表的数据集都有的同名列保留原名**
  （值仍是各自词表，需人工核对），**只出现在部分数据集的列加 `<数据集>::` 前缀**
  避免撞名互相覆盖。`dataset`/`batch`/`group` 这三个保留名不会被临床列覆盖。
- 一行只填它所属数据集的那部分列，别的数据集那些列是空——这是故意的：
  一眼能看出"这一列这个样本没有信息"。

**② `<prefix>_grouping_candidates.txt`** —— 给用户看的主材料：

- **逐数据集的候选分组列**与取值计数（取值 2..min(10, n/2) 个、每类 ≥2 个）
- 每个取值后面的 `-> 分组名`：**这是映射结果的审计**。
  `-> 混合(case/control)` 说明这个取值被映射到了多个分组（通常意味着映射写错了）；
  `-> (未进入矩阵, 已剔除)` 说明整类样本被 `--group_na=drop` 剔除了；
  `[保留 3/10]` 说明该取值只有一部分样本进了矩阵
- **跨数据集的同名临床列矩阵**，带 `*` 的是"在某个数据集里算候选分组列"的列——
  这一节大多是平台/联系方式等元数据，没有标记看不出哪几列值得看
- 下一步命令模板（含"用户自己的分组表"这一条）

**③ `<prefix>_group_template.csv`** —— 分组模板：`sample,group` 两列、全部样本已按
矩阵顺序列好，已推导出分组时预填，未分到组的格子留空（不是字面量 `NA`）。
用户填好回传、或直接换成他自己手上的表，都能当 `--group_file=` 用。

> **候选判定用的是该数据集的完整临床表**（而不是"幸存样本"）。否则某个取值整类
> 被剔除后，这列在幸存样本里只剩 1 个取值，会被当成"单一取值"直接跳过——而这
> 恰恰是用户最需要看到的情况（"我把 RA 全剔掉了，是不是连列都选不出来了"）。

> **第四种回答方式（也是最省事的）：让用户直接上传他的分组列表。**
> 用户手上常已经有分组表（样本清单、表型汇总、文章补充材料）。**两列就够**：
> 第一列样本名、第二列分组，表头可有可无，逗号或制表符都行，**不必覆盖全部样本**。
> ```bash
> 02_deg_plots.R --expr=<prefix>_merged.csv --group_file=<用户的分组表> --contrast=case-control
> # 或合并阶段就用它（ComBat 的 mod 会跟着保护该分组）:
> 03b_merge_from_geo.R --stage=merge --group_file=<用户的分组表>
> # 样本少时也可行内给: --group_spec='case=GSM1,GSM2;control=GSM3'
> ```
> 样本名要与合并矩阵列名一致（= 各 01 步产物列名，可用 `<prefix>_batch_map.csv`
> 第 1 列核对）；02 步会输出 `<prefix>_sample_group_map.csv` 供核对；未覆盖的样本
> 会被 02 步自动剔除。**不要替他改写组名**，`--contrast` 用原组名（`make.names` 后）。

拿到用户答复后（他按"哪一列 / 哪些取值 / 什么方向"回答的那三条路）：

| 路径 | 命令 | 特点 |
|---|---|---|
| A. 重建分组表 | `03b_merge_from_geo.R --stage=group --group_cols=... --group_maps=...` | **不重算 ComBat**；样本范围取自合并临床表（被剔掉的样本不会回来） |
| B. 直接用某列 | `02_deg_plots.R --clinical=<prefix>_clinical.csv --group_from=<列> [--group_keep=A,B]` | 该列取值即最终分组名；不适用于需要"取值→统一名"映射的场景 |
| C. 用户自己的表 | `--group_file=<用户的分组表>` | 最省事；两列即可，不必覆盖全部样本；合并阶段用还会影响 ComBat 的 `mod` |

> ⚠️ 路径 A 换分组**不会**改变合并矩阵里 ComBat 用的 `mod`——它保护的仍是合并时
> 那个分组。要让 ComBat 也保护新分组，得带 `--group_cols/--group_maps` 重跑
> `--stage=merge`（ComBat 很贵，所以默认不做）；路径 C 在合并阶段用就是这个效果。
>
> ⚠️ 也不要图省事把**数据集/批次**当分组：若"数据集 A 全是 case、B 全是 control"，
> `mod=~group` 与批次列共线，脚本会在调用 ComBat 前拦住（6.1 节）。

---

## 7. NA 与去重

### 7.1 01 → 03 时 NA 是从哪来的（2026-09-17 实测）

用户常在这里疑惑："原始 data 明明是干净的，怎么合并时报 NA？"

**来源是 `01` 的 log2 步骤**：文档化行为——log2 前会把 **≤0 的值置为 NA**
（见 `pitfalls.md` 3.1），而且 01 会把这件事打印出来，不是静默操作：

```
[2/6] 判断是否需要 log2 转换 ...
  已做 log2 转换 (9 个 <=0 的值置为 NA)
  分位数: min=0 25%=17.9 中位数=63.7 75%=172.7 max=18066.8
```

实测 GSE55235：原始数据是**未取对数的荧光强度**（Q75=172.7、max=18066），
01 正确地判定需要 log2；同时把 9 个 ≤0 的值置为 NA。按基因取 max 聚合后，
仍有 **4 个单元格**是 NA（2 个基因的各别样本）。

于是 03 步按默认 `--na=error` 停下来，提示可选 `--na=rowmean`。
**这不是 bug，是两件正确行为的组合**：01 如实标记出不可对数化的值，03 拒绝把 NA
喂给 `ComBat`（sva 遇到 NA 只会抛一句无信息量的
`Data contains NA. To proceed, please remove or impute missing values.`）。

**遇到时怎么处理**：

- 只有零星几个 NA → `--na=rowmean`（用该基因的行均值填补，最常用）
- NA 占比不小 → 说明数据可能已经做过背景校正、或本来就已取对数，别硬转：
  回 01 用 `--force_log=FALSE` 重跑，或 `--norm=none` 跳过归一化
- 用 `--match=union` 时 NA 必然大量出现（某数据集没测到的基因），此时必须
  配 `--na=rowmean`

### 7.2 ComBat 与去重

- **`ComBat` 不接受 NA**（`stop("Data contains NA. ...")`）。默认 `--na=error`，
  遇到 NA 会打印缺失最多的基因名并给出三个可选项，而不是让用户撞 sva 的报错。
  `--na=rowmean` 用该基因的行均值填补（最常用），`--na=zero` 置 0。
- **同名基因**：默认 `--dup=max`，先按行均值降序排再每个名字取第一条
  （等价于"保留整体表达量最高的探针"）；`--dup=first` 保留首次出现的。
- **占位基因名**（空、`---`、`NA`）会被剔除并报告个数。
- **样本名跨数据集重复**时脚本自动加后缀并**明确警告**——这会影响
  `--group_file` / `--batch_file` 的匹配，应在上游改唯一。

---

## 8. QC 判读

### 8.1 批次效应量化

`<prefix>_overlap_report.txt` 里有一张客观表格：用 `lm(PC_k ~ batch)` 的
决定系数 `R²` 与 ANOVA p 值衡量"批次能解释多少主成分方差"。

```
批次效应量化 (批次对主成分解释的方差比例, ANOVA by batch):
                矫正前                矫正后
  PC1   42.3% (p=1.2e-80)     3.1% (p=0.021)
  PC2   11.7% (p=3.4e-12)     2.0% (p=0.180)
```

- 矫正后 `R²` **大幅下降**、`p` 不再极端 → 批次效应被扣掉了
- 矫正后仍很高 → 批次没扣干净，或该数据集本身与其它数据集差异过大
- 同时看 `<prefix>_pca_after_group.pdf`（给了分组时自动输出）：
  **分组仍要分得开**。若分组也被扣没了，说明 `mod` 没起作用——
  查报告里 `ComBat 保护生物变异 : 是 / 否`

### 8.2 批次统计表

`<prefix>_batch_stats.csv` 每个批次的 `mean/median/sd` 前后对照。
矫正后各批次的 `mean`/`median` 会明显靠拢；`mean.only=TRUE` 时 `sd` 不变。

### 8.3 图

- `_boxplot_before/.after`：**每个样本一个箱体**，按批次着色。矫正后各批次的箱体
  应该落在同一高度带。❗样本数大时（如 967）箱体极窄，因此脚本**不画边框**
  （`colour = NA`）——实测带灰色边框会把批次填充色完全盖住、整张图糊成一片。
- `_boxplot_batch_before/.after`：**每个批次一个箱体**（箱体来自批内各样本的中位数）。
  样本数很大时这张图仍然清晰，是判断"各批次整体水平是否被拉齐"最直接的一张。
  矫正前批次间箱体高度差很大，矫正后应基本齐平。
- `_density_before/.after`：样本密度曲线。样本过多时按 `--density_max` 抽样
  （默认 200，`set.seed(1)` 可复现），抽样会打印提示。
- `_cor_heatmap`（`--cor_heatmap=TRUE`，默认关）：样本相关性热图，按批次抽样
  后按批次排序、不做聚类，块状结构变模糊说明批次效应被削弱。

---

## 9. 验证记录（2026-09-17）

### 9.1 数值一致性 vs `7.R`

验证方法：在同一台机器上用 `7.R` 的流程（两个代码块）产出真值并
`saveRDS()`，再用本脚本跑同一份输入，把内存中的矩阵逐元素比较
（`gt_allmerge1.rds` / `gt_combat_block1.rds` vs `03` 步 `<prefix>.RData`
里的 `merged_raw` / `merged_combat`）：

| 项目 | 结果 |
|---|---|
| 维度 | 均为 15913 × 967 |
| 未矫正合并矩阵 | `identical() = TRUE`、`max\|diff\| = 0` |
| ComBat 输出 | `identical() = TRUE`、`max\|diff\| = 0`、`mean\|diff\| = 0`、相关 = 1 |
| 全矩阵均值/中位数/极差 | 两者完全相同 |

**顺带验证了一个实现判断**：`7.R` 传的批次是整数 `c(1,1,…,2,2,…,3,3,…)`，
本脚本传的是**以数据集名命名的 factor**。两者结果**逐位相同**——
说明 `sva::ComBat` 在 `par.prior=TRUE` 下每个批次的参数是独立估计的，
结果与批次标签的取值/顺序无关，只与"哪些样本属于同一批"有关。

（`7.R` 的两个代码块本身互为显式版/循环版，都得到 15913 基因、967 样本。）

### 9.2 数据集与规模

| 数据集 | 基因 × 样本 | `auto` 判定 |
|---|---|---|
| `TCGA_LUSC_TPM.txt` | 19937 × 553 | `TCG*` → `log2(x+1)` |
| `GSE30219.txt` | 20824 × 307 | `GSM*` → 仅 `quantile` |
| `GSE74777.txt` | 30905 × 107 | `GSM*` → 仅 `quantile` |

- `intersect(1,2)` = 17138 → `∩ 3` = **15913**
- 合并矩阵：**15913 × 967**，无 NA
- `--toupper_scan` 大写收益 **+0**
- ComBat（`par.prior=TRUE`、`mod=NULL`）约 **11~15 秒**；
  sva 提示 `Found 22 genes with uniform expression within a single batch ...`

**该数据上批次矫正的效果（`mod=NULL`）**：

| 指标 | 矫正前 | 矫正后 |
|---|---|---|
| 批次对 PC1 的解释方差 | **99.3%**（p≈0） | **0.0%**（p=1） |
| 批次对 PC2 的解释方差 | 97.3% | 0.0%（p=0.998） |
| PC1 方差占比 | 79.7% | 14.7% |
| 各批次全矩阵均值 | 3.289 / 6.673 / 5.909 | **4.650 / 4.658 / 4.657** |
| 各批次样本中位数 | 3.370 / 6.769 / 5.994 | 4.795 / 4.776 / 4.766 |

PCA 图上矫正前三个数据集是完全分离的三团，矫正后大幅重叠（仅 TCGA 侧残留
一小簇，属该队列自身的亚群结构，不是矫正失败）。

### 9.3 性能与内存

| 环节 | `read.table`（7.R） | `fread`（本脚本） |
|---|---|---|
| TCGA 19937×553（79 MB） | 26.9 s | **2.4 s** |
| GSE30219 20824×307（58 MB） | 19.8 s | **1.7 s** |
| GSE74777 30905×107（28 MB） | 10.5 s | **0.8 s** |

- 15913×967 的三份矩阵（raw / combat / quantile）同时驻留内存时，
  进程峰值约 **1.6 GB**（本机 16 GB）。
- **不要同时跑多个大 R 任务**：并发时物理内存被挤占会触发换页，
  PCA 一步会从几十秒退化到数分钟（实测踩过）。

各步骤耗时（15913×967，本机实测）：

| 步骤 | 耗时 |
|---|---|
| 读 3 个矩阵（`fread`） | ~5 s |
| 3 个数据集预处理（含 2 次 `normalizeBetweenArrays`） | ~20 s |
| 取交集 + `cbind` | ~5 s |
| `ComBat(par.prior=TRUE)` | **11~15 s** |
| PCA（`svds`，前后各一次） | ~13 s |
| 箱线图（967 个箱体） | 每个约 **20 s**（PDF+PNG） |
| 密度图（抽 200 条曲线） | ~5 s |
| **导出 4 个文本文件**（raw/merged × csv/txt） | **数分钟（整步最慢）** |

- **导出往往是最慢的环节**：脚本用 `data.table::fwrite`（比 `write.csv` 快约一个数量级，
  与 `fread` vs `read.table` 同理）。只要 02 步的输入时可加
  `--save_raw=FALSE --save_txt=FALSE`，只留 `<prefix>_merged.csv`。

### 9.4 图内文字必须用 ASCII

本次实测：中文图标题/图例在 `ragg::agg_png` 与 `cairo_pdf` 下都被渲染成**乱码 ASCII**
而不是方框——`PCA - 批次矫正前` → `PCA - f 9f,lg □+f-#e □□`，图例 `批次` → `f 9f,!`。
很容易被误判成"标题字符串写错了"。**本技能约定绘图标签一律 ASCII**（02 步一向如此），
控制台 `message()` 里的中文不受影响。详见 `pitfalls.md` 7.8。

### 9.5 编排器（`03b_merge_from_geo.R`）端到端验证（2026-09-17）

数据取自真实项目（TCGA_LUSC + 两个 GSE 芯片数据集）里的
**GSE55235 + GSE55457**（同平台 GPL96、同为 OA 滑膜），并复制到中文目录
`D:/Rtmp/_测试数据/` 下，以同时验证非 ASCII 路径。

刻意构造**最坏情况**：中文的 skill 路径 + 中文数据目录 + 中文输出目录，且故意保留
环境里的 `LC_ALL=C.UTF-8`（即不经 `run_geo.sh` 启动），检验 locale 自愈。

| 阶段 | 结果 |
|---|---|
| locale 自愈 | `C` → `Chinese (Simplified)_China.utf8`，中文路径随即可用 |
| 文件发现 | `--datadir` 递归找到两个 GSE 的 matrix + soft（修 locale 前全部报"未找到"） |
| `--stage=normalize` | 两个 GSE 均成功：GSE55235 **13237 基因 × 30 样本**、GSE55457 **13237 × 33** |
| STOP 点输出 | 正确打印两数据集全部候选分组列与取值分布 |
| `--stage=merge` 分组推导 | GSE55235 `source_name_ch1` → case=10/control=10；GSE55457 `clinical status:ch1` → case=10/control=10；23 个 RA 样本留 NA 并按 `--group_na=drop` 剔除 |
| 最终合并矩阵 | **13236 基因 × 40 样本**（`--gene_case=upper`） |
| 批次矫正 | 批次对 PC1 解释方差 **99.7%（p=3.05e-50）→ 0.0%（p=0.903）** |
| 批次均值 | 6.256 / 7.836 → 7.056 / 7.035 |
| 产物 | 27 个文件（含 `per_gse/` 各 7 个 + `logs/` + `merge/` 全部 QC） |

**与用户历史手工流程的对照**：他们的 `merge处理.R` 合并 `55235_row.csv` +
`55457_row.csv`，再 `subset(sample == 'OA' | 'Normal')`，最终也是 **40 个样本**，
与编排器一致。基因数不同（手工 13433 / 13435 vs 本流程 13236）来自注释与占位符
清理口径的差异——**13237 正是 GPL96 在 `01` 下的基线值**（与 GSE55584 基线一致），
说明本流程自身一致，不是回归。

这次还顺带跑出一个 01→03 衔接的真实案例：GSE55235 的原始数据是**未取对数的荧光强度**
（Q75=172.7、max=18066），`01` 正确判定需要 log2，同时把 9 个 ≤0 的值置为 NA；
按基因取 max 聚合后仍有 4 个单元格是 NA，于是 `03` 按默认 `--na=error` 停下并提示
`--na=rowmean`。详见第 7.1 节。

### 9.6 合并临床信息与分组重建的验证（2026-09-17）

同一对数据集（GSE55235 + GSE55457，最终 13236 × 40）上：

| 项目 | 结果 |
|---|---|
| `--stage=merge` 自动传临床表 | `临床表: 2/2 个数据集有临床信息`；产出 `m3_clinical.csv`（40 样本 × 43 个临床列 + `dataset`/`batch`/`group`，其中 10 列带 `<数据集>::` 前缀）与 `m3_grouping_candidates.txt` |
| 报告里的映射审计 | GSE55235 的 `synovial tissue from rheumatoid arthritis joint` 显示 `10 [保留 0/10] -> (未进入矩阵, 已剔除)`；`gender:ch1` 显示 `-> 混合(case/control)`（因为分组用的是别的列） |
| 报告里的跨数据集列矩阵 | 33 个同名列表在 ≥2 个数据集，带 `*` 标出候选列（`source_name_ch1 *` / `characteristics_ch1 *`） |
| `--stage=group` 重建分组表（映射成 3 组） | 样本数仍是 **40**（20+20）——RA 的 23 个样本不会因为映射里多了 `ra` 就回来；分组计数 case=20 / control=20，NA=0 |
| 02 步直接吃合并临床表 | `02_deg_plots.R --clinical=m3_clinical.csv --group_from=group --contrast=case-control` 正常读入 40 样本，跑完 limma 出图 |
| 合成小数据专项（2 数据集 × 10 样本） | 列名前缀规则正确（`tissue:ch1` 为共有列用原名、数据集专属列带 `<数据集>::`）；`--group_na=drop` 后临床表与候选报告的样本范围同步收缩 |

**过程中修掉的一个真实 bug**：03b 自动定位临床表时把文件名拼成了 `<GSE>_clinical.csv`，
而 `01` 实际写的是 `clinical_<GSE>.csv`（**前缀在前**）。后果是"报告说没有任何数据集
找到临床表"，但 `per_gse/` 里明明有——因为 `ifelse(file.exists(p), p, "")` 把不存在的
路径静默转成了空串。修好后 `临床表: 2/2`。

---

## 10. 快速排错

| 现象 | 原因 / 处理 |
|---|---|
| `The covariate is confounded with batch!` | 批次与分组混杂，见 6.1 |
| `The covariates are confounded!`（分组明明没混杂） | `mod` 写成了 `~ 0 + group`，见第 6 节 |
| `Data contains NA. To proceed, ...` | 用 `--na=rowmean`，见第 7 节 |
| `reference level ref.batch is not one of the levels` | `--combat_ref` 的值必须是批次标签之一 |
| `Note: one batch has only one sample, setting mean.only=TRUE` | 正常提示，非错误 |
| `共同基因为 0 个` | 基因命名不一致，试 `--gene_case=upper` 或 `--match=union` |
| 交集比预期少很多 | 看报告的大小写敏感性一行；`upper` 可能有正收益 |
| 批次效应 PCA 矫正后没改善 | 看报告 `ComBat 保护生物变异 : 是/否`；也可能该数据集本身差异过大 |
| `--group_file` 匹配率低 | 样本名大小写/后缀不一致；或跨数据集样本名重复被加了后缀 |
| 输出 PDF 是空白 | `pheatmap` 必须用设备包裹（脚本已处理）；勿改回 `ggsave` |
| `invalid font type` / Arial not found | 默认 `pdf()` 不认 Arial，必须 `cairo_pdf`（脚本已处理） |
| 想直接拿 `_batch_map.csv` 当 02 步的 `--group_file` | **不行**，它的第 2 列是 `dataset`。用 `<prefix>_group.csv` |
| 没有 `<prefix>_clinical.csv` / `_grouping_candidates.txt` | 03 没收到 `--clinical_files`。03b 会自动填 `per_gse/clinical_<GSE>.csv`（**注意是前缀在前**）；手工跑 03 时要自己给 |
| `--stage=group` 报找不到合并临床表 | 上次 `--stage=merge` 没带临床表（或没跑过 merge）。重跑 `--stage=merge` |
| `--stage=group` 里 RA 样本"消失"了 | 正常：样本范围取自合并临床表，被 `--group_na=drop` 剔掉的样本不会回来 |
| 02 步报 `临床表中无列 'xxx'` | 合并临床表里跨数据集的列名带 `<数据集>::` 前缀，报错信息会列出全部可选列 |
