# 方法学说明：两级 KM 策略、surv_cutpoint 与最佳阈值的代价

本文件是 `SKILL.md` 的展开，讲清"为什么这样设计"和"这样设计的代价是什么"。
改动脚本前建议先读第 3 节。

---

## 1. 第一级：分位值分组（常规做法）

对应参考脚本（常规分组 KM）标准写法：

```r
group <- ifelse(rt[, "GENE"] > quantile(rt[, "GENE"], seq(0, 1, 1/2))[2], "High", "Low")
diff  <- survdiff(Surv(time, state) ~ group, data = rt)
pValue <- 1 - pchisq(diff$chisq, df = length - 1)
```

实现要点（`quantile_split()`）：

- **用 `>` 而不是 `>=`**，与参考代码逐字一致。代价是**并列值全部落到 Low 侧**。
  表达量里有大量 0（未检出基因）或大量重复值时，这会让分组严重失衡，
  甚至出现 Low=95% / High=5% 这种"分组其实不成立"的情况。
  脚本在任一组不足总量的 10% 时打印 `[警告: 分组严重不平衡]`。
- 分位点由 `--pct` 控制：`0.5` 是中位值（参考代码默认 / 注释里的 `1/2`），
  `0.3333` 对应参考代码注释掉的 `1/3`（三分位下界）。
- 这个 p 值**不需要任何校正**：切点是预先固定的，不依赖结局数据。
  所以只要 p<0.05，就直接报告，不需要额外说明。

`km_stats()` 里顺带算了 Cox 单因素：

```r
coxph(Surv(time, event) ~ group)   # group 为 factor，levels = c("Low","High")
```

**HR 的语义恒为 High vs Low**，与画图风格、与图形里 Low/High 的先后无关。
levels 顺序被显式固定为 `c("Low","High")`，因此系数就是 High 相对 Low 的风险比。

---

## 2. 第二级：最佳截断值（`surv_cutpoint`）

对应参考脚本（最佳截断值 KM）标准写法：

```r
best_threshold_surv <- surv_cutpoint(data, time="OS.time", event="OS",
                                     variables="risk_score", minprop=0.3)
data <- surv_categorize(best_threshold_surv)
```

### 2.1 它内部在做什么

`survminer::surv_cutpoint` 逐个候选切点计算**标准化 log-rank 统计量**，
取使该统计量最大的切点作为"最佳阈值"：

```
maxstat.test(Surv(time, event) ~ var, data,
             smethod = "LogRank", pmethod = "none", minprop = minprop)
```

- `minprop`：**每组**至少要占的比例（默认脚本用 0.3）。它决定了候选切点的
  搜索范围 —— `minprop=0.3` 时切点必须落在 30%~70% 分位之间。
  本机实测：`re.csv` 的 SOX9 在 `minprop=0.3` 下切点为 153.886，
  分成 Low=192 / High=390（33% / 67%），符合约束。
- `surv_categorize(cp, labels=c("low","high"))` 返回**字符向量**（不是 factor），
  取值按 `数值 ≤ 切点 → low`、`> 切点 → high` 划分。
  脚本显式传 `labels = c("Low","High")` 并自己 `factor(..., levels=c("Low","High"))`，
  保证 HR 方向可控。

### 2.2 返回结构（实测，survminer 0.5.2）

```r
names(cp)              # c(<变量名>, "data", "minprop", "cutpoint")
cp$cutpoint            # data.frame: cutpoint, statistic（标准化 log-rank 统计量 M）
cp[[<变量名>]]         # maxstat "maxtest" 对象
cp[[<变量名>]]$p.value # **恒为 NA**（内部用的 pmethod="none"）
```

**两个容易踩的点**：

1. `maxstat` **不导出** `pvalue()`（`ls("package:maxstat")` 里没有它），
   所以网上常见的 `maxstat::pvalue(x)` 写法会直接报
   `'pvalue' is not an exported object`。
2. 即使能取到，`cp[[var]]$p.value` 也是 NA —— 因为 `surv_cutpoint` 用
   `pmethod="none"` 调用，根本没算 p。

所以**校正 p 必须自己重跑一次 maxstat**（脚本已内建，见第 3 节）。

### 2.3 什么情况下会失败

`surv_cutpoint` 抛错时脚本会捕获并打印原因，最常见的是：

```
Error in maxstat.test(...) : No cutpoint of sufficient size
```

含义是**在 `minprop` 约束下找不到合法的候选切点**，典型成因：

- 变量取值过于集中（例如一列里 80% 是 0），30%~70% 分位之间根本没有不同的取值
- 样本量太小，`minprop` 要求的分组规模达不到

脚本的处理：打印原因 + 三条出路（调小 `--minprop`、改 `--pct`、
先按常规分组出图再人工定阈值），**然后自动退回第一级的分位值分组结果**，
不会静默失败、也不会卡住不出图。

---

## 3. ⚠️ 最佳阈值的代价：必须报告的"校正 p"

这一节是本技能最需要交付给用户的判断。

### 3.1 问题

最佳阈值是从**同一批数据**里选出来的。选切点时你已经把所有候选切点都看了一遍，
挑了 log-rank 统计量最大的那个 —— 等价于做了成百上千次隐式检验并只保留最好的。
于是：

- 图上报的 `log-rank p`（对切分后的两组做一次 `survdiff`）是**朴素 p**，
  它完全不知道"这个切点在被选出来之前试过多少候选项"，因此**系统性偏乐观**。
- 这就是所谓的 **optimally selected cutpoint 问题 / 择优偏倚**。
  在纯噪声数据上，只要候选切点足够多，也总能找到一组"看起来显著"的分组。

### 3.2 实测（2026-09-18，`D:/phlp/tcgay预后/re.csv`，582 例 / 119 事件）

| 基因 | 常规中位值分组 p | 最佳阈值 | 切分后**朴素** log-rank p | **maxstat 校正 p (Lau92)** |
|---|---|---|---|---|
| TIMP1 | 0.0312 | —（未回退） | — | — |
| SOX9 | 0.4304 | 153.886 | **0.0318** | **0.2998** |
| COMT | 0.6453 | 43.411 | 0.1570 | 0.7393 |

SOX9 是最典型的例子：朴素 p 是 0.0318，看起来显著；**扣除择优偏倚后 p=0.30，
完全不显著**。如果不报校正 p，这个结果会被当成阳性结论写进文章。

### 3.2b 批量场景下的放大效应（39 基因实测）

把 `re.csv` 的 39 个基因整批跑一遍（`--mode=auto`，其余用默认值）：

| 项 | 数量 |
|---|---|
| 基因总数 | 39 |
| 采用中位值分组 | 7 |
| 采用最佳阈值（回退） | **32** |
| 朴素 log-rank p < 0.05 | **23** |
| —— 其中 median 来源 | 7 |
| —— 其中 cutpoint 来源 | **16** |
| **cutpoint 来源且 maxstat 校正后仍 < 0.05** | **1（CCNB1）** |
| 校正 p 成功取到值 | 32 / 32 |

也就是说：**在 39 基因里"筛出 23 个显著基因"这个结论，把 32 个基因的
最佳阈值全部换掉后，只有 CCNB1 一个还能站住。** 校正前后对比（按朴素 p 升序前 8）：

| 基因 | 常规中位值 p | 最佳阈值朴素 p | maxstat 校正 p |
|---|---|---|---|
| CCNB1 | 0.0725 | 0.00090 | **0.0108** ✅ |
| MYC | 0.188 | 0.0087 | 0.108 |
| MMP12 | 0.0767 | 0.0093 | 0.101 |
| RRM2 | 0.564 | 0.0107 | 0.139 |
| CXCL10 | 0.255 | 0.0150 | 0.160 |
| PARP1 | 0.307 | 0.0192 | 0.218 |
| PTGS2 | 0.0858 | 0.0211 | 0.203 |
| BUB1B | 0.274 | 0.0217 | 0.214 |

**这就是"用最佳阈值把 p 做出来"的真实代价。** 遇到"批量筛预后基因"的需求时，
本技能能照做，但**交付必须同时给三列**：常规 p、阈值朴素 p、maxstat 校正 p，
并明确说清"按校正 p 只剩 N 个"。把 32 个基因的阈值都调一遍再报 23 个阳性，
在审稿和独立复现面前是站不住的。

### 3.3 脚本怎么做

`best_cutpoint()` 在拿到 `surv_cutpoint` 的结果后，**用相同的 `minprop`、
相同的数据再跑一次 maxstat**，取校正 p：

```r
maxstat::maxstat.test(Surv(time, event) ~ V, data = ...,
                      smethod = "LogRank", pmethod = pmethod, minprop = minprop)
```

实测该重跑**逐位复现** `surv_cutpoint` 的切点与统计量
（SOX9：切点 153.886、M=2.1028，两者一致）。

`pmethod` 的选择（`--maxstat_pmethod=`，默认 `Lau92`）：

| pmethod | 实测表现 | 采用 |
|---|---|---|
| `Lau92` | p=0.2998，取值合法、秒级 | ✅ 默认 |
| `Lau94` | p=**1.156**（越界 >1，数值上不可用） | 需显式指定 |
| `exactGauss` | 作为 `Lau92` 不可用时的备选 | 仅兜底 |
| `HL` / `condMC` | 明显更慢（蒙特卡洛），不适合批量 | 不采用 |

越界（>1）或非有限的结果**一律置 NA 并说明**，不猜、不截断。

### 3.4 交付话术（不能省）

用了最佳阈值时，输出里会带上：

```
4) 已使用数据驱动的最佳截断值：阈值是在同一批数据上选出来的，属探索性分析。
   报告时不要只写 log-rank p，请一并给出 maxstat 校正 p；正式结论需在独立队列复现。
```

给用户讲解时的要点：

1. 说清用的是常规分组还是最佳阈值（看 `*_km_decision.txt` 第二级那行）
2. 若用了最佳阈值，**两个 p 都给**：朴素 p 对应图上的字，校正 p 对应可信度
3. 校正后不显著就**说它不显著**。不要为了 p<0.05 继续调 `--minprop` / `--pct`
   直到凑出来 —— 那是把择优偏倚叠成两层，比不做校正更糟
4. 想要一个能站住脚的结论，正道是**独立队列验证**（用 A 队列选阈值，
   在 B 队列用这个阈值分组再算 p）

---

## 5. 出图风格：为什么默认是 unified

`生存分析-km.R` 与 `最佳阈值-km.R` 这两份参考代码**画出来的图不是一套版式**：

| 外观项 | `生存分析-km.R` | `最佳阈值-km.R` |
|---|---|---|
| 配色 | Low `MediumSeaGreen` / High `Firebrick3` | Low `#3090a1` / High `#bc5148` |
| 图例 | 标题 `<基因> expression`，标签 `Low`/`High` | 无标题，标签 `<基因>_low`/`<基因>_high` |
| 置信带 | 关 | 开 |
| 主题 | survminer 默认 | `theme_minimal(base_size=14)`，无网格、黑轴线 |
| p/HR 标注 | ggsurvplot 的 `pval`（面板内左上） | 手工 `annotate`（面板内左下、左对齐） |
| p/HR 顺序 | p 在前 | HR 在前 |
| y 轴标题 | `Survival probability`（survminer 默认） | `Survival Probability` |
| 风险表标签 | 黑色 | 跟曲线配色 |

如果按"哪种分组就用哪种参考的风格"来出图，同一批结果里就会出现两种版式 ——
**看图的人一眼就能分辨哪张是回退来的**，这在成图交付时是不能接受的。

所以默认 `--style=unified`，把冲突项一次性定死：

- **取自主参考 `生存分析-km.R`**：配色（红/绿）、图例形式（标题带基因名）、
  survminer 默认主题、`pval` 标注机制、p 在前的两行排布、置信带默认关
- **取自 `最佳阈值-km.R`**：风险表标签跟曲线配色（`risk.table.col="strata"`）
- **显式写死**：y 轴标题 `Survival probability`、图例位置 (0.8, 0.8)、
  横向刻度间隔 —— 免得随 survminer / ggplot2 版本漂移

`classic` 与 `unified` 的唯一外观差别就是风险表标签配色（前者黑色，严格复刻）。
`threshold` 完整复刻参考 2。

### 5.1 风格与分组方式必须解耦（可核验）

这一点是硬约束，不是约定：风格只由 `--style` 决定，代码里**没有任何一处**
让 `mode_used`（median / quantile / cutpoint）参与绘图分支。
`make_km_plot(fit, d, gene, style = ...)` 里 `style` 是唯一的外观开关，
`fit`/`d` 只携带数据。

可核验：同一基因分别用 `--mode=median` 与 `--mode=cutpoint` 出图，把两张图的
`plot$theme`、`plot$labels`、`table$theme`、`table$labels`、图层构成、
colour scale 调色板逐一比对，应当**逐项 identical**。实测（`re.csv` 的 SOX9，
`--style=unified`）12 项指纹全部一致：

```
plot_theme_n / plot_theme_v / table_theme_n / table_theme_v
plot_labels  / table_labels / plot_layers  / table_layers
plot_class   / table_class  / guides_n     / theme_complete    -> 全部 一致
配色  median : MediumSeaGreen, Firebrick3, #6E568C, #223D6C
配色  cutpoint: MediumSeaGreen, Firebrick3, #6E568C, #223D6C
```

批量也同样成立：`02_km_batch.R` 里 `style` 是整批共用的单个值，
混合分组的批次（如 TIMP1 走 median、SOX9 走 cutpoint）出的图版式完全一致。

> ⚠️ **想换配色请用 `--palette=色1,色2`，不要切到 `--style=threshold`。**
> 后者会连带换掉置信带、主题、图例形式、标注位置，统一性就没了。

---

## 6. 时间单位

参考代码里时间处理是隐式的，也是全流程最容易出错的一环：

- `生存分析-km.R`：`cli$time = cli$time/12`，横轴 `Time(years)`
  → 前提是输入 time 的单位是**月**
- `最佳阈值-km.R`：`OS.time/30.44`，横轴 `Time (months)`
  → 前提是输入是**天**

同一份"TCGA 生存数据"在不同下游脚本里可能是天、月或年，没有任何元数据标记它。
本机的实例：

| 文件 | 中位数 | 最大 | 实际单位 | 依据 |
|---|---|---|---|---|
| `tcgay预后/time.csv` | 730 | 4502 | day | TCGA 原始天数 |
| `tcgay预后/risk.txt` | 1.87 | 12.3 | year | = 天数/365（612/365=1.6767，逐位吻合） |
| `tcgay预后/re.csv` | 同上 | | year | 与 risk.txt 同队列 |

脚本的推断规则（`detect_time_unit()`，按最大值分档）：

| max(time) | 推断 |
|---|---|
| ≤ 25 | year |
| 25 ~ 400 | month |
| > 400 | day |

这是**启发式**，在"月 vs 年"区间会撞车。所以脚本每次都会打印推断结果，
并在决策文件里列为第 1 条必须核对项。**拿不准就问用户，不要让用户猜。**

换算一律经"天"中转：

```
day = 1      month = 30.44      year = 365
```

- `30.44` 与 `最佳阈值-km.R` 的 `/30.44` 一致
- `365/30.44 = 11.9908`，与 `生存分析-km.R` 的 `/12` 差 0.08%
- 要**逐位复刻**参考代码的除法，用 `--time_div=<数值>`：它直接 `t/time_div`，
  绕过上面的换算表（此时横轴标签仍由 `--time_out=` 决定）

横轴刻度间隔默认随输出单位走：年→2（同参考 1）、月→24（同参考 2）、天→365。
用 `--break_time=` 覆盖。

---

## 7. 样本名与样本筛选（TCGA 情形）

参考代码的顺序**不能改**：

1. **先**按完整条码的第 4 段首字符筛肿瘤样本
   `TCGA-90-A4EE-01A-11R-A24Z-07` → 第 4 段 `01A` → 首位 `0` = 原发肿瘤
2. 再转置
3. **最后**把样本名截成前 12 位 participant ID，并把 `.` 换成 `-`

第 3 步会把第 4 段抹掉，所以顺序反了筛选必然失效（脚本实测日志里
`[filter] 仅保留 sample type 首字符 ∈ {0}：701 -> 650` 就是这个环节）。

样本身份码含义：`0` 原发瘤（01-09）、`1` 正常（10-19）、`2` 复发、`6` 转移。
`--tumor_codes=0,2,6` 可放宽。**认不出身份码的样本（非 TCGA 命名、或已被
截断的名字）一律保留**，不做判定 —— 否则会把 GEO 数据整批丢掉。

规范化只对**以 `TCGA-` 开头的名字**生效（`--id_norm=auto`）。GSM 号
（`GSM1234567`）与自定义命名原样保留 —— 无条件 `substr(1,12)` 会截断
`Sample_001_rep1` 这类名字，导致交集为空。

### 一个必须留意的副作用

截到 12 位后，同一患者的多个 aliquot / 多个样本会**变成重复 ID**。
实测 `TCGA_CRC_TPM.txt`：650 列里规范化后有 **26 个重复 ID**。
脚本会打印

```
[merge] 注意: 规范化后重复 ID — 表达矩阵 26 个、临床表 0 个，重复者只取首次出现
```

处理方式是**取首次出现**（与参考代码 `data[common,]` 的行为一致）。
若用户在意"同一患者多样本"，需要先向用户确认去重策略再跑。
