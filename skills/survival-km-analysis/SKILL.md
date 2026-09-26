---
name: survival-km-analysis
description: 肿瘤/临床 Kaplan-Meier 生存分析的「两级策略」流程。第一级用常规分位值（默认中位值）分组做 log-rank 检验；当 p 值未达阈值（默认 0.05）时，自动改用 survminer::surv_cutpoint 的**最佳截断值**重新分组，并取 p 更小者作为最终方案；输出 KM 曲线（含风险表、中位生存线、HR 与 95% CI）、逐样本分组表与决策记录，另可批量跑整张基因列表。覆盖 TCGA 表达矩阵取样、条码规范化到 12 位 participant ID、肿瘤/正常样本筛选、时间单位（天/月/年）换算、Cox 单因素 HR、maxstat 最佳阈值、survdiff、ggsurvplot。当用户提到 生存分析 / 生存曲线 / KM 曲线 / Kaplan-Meier / log-rank / 预后分析 / 高表达低表达分组 / 最佳截断值 / 最佳阈值 / surv_cutpoint / maxstat / HR / 风险表 / time_LUSC / riskScore 分组 / 批量做生存 时使用。
agent_created: true
---

# KM 生存分析：常规分组 → 最佳阈值回退

## 这个技能要解决的核心问题

常规做法（按表达量中位值分 High/Low）经常做不出 p<0.05。本技能把"退而求其次"
变成一个**自动的、可追溯的**两步：

| 级 | 做法 | 等价于 |
|---|---|---|
| 第一级 | 按分位值切分（`--pct=0.5` 即中位值） | `<参考脚本>/生存分析-km.R` |
| 第二级 | 第一级 p ≥ `--p_threshold`（默认 0.05）时，用 `survminer::surv_cutpoint` 找**最佳截断值**重新分组 | `<参考脚本>/最佳阈值-km.R` |

最终**取两级里 p 更小者**，并把两级的 p 值、阈值、分组样本数**同时写进决策文件**，
所以"到底有没有回退、为什么回退"永远可查。

> ⚠️ 第二级是**数据驱动**的阈值，属于探索性分析。技能除了照做，还会**额外算一个
> maxstat 校正 p**，因为图上的 log-rank p 是"切完之后再检验一次"的朴素 p，
> 系统性偏乐观。
>
> **批量场景下这个偏差会被放大**（`re.csv` 39 基因实测）：
>
> | | 数量 |
> |---|---|
> | 回退到最佳阈值的基因 | **32 / 39** |
> | 朴素 log-rank p < 0.05 | **23** |
> | **校正 p 后仍 < 0.05** | **1**（CCNB1，0.0108） |
>
> 单基因版的例子：SOX9 常规 p=0.43 → 最佳阈值切分后朴素 p=**0.0318** →
> 校正 p=**0.2998**（不显著）。
>
> **交付时必须把校正 p 一并给用户**，见 `references/methodology.md` 第 3 节。

## 何时使用

- 给一个基因 + 表达矩阵 + 临床表，要出 KM 曲线
- 给一个风险评分（riskScore / 模型输出）+ 生存表，要出 KM 曲线
- 有一批基因，要批量出 KM 并汇总哪些显著
- 用户提到 surv_cutpoint / 最佳阈值 / maxstat / log-rank / HR / 风险表

**只做单因素 KM 与 Cox 单因素**。多因素 Cox、列线图、ROC/时间依赖 AUC、
C-index、LASSO 建模都不在本技能范围内。

## 关于"分组方式由谁定"

| 环节 | 谁决定 |
|---|---|
| 用哪个基因 / 哪个评分列 | **用户必须指定**（`--gene=` 或 `--score_col=`） |
| 时间列、状态列、时间单位 | 脚本自动探测 + **必须让用户核对一次** |
| 常规 vs 最佳阈值 | 用户可用 `--mode=` 强制；默认 `auto` 自动回退（用户要的就是这个） |

时间单位是**最容易出错**的一项：脚本会推断，但推断只是启发式，
见下文"时间单位"一节的核对要求。**不要跳过这一步。**

## 前置检查

1. **定位 R**：见上文——`run_surv.sh` 自动探测 Rscript。

   ```bash
   RSCRIPT="/path/to/Rscript"
   ```

2. **必须用包装器 `run_surv.sh` 启动** —— 它自动设好临时目录、COMSPEC 与 locale：

   ```bash
   SK="<TOOLKIT>/skills/survival-km-analysis/scripts"
   bash "$SK/run_surv.sh" "$SK/01_km_survival.R" --help
   ```

   直接调 `Rscript` 会踩两个坑：
   - **中文用户名把 `TMP/TEMP` 搞坏** → R 的 `tempdir()` 不可写，崩在临时目录上
   - **环境里带着 `LC_ALL=C.UTF-8`**，Windows 版 R 应用不了、退化到 `C` locale，
     此时**完全无法寻址含非 ASCII 字符的路径**（`file.exists` 静默为假）。
     `run_surv.sh` 已 `unset LC_ALL LANG ...`，这一步不能省。

   手动跑等价于：

   ```bash
   mkdir -p /tmp/rtmp   # Windows: mkdir C:/Rtmp
   export TMP=/tmp/rtmp TEMP=/tmp/rtmp TMPDIR=/tmp/rtmp   # Windows: 指向纯 ASCII 目录如 C:/Rtmp
   export COMSPEC='C:\WINDOWS\system32\cmd.exe'
   unset LC_ALL LANG LC_CTYPE LC_COLLATE LC_MONETARY LC_TIME
   ```

3. **依赖包**（安装命令见套件根目录 INSTALL.md）：`survival`、`survminer`（内部
   依赖 `maxstat` 做最佳截断值）、`ggplot2`、`ggpubr`、`ragg`（PNG 渲染）。
   缺 `svglite` 时矢量图一律输出 PDF。

> **已知无害现象**：R 退出时报 `Segmentation fault`（退出码 139），发生在所有
> 工作完成之后。**判断成败看产物文件，不要看退出码。**

## 两种输入形态

### A) 表达矩阵 + 临床表

```bash
SK="<TOOLKIT>/skills/survival-km-analysis/scripts"
bash "$SK/run_surv.sh" "$SK/01_km_survival.R" \
  --expr="D:/data/TCGA_CRC_TPM.txt" \
  --clin="D:/data/time.csv" \
  --gene=TIMP1 --outdir=./km --prefix=TIMP1
```

脚本按参考代码的顺序处理：**先按表头第 4 段首字符筛肿瘤样本，再转置，
再截断成 12 位 participant ID，最后取交集**。顺序不能反 —— 截断到 12 位后
第 4 段就没了，筛选必然失效。

### B) 单表（时间/状态/评分在同一张表）

```bash
bash "$SK/run_surv.sh" "$SK/01_km_survival.R" \
  --table="D:/data/H19.csv" --score_col=risk_score \
  --time_col=OS.time --event_col=OS --outdir=./km --prefix=H19
```

`risk.txt`（含 `time/state/<基因...>/riskScore/risk`）这类 TCGA 建模产物
走这条路最省事：`--table=risk.txt --score_col=riskScore --time_col=time --event_col=state`。

## 时间单位：必须核对

脚本用 `detect_time_unit()` 推断输入单位，规则是（按最大值分档）：

| 最大值 | 推断为 |
|---|---|
| ≤ 25 | year |
| 25 ~ 400 | month |
| > 400 | day |

实测：`time.csv` 中位数 730、最大 4502 → `day` ✓（TCGA 经典天数）；
`risk.txt` 中位数 1.87、最大 12.3 → `year` ✓（CTRP 里已是 /365 的年）。

**推断会在三种口径间撞车**（一个月的数据最大值 36，一年的数据最大值 12，
都可能误判）。所以每次运行的日志与决策文件都会打印：

```
[time] 自动推断输入单位 = X（中位数 a，最大 b）
[time] !! 单位推断只是启发式，请核对：肿瘤 OS 通常 1~5 年 = 12~60 月 = 365~1825 天。
```

**交付时必须把这一行指给用户确认。** 拿不准就让用户直接写 `--time_in=`。
用 `--time_out=` 选分析/画图单位（默认 `year`，与 `生存分析-km.R` 的
`Time(years)` 一致）。

换算经"天"中转：`day=1 / month=30.44 / year=365`。其中 `30.44` 与
`最佳阈值-km.R` 的 `/30.44` 一致；`365/30.44 = 11.99` 与 `生存分析-km.R`
的 `/12` 差 0.08%。要**逐位复刻参考代码的除法**（例如它除的是 12 而不是
11.99），用 `--time_div=12` 直接相除，此时按 `--time_out=` 定的横轴标签不变。

## 分组与回退

```bash
# 默认：中位值 -> 不足则最佳阈值
--mode=auto        # 默认

# 强制只看常规分组（复刻 生存分析-km.R）
--mode=median      # pct=0.5
--mode=quantile --pct=0.3333      # 三分位下界，对应参考代码注释里那行

# 强制只看最佳阈值（复刻 最佳阈值-km.R）
--mode=cutpoint --minprop=0.3
```

其它常用参数：

| 参数 | 说明 |
|---|---|
| `--p_threshold=0.05` | 第一级 p 低于它就不回退 |
| `--minprop=0.3` | `surv_cutpoint` 每组最少样本比例 |
| `--maxstat_pmethod=Lau92` | 算校正 p 用的 maxstat 近似法。**只有 Lau92 在这批数据上给出合法取值**（Lau94 会返回 1.156 这种越界值）；越界一律置 NA 并说明 |
| `--style=unified\|classic\|threshold` | 绘图风格，默认 `unified`（统一版式）。见下表 |
| `--outdir=./km` `--prefix=GENE_KM` | 产物目录与前缀 |
| `--pdf=TRUE --png=TRUE --dpi=300` | 输出格式 |
| `--tumor_only=TRUE --tumor_codes=0` | 只留原发肿瘤（0=原发瘤，1=正常，2=复发，6=转移） |
| `--id_norm=auto\|tcga\|none` | 是否把 TCGA 条码截到 12 位。**只对以 `TCGA-` 开头的名字生效**，GSM/自定义命名原样保留 |
| `--min_time=0` | 剔除时间 ≤ 该值的样本 |
| `--event_levels=` | 状态列不是 0/1 时用（如 1/2 编码写 `--event_levels=2`） |
| `--time_col=` `--event_col=` | 显式指定列名，默认自动探测 |

### 出图风格

**风格是纯出图偏好，与分组方式无关。** 同一批里混着"走了最佳阈值回退"和
"没走回退"的基因，出的图版式必须完全一样 —— 看图不该看出它用了哪种分组。
默认的 `--style=unified` 就是为此而设；另两个选项只是用来**严格复刻**你给的两份
参考代码，不是日常产出该选的。

| | `--style=unified`（默认） | `--style=classic` | `--style=threshold` |
|---|---|---|---|
| 用途 | **统一版式，日常都用它** | 严格复刻 `生存分析-km.R` | 严格复刻 `最佳阈值-km.R` |
| 配色 | Low `MediumSeaGreen` / High `Firebrick3` | 同左 | Low `#3090a1` / High `#bc5148` |
| 图例 | 标题 `<基因> expression`，标签 `Low`/`High` | 同左 | 无标题，标签 `<基因>_low`/`<基因>_high` |
| 风险表标签 | **跟曲线配色**（`risk.table.col="strata"`） | 默认黑色 | 跟曲线配色 |
| 置信带 | 关（`--conf_int=TRUE` 打开） | 关 | **开** |
| 主题 | survminer 默认 | 同左 | `theme_minimal(base_size=14)`，无网格、黑轴线 |
| p/HR 标注 | 面板内 `p = …` 换行 `HR = …` | 同左 | 面板内左下、左对齐的 `HR = …` 换行 `p = …` |
| 图例位置 | (0.8, 0.8) | 同左 | (0.85, 0.85) |
| 横轴刻度 | 年→2 / 月→24 / 天→365（自动） | 同左 | 同左 |
| 中位生存线 | `hv` | 同左 | 同左 |

`unified` 把两份参考代码里**互相冲突**的外观项定死了，不必每次去挑：配色 / 图例
标题 / 主题 / 标注机制 / 置信带取自主参考 `生存分析-km.R`；风险表标签跟曲线配色
取自 `最佳阈值-km.R`；y 轴标题显式写死，免得随 survminer 版本漂移。
`unified` 与 `classic` 的唯一外观差别就是风险表标签配色。

**改风格只改外观，不改统计口径**：HR 恒为 **High vs Low**，与风格无关。
配色也可用 `--palette=色1,色2`（Low,High）单独覆盖。

> 📌 想换配色但保持版式统一，请用 `--palette=` 而不是切到 `threshold` 风格 ——
> 后者会连置信带、主题、图例形式一起换掉，那就又不一致了。
> 真要复刻参考代码时，建议 `--style=` 只在单独重现某张已发表图时使用。

## 输出

单基因（`01_km_survival.R`）：

| 文件 | 用途 |
|---|---|
| `<prefix>_KM.pdf` / `.png` | 生存曲线（PNG 300 dpi） |
| `<prefix>_stats.txt` | 全部统计量，key=value，便于程序读 |
| **`<prefix>_km_decision.txt`** | **两级 p 值并列 + 必须核对清单，交付时先给用户看这个** |
| **`<prefix>_group_map.csv`** | **逐样本 取值/分组/时间/事件，供逐样本核验** |

批量（`02_km_batch.R`）：`<prefix>_km_summary.csv`（每基因一行，含
`p_text`/`hr_text`/`mode_used`/`fallback_used`，**显著者排前、按 p 升序**）、
`<prefix>_km_multi.pdf`（多页）、`per_gene/<GENE>_KM.pdf|png` 与
`per_gene/<GENE>_group_map.csv`。

```bash
bash "$SK/run_surv.sh" "$SK/02_km_batch.R" \
  --expr="D:/data/TCGA_CRC_TPM.txt" --clin="D:/data/time.csv" \
  --genes=TIMP1,GSTP1,SOX9 --outdir=./km --prefix=CRC

# 基因列表也可来自文件；文件是"每行一个"或"索引+基因"两列都支持
bash "$SK/run_surv.sh" "$SK/02_km_batch.R" \
  --table="D:/data/re.csv" --genes_file="D:/data/merge2.csv" \
  --time_col=time --event_col=state --outdir=./km --prefix=CRC39
```

> 批量耗时主要是 `ggsurvplot` 渲染与每次 `surv_cutpoint` 的搜索：
> 实测 5 个基因约 40 秒，**39 个基因 4 分 09 秒**（`--png=FALSE`）。
> **长批次一定要后台运行或把输出重定向到日志文件再读**，
> 用默认超时硬等会被打断（实测被 SIGTERM 杀掉后日志可能为空）。

## 交付前必须核对的三件事

1. **时间单位**：日志里 `[time] 自动推断输入单位 = ?` 是否与数据实际一致
   （看中位数与最大值两个数）；不符就用 `--time_in=` 重跑。
2. **分组样本数**：`Low=? / High=?`。用 `>` 切分时，表达量含大量并列值
   （例如 TPM 里大量 0）会把并列样本全推给 Low，出现 90%/10% 这种极端比例 ——
   脚本会打印 `[警告: 分组严重不平衡]`，此时中位值分组本身就不成立，
   应改 `--pct` 或直接用 `--mode=cutpoint`。
3. **有没有真的回退**：看 `<prefix>_km_decision.txt` 的第二级那行。
   - 若回退成功 → 图上的 p 是朴素 log-rank p，**必须同时报 `maxstat 校正 p`**
     （决策文件、`_stats.txt` 的 `p_maxstat_adjusted`、批量汇总表的
     `p_maxstat_adjusted` 列都有）
   - 若校正后不显著 → **就说它不显著**。不要为了凑 p<0.05 反复调
     `--minprop` / `--pct`，那等于把择优偏倚叠了两层
   - 若第二级也失败（`No cutpoint of sufficient size`）→ 脚本已自动退回
     第一级结果，把失败原因转告用户，别硬造阈值
4. **图看齐不看齐**：一批图里**不要**因为"这个基因走了回退"就把它的风格换掉。
   整批统一用一个 `--style`（默认 `unified` 即可），要换配色用 `--palette=`。
   交付前可肉眼扫一遍整批图，版式（配色/图例形式/风险表/标注位置/字号）应当一致。

## 关键决策点（不要替用户默默决定）

1. **用哪个基因/评分列** —— 必选项，没有默认。
2. **时间列、状态列** —— 自动探测会打印结果，让用户扫一眼；
   列名歧义（有多个候选）时必须问。
3. **时间单位** —— 见上，必核对。
4. **分组策略** —— 默认 `auto` 回退是用户明确要的；但若用户想复刻某张
   已发表的图，要问清是常规分组还是最佳阈值，再用 `--mode=` 固定。
5. **回退的代价** —— 只要用了最佳阈值，交付话术里必须带上 maxstat 校正 p
   （以及"探索性/需独立验证"），不能只报那个好看的朴素 p。
6. **出图风格** —— 默认 `unified` 就别动。**不要**因为"这个基因走了回退"就把它
   的风格或配色换掉；要换配色用 `--palette=`。`classic`/`threshold` 只用于
   单独复刻某张参考图。

## 资源

- `scripts/run_surv.sh` —— **首选入口**，包装环境变量后调用 R
- `scripts/01_km_survival.R` —— 单基因/评分：两级策略 + 出图 + 报告
- `scripts/02_km_batch.R` —— 批量基因：逐基因走同一套两级策略 + 汇总表
- `scripts/lib_surv_common.R` —— 公共库：稳健读表、样本名规范化、时间单位、
  `two_tier_km()`、`km_stats()`、`make_km_plot()`。**改脚本前先读它**
- `references/methodology.md` —— 两级策略的统计依据、`surv_cutpoint` 的返回结构、
  **最佳阈值的择优偏倚与校正 p**、时间单位约定、TCGA 样本名与筛选
- `references/pitfalls.md` —— 踩坑手册 + **回归基线表（三组实测数字，
  对不上说明中间出了错）**

两个脚本都支持 `--help`。遇到报错先查 `references/pitfalls.md`。

> 📌 **改动本技能的脚本时，图内文字一律用 ASCII。** Windows 下中文字形缺失
> 会被渲染成乱码 ASCII（不是方框），很难一眼看出来。`ascii_guard()` 会提醒。
