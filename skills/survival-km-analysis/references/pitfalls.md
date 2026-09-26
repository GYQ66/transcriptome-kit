# 踩坑手册（survival-km-analysis）

改脚本前先读这份。条目都带日期与复现现象，不是猜测。

---

## 1. 环境 / R 启动

### 1.1 必须走 `run_surv.sh`

直接 `Rscript script.R` 会踩两个坑：

- **中文用户名把 `TMP/TEMP` 搞坏** → R 启动时据此建的 `tempdir()` 是不可写路径，
  症状是 `Failed to create directory ... C:/Users/i+i:`。
  `tempdir()` 在 R 启动时就固定，**脚本内部改不掉**。
- **`LC_ALL=C.UTF-8` 让中文路径全部失效**。Windows 版 R 应用不了 `C.UTF-8`
  （只打印一行 `Setting LC_CTYPE=C.UTF-8 failed` 就继续），退化成 `C` locale。
  在 `C` locale 下 R **完全无法寻址含非 ASCII 字符的路径**：
  `dir.exists` / `file.exists` / `readLines` 一律失败且不报错，
  `list.files()` 静默返回空 —— 看起来就像"文件不存在"。

  `run_surv.sh` 已 `unset LC_ALL LANG LC_CTYPE LC_COLLATE LC_MONETARY LC_TIME`
  并设好 `TMP/TEMP/TMPDIR/COMSPEC`，所以**一律走包装器**。

### 1.2 R 退出码 139 / Segmentation fault 是无害的

加载过 `survminer` 后，R 在**完成全部工作之后**退出时会报
`Segmentation fault`（139），伴随一行 `null device`。
**判断成败看产物文件，不要看退出码。** 实测多次均有此现象，产物完整。

### 1.3 脚本路径含中文时包装器会复制脚本

`run_surv.sh` 检测到脚本路径含非 ASCII 就复制到 `D:/Rtmp/`，
于是 `SELF_DIR` 变成临时目录、找不到同目录的 `lib_surv_common.R`。
**脚本里用 `find_surv_lib()` 定位**（依次尝试脚本目录 → 环境变量
`SURV_SKILL_DIR`（包装器导出）→ 当前目录 → `./scripts`）。改动此处要小心。

> 临时脚本名必须带 PID 与本脚本名（`_surv_run_$$_<basename>`）。
> 早期写法用固定名 `_surv_run_tmp.R`，**并行跑两个脚本时互相覆盖**，
> 症状是"调用 01 却执行了 02"。2026-09-18 实测踩到。

### 1.4 Bash 工具被沙箱拦

本机 Bash 工具曾出现 `bash 脚本.sh` 被路由到被黑名单的 `wsl.exe`
（`PROGRAM BLOCKED BY SECURITY POLICY`）。`run_surv.sh` 目前可正常调用。
若被拦，备用通道是 PowerShell 工具直接调 Rscript，在**同一条命令**里先设环境变量：

```powershell
$env:TMP="D:\Rtmp"; $env:TEMP="D:\Rtmp"; $env:TMPDIR="D:\Rtmp"; $env:COMSPEC="C:\WINDOWS\system32\cmd.exe"
$log = & "C://path//to//Rscript.exe" --vanilla <脚本> <参数...> 2>&1
[IO.File]::WriteAllText("<日志路径>", (($log | ForEach-Object { "$_" }) -join "`r`n"), [Text.Encoding]::UTF8)
```

PowerShell 里 `unset LC_ALL` 的等价写法：`Remove-Item Env:LC_ALL -ErrorAction SilentlyContinue`
（`LANG`/`LC_CTYPE` 等同理；不清理则中文路径再次失效）。

---

## 2. 读表：三种真实会遇到的坑

### 2.1 `write.csv` 产物带字面引号 → 按列名取值全失败

`read.table(quote = "")` 关闭引号解析（避免基因名/ID 里的引号导致行截断），
代价是**表头与行名的引号会原样保留**：列名变成 `"time"` 而不是 `time`。

实测（2026-09-18）：读 `re.csv` 报

```
Error: 时间列 'time' 不在表中。可用列: "time", "state", "BUB1", ...
```

注意可用列是**带引号**打印的 —— 这就是判据。
修法：`strip_quotes()` 在读表后统一剥掉首尾成对引号（列名、行名、ID 都做）。

### 2.2 表头首字段为空（`,time,state`）→ `invalid 'row.names' length`

`time.csv` 的表头是 `,time,state`（首字段为空），读出来列名是
`c("", "time", "state")`。此时 `df[[""]]` 返回 `NULL`，
`length(ids)=0 != nrow(out)`，报
`Error in .rowNamesDF<- : invalid 'row.names' length`。

修法：`read_clin_auto()` **全程用位置索引**（`raw[[idx]]`、`-idx`），不用列名。
2026-09-18 实测，这是 `01_km_survival.R` 在表达矩阵路径上崩掉的直接原因。

### 2.3 基因列表文件不是"每行一个基因"

`merge2.csv` 的 `read.csv(row.names=1)` 产物是 `"1","CCNB1"` 这种
**索引 + 基因**两列结构。按行读会得到 `1","CCNB1` 这种垃圾 token。

修法：`read_gene_list()` 逐行按 `,`/TAB/`;` 切栏 → 剥引号 → 取**第一个像基因符号的
字段**（`^[A-Za-z][A-Za-z0-9_.-]*$` 且长度≥2）→ 去掉表头词。
"每行一个"的朴素列表也能正常处理。

### 2.4 分隔符判别

按第一行里 TAB 与逗号的个数取多者。第一行只有 1 个字段时（单列文件）会退化，
但那种文件本技能用不上。

---

## 3. 绘图（ggplot2 4.0 / survminer 0.5.2）

### 3.1 图例位置：写进 `ggtheme` 无效

`ggsurvplot` 的默认 `legend = "top"` 会**覆盖** `ggtheme` 里的
`legend.position`。所以 `--style=threshold` 里那句"图例放在面板内
(0.85, 0.85)"**必须通过 `legend=` 参数传**，只写在 ggtheme 里图例会跑到图外顶部。

顺带记录（ggplot2 4.0.1 实测）：`theme(legend.position = c(0.85, 0.85))`
仍可用 —— 内部会被转成 `legend.position = "inside"` +
`legend.position.inside = c(0.85, 0.85)`。可用 `legend.position.inside` 自行验证。

### 3.2 `Ignoring unknown labels: • colour : "XXX expression"`

构建（`ggsurvplot()`）与渲染（`print()`）两个阶段各出现一次，
ggplot2 4.0 用 cli 输出，既可能走 warning 也可能走 message。
**图例标题其实正常渲染**，纯噪音。`quiet_gg()` 用 `withCallingHandlers`
同时拦 warning 与 message 压掉它（两个都要拦，只拦 warning 会漏）。

### 3.3 风险表被图例标题塞进一个多余的 y 轴标签

`legend.title = "<基因> expression"` 会被 ggsurvplot 顺带写进风险表的 y 轴标题，
在图的左下角竖着排一行多余文字，与风险表行名重叠。

修法（`make_km_plot()` 末尾）：

```r
pl$table <- pl$table + ggplot2::ylab(NULL)
```

### 3.5 `--help` 的长字符串里不能出现英文双引号（踩了两次）

`cat("...")` 这个多行帮助文本是一个普通 R 字符串字面量，**里面出现 `"` 会直接
把字符串截断**，报 `unexpected symbol in: "..."` —— 而且**只有加 `--help`
时才会暴露**，正常跑参数解析根本走不到那里，极易漏掉。

同一天踩了两次：一次是决策文件里写 `"log-rank p 值"`，一次是 02 的 `--help`
里写 `所以"走了最佳阈值回退"和"没走回退"的…`。

**规矩**：长中文文本里要引用/强调，用 `【】`、`「」` 或中文引号，别用英文 `"`。
改完帮助文本或长字符串后，**务必跑一次 `--help`** 确认没被截断：

```bash
bash "$SK/run_surv.sh" "$SK/02_km_batch.R" --help
```

（便宜又管用的兜底：`Rscript -e 'parse("脚本路径")'` 只做语法分析，
能秒查出这类问题。）

### 3.6 别让"走了哪种分组"影响出图风格

**这是返工率最高的一个坑。** 参考代码是两份、版式不同（`生存分析-km.R` 一套，
`最佳阈值-km.R` 另一套），所以很容易顺手写成"中位值分组就用 classic、
最佳阈值分组就用 threshold" —— 结果同一批图里出现两种版式，看图的人一眼就能
认出哪张是回退来的。

正确做法：**风格只由 `--style` 决定，整批共用一个值**。默认 `unified`。

- 代码约束：`make_km_plot(fit, d, gene, style = ...)` 里 `style` 是唯一的外观开关；
  `mode_used` 不得参与任何绘图分支。改动绘图代码时**别把 `mode_used` 引进来**。
- 核验方法见 `methodology.md` 第 5.1 节（比对 `plot$theme` / `plot$labels` /
  `table$theme` / 图层构成 / colour scale 调色板，应逐项 identical）。
- 想换配色用 `--palette=色1,色2`，**不要**切到 `threshold` 风格 —— 那会连
  置信带、主题、图例形式、标注位置一起换掉。

### 3.4 多页 PDF 与并行

- 多页 PDF 用 `grDevices::pdf(..., onefile = TRUE)` + 循环 `print_km()`。
- 每张图 PDF 约 9 KB、PNG（300 dpi, 8×6.25 in）约 150 KB。
- **耗时**：单基因含回退约 8~15 秒；5 个基因约 40 秒；
  39 个基因超过 2 分钟 → **长批次用后台运行或输出重定向到日志文件**，
  别用默认超时硬等（实测被 SIGTERM 打断后日志可能为空）。
- 批量时 `--png=FALSE` 可省掉一半渲染时间。

---

## 4. 统计

### 4.1 HR 方向

`km_stats()` 恒用 `factor(group, levels = c("Low","High"))`，
所以 `coxph` 的系数就是 **High vs Low**，与画图风格无关。
若将来改成别的水平顺序，函数里有取倒数的兜底分支，但**别依赖它**——
新增水平顺序一定同步改 `hr_low`/`hr_hi` 的互换逻辑。

### 4.2 分组不平衡

`quantile_split()` 用 `>`（忠于参考代码），并列值全归 Low。
TPM 里大量 0 的基因会出现 90%/10%。脚本在任一组 <10% 时告警。
`km_stats()` 里还会因某组事件数为 0 导致 `coxph` 不收敛 → HR 置 NA（已 tryCatch）。

### 4.3 `maxstat` 校正 p 的三个坑

1. `maxstat::pvalue()` **不存在**（未导出），别用。
2. `cp[[<变量名>]]$p.value` **恒为 NA**（`surv_cutpoint` 用 `pmethod="none"`）。
3. `pmethod="Lau94"` 在 `re.csv` 上返回 **1.156**（越界）。默认改用 `Lau92`；
   结果越界/非有限一律置 NA，并在日志里说明。

详见 `methodology.md` 第 3 节。

### 4.4 `sprintf()` 遇到 NULL 参数会返回**零长度字符向量**，`c()` 会静默丢掉它

`km_stats()` 返回的字段是 `n_per`（一个 `table`），**没有** `n_low` / `n_high`
（那两个是 `best_cutpoint()` 的字段）。早期决策文件里写成
`sprintf("...Low=%d High=%d...", cp$cutoff, minprop, st_c$n_low, st_c$n_high, ...)`，
`st_c$n_low` 是 NULL → `sprintf` 返回 `character(0)` → `c()` 把它当成"没有元素"
直接丢掉，**整行"第二级 最佳截断值"从决策文件里消失**，不报错、不留痕。

症状极具迷惑性：`--mode=auto` 且未回退时（`st_c` 为 NULL，走另一分支）一切正常，
只有真正回退/强制 cutpoint 时那一行才不见。

修法：字段一律从确实携带它的对象取（`cp$n_low` / `cp$n_high`）。
**通用教训**：任何 `sprintf()` 的实参只要可能为 `NULL`，先确认它不会让结果变成长度 0。

---

## 5. 回归基线（数字对不上说明中间出了错）

> 📌 **2026-09-18 起默认风格由 `classic` 改为 `unified`。**
> 这只是出图版式的默认值变化（风险表标签改成跟曲线配色），
> **下面所有统计数字都不受影响** —— 分组、p、HR、阈值逐位不变。
> 要复现更早那张"黑色风险表标签"的图，加 `--style=classic`。

### 5.1 `risk.txt`（单表路径，582 例 / 119 事件，时间已是"年"）

```bash
bash run_surv.sh 01_km_survival.R --table="D:/phlp/tcgay预后/risk.txt" \
  --score_col=riskScore --time_col=time --event_col=state --prefix=riskScore
```

| 项 | 值 |
|---|---|
| 自动推断时间单位 | `year`（中位数 1.87，最大 12.3） |
| 中位值分组阈值 | 1.04198426547796 |
| 分组 | Low=291 / High=291 |
| log-rank p | 0.0014946 |
| HR (High vs Low) | 1.82（95% CI 1.25–2.64） |
| Cox p | 0.0017571 |
| 中位生存期 | Low 未达到 / High 4.77 年 |
| 决策 | `mode_used=median`，未回退 |

### 5.2 `re.csv` 5 个基因（单表 + `--genes_file=merge2.csv`）

| 基因 | 常规 p | 采用 | 阈值 | 切分后朴素 p | maxstat 校正 p | HR (H vs L) |
|---|---|---|---|---|---|---|
| TIMP1 | 0.0312 | median | 601.019 | — | — | 1.49 (1.03–2.14) |
| CDC25C | 0.00324 | median | 9.6507 | — | — | 0.578 (0.399–0.836) |
| SOX9 | 0.4304 | **cutpoint** | 153.886 | 0.0318 | **0.2998** | 0.670 (0.463–0.968) |
| GSTP1 | 0.2130 | **cutpoint** | 1377.042 | 0.0680 | — | 0.684 (0.454–1.031) |
| COMT | 0.6453 | **cutpoint** | 43.411 | 0.1570 | 0.7393 | 0.765 (0.527–1.110) |

- SOX9 的切分后分组：Low=192 / High=390（`minprop=0.3` 约束下 33%/67%）
- 汇总 CSV 顺序：显著者在前，组内按 `p_used` 升序

### 5.2b 全量 39 基因（`--genes_file=merge2.csv`，实测 4 分 09 秒）

```
共 39 个基因：显著 23，未达显著 16，失败 0
其中 32 个基因因常规分组不显著而改用最佳截断值
```

拆开看：

| 项 | 数量 |
|---|---|
| 采用 median / cutpoint | 7 / 32 |
| 朴素 p<0.05 | 23（median 7 + cutpoint 16） |
| cutpoint 来源且校正 p 仍 <0.05 | **1（CCNB1，0.0108）** |

→ 这是"最佳阈值把 p 做出来"的放大效应实证，详见 `methodology.md` 第 3.2b 节。
**耗时**：39 基因 4 分 09 秒（`--png=FALSE`）。含 PNG 会明显更久，
长批次务必后台跑。

### 5.3 `TCGA_CRC_TPM.txt` + `time.csv`（表达矩阵路径）

```bash
bash run_surv.sh 01_km_survival.R --expr="D:/phlp/tcgay预后/TCGA_CRC_TPM.txt" \
  --clin="D:/phlp/tcgay预后/time.csv" --gene=TIMP1 \
  --time_out=month --prefix=expTIMP1
```

| 项 | 值 |
|---|---|
| 表达矩阵 | 19937 行 × 701 列 |
| 肿瘤样本筛选 | 701 → 650（`0`=650，`1`=51） |
| 与临床表交集 | **557**（规范化后表达矩阵有 26 个重复 ID） |
| 事件数 | 116 |
| 自动推断时间单位 | `day`（中位数 730，最大 4502） |
| 输出单位 | month（中位数 24，最大 148） |
| TIMP1 中位值分组 | 阈值 603.312，Low=279 / High=278 |
| log-rank p | 0.04994（刚好压线） |
| HR (High vs Low) | 1.44（95% CI 1.00–2.08） |

### 5.4 风格一致性（改绘图代码后必跑）

同一基因分别 `--mode=median` 与 `--mode=cutpoint` 出图，比对两张图的
`plot$theme` / `plot$labels` / `table$theme` / `table$labels` / 图层构成 /
图层类 / colour scale 调色板，**应逐项 identical**。实测（`re.csv` 的 SOX9，
`--style=unified`）：12 项指纹全部一致、配色同为
`MediumSeaGreen, Firebrick3, #6E568C, #223D6C`。

改 `make_km_plot()` 或任何绘图分支后，**把 `mode_used` 引进来就会破坏这条**，
务必重跑这个核验（思路见 `methodology.md` 第 5.1 节）。

---

## 6. 依赖与能力边界

- 本机已装：`survival` 3.8.3、`survminer` 0.5.2、`maxstat` 0.7.26、
  `ggplot2` 4.0.1、`ggpubr` 1.0.0、`ragg` 1.5.0、`rms` 8.1.1、`data.table` 1.17.8
- **未装**：`svglite`（故矢量图只用 PDF）、`timeROC`、`survivalROC`
  → 时间依赖 ROC / AUC 相关需求本技能不覆盖，不要硬接
- PNG 优先用 `ragg::agg_png`（缺则退回 `grDevices::png`）
