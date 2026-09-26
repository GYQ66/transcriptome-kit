# 坑与实测基线

## A. 参考代码「注意事项」8 条的实测复核

| # | 参考代码的说法 | 实测结论 |
|---|---|---|
| 1 | `CoREST/LSD1/MACROH2A2/H2AX` 不是有效 symbol，会被 `removeUnmappedRows=TRUE` 静默剔除，final 图节点数 < 13 | **在本机 STRING 11.5 上不成立**：这 4 个名字**全部能直接映射**，13/13 命中，别名表一次都没触发。但「静默剔除」这个风险是真的 —— 换成 STRING 12 或别的物种就可能踩到。本 skill 的处理：两阶段映射（原名先试，失败的才换别名）+ 未映射基因告警 + `*_unmapped_genes.csv` |
| 2 | 步骤 04 的 `genes` 与步骤 01 的 `gene_list` 重复定义，可能错位 | 确认。本 skill 只有一份基因表；`--logfc` 个数与基因数不一致**直接报错退出**（不静默错位） |
| 3 | logFC 是示例数据，不是原文献真值 | 确认。报告里也写着；复现真图必须换成自己的差异表达结果 |
| 4 | `score_threshold=400` 可能过滤掉部分边，导致孤立节点/不连通 | 确认。本 skill 直接报「连通分量 N / 孤立节点 M」，并可 `--add_nodes` 补邻居 |
| 5 | `get_interactions()` 可能返回空边表，`max()` 返回 `-Inf` 报错 | 确认。本 skill 在 0 边时要么只画节点出图（`--keep_isolated=1`），要么给出带原因的中文报错 |
| 6 | 首次调用需联网下载物种数据 | 确认，而且**本机这条最要命**：aliases 20.7 MB + links 72 MB，下行约 80 KB/s，且经常断 |
| 7 | `library(tidygraph)` 未被使用，且不是常规包 | 后半句不准确：`tidygraph` 是 `ggraph` 的硬依赖，CRAN 上就有（本机 1.3.1），装了不报错。前半句对：那份代码确实用不到它。本 skill 不依赖 tidygraph |
| 8 | `repel=TRUE` 依赖 ggrepel | 确认。本机已装 ggrepel 0.9.6；缺了脚本会告警 |

## B. 本机新踩到的坑

### B1. STRING REST API 返回 0–1 标度 → 阈值 400 把边全滤光（最坑）

`version11_5.string-db.org/api/tsv/network` 返回的 `score` 是 **0–1 小数**：

```
9606.ENSP00000262241  9606.ENSP00000434024  0.427
```

按官方文档的 0–1000 口径去卡 `>= 400`，**31 条边全被滤掉**，图上剩 13 个孤立点。
脚本现在会检测 `max(score) <= 1` 并 ×1000 归一，同时请求侧只传
`required_score = score_threshold/1000`（0–1 口径）。见 `methodology.md` 第 2 节。

### B2. 只有 `version11_5` 这个域名可达

| 域名 | 结果 |
|---|---|
| `string-db.org` | ❌ `SSL connect error` |
| `version11.string-db.org` | ❌ `SSL connect error` |
| `version12_0.string-db.org` | ❌ `SSL connect error` |
| `version11_5.string-db.org` | ✅ 可达 |

→ 本机**不要**改 `--version=12.0`，改了就连不上。要换版本先自测域名连通性。

### B3. Bioconductor 主站不可达，装不了 STRINGdb

```
Error : Bioconductor version cannot be validated; no internet connection?
Warning: URL 'https://bioconductor.org/config.yaml': status was 'SSL connect error'
```

绕法（本机已用此法装成 STRINGdb 2.22.0）：直接下 Windows 二进制，不走 BiocManager。

```r
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN"))
install.packages(c("png","sqldf","plyr","igraph","httr","RColorBrewer","gplots","hash","plotrix"))
install.packages("https://mirrors.westlake.edu.cn/bioconductor/packages/3.22/bioc/bin/windows/contrib/4.5/STRINGdb_2.22.0.zip",
                 repos = NULL, type = "binary")
```

（`install.packages(repos=<bin/windows/contrib 目录>)` 这条路走不通 —— R 会再拼一层
`src/contrib` 和 `bin/windows/contrib/4.5`，直接 404。）

### B4. `STRINGdb$new()` 的版本校验要连 `string-db.org`

```
[warn] STRINGdb$new() 失败：cannot open the connection to
       'https://string-db.org/api/tsv-no-header/version'
```

该域名在本机时通时断 → `--engine=stringdb` 不稳定。脚本已加 3 次重试 + 明确提示，
失败就自动回退 `api`。**日常就用 `--engine=api`。**

### B5. `downloadAbsentFile` 只看 `file.size > 0` → 残缺缓存静默污染结果

STRINGdb 内部：

```r
downloadAbsentFile <- function(urlStr, oD = tempdir()) {
  temp <- paste(oD, "/", fileName, sep = "")
  if (!file.exists(temp) || file.info(temp)$size == 0) download.file(urlStr, temp)
  ...
}
```

下载中断留下 1.8 MB（应为 20.7 MB）的残缺 `.gz` 之后，**下次运行会直接跳过下载并使用它**，
后果不是报错而是**映射结果悄悄变成垃圾**：实测 13 个基因只映射上 1 个，
而日志里只是轻描淡写一句 `map() 失败`。

本 skill 的对策：启用 STRINGdb 前用 HTTP `Content-Length` 比对本地文件长度，
不一致就**拒绝启用该引擎**并列出该删的文件：

```
[error] STRINGdb 数据缓存残缺，拒绝启用该引擎（残缺文件会让映射静默变成垃圾）：
  9606.protein.aliases.v11.5.txt.gz  本地 1801969 字节 / 应为 21726816 字节
```

另外脚本还加了「STRINGdb 引擎映射率 < 50% 就大声告警」的兜底。

### B6. `options(timeout=)` 被 `on.exit` 提前还原

STRINGdb 的下载发生在 `STRINGdb$new()` **返回之后**（在 `map()` / `load()` 里）。
所以在 `make_engine_stringdb()` 里写 `on.exit(options(timeout = old))` 是错的：
函数一返回超时就被还原成默认 60 s，20.7 MB 的文件必然下到一半断掉。
正确做法是在引擎选定后**全局**设置 `options(timeout = max(timeout, 1800))`，不要 `on.exit`。

### B7. R 4.5.2 不接受 `$__xxx` 双下划线成员名

```r
p <- list(`__help` = TRUE); p$__help
# Error: unexpected input in " p$_"
```

写这个 skill 时踩到的：参数默认值里用 `__help` 当哨兵键，直接 parse error（不是运行时错误，
是**整个脚本解析失败**）。已改名 `help_flag`。以后写 R 参数解析别用 `__` 前缀。

### B8. 其它环境项

- `run_string.sh` 必须先补 `PortableGit/.../usr/bin` 进 PATH：本机 Bash 的 coreutils 时有时无，
  否则 `dirname`/`cygpath` 一起报 not found。
- `layout="stress"` 需要 `graphlayouts`（ggraph 的 Suggests，不是 Imports）。
  缺了脚本会告警并退回 `fr`。
- 本机没装 `svglite` → SVG 走 `grDevices::svg`；PNG 走 `ragg::agg_png`（1.5.0 已装）。
- `package 'ggraph' was built under R version 4.5.3`、R 退出码 139（Segmentation fault，
  发生在产物写完之后）都是**无害**现象，判断成败看产物文件。

### B9. `--format` 静默降级 → 扩展名与内容不符的「未知格式」文件（2026-09-25 修复）

v1.0 的 `open_dev()` 对未知格式走 `else` 分支**静默按 PDF 渲染**，而文件名却拿原始
`--format` 值拼扩展名。两种翻车方式（均已实测复现）：

| 输入 | v1.0 产物 | 后果 |
|---|---|---|
| `--format=jpg` | `*_ppi_network.jpg`，文件头 `%PDF-1.4` | 图片查看器打不开，报「未知格式」 |
| `--format=pdf,png` | `*_ppi_network.pdf,png` | 扩展名整个不认识 |

修复（v1.1，v1.2 扩展）：`normalize_formats()` 白名单过滤（pdf/png/svg）+ 非法值告警忽略 +
**pdf/png/svg 三格式强制同时输出**（v1.2 起，含 SVG）；`open_dev()` 对白名单外的值直接
`stop()`；`save_plot_formats()` 逐格式落盘并登记字节数，`print()` 出错也保证关闭设备
（否则 139 段错误会留下残缺文件）。

## C. 回归基线表（2026-09-18 实测，R 4.5.2 / STRING 11.5 / 人类 9606）

| 用例 | 命令要点 | 结果 |
|---|---|---|
| **参考图复刻** | 13 个基因（含 CoREST/LSD1/MACROH2A2/H2AX）+ 参考 logFC，`--color_limits=-2,0 --engine=api` | 映射 **13/13**；原始边 46；过阈值(400) **31**；节点 13；**连通分量 1**；边权 417~999；别名替换 0 次 |
| `--table` 输入 | 14 行 csv（含一个假基因 `NOTAGENE123`） | 映射 13/14；未映射 1 个并告警；网络与上一条一致；颜色区间 auto = -2.02~-0.35 |
| `engine=local` | 合成 info/links（3 蛋白 2 边），阈值 200 | 3 节点 2 边，出图正常 |
| 0 条边 | `--score_threshold=900 --keep_isolated=1` | 3 节点 0 边，**出图成功**（只画点，不画边） |
| 0 条边 + `--keep_isolated=0` | 同上 | 友好报错，说明「阈值太高 + 孤立节点也去掉了」并给对策 |
| `--add_nodes=5` | 3 个基因（HDAC2/KDM1A/RCOR1） | 节点 3 → **8**，边 23；自动 `prune=none` 并打 note |
| `--graphml=1 --format=svg` | 同上 | 出 `.graphml`（5.3 KB）与 `.svg`（116 KB） |
| `--genes_file` 纯文本 | 每行一个基因、无表头（5 个） | 输入 5，映射 5/5，8 条边（表头不会被当成基因） |
| `--genes_file` 单列表格 | 首行是 `gene` 表头 | 同上，5/5 |
| 拼错参数 | `--score_threhold=400` | 告警 `无法识别的参数：--score_threhold`（不再静默吞掉） |
| 残缺 STRINGdb 缓存 | aliases 本地 1.8 MB / 应为 20.7 MB | 拒绝启用 + 列出该删的文件；同一缓存下 `--engine=auto` 自动回退 `api` 并成功出图 |
| **STRINGdb 全链路** | **大肠杆菌 K-12 `--species=511145 --engine=stringdb`**（aliases 995 KB + info 250 KB + links 5.2 MB） | ✅ **6/6 映射，12 条边，边权 481~999，单分量**。重试后 `STRINGdb$new()` 通过（B4 是概率性的） |
| STRINGdb 全链路（人类） | `--engine=stringdb --species=9606` | ⚠️ 要下 20.7 MB aliases + 72 MB links，本机 ~80 KB/s 且易断；`load()` 还要把 1200 万行 links 读进内存。**人类数据建议直接用 `--engine=api`** |

## C2. v1.1 回归基线（2026-09-25 实测，R 4.5.2 / STRING 11.5 / 人类 9606 / 13 基因复刻用例）

| 用例 | 命令要点 | 结果 |
|---|---|---|
| **双格式输出** | 13 基因复刻用例，默认 `--format=pdf,png` | `*_ppi_network.pdf`（26,044 B）+ `*_ppi_network.png`（267,464 B），magic bytes 分别为 `%PDF-1.4` / `\x89PNG`，日志逐个登记字节数 |
| 非法格式值 | 同上 + `--format=jpg` | `[warn] 未知输出格式 'jpg'，已忽略`，照常出 pdf+png（不再产出伪 jpg） |
| **交互第 1 步** | 同上 + `--interactive=1` | `*_layout_editor.html`（12.4 KB，自包含）+ `*_positions_auto.csv` + `*_preview.png`（110,607 B）；不出最终图；打印两步流程指引 |
| **交互第 2 步** | 同上 + `--positions=<*_positions_auto.csv>` | `布局来源: manual（覆盖 13/13）`，最终 pdf+png 与自动布局视觉一致（无损往返） |
| 手改坐标重跑 | python 把 MACROH2A2 的 x+2.0 / y+0.5 存成新 csv 再 `--positions=` | 覆盖 13/13；成图中 MACROH2A2 精确移到新位置，其余 12 节点不动 |
| 编辑器浏览器自测 | Chrome headless 注入 harness（真实浏览器跑真实脚本） | PASS×5：CSV 表头/13 行、变换往返 maxerr=2.2e-16、模拟拖动后导出坐标正确、重置后 CSV 逐字节一致 |
| 编辑器渲染 | Chrome headless 截图 1280x960 | 13 节点 logFC 配色 + 31 边线宽 + 图例(-2~0) + 按钮 + 点阵网格全部正常；编辑器排布与 R 成图方向一致（y 翻转正确） |

## C3. v1.2 回归基线（2026-09-25 下午实测，默认参数裸跑）

| 用例 | 命令要点 | 结果 |
|---|---|---|
| **默认三格式 + 编辑器默认生成** | 13 基因复刻用例，**不带** --format/--interactive | `*_ppi_network.pdf` / `.png`（267,464 B）/ `.svg`（130,200 B，`<?xml` 头）+ `*_layout_editor.html` + `*_positions_auto.csv` 一次全出；打印 [layout] 微调指引 |
| `--interactive=0` | 5 基因用例 | pdf+png+svg（59,927 B svg），**无**编辑器/坐标快照（纯一键出图） |
| `--positions` 重跑（v1.2） | v1.1 拖拽坐标（MACROH2A2 x+2/y+0.5）+ 默认参数 | manual 覆盖 13/13；pdf+png+svg 全出；编辑器以拖后摆位为起点再生成（可继续迭代） |
| 编辑器交互自测（v1.1 遗留有效） | Chrome headless 注入 harness | PASS×5（CSV 结构/变换往返 2.2e-16/模拟拖动导出/重置一致）——v1.2 未动编辑器生成逻辑 |

## D. 报错对照表

| 现象 | 原因 | 处理 |
|---|---|---|
| `SSL connect error` / `cannot open the connection to 'https://string-db.org/...'` | 域名在本机不可达（见 B2） | 用 `--engine=api`（走 version11_5）；别改 `--version=12.0` |
| `Bioconductor version cannot be validated` | 主站不可达（B3） | 走西湖镜像装二进制，别用 BiocManager |
| 图上只剩孤立点、边数 0、阈值明明是 400 | score 标度（B1） | 已内置自动归一；若仍为 0，看日志里 `[score]` 那行有没有打出来 |
| `[error] STRINGdb 数据缓存残缺` | 下载中断留下的残缺 gz（B5） | 删掉列出的 `.txt.gz` 重跑，或用 `--engine=api` |
| `最终图一个节点都没有，无法出图` | 阈值过高 + `--keep_isolated=0` | 降 `--score_threshold` 到 200，或恢复 `--keep_isolated=1` |
| `--logfc 个数 (n) 与 --genes 个数 (m) 不一致` | 两个向量没对齐 | 一一对应重排（这是有意为之的硬报错） |
| `找不到 lib_string_common.R` | 没通过 `run_string.sh` 启动、或技能目录被移动 | 用 `bash run_string.sh 01_string_network.R ...` |
| `package 'ggraph' was built under R version 4.5.3` | 版本警告 | 忽略 |
| 退出码 139 | Windows R 收尾时的段错误 | 看产物文件有没有生成，有就是成功 |
