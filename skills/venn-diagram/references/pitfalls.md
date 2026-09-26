# Venn 图踩坑手册（VennDiagram 1.8.2 / R 4.5.2 / Windows）

**改动 `01_venn.R` 之前先把这份读完。** 下面每一条都是在本机实测出来的，不是猜的。

---

## 0. ★ 别把"退出码 3 / 4"当成报错

出图前有**两道组名闸门**（都在 `01_venn.R`），组名一律由用户手动输入，脚本不自动识别。

**第 1 道：命名闸门（退出码 4）** —— 名字没给齐时触发：

| 情形 | 结果 |
|---|---|
| `--set=` 缺标签（`--set=a.txt` 这种写法） | 收集全部缺名文件（按 `--set` 出现顺序）→ 打印"待命名清单" → `quit(status = 4)`，**不读文件、不画图** |
| `--table` 模式没给 `--by_levels` | 打印分组列的全部取值 → `quit(status = 4)` |

（2026-09-25 之前的行为：`--set` 缺标签会用**文件名去后缀**当组名、`--table` 会按
字母序自动取分组值当组名 —— 已废除，组名必须用户手动输入。）

**第 2 道：确认闸门（退出码 3）**（`01_venn.R` 第 4b 节）。名字齐了之后：

| 情形 | 结果 |
|---|---|
| 没给 `--labels_confirmed=TRUE` | 打印组名清单 → `quit(status = 3)`，**不开任何设备、不画图** |
| 给了标志，但 `<outdir>/<前缀>_labels.txt` 不存在 | 同上（第一次跑，还没有清单可比对） |
| 给了标志，但本次组名和清单里的**逐字不一致** | 打印差异 → 退出 3（名字改过 = 用户还没看过新名字） |
| 给了标志且组名完全一致 | 继续出图，退出码 0 |

两个容易踩的点：

1. **每次运行都会重写 `<前缀>_labels.txt`**（包括停在确认的那次）。所以"改名 + 带标志"
   那一次本身也会把新名字写进去 —— 于是**必须再跑一次**才能出图。别以为第一次就画了。
   顺序恒为：记录/展示 → 确认 → 再跑 → 出图。
2. **比对只看第 2 列（显示名），不看基因数**。基因列表换了但名字没变 → 直接出图
   （这是有意的：确认的是"图上叫什么"，不是"圈里有多少个"）。

比对用的 `lab_names` 是 `vapply(..., USE.NAMES = FALSE)` 的结果。**漏掉
`USE.NAMES = FALSE` 会让 `identical()` 因为"带名字的向量 vs 不带名字的向量"而恒为 FALSE**，
表现为"明明没改名却一直说对不上"（实测踩过一次，已修）。

## 1. 参考代码里那批"画布参数"其实全是空的

```r
p1 <- venn.diagram(..., filename = NULL, imagetype = "pdf",
                   height = 180, width = 180, resolution = 200, margin = 0.25)
```

`venn.diagram()` 的源码里，开绘图设备那段是包在 `if (!is.null(filename))` 里的：

```r
if (!is.null(filename)) {
  if ("tiff" == imagetype) tiff(filename = filename, height = height, width = width,
                                units = units, res = resolution, ...)
  else if ("png" == imagetype) png(filename = filename, height = height, width = width,
                                   units = units, res = resolution)
  ...
}
```

`filename = NULL` 时**一个设备都不会开**。所以 `height/width/resolution/imagetype`
四个参数在这个调用里完全没有作用 —— 参考代码实际是把 grob 拿回来自己
`grid.draw()`，画布大小由当时的设备决定。

另外 `venn.diagram` 的默认值是 `height = 3000, width = 3000, units = "px",
resolution = 500` —— 注意 **`units` 默认是 `"px"`**。真要走 `filename=` 那条路时，
`height = 180, width = 180, resolution = 200` 得到的是 180x180 像素、0.9x0.9 英寸的图。

`margin` 更彻底：它不是 `venn.diagram` 的形参，会被塞进 `...` 转发给
`draw.pairwise.venn` 等函数；而 **1.8.2 里那些 draw 函数已经没有 `margin` 形参了**
（实测：`"margin" %in% names(formals(draw.pairwise.venn))` 为 FALSE，函数体里也
完全没出现过 `margin`）。它最终落在 draw 函数的 `...` 里，悄无声息。

→ 本技能的做法：**自己开设备**（`--width/--height/--dpi/--format` 真生效），
**自己用视口实现留白**（`--margin` 真生效）。

## 2. `disable.logging` 这个名字是反的

```r
if (disable.logging) {
  flog.appender(appender.console(), name = "VennDiagramLogger")
} else {
  flog.appender(appender.file(paste0(if (!is.null(filename)) filename else "VennDiagram",
                                     ".", time.string, ".log")), name = "VennDiagramLogger")
}
out.list <- as.list(sys.call()); out.list[[1]] <- NULL
out.string <- capture.output(out.list)
flog.info(out.string, name = "VennDiagramLogger")
```

- `disable.logging = FALSE` + `filename = NULL` → 在**当前工作目录**写一个
  `VennDiagram.<时间戳>.log`。参考代码就是这个写法，所以每跑一次就多一个日志文件。
- `disable.logging = TRUE` → 不是"不记录"，而是**改打到控制台**，而且记的是
  **整个调用参数列表**：`x` 里全部基因名会被逐行 `print` 出来，几千行噪声。

→ 本技能恒用 `TRUE`，并在调用前把那个 logger 的门槛抬高：

```r
if (requireNamespace("futile.logger", quietly = TRUE)) {
  try(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"),
      silent = TRUE)
}
```

（`futile.logger` 是 VennDiagram 的依赖，一定在。加了这行之后日志从 4000+ 行降到 30 行。）

## 3. "不传" ≠ "传 NULL"

```r
venn.diagram(..., cat.pos = NULL)
# Error: Unexpected parameter length for "cat.pos"

venn.diagram(..., cat.dist = NULL)          # 4 组图时
# Error: Unexpected parameter length for "cat.dist"
```

这些参数在 draw 函数里是**按长度校验**的（`if (length(cat.pos) != n) stop(...)`），
`NULL` 的长度是 0，会直接报错。要在 4/5 组图里用包自带默认值，只能**根本不传**。

→ 本技能在组装参数时先把 `NULL` 项整个丢掉：

```r
extra <- extra[!vapply(extra, is.null, TRUE)]
```

## 4. 多传的参数会被静默吞掉

`venn.diagram` 把 `...` 一路转发给 `draw.pairwise.venn` / `draw.triple.venn` /
`draw.quad.venn` / `draw.quintuple.venn`，而这些函数签名末尾都有 `...`，
**不认识的参数既不报错也不生效**。实测哪些参数对哪些 n 无效：

| 参数 | 2 组 | 3 组 | 4 组 | 5 组 |
|---|---|---|---|---|
| `scaled` / `euler.d` | ✓ | ✓ | ✗（无此形参，被吞） | ✗ |
| `sep.dist` | ✓ | ✓ | ✗ | ✗ |
| `margin` | ✗ | ✗ | ✗ | ✗（1.8.2 全都没有） |
| `print.mode` / `sigdigs` | ✓ | ✓ | ✓ | ✓（venn.diagram 会显式转发） |

→ 本技能按目标 draw 函数的 `formals()` 过滤一遍，并把被丢弃的参数**打印出来**，
避免"我明明设了却没效果"。

## 5. ★ 圈外标签被画布边缘裁掉（最常见的翻车点）

现象：图里左右两侧的标签只剩半截（`hdWGCNA \Fibroblasts M2` 变成 `GCNA \Fibroblasts M2`）。

原因：`cat.pos = c(-90, 90)` 把标签放在圈的正左/正右，而 VennDiagram 计算出的
锚点是 **npc 0.016 / 0.984** —— 几乎就贴在画布边上（实测：`cat.dist` 从 0 到 0.3
变化时，锚点只在 0.021 → 0.013 之间晃）。标签是**以锚点为中心**画的，
于是半个标签必然甩到画布外。

为什么"放大画布"没用：锚点到画布边缘的距离只占画布宽度的 1.6%，画布从 4 英寸涨到
11 英寸，这点距离只从 0.06 英寸涨到 0.18 英寸，而标签半宽是**固定的英寸数**
（实测 `hdWGCNA \Fibroblasts M2` 在 cat.cex=1.3 时宽 2.52 英寸，半宽 1.26 英寸）。
实测 4/5.2/6.8/8.8/11.4 英寸的画布，标签**始终**是裁的。

正确的解法：**加大留白** —— 图形整体往内缩，锚点跟着往内走（每加 1 份留白，
锚点的英寸位置就往外挪 0.97 份），半个标签就有地方了。这也正是当年 VennDiagram
的 `margin` 参数干的事（参考代码写 0.25 就是出于这个目的，可惜在 1.8.2 里失效了）。

→ 本技能的 `--autofit=TRUE`（默认）把这件事算准：

1. 从 grob 里把每个 `text` 元素捞出来（gList / gTree 递归）；
2. 用 `convertWidth(grobWidth(tx), "in")` 量出它的真实英寸尺寸；
3. 解不等式：锚点位置 `P(m) = T*u + m*T*(1-2u)`，要求 `P(m) ± size/2` 落在 `[pad, T-pad]`；
4. 取所有约束里最大的下界作为留白；若超过 0.36 上限就把画布放大 1.25 倍重算；
5. 出图前再用**低分辨率试渲 + 读回像素找内容 bbox** 复核一遍。

一个坑中坑：`u = 0.5` 的元素（圈里的交集数字）会让 `(1-2u) = 0`，除零得 `Inf`，
整个结果变成 `NA` —— 这类元素跟留白无关，**必须跳过**。

另一个坑：**`dev.capture()` 在 Windows 的 png 设备上不可用**，实测报
`raster capture is not available for this device`。核裁切只能走
"渲成 png 文件 → `png::readPNG()` 读回来"。

## 6. 2~3 组图的圈只在正方形画布上才是正圆

VennDiagram 用 npc 坐标画图，画布被拉长时圈会跟着变成椭圆。实测 5x5 时
长短轴比 = 1.000，5x10 时 = 2.0。→ 想让圈是正圆，`--width = --height`。

4~5 组图那种"椭圆"**不是被拉的**：VennDiagram 的 4/5 组版式本来就用旋转椭圆拼
（实测 5x5 画布下 n=4 的长短轴比是 1.75~2.33、n=5 是 1.92），跟画布比例无关。

**量圆的宽高比时不要用"某个颜色的像素范围"**：左圈单独露出来的那块是**月牙形**，
量出来必然比圆窄（我第一次就被这个骗了，以为圈是椭的）。要量就用
`convertX(unit(polygon$x, "npc"), "in")` 把多边形顶点换算成英寸再做 PCA，
或者干脆量圆的完整高度。

## 7. 5 组图必然很难看

31 个区、几十个数字堆在 5 个椭圆里，实测数字会互相压、看不清。这不是参数问题。
需要 5 组以上时，优先建议 UpSet 图 / 花瓣图，别硬出 Venn。

## 8. 其它小的

- `category.names` 的长度必须**恰好等于**集合数，多了少了都报
  `Unexpected parameter length`。
- `\n` 在标签里是有效的换行（R 4.x 的 grid 支持 `textGrob` 多行文本，实测 3 行正常）。
  bash 里写 `--set="A\nB=file"`，双引号里 `\n` 会原样传给 R，脚本再做 `gsub("\\\\n", "\n")`。
- **图内文字一律 ASCII**：Windows 下中文字形缺失，会渲染成乱码方块/问号。
- 集合名里带 `/` `&` 不影响出图，但导出 CSV 时要保证 `_regions.csv` 的 `region` 列
  不含换行（脚本已把标签里的 `\n` 换成 `/`）。
- 重复基因：`force.unique` 默认 TRUE（去重）；`na` 默认在本技能里是 `"remove"`
  （包默认是 `"stop"`，遇到 NA 直接中止）。
- `--col` 是**圈线颜色**，不是"列名"；列名参数是 `--gene_col`。这两个名字撞过车，
  实测报的错很有迷惑性：`invalid color name 'gene'`。
- **标签里不能写 `=`**：`--set=` 用第一个 `=` 分隔标签和文件，写
  `--set=A B=x=file` 会把 `x=file` 整段当路径，报"基因列表文件不存在"（实测踩过）。
- 顶层 `if (...) { try(...) }` 在 Rscript 里会把返回值 `NULL` 打到 stdout，
  日志里会多一行莫名其妙的 `NULL` —— 用 `invisible(...)` 包住。
