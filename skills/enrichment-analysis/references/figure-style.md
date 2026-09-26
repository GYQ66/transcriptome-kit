# 美化主图的版式与配色（逐项对照参考代码）

参考代码：`<参考脚本>/kegg美化.R`（`gground` + `ggprism`）。
本文档记录**每一个版式参数的来源**，以及本技能**故意改掉**的地方。

## 1. 图形结构（从外到内）

```
x 轴范围：[-1.6, xaxis_max]，expand = c(0,0)，breaks = seq(0, xaxis_max, 2)
  其中 xaxis_max = max(-log10(p.adjust)) + 1

  ┌─ 分类色块区 ─┬─ 基因数圆点 ─┬─ 数据区（圆角柱 + 文字）────────────┐
  │ x ∈ [-1.5,-1]│ x = -0.5     │ 柱：x ∈ [0, -log10(p.adjust)]      │
  │ 圆角矩形      │ shape=21 圆  │ 通路名：x = 0.05, hjust = 0        │
  │ 内写分类名    │ 中写 Count   │ 基因名：x = 0.10, vjust = 2.6,     │
  │              │              │         斜体、分类色              │
  └──────────────┴──────────────┴────────────────────────────────────┘
  底部另有 annotate("segment", 0 → xaxis_max, y = 0, linewidth = 1.5) 作 x 轴线
```

| 元素 | 参考代码写法 | 本技能写法 | 是否改动 |
|---|---|---|---|
| 圆角柱 | `geom_round_col(aes(y = Description), width = 0.6, alpha = 0.8)` | 同 | 一致 |
| 通路名 | `geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 5)` | 同 | 一致 |
| 基因名 | `geom_text(aes(x = 0.1, label = geneID, colour = ONTOLOGY), hjust = 0, vjust = 2.6, size = 3.5, fontface = "italic", show.legend = FALSE)` | 同，但默认改成 `y = index - gene_dy`（**贴在柱下方**，见 §4.5）；取 `geneName`（symbol）列 | **改了** |
| 基因数圆点 | `geom_point(aes(x = -width, size = Count), shape = 21)` + `geom_text(aes(x = -width, label = Count))` | 同 | 一致 |
| 圆点尺寸 | `scale_size_continuous(name = "Count", range = c(5, 12))` | 同 + 图例最多 4 档（见 §4.6） | 略改 |
| 分类色块 | `geom_round_rect(aes(xmin = -3*width, xmax = -2*width, ymin = lag(cumsum(n), default=0)+0.6, ymax = cumsum(n)+0.4, fill = ONTOLOGY), radius = unit(2,"mm"), inherit.aes = FALSE)` | 同，但 `xmin` 会**往左自适应加宽**（见 §4.4） | **改了** |
| 色块内文字 | `geom_text(aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label = ONTOLOGY))` | 同，但 label 是**折行后**的文本（`KEGG` → `KE`/`GG`） | **改了**，见 §4.4 |
| x 轴线 | `geom_segment(aes(x=0,y=0,xend=xaxis_max,yend=0), linewidth = 1.5, inherit.aes = FALSE)` | `annotate("segment", ...)` | **改了**，见 §4 |
| 主题 | `theme_prism()` + 隐藏 y 轴文字/轴线/刻度，`legend.title = element_text()` | 同 + `coord_cartesian(clip = "off")` | **改了**，见 §4 |
| 图例 | `scale_fill_manual(name = "Category", values = pal)` + `scale_colour_manual(values = pal)` | 同，位置可用 `--legend=` 调 | 一致 |
| x 轴标签 | `labs(y = NULL)`，无标题 | 另加可选 `--title` | 一致（默认无标题） |

`width <- 0.5`（参考代码里的局部变量）→ 本技能的 `--bar_width=0.5`。

## 2. 挑条目的三种方式

### 2.1 `--mode=per_ontology`（默认，对齐参考代码）

参考代码分两步：

```r
group_by(GO, ONTOLOGY) %>% top_n(5, wt = -p.adjust) %>%   # 每个 GO 本体取 p 最小的 5 条
  group_by(p.adjust) %>% top_n(1, wt = Count) %>%          # 同一 p.adjust 只留 Count 最多的
  rbind( top_n(KEGG, 5, -p.adjust) %>% group_by(p.adjust) %>% top_n(1, wt = Count) %>%
           mutate(ONTOLOGY = 'KEGG') )
```

本技能用 `order(p.adjust, -Count)` + `head(top)` 达到同样效果（`--top` 默认 5），
去重由 `--dedup_padj=TRUE` 控制。**注意去重会让每类实际条目少于 5 条**（实测
GO_up 20 条 → 15 条），这是参考代码的既有行为，不是 bug。

### 2.2 `--mode=global`（"只看排名前 10"）

全库按 `p.adjust` 排名取前 `--top` 条，不再限制每类 5 条。

> ⚠️ `--dedup_padj=TRUE`（默认）会把 **p.adjust 完全相同**的通路合成一条。
> GO 里常有基因集完全重叠的通路（`extracellular matrix organization` 与
> `extracellular structure organization` 的 p 值、基因集都一样），
> 实测 `--top=10` 最终只出 **6 条**。要足额就 `--dedup_padj=FALSE`。

### 2.3 `--category`（"只看某个功能方向"）

配合 `04_enrich_summary.R` 的 `<前缀>_func_pathways_<方向>.csv`：

- 逐条通路先被归到功能大类（`func_cat`），再按大类内部**基因重叠**去冗余，
  留下 `representative = TRUE` 的代表通路
- `--top=0` = 不限制条数，`--use_representative=TRUE`（默认）只画代表
- 与 `--dedup_padj` 的区别：**`--dedup_padj` 按 p 值去重**（同一 p 值只留一条，
  可能误删不同的通路）；**04 的代表筛选按基因重叠去重**（p 值不同但基因集几乎相同的
  通路才会被折叠）。用了 `--func_csv` 时脚本默认关掉 `--dedup_padj`。

三种方式画的都是同一套版式，只是"选中哪些条目"不同。

## 3. 分类顺序与配色

参考代码：

```r
pal <- c('#c3e1e6', '#f3dfb7', '#dcc6dc', '#96c38e')            # 最终采用的一套
ONTOLOGY = factor(ONTOLOGY, levels = rev(c('BP','CC','MF','KEGG')))
```

`levels = rev(...)` → `c('KEGG','MF','CC','BP')`，配合 `scale_fill_manual(values = pal)`
得到：

| 分类 | 色号 | 说明 |
|---|---|---|
| KEGG | `#c3e1e6` | 浅青 |
| MF | `#f3dfb7` | 浅黄 |
| CC | `#dcc6dc` | 浅紫 |
| BP | `#96c38e` | 浅绿 |

**图里 KEGG 在最下、BP 在最上**（按 levels 顺序 `arrange` 后再取 `index`，
`y = index` 自下而上）。

另外两套（参考代码里被注释掉的候选，可用 `--pal_idx=1|2` 切换）：

| 套 | KEGG | MF | CC | BP |
|---|---|---|---|---|
| `pal1` | `#7bc4e2` | `#acd372` | `#fbb05b` | `#ed6ca4` |
| `pal2` | `#eaa052` | `#b74147` | `#90ad5b` | `#23929c` |
| `pal3`（默认） | `#c3e1e6` | `#f3dfb7` | `#dcc6dc` | `#96c38e` |

### ★ 与参考代码的**唯一实质性差异**

参考代码传的是**裸向量** `values = pal`。ggplot 按 levels 顺序逐一取色，所以：

- 四类齐全 → 与上表一致；
- 少一类（比如 CC 一条都没进 TOP5）→ 颜色**整体平移一格**，
  `MF` 会拿到 `#c3e1e6`（本来是 KEGG 的颜色），BP 拿到 MF 的颜色……

本技能改成**按分类名绑定**（`build_pal()`，见 `lib_enrich_common.R`）：

```r
CANON_CATS   <- c("KEGG", "MF", "CC", "BP")          # 与 levels = rev(c('BP','CC','MF','KEGG')) 对应
base_named   <- setNames(pal[seq_along(CANON_CATS)], CANON_CATS)
```

于是 `all` / `up` / `down` 三张图里同一个分类永远是同一个颜色，便于并排比较；
四类齐全时与参考代码**完全一致**。超出的分类（例如叠加了 MSigDB 的 `H`、`C2`）
用 `EXTRA_CATS_PAL` 依次取色，或用 `--pal=` 显式给 4 个色号。

**按功能方向出图时这条尤其重要**：一个功能方向（如「细胞外基质与黏附」）可能只落在
CC 和 KEGG 上，若用裸向量配色，CC 会拿到本该属于 KEGG 的浅青；名称绑定后 CC 永远是浅紫、
KEGG 永远是浅青，和四分类的图对照着看不会看错。

## 4. 其他稳健性改动

| # | 参考代码 | 问题 | 本技能 |
|---|---|---|---|
| 1 | `mutate(Description = factor(Description, levels = Description))` | 同一 `Description` 在 BP/CC/MF 重复出现时，R 4.x 直接报 `duplicated levels in factors are deprecated` 并中断 | 改成 `levels = unique(Description)` |
| 2 | `geom_segment(aes(x=0,y=0,xend=xaxis_max,yend=0), inherit.aes = FALSE)` | 标量美学按数据行数广播，抛 `All aesthetics have length 1, but the data has N rows` | 改用 `annotate("segment", ...)` |
| 3 | 一次只画一套、只叠 GO+KEGG | 分类被筛空 → 整个脚本报错或出错图 | 逐类检查，某类为空时打印告警并退化为"该类全部条目"；全空时退化为不筛选并明确告警 |
| 4 | `ggsave`/`pdf` 固定 12×8 | 条目多时文字挤压 | 高度按条目数自适应：`max(--height, 0.34 * n + 2.2)` |
| 5 | 基因名全列 | 基因多时冲出画布 | 见 §4.1 |
| 6 | 无 `clip` 设置（默认裁切到面板） | **左侧分类色块里的 `KEGG` 被裁成 `EGG`**（色块只有 0.5 个数据单位宽，装不下 4 个字符） | `coord_cartesian(clip = "off")`，并把基因名字体按最长字符串自动缩小 |
| 7 | `geom_text(..., size = 3.5)` 固定 | 基因名长时会冲进右侧图例区 | 自动缩字号（见 §4.2），或 `--auto_width=TRUE` / `--max_genes=N` |

### 4.1 x 轴 / 面板裁切的坑

参考代码里色块在 `x ∈ [-3w, -2w] = [-1.5, -1]`，文字画在中心 `x = -1.25`。
面板左边界正好在 `-1.5`（`expand = c(0, 0)`），于是 **`KEGG` 的左半边被裁掉**，
出图上看到的是 `EGG`。实测本机 12 英寸、`xaxis_max ≈ 16` 时必然发生。

两种修法：
- **扩大左侧范围**（改 `scale_x_continuous(limits=...)`）：会改变与参考代码的
  版面比例，且需要的左边距随 `xaxis_max` 变化（柱子越长，1 个数据单位对应的
  像素越少，同样宽的文字反而占更多数据单位）；
- **`clip = "off"`**（本技能采用）：文字可以画到面板外。副作用是基因名也不再被裁，
  所以要配合 §4.2 的自动缩字号。

### 4.2 基因名长度是完全由数据决定的

一条通路富集到 40 多个基因时，`geneID` 拼出来的字符串可以到 **260+ 字符**
（实测 `GSE62452_T_GO_all.csv` 的 `external encapsulating structure`）。
按参考代码的固定 `size = 3.5`，这段文字需要约 **14 英寸**，而画布只有 12 英寸。

三种策略，本技能默认第一种：

| 策略 | 参数 | 效果 |
|---|---|---|
| **自动缩字号**（默认） | `--auto_gene_size=TRUE` | 261 字符 → `3.5` 缩到 `2.27`；内容一个不丢 |
| 自动加宽画布 | `--auto_width=TRUE` | 保持 3.5 字号，画布加宽到约 17.6 英寸 |
| 截断并标注 | `--max_genes=20` | `A/B/C ... (+23)` |

字号换算（脚本里用的经验公式）：`geom_text` 的 `size` 单位是 **mm**，
实际 pt = `size * 2.845`；基因名以大写字母、数字、斜杠为主，平均字宽约 `0.4 em`
→ 每字符约 `0.0158 * size` 英寸。按 `0.5 em` 估会高估 25%，
把字号压得过小（实测从 2.27 压到 1.9），别用 0.5。

可用横向空间取 `0.78 * width` —— 右侧要留 `Category` / `Count` 两个图例的位置；
用 0.82 时最长的那几条会正好贴到图例上（实测）。

### 4.3 临时文件

`<前缀>_enrich_<方向>_selected_terms.csv` 是**入选清单**（含 `index`/分类/p.adjust/Count/基因名）。
图出来后想改 `--top` 或配色，直接改 CSV 重画也行，不用重跑富集。

### 4.4 ★ 分类标签折行 + 色块自适应加宽

参考代码把色块写死成 x 属于 [-3w, -2w] = [-1.5, -1.0]，宽度 0.5 个**数据单位**。
问题是"数据单位 → 英寸"的换算随 x 轴范围变化：

| xaxis_max | 面板数据范围 | 色块实际宽度 | KEGG 需要 | 结果 |
|---|---|---|---|---|
| 16 | 17.5 单位 / 9.4 in | 0.27 in | 约 0.37 in | **装不下，被裁成 EGG** |
| 60 | 61.5 单位 / 9.4 in | 0.076 in | 约 0.37 in | 连一个字符都放不下 |

所以参考代码的 KEGG 在多数图上都会被面板左边界裁掉一半；加上 `clip="off"` 之后
虽然画全了，却变成"文字冲出彩色框"，同样难看。

本技能分两步做：

1. **折行**：按"这块宽度能放几个字符"把标签切开（KEGG → KE + GG 两行）。
   字符宽度用 `CHAR_IN_CAT = 0.0245 * size` 英寸 —— 全大写标签比混合大小写的
   基因名宽，若沿用基因名的 0.0158 会低估 35%，算出来的框仍然装不下。
2. **色块往左加宽**到刚好包住折行后的文字：
   `box_units = max(w, 最宽一行字符数 * 字符宽 * 1.15 / 每单位英寸数)`。
   **右边界固定不动**（仍是 -2w），所以：左侧基因数圆点（x = -w）与 x = 0 的位置
   都不变，只是左边距占多一点，版面右半边完全不受影响。

因为"块宽 ← 需要多宽 ← 块宽决定的数据单位换算"是循环依赖，实现里迭代 3 次收敛。
`--cat_label_chars=N` 可强制每行字数（想固定版式时用），`--cat_size` 调字号。

实测：xaxis_max=16 时色块 0.51 单位、xaxis_max=24 时 0.62 单位
（参考代码都是 0.50），KEGG 都折成 KE/GG 且被框完整包住。

### 4.5 ★ 基因名的位置：贴在彩色柱下方（跟着柱高走）

参考代码用 `vjust = 2.6` 把基因名往下推。**`vjust` 的单位是文字高度**，
所以它的位置只跟字号有关、**和柱子多高完全无关** —— 结果基因名总是压在彩色柱
里面（柱高 0.6 个 y 单位，文字下移 2.6 个字高后仍落在柱子范围内）。

本技能改成显式给 y：

```r
line_in      <- size * 2.845 / 72 * 1.28   # 一行文字高度（含行距），英寸
text_h_units <- line_in / row_in           # row_in = 画布高 / 条目数
gene_dy      <- bar_height / 2             # 柱下沿
                + text_h_units / 2         # 再往下半个字高（文字顶边贴住柱下沿）
                + 0.05                     # 一点空隙
y = index - gene_dy
```

于是 `--bar_height` 一改，基因名**自动跟着柱子移动**，不用手调。
`--gene_pos=inside` 可切回参考代码那种"压在柱子里"的效果；
`--gene_dy=<数值>` 可手动指定偏移。

实测（12 英寸宽、16 行、高 7.64 英寸、柱高 0.6、字号 2.27）：
gene_dy = 0.470 行 = 0.225 英寸 → 文字落在柱下沿之外的空白里，
与下一行柱顶之间还有约 0.11 行余量，不会互相压。

### 4.6 高度双向自适应 + 图例档位

参考代码固定 12×8 英寸。条目多时文字会挤，条目少时柱子会胖得难看
（实测 4 条时每行 2 英寸、柱高 1.2 英寸，而且两个图例竖着排不下、
`Count` 标题被裁掉）。

- **高度**：不给 `--height` 时用 `max(4.6, 0.34 * 条目数 + 2.2)`；
  显式给了 `--height` 就完全按你给的来（不再自适应）。
  实测 16 条 → 7.64 in、15 条 → 7.30 in、4 条 → 4.60 in。
- **图例**：`Count` 的档位压到最多 4 个（`seq(min, max, length.out = 4)` 取整去重），
  否则矮画布下图例竖排会超出设备高度被裁。

## 5. 出图前先核对的表

`03_enrich_plot.R` 会把最终入选的条目落成
`<前缀>_enrich_<方向>_selected_terms.csv`（`index/ONTOLOGY/ID/Description/p.adjust/Count/geneName`）。
**先看这张表**：条目数、每类几条、分类顺序对不对，再去看图。

## 6. 数值口径提醒

- 柱长 = `-log10(p.adjust)`，**用的是 `p.adjust` 不是 `pvalue`**。
  想换成 `pvalue` 就改 `--fig_filter`（过滤）配合 `01` 的输出表，
  柱长口径目前固定为 `p.adjust`（与参考代码一致）。
- GSEA 的条形图用另一套约定：横轴 `NES`，**红 `#B2182B` = 激活（NES>0）、
  蓝 `#2166AC` = 抑制（NES<0）**，与 GEO 技能的期刊风格配色保持一致。
