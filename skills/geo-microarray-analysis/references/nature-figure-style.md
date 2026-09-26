# Nature 期刊风格绘图（火山图 / 热图）

`02_deg_plots.R` 内置了 `--fig_style=nature`，把火山图与热图切换为
高影响因子期刊（Nature 系）常用的版式与配色约定。默认 `--fig_style=classic`
保持原有输出完全不变。

> **来源与署名**：本风格改编自开源项目 **nature-skills**
> <https://github.com/Yuan1z0825/nature-skills>（Apache-2.0）中的 `nature-figure`
> 技能，取其已发布的火山图模板、期刊主题、配色体系与导出策略，改写为自包含的
> R 实现。原实现为 Python/matplotlib，本处为 R/ggplot2 + pheatmap 等价移植。

## 启用方式

```bash
bash "$SK/run_geo.sh" "$SK/02_deg_plots.R" \
  --expr=./out/GSE15471.csv \
  --group_file=group.csv \
  --contrast=case-control \
  --prefix=GSE15471_cc --fig_style=nature --label_top=15 --top=40
```

## 与 classic 的差异

| 项目 | classic（默认） | nature |
|---|---|---|
| 主题 | `theme_bw(base_size=14)` | `theme_nature()`：`theme_classic` + Arial(缺省回退无衬线) + 0.35pt 细轴线 + 无网格 + 5~7pt 字号 |
| 火山图配色 | 红 `#C31E1F` / 蓝 `#1F6FC3` / 灰 `#898989` | 红 `#B2182B` / 蓝 `#2166AC` / 灰 `#B3B3B3` |
| 点大小 / 透明度 | 1.75 / 0.4 | 0.7 / 0.55 |
| 图例 | `Up`/`Down`/`NOT` | **带计数** `Up (n)` / `Down (n)` / `Not significant (n)`，置于图下方 |
| 标题 | 有（阈值摘要两行） | 默认无（期刊图不带标题）；`--fig_title=` 可强制加 |
| 基因标注 | 仅 `--label_genes` 手动 | `--label_genes` ∪ `--label_top=N`（自动选 p 值最小的前 N 个显著基因） |
| 热图色带 | `#1F6FC3`-白-`#C31E1F` | 发散 `#2166AC`-`#F7F7F7`-`#B2182B` |
| 热图边框 | `grey60` | 无（`border_color=NA`） |
| 分组注释色 | Set2 系 | 蓝/青/紫/红期刊色系 |
| 导出格式 | PDF + PNG(300dpi) | PDF(cairo) + PNG(300dpi) + **600dpi TIFF** + SVG（需 svglite，缺失则跳过） |

## 配色体系（`NATURE_PALETTE`）

```
blue_main      #0F4D92   # 深蓝 —— 主方法 / 关键对象
blue_secondary #3775BA   # 中蓝
red_strong     #B64342   # 强调红
neutral_light  #CFCECE   # 中性浅灰
neutral_mid    #767676   # 中性中灰
neutral_dark   #4D4D4D   # 中性深灰
```

火山图三态用已验证的 `NATURE_VOLCANO_COLS`：`UP=#B2182B`、`DOWN=#2166AC`、
`NOT=#B3B3B3`（蓝=下调、红=上调、灰=不显著）。

## 尺寸

- 火山图：期刊单栏约 **89 × 82 mm**，可用
  `--vol_width_mm=89 --vol_height_mm=82` 精确指定。
- 热图：单栏约 **89–120 mm**，或按基因数调整；用 `--hm_width_mm/--hm_height_mm`。
- 仍兼容原英寸参数 `--vol_width/--vol_height`、`--hm_width/--hm_height`（mm 优先）。

## 导出与依赖

- 导出设备（本机实测）：
  - PDF 必须用 **`cairo_pdf`**。默认 `pdf()` 在未注册 Arial 时会抛
    `invalid font type` 并使脚本中断——脚本已自动改用 `cairo_pdf`。
  - PNG/TIFF 用 **`ragg`**（能解析系统字体），600 dpi 出 TIFF。
  - SVG 需要 **`svglite`**；本机未安装时会**自动跳过**并打印提示，
    **矢量图可用 PDF 替代**（`cairo_pdf` 输出的 PDF 文字同样可编辑）。
    需要 SVG 时安装：`install.packages("svglite")`。
- 期刊投稿通常要求 **矢量（SVG 或 PDF）+ 600dpi TIFF**；本流程在
  nature 风格下默认产出 PDF + PNG + TIFF（有 svglite 时再加 SVG）。

## 期刊版式约定（取自 nature-figure 设计要点）

- 只保留左/下轴线，无网格；用稀疏的刻度引导视线。
- 每个统计比较在图例或源数据注中标明 n、中心、离散度、检验与校正。
- 颜色语义一致：同一含义在全图用同一色相；连续量用浅→深顺序色带。
- 发散色带（红-白-蓝）仅用于有正负方向的量（如 z-score 行标准化）。

## 组合示例

```bash
# 期刊风格 + 自动标注 top 15 基因 + 热图展示 top 40
bash "$SK/run_geo.sh" "$SK/02_deg_plots.R" \
  --expr=./out/GSE15471.csv --group_file=group.csv \
  --contrast=case-control --prefix=GSE15471_cc \
  --fig_style=nature --label_top=15 --top=40 \
  --vol_width_mm=89 --vol_height_mm=82 --hm_width_mm=110 --hm_height_mm=140 \
  --p_type=adj.P.Val

# 同时手动指定几个关键基因
  # ... --label_genes=POSTN,INHBA,ALB
```
