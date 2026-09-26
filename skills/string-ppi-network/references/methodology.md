# 方法论：口径、参数与与参考代码的差异

## 1. 三步取数的口径

参考代码的链路是：

```
gene_list → data.frame → string_db$map() → mapped(STRING_id)
        → string_db$get_interactions() → ppi_edge(from/to/combined_score)
        → graph_from_data_frame() → igraph → ggraph
```

本 skill 三个引擎在同一条链路上替换「取数」这一段，其余完全一致：

| 环节 | stringdb | api | local |
|---|---|---|---|
| symbol → STRING_id | `sdb$map()`（查 aliases 表 + protein.info） | `GET /tsv/get_string_ids?limit=1&echo_query=1` | 查 `protein.info` 的 `preferred_name` 列 |
| STRING_id → 边 | `sdb$get_interactions()`（`load()` 全量 links 表 + `graph.data.frame`） | `GET /tsv/network` | 读本地 `protein.links` |
| 服务端过滤 | `score_threshold` 在 `load()` 里生效 | `required_score`（**0–1 口径**） | 本地过滤 |
| 本地过滤 | `combined_score >= score_threshold` | 同左（**先做标度归一**） | 同左 |

三个引擎查的都是 STRING 11.5（默认），**结果同源**。差异只在「谁做过滤」和「有没有本地缓存」。

## 2. combined_score 的标度陷阱（最容易出事的一处）

STRING 的 `combined_score` 官方口径是 **0–1000 的整数**（400 = medium confidence，
700 = high，900 = highest）。但**官方 REST API 在不同版本域名上返回过 0–1 的小数**：

```
# version11_5.string-db.org/api/tsv/network 实测返回
stringId_A  stringId_B  ...  score
9606.ENSP00000262241  9606.ENSP00000434024  0.427      ← 0-1 标度！
```

直接按 0–1000 口径卡 `>= 400`，**所有边都会被滤光**（实测：31 条边 → 0 条边，
然后图里 13 个孤立点）。脚本的处理：

1. 拿到边表后先看 `max(combined_score)`；
2. `max <= 1.0001` → 判定为 0–1 标度，**×1000 归一**，并打印
   `[score] 引擎返回 0-1 标度，已 x1000 归一到 STRING 0-1000 口径`；
3. 之后所有过滤、`score_norm`、报告里的「边权范围」都按 0–1000 处理。

请求侧也做了对应处理：`required_score` 一律传 `score_threshold/1000`（0–1 口径），
这样即使服务端理解成 0–1 也不会把边全砍掉，最终以本地过滤为准。

`score_norm` 的定义与参考代码一致：`score_norm = combined_score / max(combined_score)`，
所以**线宽是相对于这张图的**，不同图之间的线宽不能直接比。

## 3. 置信度阈值怎么选

| 值 | STRING 官方语义 | 什么时候用 |
|---|---|---|
| 150 | low | 只想看有没有连接 |
| **400** | **medium（参考代码的取值）** | 默认，Nature 这类图的常规起点 |
| 700 | high | 要「确定相关」的子网络 |
| 900 | highest | 极严，通常只剩核心复合体 |

同一批基因降阈值会让网络变密、连通分量变少。参考代码注意事项第 4 条提到
`score_threshold=400` 可能让网络出现孤立点 —— 本 skill 会直接报出
`连通分量` 与 `孤立节点` 数量，不用自己数。

## 4. network_type

- `functional`（默认）：包含功能关联（共表达、共现、文本挖掘等），边多，适合看通路/复合体。
- `physical`：只保留物理结合证据（实验、数据库、融合），边少但更硬。
  `engine=api` 时对应 `network_type=physical`；`engine=stringdb` 时对应
  `STRINGdb$new(network_type="physical")`（加载的是 `protein.physical.links.*` 文件）。
  两者数据文件不同，**不能混用缓存**。

## 5. add_nodes：要不要让 STRING 补邻居

参考代码没有这一步（`get_interactions` 只返回输入集合内部的边）。本 skill 的
`--add_nodes=N` 让 STRING 额外带 N 个最相关的外部蛋白进来，用途是
「给定几个种子蛋白，看它们的第一圈邻居」。

两个必须知道的后果：

1. 补进来的蛋白**没有 logFC**，颜色是灰色（`na.value`），正文里要说明；
2. 脚本会自动把 `--prune` 改成 `none`（否则这些邻居会被剪掉，白补），
   并在日志里打 `[note] --add_nodes>0 且未显式给 --prune，自动改为 prune=none`。

## 6. prune 与 keep_isolated：两处跟参考代码不同的默认

| | 参考代码 | 本 skill 默认 | 为什么 |
|---|---|---|---|
| 无 logFC 的节点 | `induced_subgraph(which(!is.na(logFC)))` 删掉 | `--prune=logfc`（**同参考代码**） | 保持一致；不给 logFC 时自动变 `none`，避免把图剪空 |
| 没有任何边的输入基因 | `graph_from_data_frame(边表)` 根本不会创建这些顶点 → **静默消失** | `--keep_isolated=1` 显式补回 | 少了节点在图上看不出来，属于静默错误 |

`--keep_isolated=1` 的代价：如果阈值过高导致一条边都没有，图里会出现一堆孤立点
（这时脚本会打 `[warn] 网络不连通，有 N 个孤立节点`）。

## 7. 配色与颜色区间

参考代码的配色是 7 段蓝→浅黄，并**硬编码** `limits = c(-2, 0)`：

```r
scale_colour_gradientn(colours = c("#08306B","#2171B5","#1F9BCD","#41B6C4",
                                   "#7FCDBB","#C7E9B4","#FFFFCC"),
                       limits = c(-2, 0))
```

本 skill 的 `--palette=nature` 就是这 7 个色（默认），但 `--color_limits` 默认是
`auto`（取实际范围）。**要跟参考图逐格对色，必须显式写 `--color_limits=-2,0`。**

另外参考代码的 `limits` 遇到超范围值会被 ggplot 默认的 `oob`（censoring）变成灰色 ——
本 skill 换成了 `squish_oob`（压到边界色），避免出现莫名的灰点。
`--add_nodes>0` 时邻居蛋白是真正的 `NA`，仍然是灰色，这是有意为之。

## 8. 布局

`layout = "stress"` 是 `graphlayouts` 包的 stress-majorization 布局，也是 Nature 那张图的观感来源。
参考代码有 `set.seed(42)`，本 skill 用 `--seed=42` 固定 —— **同一个图换个种子，
点和边的相对位置会变，拓扑不变**。想让图更"开"可以试 `fr` 或 `graphopt`；
小网络（<20 节点）`circle` 反而更清楚。

## 9. 别名修正的两阶段策略

参考代码的坑是「名字写错了 → `removeUnmappedRows=TRUE` → 静默少节点」。
本 skill 的做法不是「一上来就改名」，而是：

1. **先用原样名映射一遍**（避免把本来合法的 symbol 改坏）；
2. 只对**第一轮没映射上**的基因查内置别名表，换名后**再映射一轮**；
3. 两轮都失败的才进 `*_unmapped_genes.csv` 并大声告警。

所以 `--alias_fix=auto`（默认）是安全的：它只救失败的，不会动成功的。
`--alias_fix=all` 才是无条件替换（只在明确知道自己那批名字都过时时用），
`--alias_fix=none` 完全交给 STRING。

> 实测结论与参考代码笔记相反：在 **STRING 11.5** 里 `CoREST / LSD1 / MACROH2A2 / H2AX`
> 都能直接映射成功（13/13），别名表不会触发。换到 STRING 12 或别的物种时可能就不一样了 ——
> 以脚本的 `[alias]` 日志为准。

## 10. 与参考代码的逐条差异

| 位置 | 参考代码 | 本 skill | 影响 |
|---|---|---|---|
| 取数 | 只有 STRINGdb 包 | `--engine=api` 可绕过包 | 无网络障碍时结果一致 |
| score 标度 | 不需要处理（包内固定 0–1000） | 自动归一，防 API 0–1 | 不处理会导致边全被滤光 |
| 孤立节点 | 静默丢弃 | 默认补回 | 节点数会**多于**参考代码 |
| 无 logFC 节点 | 剪掉 | 同参考代码 | — |
| 颜色区间 | 硬编码 -2~0 | `auto`，可显式指定 | 配色深浅会不同 |
| 超范围取值 | ggplot 默认 → 灰点 | 压到边界色 | 出图更好看 |
| `library(tidygraph)` | 写了但没用 | 不依赖 | 少一个装不上的包 |
| `library(tidyverse)` | 写了 | 不依赖 | 只用到 igraph/ggraph/ggplot2 |
| logFC 与基因列表 | 两份手写列表，易错位 | 一份；个数不匹配直接报错 | 防静默错位 |
| 未映射基因 | 静默剔除 | 告警 + 落 CSV | — |
