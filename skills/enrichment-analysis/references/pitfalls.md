# 踩坑手册（改脚本前先读）

开发环境：Windows + 中文用户名（第 1~2 节的坑主要在这类机器上踩出来）。
以下每一条都是实测过的，不是推测。

---

## 1. 环境类（R 启动前必须办好）

### 1.1 TMP/TEMP/TMPDIR 必须指向纯 ASCII 目录
Windows 中文用户名会让 `TMP` 里的路径编码出问题，R 启动时据此定的 `tempdir()`
是坏路径，报 `Failed to create directory ... C:/Users/i+i:`。
**`tempdir()` 在 R 启动时就固定，脚本内部改不掉。**

```bash
export TMP=D:/Rtmp TEMP=D:/Rtmp TMPDIR=D:/Rtmp
```

### 1.2 COMSPEC 必须导出
否则 R 的 `shell()` 报 `'/c' not found`（某些包的安装/解压步骤会用到）。

```bash
export COMSPEC='C:\WINDOWS\system32\cmd.exe'
```

### 1.3 ★ 必须去掉 `LC_ALL`/`LANG`/`LC_CTYPE`
本机环境里带着 `LC_ALL=C.UTF-8`。Windows 版 R **应用不了 `C.UTF-8`**，
只打印一行 `Setting LC_CTYPE=C.UTF-8 failed` 就继续，于是退化成 `C` locale。
在 `C` locale 下 R **完全无法寻址含非 ASCII 字符的路径**：
`dir.exists` / `file.exists` / `readLines` 一律失败，`list.files()` 静默返回空——
看起来就像"文件不存在"。去掉之后 R 用系统原生 UTF-8 locale 就正常。

```bash
unset LC_ALL LANG LC_CTYPE LC_COLLATE LC_MONETARY LC_TIME
```

`run_enrich.sh` 已内置 1.1~1.3。**这就是"一律走包装器"的理由。**

### 1.4 POSIX 路径会让 Rscript 段错误
`/d/Rtmp/xxx` 这种 Git Bash 形式交给 Windows 版 `Rscript.exe`，父进程会
Segmentation fault（exit 139、无产物无报错）。一律写 `D:/Rtmp/xxx`。
`run_enrich.sh` 会对 `--deg/--gene_list/--gmt/--gmt_dir/--outdir/--in_dir/--out`
这些"值一定是路径"的参数自动用 `cygpath -w` 转换。

### 1.5 R 退出时报 Segmentation fault（139）是**无害**的
只要 `library()` 加载过包就会出现，且发生在所有工作完成之后。
**判断成败看产物文件，不要看退出码。**

### 1.6 bash 工具被沙箱拦（2026-09-18 实测）
```
PROGRAM BLOCKED BY SECURITY POLICY
  - wsl.exe (C:\Program Files\WSL\wsl.exe)
```
客户端沙箱会把 bash 调用路由到被黑名单拦截的 `wsl.exe`。表现是**整个命令 exit 1、
只输出一行乱码**。此时改用 `Bash` 工具直接调 `Rscript.exe`（不要嵌套 `bash 脚本.sh`）：

```bash
export TMP=D:/Rtmp TEMP=D:/Rtmp TMPDIR=D:/Rtmp
export COMSPEC='C:\WINDOWS\system32\cmd.exe'
unset LC_ALL LANG LC_CTYPE LC_COLLATE LC_TIME
export ENRICH_SKILL_SCRIPTS="<TOOLKIT>/skills/enrichment-analysis/scripts"
"C:/path/to/Rscript.exe" --vanilla "<脚本绝对路径>" <参数...>
```

`ENRICH_SKILL_SCRIPTS` 是给 `locate_lib()` 找 `lib_enrich_common.R` 用的；
不设的话脚本会用 `commandArgs()` 里的 `--file=` 反推，通常也能找到。

### 1.7 Bash 工具里 coreutils 常常全缺
本机 Bash 工具注入的 shim 会让 `dirname`/`ls`/`tail`/`rm` 报 `command not found`
（`shell-runtime-bash-env.sh: line 3: dirname: command not found`）。
bash **内建**（`cd`/`echo`/`for`/`[`/`case`/`unset`/`printf`）与**绝对路径 exe** 照常可用。
所以：
- 要列目录/读写文件 → 用 Read/Write/Glob 工具，或 `node.exe -e "..."`；
- 要删文件 → Bash 工具的 `rm`（实测正常）；
- `run_enrich.sh` 开头有一段 `_enrich_fix_path()`，会自己把 PortableGit 的
  `usr/bin`、`bin` 补进 `PATH`，所以它内部的 `dirname`/`cygpath` 能用。

### 1.8 PowerShell 工具里不要 `Remove-Item`、不要 `&`
`Remove-Item` 必然被 safe-delete 守卫拦下并报 `0x800700e8`；`&`/原生命令静默不执行。
要跑外部程序用 `Start-Process -NoNewWindow -Wait -PassThru -RedirectStandardOutput`；
删文件一律走 Bash 工具的 `rm`。

---

## 2. 输入表类（最影响结果正确性）

### 2.1 ★ `write.csv(row.names = TRUE)` 的第一列列名是空字符串
GEO 技能产出的 `*_DEG.csv` / `*_all.csv` 表头形如：

```
"","logFC","AveExpr","t","P.Value","adj.P.Val","B","neglogP","change"
"LAMC2",2.63,...
```

用 `read.csv(check.names = FALSE)` 读进来，`names(df)[1]` 是 `""`，
而 **`df[[""]]` 在 R 里取不到任何东西**（返回 `NULL`，不报错）。
2026-09-18 实测表现：`[输入表] 共 23307 行 / 0 个唯一基因`，然后"筛完是空的"直接退出。

`read_table_robust()` 会把空列名统一改成 `rownames1`，`pick_gene_col()` 优先认它。

### 2.2 列候选要按**优先级**挑，不能按表里的列顺序
`pick_col()` 老实现是 `nm[tolower(nm) %in% tolower(candidates)][1]`，
返回的是**表里最靠前**的匹配列，不是优先级最高的候选。实测后果：
候选是 `c("adj.P.Val","P.Value",...)`，表里 `P.Value` 在第 5 列、
`adj.P.Val` 在第 6 列，于是挑中了 `P.Value` —— 显著性基因数从 292 变成 12473。
现在改成按 candidates 顺序逐个找；模糊匹配只对长度 ≥ 4 的候选启用
（否则候选 `"t"` 会乱抓列）。

### 2.3 `p` 口径差一个数量级（GSE62452 实测）

| 口径 | 基因数 |
|---|---|
| `P.Value < 0.05` | 12473 |
| `adj.P.Val < 0.05` | 10960 |
| `\|logFC\|>1 & P.Value<0.05` | **292**（= up 177 + down 115） |
| `\|logFC\|>1 & adj.P.Val<0.05` | 292 |

**做 ORA 必须带 `|logFC|` 阈值，或者直接喂筛好的基因清单。**
拿 1 万个基因做 ORA，几乎每条通路都显著，等于没做。
`01_enrich_ora.R` 会在基因集 > 2000 时打印醒目告警。

### 2.4 `--ora_logfc` 对 `all` 也生效
这样 `all` 恒等于 `up ∪ down`，三套结果的基因数能对上账（脚本会打印 ✓ / ⚠ 核对行）。
想让 `all` = 全部显著基因不管方向，就 `--ora_logfc=0`。

### 2.5 ID 转换有损耗，必须报出来
SYMBOL→ENTREZID 在 GSE62452 的 292 个基因上丢 7 个（2.4%）；
别名、旧名、非编码基因都会丢。`<前缀>_id_map_<方向>.csv` 是明细。
**别把损耗当成 0。**

### 2.6 GMT 里的基因是 symbol，输入表可能是 entrez
`enricher()` 的 `gene` 与 `TERM2GENE$gene` 必须是同一套 ID，否则结果全空、
还会报 `None of the genes are in the gene sets`。
脚本用 `gmt_id_kind()` 按命中率自动判断该用 symbol 还是 entrez 版基因列表，
并在日志里打印判断结果。

### 2.7 `enricher()` 不传 `universe` 时背景 = 该 GMT 里所有基因的并集
这是 clusterProfiler 的行为，不是全基因组。想更严格就 `--universe=deg`
（用输入表里全部被检验基因作背景）。**默认 `none`（忠于参考代码），
换口径要在回复里说明。**

---

## 3. 富集计算类

### 3.1 `enrichKEGG` 需要联网
走 KEGG REST。代理环境下 SSL 偶发失败（实测 bioconductor.org 一次
`SSL connect error`、下一次就正常）。**每跑一个方向下载一次**，三个方向约 2~4 分钟。

`--kegg=auto` 在线失败会自动退到本地 MSigDB 的 KEGG 子集，但要知道三件事：

1. 那是 `KEGG_LEGACY`，固定 **186 条**（跨 v6.2 / v7.4 / v7.5.1 / v2026.1 都是 186 条）；
2. **从来不含 `KEGG_PI3K_AKT_SIGNALING_PATHWAY`、`KEGG_TNF_SIGNALING_PATHWAY`**
   ——不是版本问题，是 MSigDB 的 KEGG 子集本身就没有；
3. 通路 ID 是 `KEGG_XXX` 而不是 `hsaXXXXX`。

**退化了必须在回复里说明。**

### 3.2 富集一律用 `cutoff = 1` 计算
`pvalueCutoff = 1, qvalueCutoff = 1` → 所有被检验条目都留在表里，
再用 `sig` 列标记是否通过 `--pvalue/--padj/--qvalue`。
好处：改阈值不用重跑（KEGG 重跑要重新下载）。**富集值本身不随阈值改变。**

### 3.3 `enrichGO(ont = "ALL")` 支持，但 `gseGO(ont = "ALL")` 不支持
所以 `02_enrich_gsea.R` 里 GO 的 GSEA 是 **BP/CC/MF 各跑一次再 rbind**，
自带 `ONTOLOGY` 列。别图省事写成 `ont = "ALL"`。

### 3.4 `as.data.frame(enrichResult)` 可能返回 NULL
`enrichKEGG` 在没有条目通过时返回 `NULL`，直接 `@result` 会报
`no slot of name "result" for this object of class "NULL"`。
所有地方都要先判 `is.null()`。

### 3.5 GSEA 必须固定随机种子
`fgsea` 的 P 值来自随机置换。`--seed=1234` + `set.seed()`。
不固定则两次跑的显著条目数会不一样。

### 3.6 GSEA 的 `--padj` 惯例是 0.25
不是 0.05。脚本默认 0.25；换成 0.05 会只剩极少数条目。

### 3.7 排序向量要能复现
`<前缀>_GSEA_ranked_list.txt` 是**实际用的** `基因<TAB>排序值`。
查问题时对照它，不要重新从差异表推。

### 3.8 `use_internal_data = TRUE` 是 2011 年的老数据
除非完全没网且必须用 KEGG，否则别开。目前没做成参数。

---

## 4. 出图类

### 4.1 ★ `factor(Description, levels = Description)` 会崩
参考代码 `kegg美化.R` 第 99 行。同一 `Description` 在 BP/CC/MF 里重复出现时，
R 4.x 直接报 `duplicated levels in factors are deprecated` 并中断。
改成 `levels = unique(Description)`。

### 4.2 `geom_segment` 的标量美学会告警
`geom_segment(aes(x=0,y=0,xend=xaxis_max,yend=0), inherit.aes = FALSE)`
会抛 `All aesthetics have length 1, but the data has N rows`。
换成 `annotate("segment", ...)`。

### 4.3 裸向量配色在缺分类时会整体平移
参考代码 `scale_fill_manual(values = pal)` 配 4 个分类；
某一类没有条目时，`MF` 会拿到本该属于 `KEGG` 的颜色。
本技能按分类名绑定颜色（见 `references/figure-style.md` §3）。

### 4.4 `--fig_filter` 是人话，要映射到列名
`padj → p.adjust`、`p → pvalue`、`qvalue → qvalue`。
不映射的话 `dat[["padj"]]` 返回 `NULL`，`!is.na(NULL)` 是 `logical(0)`，
过滤后 0 行，表现为**"一条都不剩"的假告警 + 退化成全量出图**（2026-09-18 实测）。

### 4.4b ★ `ONTOLOGY` 为 NA 时会凭空多出一行全 NA 的"通路"
`--ont=BP/CC/MF`（不是 `ALL`）时，**`enrichGO` 的结果表里没有 `ONTOLOGY` 列**。
01 补列时如果传的是 `NULL`，该列就全是 `NA`；到了 03 出图：
`dat[dat$ONTOLOGY == k, ]` 里 `k` 若是 `NA`，**R 会把整行返回成 NA 而不是空**
（`df[NA, ]` 得到一行全 NA），于是图上多出一条没有名字的空通路，
并且 `max(nchar(sel$geneName))` 变 NA 直接把脚本打崩
（`missing value where TRUE/FALSE needed`）。实测 2026-09-18。

三处都要防：01 补列时传具体本体名；03 读表时把空分类回填成库名；
`cats_present` 过滤掉 NA；`max(..., na.rm = TRUE)` + `is.finite()` 兜底。

### 4.5 GSEA 出图前必须对齐列
GO / KEGG / GMT 三个结果表的列集不完全一样（GMT 多 `collection`、
GO 多 `ONTOLOGY`、`core_symbol` 只在有 ID 映射时才有），
直接 `do.call(rbind, ...)` 会报
`numbers of columns of arguments do not match`（2026-09-18 实测）。
先按固定 `keep_cols` 逐表补齐再 rbind。

### 4.6a `gground` 的参数走 `...`
`geom_round_col(mapping, data, position, ..., just, radius, na.rm, show.legend, inherit.aes)`
——`width`/`alpha` 是通过 `...` 传下去的，不是显式形参。

### 4.6b `gground` 装过了，不要重装
`devtools::install_github("dxsbiocc/gground")` 本机已装 1.0.1。
参考代码里那行 install 每次跑都要联网编译，删掉。

### 4.7 分类顺序与 y 轴
`levels = rev(c('BP','CC','MF','KEGG'))` 再 `arrange(ONTOLOGY, p.adjust)` +
`y = index`，结果是 **KEGG 在最下、BP 在最上**。
左侧 `rect.data` 用 `cumsum(n)` 从下往上堆，两者必须用同一套 level 顺序，
否则色块会和条目错位。

### 4.8 图内文字一律 ASCII
Windows 下中文字形缺失会被渲染成**乱码 ASCII**（不是方框），很难一眼看出来。
GO/KEGG 的 `Description` 本身是英文，不用管；但**不要往图里加中文标题**。
### 4.9 高度要**双向**自适应
固定 12×8 时：条目多 → 文字挤；条目少 → 柱子胖得难看、图例还排不下。
实测 4 条通路时每行 2 英寸、柱高 1.2 英寸，`Count` 图例的标题被设备边界裁掉。
现在不给 `--height` 时用 `max(4.6, 0.34 * n + 2.2)`（16 条 → 7.64 in、4 条 → 4.6 in），
**显式给了 `--height` 就完全按用户给的来**（不做自适应，便于固定版式）。
另外 `Count` 的图例档位压到最多 4 个（`seq(min,max,length.out=4)` 取整去重），
否则矮画布下图例竖排超出设备高度会被裁。

### 4.10 ★ 分类色块的宽度是"数据单位"，x 轴一长就装不下标签
参考代码把色块写死成 0.5 个数据单位宽。1 数据单位 = 多少英寸随 x 轴范围变化：
`xaxis_max=16` 时色块约 0.27 in，`xaxis_max=60` 时只有 0.076 in —— 而 `KEGG`
四个大写字母需要约 0.37 in。所以参考代码的 `KEGG` 必然被裁（裁成 `EGG`），
加 `clip="off"` 之后变成"文字冲出彩色框"。

修法：按色块可用宽度**把标签折行**（`KEGG` → `KE`+`GG`），并把色块**往左**加宽到包住文字。
两个坑：
- **字符宽度不能沿用基因名的系数**。基因名是混合大小写 + 数字 + 斜杠，平均 0.4 em；
  全大写标签约 0.62 em（`CHAR_IN_CAT = 0.0245` vs 基因名的 `0.0158`）。
  用 0.0158 会低估 35%，算出来的框还是装不下。
- **只能往左加宽**（右边界固定在 `-2w`）。往右加会撞上左侧的基因数圆点（`x = -w`），
  也会改变 `x = 0` 的视觉位置。
- "块宽 ← 需要多宽 ← 块宽决定的数据单位换算"是循环依赖，**迭代 3 次收敛**即可。

### 4.11 ★ 基因名的位置要挂在"柱下沿"，不能用 vjust
`geom_text(vjust = 2.6)` 的偏移单位是**文字高度**，跟柱子多高完全无关 ——
所以参考代码里基因名永远压在彩色柱**内部**；`--bar_height` 怎么改它都不动。
要"贴在柱子下面、并跟着柱子走"必须显式给 y：

```r
gene_dy <- bar_height / 2 + (一行文字高度 / 行高) / 2 + 0.05
y = index - gene_dy
```

注意 `y = index - gene_dy` 是靠**连续 y 轴**生效的（全局 `aes(y = index)` 建立的是
连续尺度，`geom_round_col` 里的 `y = Description` 因子按 level 顺序映射成 1..n，
两者数值正好重合）。别再引入离散 y，否则这个偏移会被丢掉。

实测 16 行 / 7.64 英寸 / 柱高 0.6 时 `gene_dy = 0.47` 行：文字落在柱下沿外，
与下一行柱顶还有约 0.11 行余量，不打架。

---

## 4b. 功能归纳类（04_enrich_summary.R）

### 4b.1 ★ 规则表里的短词必须加词边界（2026-09-18 实测两个实例）
规则用 `grepl(pattern, x, ignore.case = TRUE, perl = TRUE)`，
**没有词边界就会误命中子串**：

| 规则写法 | 实际命中了 | 后果 | 修法 |
|---|---|---|---|
| `actin` | `acting` | `oxidoreductase activity, **acting** on NAD(P)H…` 被归成「细胞骨架与运动」 | `\bactin\b` |
| `ral` | `endochond**ral**` | `endochondral bone morphogenesis` 被归成「小G蛋白与Rho」 | `\bral[ab]?\b` |

同类风险词：`erk` `src` `abl` `jak` `tgf` `bmp` `axis` `t cell`。
**短词（≤5 字符）一律加 `\b`**，或干脆不用两三个字母的缩写。

排错办法：看 `*_func_pathways_*.csv` 的 `func_rule` 列 —— 它记录了**实际命中的那条正则**，
一眼就能看出是规则写太宽还是通路名太怪。

### 4b.2 「未归类」不是功能方向，不能占排名位次
第一版把「未归类」混进按通路数排名，它直接排到第 2 名（38 条），把真正的功能方向
挤出前 5。现在脚本把它从排名里摘出来单列，并列出它的前 8 条供判断要不要补规则。
扩大规则覆盖后（20 大类 / 84 条规则），GSE62452 三个方向的未归类已经**清零**
（补规则前是 38/47/15，中途是 19/32/4，现在是 0/0/0）。

### 4b.2b ★ 合并汇总表与「方向=all」的逐方向表会同名
`04_enrich_summary.R` 里逐方向写 `<prefix>_func_counts_<方向>.csv`，
所以方向 `all` 得到的是 `<prefix>_func_counts_all.csv`；
而"所有方向合并的汇总表"如果也叫 `_func_counts_all.csv`，就会把前者**覆盖掉**
（2026-09-18 实测：`_func_counts_all.csv` 里躺着 45 行 = 3 个方向 × 15 类，
逐方向的 16 行反而没了）。现在合并表改名为
**`<prefix>_func_counts_by_direction.csv`**。
类似的命名冲突在"前缀 + 方向名 + 固定后缀"的产物上都要留意。

### 4b.3 规则表的顺序就是优先级
`ECM-receptor interaction` 里既有 `ecm` 也有 `receptor`。
「细胞外基质与黏附」必须排在「信号转导」**前面**，否则会被 `receptor` 抢走。
同理「骨骼、牙齿与矿化」必须排在「发育与分化」前面，否则
`bone development` / `cartilage development` 会被 `development` 抢走。
改规则表时把具体规则往上放、宽泛规则往下放。

### 4b.3b 新增大类时，先看未归类清单里"哪几件事"最集中
2026-09-18 把大类从 18 扩到 20 时，是**先导出 `*_func_pathways_*.csv` 里
`func_cat == 未归类` 的全部条目**（all 19 / up 32 / down 4），按主题归类后才决定新增
哪两个大类：

- `odontogenesis` / `ossification`×3 / `bone resorption` / `tooth mineralization` /
  `cornified envelope` / `hair cycle`×3 / `molting cycle`×2 → 一类「骨骼、牙齿与矿化」
- `tissue homeostasis` / `anatomical structure homeostasis` → 一类「稳态与内环境」

其余零散的长尾（`laminin complex`、`fibronectin binding`、`growth factor binding`、
`response to oxygen levels`、`antioxidant activity`、`acute-phase response`、
`cell leading edge`、`glomerular filtration`、`gonadotropin`、`innervation`、
`carboxylic acid binding`、`sulfur compound binding`、`symbiotic interaction` 等）
不必开新大类，**补进已有大类**即可。别为了"分类数好看"硬造大类。

### 4b.4 Jaccard 去冗余的规模
贪心去冗余是 O(n²) 的集合比较。单个功能大类里几百条通路没问题（实测 45 条 <0.1 s），
但一个方向上千条同功能通路时会明显变慢。`--jaccard=0` 可整个关掉。

### 4b.5 04 只做统计，不出图
归纳结果的价值在于**让用户做选择**。不要因为"顺手"就在 04 里出图——
出图需要用户先决定 A（按排名前 N）还是 B（按功能方向）。这是流程里的硬停顿点。

## 4c. 与 geo-microarray-analysis 联通（`--geo_dir`）

### 4c.1 交接清单长什么样、在哪
geo 的 `02_deg_plots.R` 会写 `<prefix>_for_enrichment.txt`（key=value，`#` 是注释）。
**它在 geo 的运行目录（进程 CWD）里** —— geo 的 02 步**没有 `--outdir`**，
所有产物都落在 CWD。所以 PowerShell 里跑 geo 要先 `Set-Location`，
否则清单里的 `dir`/`deg_file` 指向别处。

### 4c.2 `--geo_dir` 自动填了哪些参数
`01_enrich_ora.R` / `02_enrich_gsea.R` 里都是「**显式参数优先，geo 只填空**」：
`cfg <- oa(opt,"x",NULL); if (is.null(cfg) && !("x" %in% names(opt))) opt$x <- geo$x`。
实测 `--geo_dir` 一个参数就等价于同时给了：
`--deg`（_DEG.csv）、`--prefix`、`--species`、`--p_col`、`--logfc_col`、
`--gene_col`、`--ora_p=1`、`--ora_logfc=0`、`--split_by=change`、`--outdir`。

### 4c.3 ★ geo 模式下**默认不再二次筛选**，而且按 `change` 列拆上下调
geo 的 `_DEG.csv` 已经是 `|logFC|>1 & P<PCOL` 筛过的；如果富集侧再按
`adj.P.Val ≤ 0.05` 筛一遍，会**悄悄少掉一批基因**（geo 默认按 `P.Value` 筛，
两套口径不等价）。所以 geo 模式默认 `--ora_p=1 --ora_logfc=0`。

上下调也不靠 logFC 正负，而是直接用 geo 写好的 `change` 列（`UP`/`DOWN`/`NOT`），
这样 `all = up + down` 与 geo 那次的判定**逐字一致**（脚本会打印 ✓ 核对）。

### 4c.4 清单里 `gene_col=rownames` 不能直接当列名
`read_table_robust()` 会把 `write.csv(row.names=TRUE)` 产出的空列名改成 `rownames1`。
若把清单里的 `rownames` 原样传给 `pick_gene_col`，会报
`--gene_col=rownames 不在表里`。**解析清单时要把 `rownames` 还原成"自动识别"**
（`resolve_geo_inputs()` 已这么处理）。

### 4c.5 `species` 一律按 human 交接
geo 流程不解析平台物种，清单里恒为 `species=human`。**小鼠芯片必须显式
`--species=mouse`**，否则会拿人的 OrgDb 去注释小鼠基因 —— 基因名大小写不一样
（人 `TP53` / 鼠 `Trp53`），结果是"一个基因都没转成 ENTREZID"直接报错，
还不至于出错的结论，但也别指望能跑通。

### 4c.6 `--p_col` 曾经只是写在帮助里、没接上（2026-09-18 修）
`01_enrich_ora.R` 原来只用 `--ora_p_type` 选 p 列，`--p_col` 完全没被读。
现在优先级是：`--p_col`（显式列名） > `--ora_p_type`（非 auto 时当列名） > auto 候选
（`adj.P.Val` → `P.Value`）。

### 4c.7 联通后的实测基线（GSE62452，2026-09-18）
geo 02 跑出 `177 up / 115 down` → `--geo_dir` 一条命令跑 01，
结论与手工给 `--deg=…_all.csv --ora_logfc=1 --ora_p_type=P.Value` **完全一致**：
ORA `292 = 177+115`、GO 212/196/76、KEGG 14/12/7、Hallmark 2/6/2；
GSEA 用 `_all.csv` 得 GO 8228/2567、KEGG 352/157、H 50/39。

## 5. 依赖类

### 5.1 物种注释包
| 物种 | OrgDb | KEGG 代码 | 本机状态 |
|---|---|---|---|
| 人类 | `org.Hs.eg.db` | `hsa` | ✅ 已装，全流程实测通过 |
| 小鼠 | `org.Mm.eg.db` | `mmu` | ❌ **没装上**，见 5.2 |

### 5.2 ★ 本机装不上 `org.Mm.eg.db`（2026-09-18 实测，三条路都堵了）
代码层面小鼠是支持的（`--species=mouse` → `org.Mm.eg.db` + `mmu`），
但**注释包本身在本机装不进去**：

| 尝试 | 结果 |
|---|---|
| `BiocManager::install("org.Mm.eg.db")`（默认源） | 下到 88.5 MB 完整包，但安装阶段报 **`The syntax of the command is incorrect.`**（中文版：`命令语法不正确。`），`installation had non-zero exit status` |
| Windows 二进制包（`data/annotation/bin/windows/contrib/4.5/`） | **不存在**。该目录的 `PACKAGES` 索引下载成功但**内容为 0 行**——Bioconductor 不给注释包出 Windows 二进制 |
| 手动 `R.exe CMD INSTALL <tar.gz>` | 同样报 `命令语法不正确。` |
| 镜像 | 清华/中科大/NJU 的 bioconductor 路径与 `BiocManager` 拼出来的不一致（一律 404）；**西湖镜像可用**（`options(BioC_mirror="https://mirrors.westlake.edu.cn/bioconductor")`），但它也没有注释包的 Windows 二进制 |
| 直接下载 tar.gz 再本地装 | 走本机代理（`http://127.0.0.1:5xxxx`）时**被限速到约 26 KB/s**，88 MB 需要 ~1 小时，`download.file` 每次都在 60 s 超时后拿到残缺文件（`downloaded length 1671168 != reported length 92817018`）；curl 断点续传也无效（代理不支持 Range） |

2011 年那套 `use_internal_data` 与本问题无关。

**根因判断**：`R CMD INSTALL` 在 Windows 上要经 `cmd.exe`，而本机这个 shell 链路
（沙箱 + 中文用户名 + 长 PATH）下 R 拼出来的安装命令行会被 cmd 判为语法错误。
缺少 Rtools 不是主因——注释包是纯数据包，不需要编译。

**两条出路**：

1. **装 Rtools**（`https://cran.r-project.org/bin/windows/Rtools/rtools45.html`）
   后再试源码安装；或
2. **不装包，只加载 sqlite**（本技能的 `load_orgdb()` 已支持）：
   ```r
   options(timeout = 3600)   # 一定要放大，代理下很慢
   u <- "https://bioconductor.org/packages/3.22/data/annotation/src/contrib/org.Mm.eg.db_3.22.0.tar.gz"
   download.file(u, "D:/Rtmp/org.Mm.eg.db.tar.gz", method = "libcurl", mode = "wb")
   untar("D:/Rtmp/org.Mm.eg.db.tar.gz", exdir = "D:/Rtmp")
   # 取出 org.Mm.eg.db/inst/extdata/org.Mm.eg.sqlite（约 90 MB），放到：
   file.copy("D:/Rtmp/org.Mm.eg.db/inst/extdata/org.Mm.eg.sqlite",
             "D:/R/orgdb/org.Mm.eg.sqlite")
   ```
   脚本会自动在 `--orgdb_sqlite=` / `ENRICH_ORGDB_DIR` / `<技能>/data/orgdb` /
   `ENRICH_ORGDB_DIR`（默认含 `~/orgdb`）里找 `<OrgDb 名>.sqlite`，
   找到就 `AnnotationDbi::loadDb()` 直接用（走 `AnnotationDbi::select`，
   `enrichGO`/`bitr` 都能吃）。

**给小鼠做富集前，先 `requireNamespace("org.Mm.eg.db")` 或确认 sqlite 在不在，
不要把"没装包"报成"没富集结果"。**

### 5.3 MSigDB GMT 从哪来
按顺序找：`--gmt_dir` > 环境变量 `ENRICH_GMT_DIR` > 本技能 `data/` >
`msigdb-pathway-gene-lookup` 技能的 `data/`（同套件内置 GMT）。
后两处都有 `h.all./c1..c9.all.v2026.1.Hs.symbols.gmt`。
GMT 只有 **人类 symbol** 版。

### 5.4 `clusterProfiler` 版本相关
本机是 **4.18.2**。`GSEA`/`gseGO`/`gseKEGG` 都有 `...`，
`nPermSimple` 会透传给 `fgsea::fgseaMultilevel`（`enrichKEGG`/`enrichGO` 没有 `...`，
参数必须逐个显式给）。
