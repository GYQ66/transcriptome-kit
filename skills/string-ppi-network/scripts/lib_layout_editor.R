# ===========================================================================
# lib_layout_editor.R —— 交互式 PPI 布局编辑器（HTML）生成器
#
# 被 01_string_network.R source。生成一个自包含 HTML：
#   - 浏览器里按住节点拖动，连线实时跟随（原生 SVG + 原生 JS，零依赖、离线可用）
#   - 「导出坐标 CSV」下载 <prefix>_positions.csv（symbol,STRING_id,x,y）
#   - 把该 csv 用 --positions= 传回 01_string_network.R 重跑，即出最终 PDF+PNG
#
# 坐标体系（关键，不要动）：
#   R 侧 layout 坐标 y 向上；SVG y 向下。HTML 内部做了翻转：
#     X = (dx - BX) * S ;  Y = -(dy - BYT) * S
#   导出时原样反演回 R 坐标：  dx = X/S + BX ;  dy = BYT - Y/S
#   因此导出的 csv 与 R 的 create_layout("manual") 坐标完全同一体系，
#   「导出 -> 回读」是无损往返。
#
# 注意：模板全部用 R 原始字符串 r"---( )---"，内容里不允许出现 )---" 。
# ===========================================================================

esc_json <- function(s) {
  s <- as.character(s)
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\r", "", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  s
}

num_json_vec <- function(v) paste(trimws(formatC(v, format = "g", digits = 10)), collapse = ",")

write_layout_editor_html <- function(path, prefix, nodes, edges, x, y,
                                     cols, lims, meta, has_fc = TRUE) {
  stopifnot(nrow(nodes) == length(x), length(x) == length(y))

  # ---- 组 DATA JSON --------------------------------------------------------
  nj <- vapply(seq_len(nrow(nodes)), function(i) sprintf(
    '{"symbol":"%s","id":"%s","color":"%s","fc":%s}',
    esc_json(nodes$symbol[i]), esc_json(nodes$STRING_id[i]),
    esc_json(nodes$color[i]),
    if (is.finite(nodes$logFC[i])) sprintf("%.4f", nodes$logFC[i]) else "null"), "")
  ej <- vapply(seq_len(nrow(edges)), function(i)
    sprintf("[%d,%d,%.4f]", as.integer(edges$a[i]) - 1L, as.integer(edges$b[i]) - 1L,
            as.numeric(edges$w[i])), "")
  data_json <- sprintf(
    '{"prefix":"%s","hasfc":%s,"cols":[%s],"lims":[%.4f,%.4f],"meta":"%s","nodes":[%s],"edges":[%s],"x":[%s],"y":[%s]}',
    esc_json(prefix), if (isTRUE(has_fc)) "true" else "false",
    paste0("\"", esc_json(cols), "\"", collapse = ","),
    as.numeric(lims[1]), as.numeric(lims[2]), esc_json(meta),
    paste(nj, collapse = ","), paste(ej, collapse = ","),
    num_json_vec(x), num_json_vec(y))

  # ---- HTML 模板 ------------------------------------------------------------
  html <- r"---(<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>PPI 布局编辑器</title>
<style>
  :root { --line:#d9e0e8; --ink:#1f2937; --sub:#6b7280; --accent:#1d6fd1; }
  * { box-sizing: border-box; }
  body { margin:0; font-family: system-ui, "Microsoft YaHei", "Segoe UI", sans-serif;
         color:var(--ink); background:#f4f6f9; }
  header { background:#fff; border-bottom:1px solid var(--line); padding:14px 22px 12px; }
  h1 { font-size:17px; margin:0 0 4px; }
  .meta { color:var(--sub); font-size:12px; margin-bottom:8px; }
  .steps { font-size:13px; line-height:1.7; margin:0; color:#374151; }
  .steps b { color:var(--accent); }
  #legend { display:inline-flex; align-items:center; gap:6px; font-size:11px;
            color:var(--sub); margin-top:8px; }
  #legend .bar { width:140px; height:10px; border-radius:5px; border:1px solid var(--line); }
  #toolbar { padding:10px 22px; display:flex; gap:8px; flex-wrap:wrap; align-items:center;
             background:#fbfcfe; border-bottom:1px solid var(--line); }
  button { font:13px/1 system-ui, "Microsoft YaHei", sans-serif; padding:7px 14px;
           border:1px solid #c3ccd8; border-radius:6px; background:#fff; color:var(--ink);
           cursor:pointer; }
  button:hover { border-color:var(--accent); color:var(--accent); }
  button.primary { background:var(--accent); border-color:var(--accent); color:#fff; }
  button.primary:hover { opacity:.9; color:#fff; }
  #status { font-size:12px; color:var(--sub); margin-left:6px; }
  #wrap { padding:12px 22px 6px; }
  #net { width:100%; height:78vh; min-height:420px; background:#fff; border:1px solid var(--line);
         border-radius:8px; touch-action:none; cursor:default; }
  #net .nd { stroke:#2b2b2b; stroke-width:1.5; cursor:grab; }
  #net .nd.dragging { stroke:#e04f4f; stroke-width:3; cursor:grabbing; }
  #net .lb { font-size:24px; fill:#1f2937; pointer-events:none;
             font-family: system-ui, Arial, sans-serif; }
  #net .ed { stroke:#b7c0cb; stroke-linecap:round; pointer-events:none; }
  details { padding:4px 22px 22px; }
  summary { cursor:pointer; font-size:13px; color:var(--sub); }
  #csvbox { width:100%; height:120px; margin-top:6px; font:12px/1.5 Consolas, monospace;
            border:1px solid var(--line); border-radius:6px; padding:8px; }
</style>
</head>
<body>
<header>
  <h1>PPI 布局编辑器 — 拖动节点摆位，导出坐标回 R 出最终图</h1>
  <div class="meta" id="meta"></div>
  <p class="steps">
    <b>1.</b> 按住节点拖动（连线实时跟随），把网络摆到满意为止；
    <b>2.</b> 点「导出坐标 CSV」，下载 <span id="posname"></span>；
    <b>3.</b> 原命令追加 <code>--positions=&lt;该 csv 路径&gt;</code> 重跑，即出最终 PDF + PNG。
    坐标与 R 内部体系一致（y 方向已处理），导出即可直接使用；「重置布局」回到初始自动布局。
  </p>
  <div id="legend"><span>logFC</span><span id="lo"></span>
    <span class="bar" id="bar"></span><span id="hi"></span></div>
</header>
<div id="toolbar">
  <button class="primary" id="btn-export">导出坐标 CSV</button>
  <button id="btn-copy">复制到剪贴板</button>
  <button id="btn-reset">重置布局</button>
  <button id="btn-fit">适应视图</button>
  <span id="status">拖动节点开始调整</span>
</div>
<div id="wrap"><svg id="net" xmlns="http://www.w3.org/2000/svg">
  <defs>
    <pattern id="dots" width="50" height="50" patternUnits="userSpaceOnUse">
      <circle cx="1.4" cy="1.4" r="1.4" fill="#e4e9f0"></circle>
    </pattern>
  </defs>
</svg></div>
<details><summary>坐标 CSV（浏览器拦截下载时从这里手动复制）</summary>
  <textarea id="csvbox" readonly></textarea>
</details>
<script>
"use strict";
var DATA = @@DATA@@;
var N = DATA.nodes, E = DATA.edges, n = N.length;
var CANVAS = 1000, PAD = 100, R = 13, FS = 24;
var S, BX, BYT;
(function () {
  var xs = DATA.x, ys = DATA.y;
  var x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
  var y0 = Math.min.apply(null, ys), y1 = Math.max.apply(null, ys);
  if (x1 - x0 < 1e-9) x1 = x0 + 1;
  if (y1 - y0 < 1e-9) y1 = y0 + 1;
  S = CANVAS / (x1 - x0);
  BX = x0; BYT = (y0 + y1) / 2;
})();
function toSvgX(dx) { return (dx - BX) * S; }
function toSvgY(dy) { return -(dy - BYT) * S; }
function toDataX(px) { return px / S + BX; }
function toDataY(py) { return BYT - py / S; }
var P0 = [], P = [];
for (var i = 0; i < n; i++) {
  var q = { x: toSvgX(DATA.x[i]), y: toSvgY(DATA.y[i]) };
  P0.push({ x: q.x, y: q.y }); P.push(q);
}
var svg = document.getElementById("net");
var NS = "http://www.w3.org/2000/svg";
var grid = document.createElementNS(NS, "rect");
grid.setAttribute("fill", "url(#dots)");
var gEdges = document.createElementNS(NS, "g");
var gNodes = document.createElementNS(NS, "g");
svg.appendChild(grid); svg.appendChild(gEdges); svg.appendChild(gNodes);
var nodeEls = [], edgeEls = [], adj = [];
for (var i2 = 0; i2 < n; i2++) adj.push([]);
for (var k = 0; k < E.length; k++) {
  var ln = document.createElementNS(NS, "line");
  ln.setAttribute("class", "ed");
  ln.setAttribute("stroke-width", (2.5 + E[k][2] * 7).toFixed(2));
  edgeEls.push(ln); gEdges.appendChild(ln);
  adj[E[k][0]].push(k); adj[E[k][1]].push(k);
}
for (var i3 = 0; i3 < n; i3++) {
  var c = document.createElementNS(NS, "circle");
  c.setAttribute("class", "nd"); c.setAttribute("r", R);
  c.setAttribute("fill", N[i3].color); c.setAttribute("data-i", i3);
  var t = document.createElementNS(NS, "text");
  t.setAttribute("class", "lb"); t.setAttribute("text-anchor", "middle");
  t.textContent = N[i3].symbol;
  gNodes.appendChild(c); gNodes.appendChild(t);
  nodeEls.push({ c: c, t: t });
}
function redraw(i) {
  var p = P[i], el = nodeEls[i];
  el.c.setAttribute("cx", p.x); el.c.setAttribute("cy", p.y);
  el.t.setAttribute("x", p.x); el.t.setAttribute("y", p.y + R + FS * 0.95);
}
function redrawEdge(k) {
  var e = E[k], ln = edgeEls[k];
  ln.setAttribute("x1", P[e[0]].x); ln.setAttribute("y1", P[e[0]].y);
  ln.setAttribute("x2", P[e[1]].x); ln.setAttribute("y2", P[e[1]].y);
}
for (var i4 = 0; i4 < n; i4++) redraw(i4);
for (var k2 = 0; k2 < E.length; k2++) redrawEdge(k2);
function fit() {
  var a = Infinity, b = -Infinity, c2 = Infinity, d = -Infinity;
  for (var i = 0; i < n; i++) {
    a = Math.min(a, P[i].x - R - FS * 2.2); b = Math.max(b, P[i].x + R + FS * 2.2);
    c2 = Math.min(c2, P[i].y - R - FS * 1.2); d = Math.max(d, P[i].y + R + FS * 2.4);
  }
  if (!isFinite(a)) { a = 0; b = CANVAS; c2 = 0; d = CANVAS; }
  var w = (b - a) + 2 * PAD, h = (d - c2) + 2 * PAD;
  svg.setAttribute("viewBox", (a - PAD) + " " + (c2 - PAD) + " " + w + " " + h);
  grid.setAttribute("x", a - PAD); grid.setAttribute("y", c2 - PAD);
  grid.setAttribute("width", w); grid.setAttribute("height", h);
}
fit();
var dragging = -1, offX = 0, offY = 0;
function svgPt(e) {
  var pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY;
  return pt.matrixTransform(svg.getScreenCTM().inverse());
}
svg.addEventListener("pointerdown", function (e) {
  var t = e.target;
  if (t.classList && t.classList.contains("nd")) {
    dragging = +t.getAttribute("data-i");
    var p = svgPt(e);
    offX = p.x - P[dragging].x; offY = p.y - P[dragging].y;
    t.classList.add("dragging");
    if (svg.setPointerCapture) { try { svg.setPointerCapture(e.pointerId); } catch (err) {} }
    e.preventDefault();
  }
});
svg.addEventListener("pointermove", function (e) {
  if (dragging < 0) return;
  var p = svgPt(e);
  P[dragging].x = p.x - offX; P[dragging].y = p.y - offY;
  redraw(dragging);
  for (var j = 0; j < adj[dragging].length; j++) redrawEdge(adj[dragging][j]);
  e.preventDefault();
});
function endDrag() {
  if (dragging >= 0) {
    nodeEls[dragging].c.classList.remove("dragging");
    dragging = -1;
  }
}
svg.addEventListener("pointerup", endDrag);
svg.addEventListener("pointercancel", endDrag);
svg.addEventListener("pointerleave", endDrag);
function csvEscape(v) {
  v = String(v);
  if (/[",\r\n]/.test(v)) { v = "\u0022" + v.replace(/\u0022/g, "\u0022\u0022") + "\u0022"; }
  return v;
}
function buildCsv() {
  var rows = ["symbol,STRING_id,x,y"];
  for (var i = 0; i < n; i++) {
    rows.push(csvEscape(N[i].symbol) + "," + csvEscape(N[i].id) + "," +
      toDataX(P[i].x).toFixed(6) + "," + toDataY(P[i].y).toFixed(6));
  }
  return rows.join("\r\n");
}
function showCsv(s) { document.getElementById("csvbox").value = s; }
function setStatus(s) { document.getElementById("status").textContent = s; }
document.getElementById("meta").textContent = DATA.meta;
document.getElementById("posname").textContent = DATA.prefix + "_positions.csv";
if (!DATA.hasfc) { document.getElementById("legend").style.display = "none"; }
document.getElementById("bar").style.background =
  "linear-gradient(90deg," + DATA.cols.join(",") + ")";
document.getElementById("lo").textContent = DATA.lims[0];
document.getElementById("hi").textContent = DATA.lims[1];
showCsv(buildCsv());
document.getElementById("btn-export").addEventListener("click", function () {
  var csv = buildCsv(); showCsv(csv);
  var blob = new Blob([csv], { type: "text/csv;charset=utf-8" });
  var a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = DATA.prefix + "_positions.csv";
  document.body.appendChild(a); a.click();
  setTimeout(function () { URL.revokeObjectURL(a.href); a.remove(); }, 800);
  setStatus("已导出 " + DATA.prefix + "_positions.csv（若被拦截，从下方文本框复制）");
});
document.getElementById("btn-copy").addEventListener("click", function () {
  var csv = buildCsv(); showCsv(csv);
  var ta = document.getElementById("csvbox");
  ta.focus(); ta.select();
  var okv = false;
  try { okv = document.execCommand("copy"); } catch (err) { okv = false; }
  setStatus(okv ? "已复制到剪贴板" : "复制失败，请在文本框手动 Ctrl+C");
});
document.getElementById("btn-reset").addEventListener("click", function () {
  for (var i = 0; i < n; i++) { P[i].x = P0[i].x; P[i].y = P0[i].y; redraw(i); }
  for (var k3 = 0; k3 < E.length; k3++) redrawEdge(k3);
  fit(); setStatus("已重置为初始自动布局");
});
document.getElementById("btn-fit").addEventListener("click", function () {
  fit(); setStatus("已适应视图");
});
</script>
</body>
</html>)---"

  html <- gsub("@@DATA@@", data_json, html, fixed = TRUE)
  con <- file(path, open = "wb")
  on.exit(close(con), add = TRUE)
  writeLines(html, con, useBytes = TRUE)
  invisible(path)
}