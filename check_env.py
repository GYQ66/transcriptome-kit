#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
transcriptome-kit 环境检查与 R 安装助手（纯标准库，Python >= 3.8）

用法：
    python check_env.py                     # 检查 R / Python / GMT 数据
    python check_env.py --check-packages    # 顺带用 Rscript 检查各技能的 R 包
    python check_env.py --simulate-missing  # 测试"未装 R"分支（不改动本机）
    python check_env.py --print-download-url  # 只解析最新 R 安装包直链，不下载

R 未安装时：
    1. 打印官方与国内镜像的下载地址；
    2. 询问是否自动下载（**默认否**；代理/Agent 必须把问题转给用户，不得代答）；
    3. 用户确认后下载安装包到用户目录，Windows/macOS 可再确认后静默/打开安装。
"""

from __future__ import annotations

import argparse
import glob
import io
import os
import platform
import re
import shutil
import subprocess
import sys
import urllib.request

try:  # Windows 控制台/重定向时的编码兜底
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

CRAN = "https://cran.r-project.org"
MIRRORS = {
    "cran": CRAN,
    "tuna": "https://mirrors.tuna.tsinghua.edu.cn/CRAN",
}

R_PAGES = {
    "Windows": "/bin/windows/base/",
    "Darwin": "/bin/macosx/",
}

CRAN_PKGS = [
    "ggplot2", "pheatmap", "ggrepel", "data.table", "RSpectra",
    "ggprism", "survival", "survminer", "ggpubr", "ragg",
    "igraph", "ggraph", "tidygraph", "graphlayouts", "svglite",
    "VennDiagram", "png", "patchwork", "circlize", "devtools", "BiocManager",
]
BIOC_PKGS = [
    "GEOquery", "limma", "Biobase", "sva",
    "clusterProfiler", "enrichplot", "DOSE", "fgsea",
    "org.Hs.eg.db", "ComplexHeatmap",
]


def toolkit_root() -> str:
    return os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------- R 探测

def detect_r() -> str | None:
    exe = "Rscript.exe" if os.name == "nt" else "Rscript"
    hit = shutil.which(exe)
    if hit:
        return hit
    pats = []
    if os.name == "nt":
        for drive in "CDEFGH":
            pats += [
                f"{drive}:/Program Files/R/R-*/bin/x64/Rscript.exe",
                f"{drive}:/Program Files/R/R-*/bin/Rscript.exe",
                f"{drive}:/R/*/R/bin/x64/Rscript.exe",
                f"{drive}:/R/R_down/R/bin/x64/Rscript.exe",
                f"{drive}:/Users/*/AppData/Local/Programs/R/R-*/bin/x64/Rscript.exe",
            ]
    elif platform.system() == "Darwin":
        pats = [
            "/usr/local/bin/Rscript",
            "/opt/homebrew/bin/Rscript",
            "/Library/Frameworks/R.framework/Resources/bin/Rscript",
        ]
    else:
        pats = ["/usr/bin/Rscript", "/usr/local/bin/Rscript", "/opt/R/*/bin/Rscript"]
    for pat in pats:
        for hit in sorted(glob.glob(pat)):
            if os.path.isfile(hit):
                return hit
    return None


def r_version(rscript: str) -> str:
    try:
        out = subprocess.run(
            [rscript, "--version"], capture_output=True, text=True, timeout=30
        )
        m = re.search(r"(\d+\.\d+\.\d+)", (out.stdout or "") + (out.stderr or ""))
        return m.group(1) if m else "未知版本"
    except Exception:
        return "未知版本"


# ---------------------------------------------------------------- R 包检查

CHECK_PKGS_R = r"""
need_cran <- c({cran})
need_bioc <- c({bioc})
miss <- function(p) !requireNamespace(p, quietly = TRUE)
m1 <- need_cran[sapply(need_cran, miss)]
m2 <- need_bioc[sapply(need_bioc, miss)]
cat('CRAN_MISS:', paste(m1, collapse = ','), '\n', sep = '')
cat('BIOC_MISS:', paste(m2, collapse = ','), '\n', sep = '')
cat('GH_MISS:', if (requireNamespace('gground', quietly = TRUE)) '' else 'gground', '\n', sep = '')
""".format(cran=", ".join(f"'{p}'" for p in CRAN_PKGS),
           bioc=", ".join(f"'{p}'" for p in BIOC_PKGS))


def check_packages(rscript: str) -> bool:
    print("\n[3/3] 检查 R 包依赖（首次运行会加载命名空间，约 10~60 秒）…")
    try:
        out = subprocess.run(
            [rscript, "--vanilla", "-e", CHECK_PKGS_R],
            capture_output=True, text=True, timeout=600,
        )
    except Exception as e:
        print(f"  ⚠ 无法运行 Rscript 检查包：{e}")
        return False
    text = (out.stdout or "")
    m_cran = re.search(r"CRAN_MISS:[ \t]*(.*)", text)
    m_bioc = re.search(r"BIOC_MISS:[ \t]*(.*)", text)
    m_gh = re.search(r"GH_MISS:[ \t]*(.*)", text)
    ok = True
    for label, m, install in [
        ("CRAN", m_cran, "install.packages"),
        ("Bioconductor", m_bioc, "BiocManager::install"),
    ]:
        miss = [p.strip() for p in (m.group(1).strip().strip("'\"") if m else "").split(",") if p.strip()]
        if miss:
            ok = False
            print(f"  ✗ 缺 {label} 包：{', '.join(miss)}")
            print(f"    安装：{install}(c({', '.join(repr(p) for p in miss)}))")
        else:
            print(f"  ✓ {label} 包齐全")
    gh = (m_gh.group(1).strip() if m_gh else "")
    if gh:
        ok = False
        print("  ✗ 缺 GitHub 包：gground")
        print("    安装：devtools::install_github('dxsbiocc/gground')")
    else:
        print("  ✓ GitHub 包（gground）齐全")
    return ok


# ---------------------------------------------------------------- 数据检查

def check_gmt() -> bool:
    data = os.path.join(
        toolkit_root(), "skills", "msigdb-pathway-gene-lookup", "data"
    )
    gmts = glob.glob(os.path.join(data, "*.gmt"))
    if gmts:
        print(f"[2/3] MSigDB GMT 数据：✓ {len(gmts)} 个集合文件（{data}）")
        return True
    print(f"[2/3] MSigDB GMT 数据：✗ 未找到（期望目录 {data}）")
    print("      富集/检索功能需要这批数据；请重新下载完整套件。")
    return False


# ---------------------------------------------------------------- 下载

def fetch(url: str, timeout: int = 30) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


def resolve_installer(osname: str, base: str) -> str | None:
    """从 CRAN(镜像) 页面解析最新安装包直链。"""
    page_url = base + R_PAGES[osname]
    try:
        html = fetch(page_url)
    except Exception as e:
        print(f"  ✗ 无法访问 {page_url}：{e}")
        return None
    if osname == "Windows":
        m = re.search(r"R-\d+(?:\.\d+)+-win(?:\.32)?\.exe", html)
    else:
        want = "-arm64" if platform.machine() == "arm64" else "-x86_64"
        m = (re.search(rf"R-\d+(?:\.\d+)+{want}\.pkg", html)
             or re.search(r"R-\d+(?:\.\d+)+\.pkg", html))
    if not m:
        print(f"  ✗ 页面上没解析到安装包文件名（{page_url}）")
        return None
    return page_url + m.group(0)


def download(url: str, dest_dir: str) -> str | None:
    fn = url.rsplit("/", 1)[-1]
    dest = os.path.join(dest_dir, fn)
    print(f"  下载 {url}")
    print(f"  ->   {dest}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r, open(dest, "wb") as f:
            total = int(r.headers.get("Content-Length") or 0)
            done, last_mb = 0, -1
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                f.write(chunk)
                done += len(chunk)
                mb = done >> 20
                if mb != last_mb:
                    last_mb = mb
                    if total:
                        print(f"\r  {done/1e6:8.1f} / {total/1e6:.1f} MB "
                              f"({done*100//total}%)", end="", flush=True)
                    else:
                        print(f"\r  {done/1e6:8.1f} MB", end="", flush=True)
        print()
        print(f"  ✓ 下载完成（{os.path.getsize(dest)/1e6:.1f} MB）")
        return dest
    except Exception as e:
        print(f"\n  ✗ 下载失败：{e}")
        if os.path.exists(dest):
            try:
                os.remove(dest)
            except OSError:
                pass
        return None


def ask(prompt: str, default_no: bool = True) -> bool:
    try:
        ans = input(f"{prompt} ").strip().lower()
    except EOFError:
        return False
    if not ans:
        return not default_no
    return ans in ("y", "yes")


# ---------------------------------------------------------------- 主流程

def offer_r_install(mirror_key: str, auto_yes: bool, print_url_only: bool) -> int:
    osname = platform.system()
    base = MIRRORS.get(mirror_key, mirror_key.rstrip("/"))
    if osname == "Linux":
        print("\nLinux 安装 R（用系统包管理器，无需下载安装器）：")
        print("  Debian/Ubuntu : sudo apt install r-base r-base-dev")
        print("  Fedora/RHEL   : sudo dnf install R")
        print("  Arch          : sudo pacman -S r")
        print(f"  更多：{base}/bin/linux/")
        return 1

    page = base + R_PAGES[osname]
    print(f"\nR 下载地址（{mirror_key}）：{page}")
    if mirror_key != "tuna":
        print(f"国内镜像（清华）        ：{MIRRORS['tuna'] + R_PAGES[osname]}")

    if print_url_only:
        url = resolve_installer(osname, base)
        print(f"最新安装包直链：{url or '（解析失败，请打开上方页面手动下载）'}")
        return 1

    prompt = "是否现在自动下载 R 安装包？[y/N] "
    if auto_yes or ask(prompt):
        dest_dir = os.path.join(os.path.expanduser("~"), "Downloads")
        if not os.path.isdir(dest_dir):
            dest_dir = os.getcwd()
        url = resolve_installer(osname, base)
        if not url:
            return 1
        dest = download(url, dest_dir)
        if not dest:
            return 1
        if osname == "Windows":
            print("\n下一步：双击运行安装包，或让我静默安装（需管理员授权弹窗）。")
            if auto_yes or ask("是否现在启动静默安装？[y/N] "):
                try:
                    subprocess.Popen([dest, "/VERYSILENT", "/NORESTART"])
                    print("  安装向导已启动；完成后重开终端，再跑一次本检查确认。")
                except Exception as e:
                    print(f"  ✗ 启动安装器失败：{e}\n  请手动运行 {dest}")
            else:
                print(f"  已取消安装。稍后手动运行：{dest}")
        else:  # macOS
            print("\n下一步：")
            if auto_yes or ask("是否现在打开安装向导？[y/N] "):
                try:
                    subprocess.run(["open", dest], check=False)
                except Exception as e:
                    print(f"  ✗ 打开失败：{e}\n  请手动双击 {dest}")
            else:
                print(f"  已取消。稍后双击 {dest} 安装。")
    else:
        print("  已跳过下载。安装 R 后重新运行本检查即可。")
    return 1


def main() -> int:
    ap = argparse.ArgumentParser(
        description="transcriptome-kit 环境检查与 R 安装助手")
    ap.add_argument("--check-packages", action="store_true",
                    help="检测到 R 后顺带检查各技能的 R 包依赖")
    ap.add_argument("--yes", "-y", action="store_true",
                    help='所有确认自动选"是"（仅代理场景、用户明确同意后使用）')
    ap.add_argument("--no-download", action="store_true",
                    help="R 缺失时只打印地址，不询问下载")
    ap.add_argument("--mirror", default="cran", help="cran | tuna | 自定义镜像根 URL")
    ap.add_argument("--print-download-url", action="store_true",
                    help="只解析最新安装包直链后退出（不下载）")
    ap.add_argument("--simulate-missing", action="store_true",
                    help="假装没装 R（用于测试缺 R 分支，不改本机）")
    args = ap.parse_args()

    print("=" * 64)
    print("transcriptome-kit 环境检查")
    print("=" * 64)

    ok = True

    # [1/3] R
    print("[1/3] 检查 R …")
    rscript = None if args.simulate_missing else detect_r()
    if rscript:
        ver = r_version(rscript)
        print(f"  ✓ Rscript：{rscript}（{ver}）")
        if args.print_download_url:
            # 已装 R 但用户只想看最新安装包直链
            return offer_r_install(args.mirror, True, True)
    else:
        print("  ✗ 未检测到 R（本套件 6 个技能中有 5 个依赖 R）")
        if args.print_download_url or not args.no_download:
            return offer_r_install(args.mirror, args.yes, args.print_download_url)
        return 1

    # [2/3] GMT 数据
    ok = check_gmt() and ok

    # [3/3] R 包
    if args.check_packages:
        ok = check_packages(rscript) and ok
    else:
        print("[3/3] R 包依赖：跳过（需要时加 --check-packages）")

    print()
    if ok:
        print("✓ 环境就绪，套件可用。")
        return 0
    print("⚠ 环境不完整：按上方提示补齐后重跑本检查。")
    return 1


if __name__ == "__main__":
    sys.exit(main())
