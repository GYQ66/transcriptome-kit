#!/usr/bin/env bash
# ===========================================================================
# _lib_r_env.sh —— 本套件所有 run_*.sh 包装器共用的环境初始化库
#
# 它解决四件事，让 Windows / macOS / Linux 都能直接跑：
#   1. 找到 Rscript：顺序为
#        环境变量 RSCRIPT > PATH 里的 Rscript > Windows 常见安装位置
#      找不到就报清晰错误并退出。
#   2. Windows + 中文用户名/含空格路径：
#      - TMP/TEMP/TMPDIR 指向纯 ASCII 临时目录（R 的 tempdir() 在启动时固定，
#        指向含非 ASCII 的用户目录会变成坏路径）
#      - COMSPEC 显式设为 cmd.exe（否则 R 的 shell() 报 "'/c' not found"）
#      - unset LC_ALL/LANG/LC_*（LC_ALL=C.UTF-8 会让 Windows R 退化成 C locale，
#        含非 ASCII 字符的路径完全无法寻址，file.exists 静默失败）
#      - POSIX 路径（/d/xxx）统一转成 Windows 形式（D:/xxx），否则 Rscript 段错误
#   3. macOS / Linux：只清掉会干扰 R 的 locale 变量，其他不动。
#
# 各技能包装器在 source 本文件后需要：
#   - 定义 <PREFIX>_SKILL_SCRIPTS 环境变量（把 scripts 目录的 Windows 形式告诉 R 端）
#   - 定义本技能路径参数里"值一定是路径"的键（PATH_KEYS / PATH_LIST_KEYS）
# ===========================================================================

# --- 自愈 PATH：保证 dirname/cygpath/sed 可用（精简 Git 环境有时缺 coreutils）---
_fix_toolkit_path() {
  local d
  for d in "/usr/bin" "/bin" \
           "/c/Program Files/Git/usr/bin" "/c/Program Files/Git/bin" \
           /c/Users/*/.workbuddy/binaries/PortableGit/versions/*/usr/bin \
           /c/Users/*/.workbuddy/binaries/PortableGit/versions/*/bin ; do
    if [ -d "$d" ]; then
      case ":$PATH:" in
        *":$d:"*) ;;
        *) PATH="$d:$PATH" ;;
      esac
    fi
  done
  export PATH
}
_fix_toolkit_path

# --- 找 Rscript ------------------------------------------------------------
find_rscript() {
  # 1) 显式环境变量
  if [ -n "${RSCRIPT:-}" ] && [ -f "$RSCRIPT" ]; then
    printf '%s' "$RSCRIPT"; return 0
  fi
  # 2) PATH 里已有（macOS: /usr/local/bin/Rscript；Linux 通常在 PATH）
  if command -v Rscript >/dev/null 2>&1; then
    command -v Rscript; return 0
  fi
  # 3) Windows 常见安装位置（unquoted glob 经 word splitting 展开目录名）
  local pattern found
  for pattern in \
      "/c/Program[ ]Files/R/R-*/bin/x64/Rscript.exe" \
      "/c/Program[ ]Files/R/R-*/bin/Rscript.exe" \
      "/c/Program[ ]Files[ (x86)]/R/R-*/bin/x64/Rscript.exe" \
      "/d/R/*/R/bin/x64/Rscript.exe" \
      "/c/Users/*/AppData/Local/Programs/R/R-*/bin/x64/Rscript.exe" ; do
    for found in $pattern; do
      if [ -f "$found" ]; then printf '%s' "$found"; return 0; fi
    done
  done
  return 1
}

RSCRIPT_CAND="$(find_rscript)" || {
  echo "[run] 找不到 Rscript。" >&2
  echo "  请先安装 R（https://cran.r-project.org），或用环境变量指定：" >&2
  echo '    RSCRIPT="/path/to/Rscript" bash run_xxx.sh ...        # macOS/Linux' >&2
  echo '    RSCRIPT="C:/Program Files/R/R-4.4.1/bin/x64/Rscript.exe" ...  # Windows' >&2
  exit 1
}
RSCRIPT="$RSCRIPT_CAND"

# --- POSIX -> Windows 路径转换 ----------------------------------------------
_slash2win() { printf '%s' "$1" | sed -e 's|^/\([a-zA-Z]\)/|\1:/|' ; }
if command -v cygpath >/dev/null 2>&1; then
  winpath() { cygpath -w "$1" 2>/dev/null || printf '%s' "$1"; }
else
  winpath() { _slash2win "$1"; }
fi

# --- 临时目录 / COMSPEC / locale（仅 Windows 需要）---------------------------
RTMP="${RTMP:-}"
_case_uname="$(uname 2>/dev/null || echo unknown)"
if printf '%s' "$_case_uname" | grep -qiE 'mingw|msys|cygwin'; then
  # Windows：临时目录必须是纯 ASCII（中文用户名会把 tempdir() 搞坏）
  if [ -z "$RTMP" ]; then
    for _t in "C:/Rtmp" "D:/Rtmp" "$HOME/Rtmp"; do
      if [ ! -d "$_t" ]; then mkdir -p "$_t" 2>/dev/null || true; fi
      if [ -d "$_t" ]; then RTMP="$_t"; break; fi
    done
    # 都不行就退回系统临时目录（英文用户名的机器上直接可用）
    [ -z "$RTMP" ] && RTMP="${TMPDIR:-/tmp}"
  fi
  [ ! -d "$RTMP" ] && mkdir -p "$RTMP" 2>/dev/null || true
  export TMP="$RTMP" TEMP="$RTMP" TMPDIR="$RTMP"
  export COMSPEC="${COMSPEC:-C:\\WINDOWS\\system32\\cmd.exe}"
  unset LC_ALL LANG LC_CTYPE LC_COLLATE LC_MONETARY LC_TIME
else
  # macOS / Linux：只清 locale 干扰，临时目录交给系统
  unset LC_ALL LC_CTYPE LC_COLLATE LC_MONETARY LC_TIME 2>/dev/null || true
  : "${LANG:=en_US.UTF-8}"; export LANG
fi

# --- 把路径值参数统一转成 Windows 形式（仅 Windows；其他平台原样透传）--------
 toolkit_conv_arg() {
  local a="$1"
  if printf '%s' "$_toolkit_is_win" | grep -q 1; then :; fi
  if [ -n "$_toolkit_is_win" ]; then
    local opt key v
    if [[ "$a" == --*=* ]]; then
      opt="${a%%=*}"; key="${a#--}"; key="${key%%=*}"; v="${a#*=}"
      if [[ " ${TOOLKIT_PATH_LIST_KEYS:-} " == *" $key "* ]]; then
        local old_ifs="$IFS" p plist=()
        IFS=','
        for p in $v; do
          if [ -z "$p" ]; then plist+=(""); else plist+=("$(winpath "$p")"); fi
        done
        IFS="$old_ifs"
        printf '%s=%s' "$opt" "$(IFS=,; printf '%s' "${plist[*]}")"
        return
      fi
      if [[ " ${TOOLKIT_PATH_KEYS:-} " == *" $key "* ]]; then
        printf '%s=%s' "$opt" "$(winpath "$v")"; return
      fi
    fi
    if [ -e "$a" ]; then printf '%s' "$(winpath "$a")"; return; fi
  fi
  printf '%s' "$a"
}

toolkit_run() {
  local script="$1"; shift
  # 解析成物理绝对路径
  if [ -n "$_toolkit_is_win" ]; then
    # Windows：相对名先按"包装器所在目录"补全，再转 Windows 形式
    if [ ! -f "$script" ] && [ -f "$_toolkit_scripts_dir/$script" ]; then
      script="$_toolkit_scripts_dir/$script"
    fi
    script="$(winpath "$script")"
  else
    if [ ! -f "$script" ] && [ -f "$_toolkit_scripts_dir/$script" ]; then
      script="$_toolkit_scripts_dir/$script"
    fi
    if [ -f "$script" ]; then
      script="$(cd "$(dirname "$script")" && pwd -P)/$(basename "$script")"
    fi
  fi
  local args=() a
  for a in "$@"; do args+=("$(toolkit_conv_arg "$a")"); done
  "$RSCRIPT" --vanilla "$script" "${args[@]}"
  exit $?
}

# --- 平台标记（在 source 本文件之后由各包装器设置）---------------------------
case "$(uname 2>/dev/null || echo unknown)" in
  *MINGW*|*MSYS*|*CYGWIN*) _toolkit_is_win=1 ;;
  *) _toolkit_is_win="" ;;
esac

# --- 记录调用方 scripts 目录（toolkit_run 用它解析相对脚本名）----------------
_toolkit_scripts_dir="${_TOOLKIT_SCRIPTS_DIR:-}"

# --- 包装器收尾打印 ----------------------------------------------------------
toolkit_banner() {
  echo "[run] Rscript : $RSCRIPT"
  [ -n "${TOOLKIT_ENV_NOTE:-}" ] && echo "[run] env     : $TOOLKIT_ENV_NOTE"
}
