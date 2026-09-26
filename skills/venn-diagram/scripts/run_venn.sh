#!/usr/bin/env bash
# ===========================================================================
# run_venn.sh —— venn-diagram 技能的环境包装器
#
# 用法：
#   bash run_venn.sh --set="组名A=a.txt" --set="组名B=b.txt" --labels_confirmed=TRUE
#   bash run_venn.sh --table=deg.csv --by=change --by_levels=UP,DOWN
#
# 注意：退出码 3 = 等待确认组名、退出码 4 = 缺组名，两者都不是报错。
# 环境变量：RSCRIPT（R 路径）、RTMP（Windows 临时目录）
# ===========================================================================
set -uo pipefail

SELF="${BASH_SOURCE[0]:-$0}"
case "$SELF" in */*) SELF_DIR="${SELF%/*}" ;; *) SELF_DIR="." ;; esac
SCRIPT_DIR="$(cd "$SELF_DIR" && pwd -P)"
# shellcheck disable=SC1091
_TOOLKIT_SCRIPTS_DIR="$SCRIPT_DIR"
source "$SCRIPT_DIR/../../_lib_r_env.sh"

TOOLKIT_PATH_KEYS=" outdir table log_file "
TOOLKIT_PATH_LIST_KEYS=""
export TOOLKIT_PATH_KEYS TOOLKIT_PATH_LIST_KEYS

if [ -n "$_toolkit_is_win" ]; then
  export VENN_SKILL_SCRIPTS="$(winpath "$SCRIPT_DIR")"
else
  export VENN_SKILL_SCRIPTS="$SCRIPT_DIR"
fi
export TOOLKIT_ENV_NOTE="TMP=$TMP"

toolkit_banner
toolkit_run "01_venn.R" "$@"
