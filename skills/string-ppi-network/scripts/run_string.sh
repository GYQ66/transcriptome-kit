#!/usr/bin/env bash
# ===========================================================================
# run_string.sh —— string-ppi-network 技能的环境包装器
#
# 用法：
#   bash run_string.sh 01_string_network.R --genes=A,B,C --outdir=./string_out
#   bash run_string.sh 01_string_network.R --table=deg.csv --engine=api
#
# 环境变量：RSCRIPT（R 路径）、RTMP（Windows 临时目录）
# ===========================================================================
set -uo pipefail

SELF="${BASH_SOURCE[0]:-$0}"
case "$SELF" in */*) SELF_DIR="${SELF%/*}" ;; *) SELF_DIR="." ;; esac
SCRIPT_DIR="$(cd "$SELF_DIR" && pwd -P)"
# shellcheck disable=SC1091
_TOOLKIT_SCRIPTS_DIR="$SCRIPT_DIR"
source "$SCRIPT_DIR/../../_lib_r_env.sh"

TOOLKIT_PATH_KEYS=" table genes_file logfc_file positions outdir info_file links_file cache_dir log_file "
TOOLKIT_PATH_LIST_KEYS=" info_files links_files "
export TOOLKIT_PATH_KEYS TOOLKIT_PATH_LIST_KEYS

if [ -n "$_toolkit_is_win" ]; then
  export STRING_SKILL_DIR="$(winpath "$SCRIPT_DIR/..")"
else
  export STRING_SKILL_DIR="$(cd "$SCRIPT_DIR/.." && pwd -P)"
fi
export TOOLKIT_ENV_NOTE="TMP=$TMP"

toolkit_banner
toolkit_run "${1:-}" "${@:2}"
