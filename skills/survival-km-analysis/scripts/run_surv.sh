#!/usr/bin/env bash
# ===========================================================================
# run_surv.sh —— survival-km-analysis 技能的环境包装器
#
# 用法：
#   bash run_surv.sh 01_km_survival.R --expr=... --clin=... --gene=TIMP1 --outdir=./km
#   bash run_surv.sh 02_km_batch.R --table=... --genes_file=... --outdir=./km
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

TOOLKIT_PATH_KEYS=" expr clin table outdir genes_file palette_file "
TOOLKIT_PATH_LIST_KEYS=" expr_files clin_files "
export TOOLKIT_PATH_KEYS TOOLKIT_PATH_LIST_KEYS

if [ -n "$_toolkit_is_win" ]; then
  export SURV_SKILL_DIR="$(winpath "$SCRIPT_DIR")"
else
  export SURV_SKILL_DIR="$SCRIPT_DIR"
fi
export TOOLKIT_ENV_NOTE="TMP=$TMP"

toolkit_banner
toolkit_run "${1:-}" "${@:2}"
