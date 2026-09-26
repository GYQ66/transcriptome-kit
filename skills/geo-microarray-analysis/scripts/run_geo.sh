#!/usr/bin/env bash
# ===========================================================================
# run_geo.sh —— geo-microarray-analysis 技能的环境包装器
#
# 用法：
#   bash run_geo.sh 01_geo_normalize.R --gse=GSE62452 --matrix=... --soft=... --dir=./out
#   bash run_geo.sh 02_deg_plots.R --expr=... --group_file=... --contrast=case-control
#   bash run_geo.sh 03b_merge_from_geo.R --gse=GSE1,GSE2 --datadir=... --stage=normalize
#
# 环境变量：RSCRIPT（R 路径）、RTMP（Windows 临时目录，默认 C:/Rtmp）
# ===========================================================================
set -uo pipefail

SELF="${BASH_SOURCE[0]:-$0}"
case "$SELF" in */*) SELF_DIR="${SELF%/*}" ;; *) SELF_DIR="." ;; esac
SCRIPT_DIR="$(cd "$SELF_DIR" && pwd -P)"
# shellcheck disable=SC1091
_TOOLKIT_SCRIPTS_DIR="$SCRIPT_DIR"
source "$SCRIPT_DIR/../../_lib_r_env.sh"

TOOLKIT_PATH_KEYS=" matrix soft datadir outdir dir annot expr group_file clinical geo_dir geo_prefix geo_manifest orgdb_sqlite extra_inputs "
TOOLKIT_PATH_LIST_KEYS=" matrices softs clinical_files extra_inputs "
export TOOLKIT_PATH_KEYS TOOLKIT_PATH_LIST_KEYS

# 把 scripts 目录告诉 R 端（R 脚本靠它找公共库）
if [ -n "$_toolkit_is_win" ]; then
  export GEO_SKILL_SCRIPTS="$(winpath "$SCRIPT_DIR")"
else
  export GEO_SKILL_SCRIPTS="$SCRIPT_DIR"
fi
export TOOLKIT_ENV_NOTE="TMP=$TMP"

toolkit_banner
toolkit_run "${1:-}" "${@:2}"
