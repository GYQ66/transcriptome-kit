#!/usr/bin/env bash
# ===========================================================================
# run_enrich.sh —— enrichment-analysis 技能的环境包装器
#
# 用法：
#   bash run_enrich.sh 01_enrich_ora.R --deg=... --species=human --outdir=./enrich
#   bash run_enrich.sh 02_enrich_gsea.R --deg=... --gmt_sets=H --plot=TRUE
#   bash run_enrich.sh 04_enrich_summary.R --in_dir=./enrich --prefix=X
#   bash run_enrich.sh 03_enrich_plot.R --in_dir=./enrich --prefix=X --mode=even
#
# 环境变量：RSCRIPT（R 路径）、RTMP（Windows 临时目录）、ENRICH_GMT_DIR（GMT 目录）
# ===========================================================================
set -uo pipefail

SELF="${BASH_SOURCE[0]:-$0}"
case "$SELF" in */*) SELF_DIR="${SELF%/*}" ;; *) SELF_DIR="." ;; esac
SCRIPT_DIR="$(cd "$SELF_DIR" && pwd -P)"
# shellcheck disable=SC1091
_TOOLKIT_SCRIPTS_DIR="$SCRIPT_DIR"
source "$SCRIPT_DIR/../../_lib_r_env.sh"

TOOLKIT_PATH_KEYS=" deg gene_list gmt gmt_dir outdir in_dir output rules func_csv geo_dir geo_prefix geo_manifest orgdb_sqlite "
TOOLKIT_PATH_LIST_KEYS=" gmt gmt_files "
export TOOLKIT_PATH_KEYS TOOLKIT_PATH_LIST_KEYS

if [ -n "$_toolkit_is_win" ]; then
  export ENRICH_SKILL_SCRIPTS="$(winpath "$SCRIPT_DIR")"
else
  export ENRICH_SKILL_SCRIPTS="$SCRIPT_DIR"
fi
export TOOLKIT_ENV_NOTE="TMP=$TMP"

toolkit_banner
toolkit_run "${1:-}" "${@:2}"
