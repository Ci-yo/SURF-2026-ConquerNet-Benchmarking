#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/xiyao/15data/W8
BASE=${ROOT}/results/02_fp_tp_meta_distributions/Total_gt_10_no_cali_IVT_matched

Rscript "${ROOT}/code/05_summarize_week8_fp_tp_meta_effect_sizes.R" \
  --reliability-min-per-class=51 \
  --output-dir="${BASE}/effect_sizes_min51"

Rscript "${ROOT}/code/06_plot_fp_tp_standard_metagene.R" \
  --min-class=51 \
  --output-dir="${BASE}/standard_metagene_min51"
