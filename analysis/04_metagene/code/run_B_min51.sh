#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUTPUT_ROOT=${OUTPUT_ROOT:-/home/xiyao/15data/W8/results/reproduction_20260827_metagene}
CLASS_MATRIX=${CLASS_MATRIX:-/home/xiyao/15data/W8/results/01_four_class_global_fdr/four_class_matrix_Total_gt_10.rds}
BENCHMARK=/home/xichen/ConquerNet/ConquerBN_basic/m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali_IVT_matched.rds
FEATURES=/home/xiyao/15data/W8/input_data/features/m6Aconquer_omicsFeatures_hg38.rds
BASE=${OUTPUT_ROOT}

Rscript "${SCRIPT_DIR}/05_summarize_week8_fp_tp_meta_effect_sizes.R" \
  --class-matrix="${CLASS_MATRIX}" \
  --benchmark="${BENCHMARK}" \
  --features="${FEATURES}" \
  --reliability-min-per-class=51 \
  --output-dir="${BASE}/effect_sizes_min51"

Rscript "${SCRIPT_DIR}/06_plot_fp_tp_standard_metagene.R" \
  --class-matrix="${CLASS_MATRIX}" \
  --benchmark="${BENCHMARK}" \
  --features="${FEATURES}" \
  --min-class=51 \
  --output-dir="${BASE}/standard_metagene_min51"
