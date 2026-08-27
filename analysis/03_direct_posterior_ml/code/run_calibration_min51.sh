#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# Override OUTPUT_ROOT to place a reproduction elsewhere. The default is a new
# directory and does not overwrite the archived W8 results.
OUTPUT_ROOT=${OUTPUT_ROOT:-/home/xiyao/15data/W8/results/reproduction_20260827_direct_posterior}

MODEL=/home/xiyao/15data/W7/06_no_cali_IVT_matched_full_model/Global_no_additional_coverage/full_no_cali_global.conquernet.rds
DATASET=/home/xichen/ConquerNet/ConquerBN_basic/m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali_IVT_matched.rds
FEATURES=/home/xiyao/15data/W8/input_data/features/m6Aconquer_omicsFeatures_hg38.rds
RF_DESIGN=/home/xiyao/15data/W8/input_data/rf_design/tp_fp_tp_tn_metadata.rds

PREPARED=${OUTPUT_ROOT}/01_four_class_global_fdr
RUNS=${OUTPUT_ROOT}/03_rf_four_class_global_fdr
SUMMARY=${OUTPUT_ROOT}/04_rf_four_class_summary_min51
FEATURE_EXPORTS=${OUTPUT_ROOT}/05_rf_feature_exports_min51

mkdir -p "${PREPARED}" "${RUNS}" "${SUMMARY}" "${FEATURE_EXPORTS}"

Rscript "${SCRIPT_DIR}/01_generate_four_class_global_fdr.R" \
  "${MODEL}" "${PREPARED}" 0.05 10 0

CLASS_MATRIX=${PREPARED}/four_class_matrix_Total_gt_10.rds

mapfile -t TECHNOLOGIES < <(
  Rscript -e 'x <- readRDS(commandArgs(TRUE)[1]); cat(colnames(x), sep="\n")' \
    "${CLASS_MATRIX}"
)

TASKS=(calibration_tp_fp error_fp_tn artifact_tp_tn)

for technology in "${TECHNOLOGIES[@]}"; do
  for task in "${TASKS[@]}"; do
    Rscript "${SCRIPT_DIR}/02_run_rf_fourclass.R" \
      "${technology}" "${task}" \
      --class-matrix="${CLASS_MATRIX}" \
      --dataset="${DATASET}" \
      --features="${FEATURES}" \
      --rf-design="${RF_DESIGN}" \
      --output-root="${RUNS}" \
      --min-reliable=51 \
      --include-exploratory=false \
      --ntree=300 \
      --max-per-class=100000 \
      --nthread=8 \
      --smoke-n=0
  done
done

Rscript "${SCRIPT_DIR}/03_summarize_rf_fourclass.R" "${RUNS}" "${SUMMARY}"

Rscript "${SCRIPT_DIR}/04_export_rf_top20.R" \
  --summary-dir="${SUMMARY}" \
  --output-dir="${FEATURE_EXPORTS}" \
  --overwrite=true

echo "Completed direct-posterior calibration screening with min class >50."
