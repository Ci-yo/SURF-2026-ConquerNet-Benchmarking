#!/usr/bin/env bash
set -euo pipefail

W8_ROOT="${W8_ROOT:-/home/xiyao/15data/W8}"
SCRIPT="${W8_ROOT}/code/03_run_week8_rf_fourclass.R"
RESULT_ROOT="${W8_ROOT}/results/03_rf_four_class_global_fdr"
LOG_ROOT="${W8_ROOT}/logs/03_rf_four_class_global_fdr_min51"
SUMMARY_DIR="${W8_ROOT}/results/04_rf_four_class_summary_min51"
FEATURE_DIR="${W8_ROOT}/results/05_rf_feature_exports_min51"

mkdir -p "${LOG_ROOT}" "${SUMMARY_DIR}"

TECHNOLOGIES=(m6A_SAC_seq Oxford_nanopore_rna004_DeepRM)
TASKS=(calibration_tp_fp error_fp_tn)

for technology in "${TECHNOLOGIES[@]}"; do
  for task in "${TASKS[@]}"; do
    log="${LOG_ROOT}/${task}.${technology}.log"
    echo "Running ${technology} / ${task}"
    Rscript "${SCRIPT}" "${technology}" "${task}" \
      --min-reliable=51 \
      --include-exploratory=false \
      --ntree=300 \
      --max-per-class=100000 \
      --nthread=8 \
      --smoke-n=0 \
      --output-root="${RESULT_ROOT}" \
      >"${log}" 2>&1
    tail -n 3 "${log}"
  done
done

Rscript "${W8_ROOT}/code/04_summarize_week8_rf_fourclass.R" \
  "${RESULT_ROOT}" "${SUMMARY_DIR}" \
  >"${LOG_ROOT}/04_summarize.log" 2>&1

Rscript "${W8_ROOT}/code/05_export_week8_rf_features.R" \
  --summary-dir="${SUMMARY_DIR}" \
  --output-dir="${FEATURE_DIR}" \
  --overwrite=true \
  >"${LOG_ROOT}/05_export.log" 2>&1

echo "Completed min51 additions and summaries."
