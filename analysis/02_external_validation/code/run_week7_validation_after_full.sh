#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/xiyao/15data/W7"
CODE="${ROOT}/code/no_cali_validation"
MODEL_DIR="${ROOT}/03_no_cali_full_model/Total_gt_10"
LABEL_DIR="${ROOT}/04_no_cali_tp_fp_labels/Total_gt_10"
IVT_DIR="${ROOT}/05_MeRIP_IVT_no_cali_TP_FP"
PERTURB_DIR="${ROOT}/06_perturbation_no_cali_TP_FP"
MODEL_FILE="${MODEL_DIR}/full_no_cali_Total_gt_10.conquernet.rds"
POSTERIOR_FILE="${MODEL_DIR}/node_state1_posterior.rds"
VITERBI_FILE="${MODEL_DIR}/viterbi_class_matrix.rds"

echo "Waiting for the formal no-cali full model."
while pgrep -f '[2]0_train_full_no_cali.R.*03_no_cali_full_model/Total_gt_10' >/dev/null
do
  sleep 60
done

if [[ ! -s "${MODEL_FILE}" || ! -s "${POSTERIOR_FILE}" || ! -s "${VITERBI_FILE}" ]]
then
  echo "Full-model process ended without all required outputs." >&2
  exit 1
fi

echo "Building formal no-cali TP/FP labels."
Rscript "${CODE}/21_build_no_cali_tp_fp_labels.R" \
  --model-dir="${MODEL_DIR}" \
  --output-dir="${LABEL_DIR}"

echo "Running Week7 B MeRIP IVT validation."
Rscript "${CODE}/22_analyze_merip_ivt_no_cali_tp_fp.R" \
  --labels="${LABEL_DIR}/selected_tp_fp_class_matrix.rds" \
  --output-dir="${IVT_DIR}"

echo "Running Week7 C perturbation response validation."
Rscript "${CODE}/23_analyze_perturbation_no_cali_tp_fp.R" \
  --labels="${LABEL_DIR}/selected_tp_fp_class_matrix.rds" \
  --output-dir="${PERTURB_DIR}"

touch "${ROOT}/WEEK7_B_C_NO_CALI_TOTAL_GT_10_COMPLETE"
echo "Week7 B/C completed."
