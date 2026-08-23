#!/usr/bin/env bash
set -euo pipefail

W8_ROOT="${W8_ROOT:-/home/xiyao/15data/W8}"
SCRIPT="${W8_ROOT}/code/03_run_week8_rf_fourclass.R"
CLASS_MATRIX="${W8_ROOT}/results/01_four_class_global_fdr/four_class_matrix_Total_gt_10.rds"
OUTPUT_ROOT="${OUTPUT_ROOT:-${W8_ROOT}/results/03_rf_four_class_global_fdr}"
LOG_ROOT="${LOG_ROOT:-${W8_ROOT}/logs/03_rf_four_class_global_fdr}"
TASKS="${TASKS:-calibration_tp_fp,error_fp_tn,artifact_tp_tn}"
ONLY_TECH="${ONLY_TECH:-}"
INCLUDE_EXPLORATORY="${INCLUDE_EXPLORATORY:-false}"
NTREE="${NTREE:-300}"
MAX_PER_CLASS="${MAX_PER_CLASS:-100000}"
NTHREAD="${NTHREAD:-8}"
SMOKE_N="${SMOKE_N:-0}"
DRY_RUN="${DRY_RUN:-0}"

mkdir -p "${OUTPUT_ROOT}" "${LOG_ROOT}"
if [[ ! -f "${SCRIPT}" ]]; then
  echo "Missing RF script: ${SCRIPT}" >&2
  exit 1
fi
if [[ ! -f "${CLASS_MATRIX}" ]]; then
  echo "Missing Week8 class matrix: ${CLASS_MATRIX}" >&2
  exit 1
fi

mapfile -t TECHNOLOGIES < <(
  Rscript -e 'x <- readRDS(commandArgs(TRUE)[1]); cat(colnames(x), sep="\n")' \
    "${CLASS_MATRIX}"
)
IFS=',' read -r -a TASK_ARRAY <<< "${TASKS}"

for task in "${TASK_ARRAY[@]}"; do
  for technology in "${TECHNOLOGIES[@]}"; do
    if [[ -n "${ONLY_TECH}" && "${technology}" != "${ONLY_TECH}" ]]; then
      continue
    fi
    log="${LOG_ROOT}/${task}.${technology}.log"
    cmd=(
      Rscript "${SCRIPT}" "${technology}" "${task}"
      "--include-exploratory=${INCLUDE_EXPLORATORY}"
      "--ntree=${NTREE}"
      "--max-per-class=${MAX_PER_CLASS}"
      "--nthread=${NTHREAD}"
      "--smoke-n=${SMOKE_N}"
      "--output-root=${OUTPUT_ROOT}"
    )
    printf 'COMMAND:'
    printf ' %q' "${cmd[@]}"
    printf '\n'
    if [[ "${DRY_RUN}" == "1" ]]; then
      continue
    fi
    "${cmd[@]}" >"${log}" 2>&1
    tail -n 2 "${log}"
  done
done

echo "Week8 RF task loop finished. Results: ${OUTPUT_ROOT}"
