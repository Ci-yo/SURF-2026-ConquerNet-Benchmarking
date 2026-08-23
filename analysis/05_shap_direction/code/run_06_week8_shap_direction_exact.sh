#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/xiyao/15data/W8
CODE="$ROOT/code"
OUT="$ROOT/results/06_shap_direction_exact"
LOG="$ROOT/logs"
mkdir -p "$OUT" "$LOG"

RUNNER="$CODE/06_run_week8_shap_direction_exact.R"
SUMMARIZER="$CODE/07_summarize_week8_shap_direction_exact.R"

run_model() {
  local technology="$1"
  local task="$2"
  local complete="$OUT/$task/$technology/COMPLETE.ok"
  local model_log="$LOG/06_shap_exact__${task}__${technology}.log"
  if [[ -s "$complete" ]]; then
    echo "SKIP complete: $technology / $task"
    return 0
  fi
  echo "RUN: $technology / $task"
  Rscript "$RUNNER" "$technology" "$task" \
    --output-root="$OUT" \
    --nsim=16 \
    --per-class=100 \
    --background-per-class=100 \
    --nthread=8 \
    2>&1 | tee "$model_log"
}

## Mandatory full-five-fold exact smoke. Its valid output is already a final
## model result and is therefore not recomputed in the subsequent list.
echo "===== FULL 5-FOLD EXACT SMOKE ====="
run_model eTAM_seq calibration_tp_fp
test -s "$OUT/calibration_tp_fp/eTAM_seq/COMPLETE.ok"
echo "===== SMOKE PASSED; RUNNING 11-MODEL MANIFEST ====="

run_model eTAM_seq calibration_tp_fp
run_model m6A_REF_seq calibration_tp_fp
run_model m6ACE_seq calibration_tp_fp
run_model MAZTER_seq calibration_tp_fp
run_model MeRIP_seq calibration_tp_fp
run_model scDART_seq calibration_tp_fp

run_model eTAM_seq error_fp_tn
run_model m6ACE_seq error_fp_tn
run_model MAZTER_seq error_fp_tn
run_model MeRIP_seq error_fp_tn
run_model scDART_seq error_fp_tn

Rscript "$SUMMARIZER" "$OUT" 2>&1 | tee "$LOG/07_summarize_week8_shap_exact.log"
echo "Week8 exact SHAP workflow complete: $OUT"
