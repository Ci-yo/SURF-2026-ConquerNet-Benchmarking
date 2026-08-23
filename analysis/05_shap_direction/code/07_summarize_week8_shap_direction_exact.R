options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args)) args[[1L]] else {
  "/home/xiyao/15data/W8/results/06_shap_direction_exact"
}

models <- data.frame(
  task = c(
    rep("calibration_tp_fp", 6L),
    rep("error_fp_tn", 5L)
  ),
  technology = c(
    "eTAM_seq", "m6A_REF_seq", "m6ACE_seq", "MAZTER_seq",
    "MeRIP_seq", "scDART_seq",
    "eTAM_seq", "m6ACE_seq", "MAZTER_seq", "MeRIP_seq", "scDART_seq"
  ),
  stringsAsFactors = FALSE
)

read_one <- function(index) {
  task <- models$task[[index]]
  technology <- models$technology[[index]]
  directory <- file.path(output_root, task, technology)
  required <- c(
    "COMPLETE.ok", "shap_fold_directions.csv",
    "shap_direction_summary.csv", "retraining_qc.csv", "protocol.csv",
    "sampled_shap.rds"
  )
  missing <- required[!file.exists(file.path(directory, required))]
  data.frame(
    task = task,
    technology = technology,
    complete = !length(missing),
    missing_files = paste(missing, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
completion <- do.call(rbind, lapply(seq_len(nrow(models)), read_one))
write.csv(
  completion, file.path(output_root, "model_completion_audit.csv"),
  row.names = FALSE
)
if (!all(completion$complete)) {
  stop("Not all 11 exact SHAP models are complete. See model_completion_audit.csv")
}

read_table <- function(filename) {
  do.call(rbind, lapply(seq_len(nrow(models)), function(index) {
    read.csv(
      file.path(
        output_root, models$task[[index]], models$technology[[index]], filename
      ),
      check.names = FALSE
    )
  }))
}

direction_summary <- read_table("shap_direction_summary.csv")
fold_directions <- read_table("shap_fold_directions.csv")
retraining_qc <- read_table("retraining_qc.csv")
protocol <- read_table("protocol.csv")

if (nrow(direction_summary) != 110L) {
  stop("Expected 11 x 10 = 110 model-feature summary rows.")
}
if (nrow(retraining_qc) != 55L ||
    !all(retraining_qc$exact_reproduction_pass) ||
    max(retraining_qc$raw_probability_max_abs_error) > 1e-12 ||
    max(retraining_qc$prior_corrected_probability_max_abs_error) > 1e-12) {
  stop("Global exact-retraining QC failed.")
}

write.csv(
  direction_summary,
  file.path(output_root, "all_models_shap_direction_summary.csv"),
  row.names = FALSE
)
write.csv(
  fold_directions,
  file.path(output_root, "all_models_shap_fold_directions.csv"),
  row.names = FALSE
)
write.csv(
  retraining_qc,
  file.path(output_root, "all_models_retraining_qc.csv"),
  row.names = FALSE
)
write.csv(
  protocol,
  file.path(output_root, "all_models_protocol.csv"),
  row.names = FALSE
)

legend <- data.frame(
  task = c("calibration_tp_fp", "error_fp_tn"),
  negative_class = c("FP", "TN"),
  positive_class = c("TP", "FP"),
  positive_shap_meaning = c(
    "positive SHAP pushes prediction toward TP",
    "positive SHAP pushes prediction toward FP"
  ),
  negative_shap_meaning = c(
    "negative SHAP pushes prediction toward FP",
    "negative SHAP pushes prediction toward TN"
  ),
  stable_direction_rule = rep(
    "same non-zero direction in at least 4 of 5 chromosome folds", 2L
  ),
  stringsAsFactors = FALSE
)
write.csv(
  legend, file.path(output_root, "shap_direction_legend.csv"),
  row.names = FALSE
)

cat("Exact SHAP summary complete.\n")
cat("Models:", nrow(models), "\n")
cat("Stable feature directions:", sum(direction_summary$stable_direction_4of5),
    "/", nrow(direction_summary), "\n")
cat("Maximum raw reproduction error:",
    max(retraining_qc$raw_probability_max_abs_error), "\n")
cat("Maximum corrected reproduction error:",
    max(retraining_qc$prior_corrected_probability_max_abs_error), "\n")
