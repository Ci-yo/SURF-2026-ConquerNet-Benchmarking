## Leave-one-technology-out benchmark
## Each fold removes one technology, retrains ConquerBN with the remaining
## 13 technologies, and evaluates against the held-out BBmix MAP call.

# Locate script ----------------------------------------------------------

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Please run this file with Rscript.")
SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", script_arg[[1]]), winslash = "/"
))

# Load settings and functions -------------------------------------------

source(file.path(SCRIPT_DIR, "config.R"))
source(file.path(SCRIPT_DIR, "common_functions", "common_functions.R"))
source(file.path(SCRIPT_DIR, "common_functions", "conquernet_model_functions.R"))
source(file.path(SCRIPT_DIR, "common_functions", "benchmark_functions.R"))

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(S4Vectors)
  library(bnlearn)
  library(extraDistr)
  library(matrixStats)
  library(pROC)
  library(PRROC)
})

# Output folders ---------------------------------------------------------

output_dir <- ensure_directory(file.path(result_dir, "benchmark", "LOTO"))
model_dir <- ensure_directory(file.path(output_dir, "models"))
qc_dir <- ensure_directory(file.path(result_dir, "benchmark", "QC"))

# Load benchmark data ----------------------------------------------------

require_files(c(benchmark_file, loto_idr_dir))
dataset <- load_benchmark_dataset(CONFIG)

technologies <- colnames(dataset)
n_sites <- nrow(dataset)
total <- assay(dataset, "Total")
probability <- assay(dataset, "m6ASiteProb")

# Apply technology-site Total>10 ----------------------------------------

filtered_probability <- probability
filtered_probability[!is.finite(total) | total <= coverage_cutoff] <- NA_real_
rank_matrix <- normalized_rank_matrix(filtered_probability)

metric_rows <- list()
audit_rows <- list()
qc_rows <- list()

# Run the 14 LOTO folds --------------------------------------------------

for (heldout in technologies) {
  message("\nHeld-out technology: ", heldout)

  heldout_index <- match(heldout, technologies)
  kept_index <- setdiff(seq_along(technologies), heldout_index)
  kept_technologies <- technologies[kept_index]

  # Retrain ConquerBN with the remaining technologies -------------------

  model <- fit_conquernet_holdout(
    dataset = dataset,
    kept_technologies = kept_technologies,
    cutoff = coverage_cutoff,
    max_iter = em_max_iteration,
    epsilon = em_convergence_epsilon,
    seed = em_seed,
    analysis_version = "Total_gt10_LOTO_v1",
    holdout_name = heldout
  )

  saveRDS(
    model,
    file.path(model_dir, paste0("Leave.", heldout, ".rds")),
    compress = TRUE
  )

  conquernet_score <- rep(NA_real_, n_sites)
  conquernet_score[model$row_index] <- model$score

  # Held-out ground truth: BBmix MAP under Total>10 ----------------------

  heldout_observed <- is.finite(filtered_probability[, heldout_index])
  heldout_label <- rep(NA_integer_, n_sites)
  heldout_label[heldout_observed] <- as.integer(
    filtered_probability[heldout_observed, heldout_index] > posterior_cutoff
  )

  # RankProduct and RRA from the remaining technologies -----------------

  kept_probability <- filtered_probability[, kept_index, drop = FALSE]
  n_available <- rowSums(is.finite(kept_probability))
  kept_ranks <- rank_matrix[, kept_index, drop = FALSE]

  scores_without_idr <- list(
    ConquerBN = conquernet_score,
    RankProduct = rank_product_score(kept_ranks),
    RRA = rra_score(kept_ranks)
  )

  # Cutoff-matched IDR continuous score ---------------------------------

  idr <- read_idr_baseline(
    directory = loto_idr_dir,
    unit = heldout,
    cutoff = coverage_cutoff,
    n_sites = n_sites
  )
  scores_with_idr <- c(scores_without_idr, list(IDR = idr$score))

  # Common evaluation domains -------------------------------------------

  common_without_idr <- !is.na(heldout_label) & n_available > 0L &
    Reduce("&", lapply(scores_without_idr, is.finite))

  common_with_idr <- common_without_idr & idr$coverage &
    Reduce("&", lapply(scores_with_idr, is.finite))

  # Three-method sensitivity analysis -----------------------------------

  for (method in names(scores_without_idr)) {
    current <- ranking_metrics(
      heldout_label[common_without_idr],
      scores_without_idr[[method]][common_without_idr],
      method,
      heldout,
      coverage_cutoff
    )
    current$evaluation_domain <- "three_method_no_IDR_common_sensitivity"
    metric_rows[[length(metric_rows) + 1L]] <- current
  }

  # Four-method IDR-common primary analysis ------------------------------

  for (method in names(scores_with_idr)) {
    current <- ranking_metrics(
      heldout_label[common_with_idr],
      scores_with_idr[[method]][common_with_idr],
      method,
      heldout,
      coverage_cutoff
    )
    current$evaluation_domain <- "four_method_IDR_common_primary"
    metric_rows[[length(metric_rows) + 1L]] <- current
  }

  # Evaluation-domain audit ---------------------------------------------

  audit_rows[[length(audit_rows) + 1L]] <- data.frame(
    heldout = heldout,
    coverage_rule = paste0("Total>", coverage_cutoff),
    coverage_cutoff = coverage_cutoff,
    ground_truth = "heldout_BBmix_MAP_m6ASiteProb_gt_0.5",
    n_total_sites = n_sites,
    n_heldout_total_passing_coverage = sum(
      is.finite(total[, heldout_index]) & total[, heldout_index] > coverage_cutoff
    ),
    n_heldout_probability_observed = sum(heldout_observed),
    n_heldout_positive = sum(heldout_label == 1L, na.rm = TRUE),
    n_heldout_negative = sum(heldout_label == 0L, na.rm = TRUE),
    n_training_sites_in_model = length(model$row_index),
    n_remaining_probability_available = sum(n_available > 0L),
    n_three_method_no_IDR_common = sum(common_without_idr),
    n_positive_three_method_no_IDR_common = sum(
      heldout_label[common_without_idr] == 1L
    ),
    n_negative_three_method_no_IDR_common = sum(
      heldout_label[common_without_idr] == 0L
    ),
    three_method_no_IDR_valid_binary =
      length(unique(heldout_label[common_without_idr])) == 2L,
    three_method_no_IDR_identifiable_at_1000 =
      sum(heldout_label[common_without_idr] == 1L) >= 1000L &&
      sum(heldout_label[common_without_idr] == 0L) >= 1000L,
    n_idr_coverage = sum(idr$coverage),
    n_four_method_IDR_common = sum(common_with_idr),
    n_positive_four_method_IDR_common = sum(
      heldout_label[common_with_idr] == 1L
    ),
    n_negative_four_method_IDR_common = sum(
      heldout_label[common_with_idr] == 0L
    ),
    four_method_IDR_valid_binary =
      length(unique(heldout_label[common_with_idr])) == 2L,
    four_method_IDR_identifiable_at_1000 =
      sum(heldout_label[common_with_idr] == 1L) >= 1000L &&
      sum(heldout_label[common_with_idr] == 0L) >= 1000L,
    idr_cutoff_matched = TRUE,
    idr_rowmin_file = idr$files[["rowmin"]],
    idr_label_file = idr$files[["label"]],
    idr_coverage_file = idr$files[["coverage"]],
    idr_score_definition = "-rowmin (smaller rowmin means stronger m6A rank)",
    stringsAsFactors = FALSE
  )

  # Model QC -------------------------------------------------------------

  qc_rows[[length(qc_rows) + 1L]] <- data.frame(
    coverage_cutoff = coverage_cutoff,
    heldout = heldout,
    n_training_sites = model$n_training_sites,
    n_kept_technologies = length(kept_technologies),
    converged = model$converged,
    iterations = model$iterations,
    runtime_minutes = model$runtime_minutes,
    score_min = min(model$score),
    score_max = max(model$score),
    passed = isTRUE(model$converged) && all(is.finite(model$score)),
    stringsAsFactors = FALSE
  )

  rm(model, conquernet_score, idr)
  gc(verbose = FALSE)
}

# Save fold-level results ------------------------------------------------

metrics <- do.call(rbind, metric_rows)
audit <- do.call(rbind, audit_rows)
qc <- do.call(rbind, qc_rows)

metrics <- metrics[
  order(metrics$fold, metrics$evaluation_domain, metrics$method),
  , drop = FALSE
]
audit <- audit[order(audit$heldout), , drop = FALSE]
qc <- qc[order(qc$heldout), , drop = FALSE]
write.csv(metrics, file.path(output_dir, "per_fold_metrics.csv"), row.names = FALSE)
write.csv(audit, file.path(output_dir, "evaluation_domain_audit.csv"), row.names = FALSE)
write.csv(qc, file.path(qc_dir, "loto_model_qc.csv"), row.names = FALSE)

# Generate summary tables ------------------------------------------------

summarize_loto(metrics, output_dir)
message("\nLOTO completed.")
message("Fold-level results: ", output_dir)
message("Model QC: ", file.path(qc_dir, "loto_model_qc.csv"))
