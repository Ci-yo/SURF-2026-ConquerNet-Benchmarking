## Leave-one-family-out benchmark
## Each fold removes one complete technology family, retrains ConquerBN with
## the remaining families, and evaluates each held-out technology against its
## own BBmix MAP call.

# Locate script -----------------------------------------------------------

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Please run this file with Rscript.")
SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", script_arg[[1]]),
  winslash = "/"
))


# Load settings and functions --------------------------------------------

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


# Output folders ----------------------------------------------------------

output_dir <- ensure_directory(file.path(result_dir, "benchmark", "LOFO"))
model_dir <- ensure_directory(file.path(output_dir, "models"))
qc_dir <- ensure_directory(file.path(result_dir, "benchmark", "QC"))


# Load benchmark data -----------------------------------------------------

require_files(c(benchmark_file, lofo_idr_dir))
dataset <- load_benchmark_dataset(CONFIG)

technologies <- colnames(dataset)
strategy_map <- conquernet_strategy_map()[technologies]
families <- unique(unname(strategy_map))
n_sites <- nrow(dataset)

total <- assay(dataset, "Total")
probability <- assay(dataset, "m6ASiteProb")


# Apply technology-site Total>10 -----------------------------------------

filtered_probability <- probability
filtered_probability[!is.finite(total) | total <= coverage_cutoff] <- NA_real_
rank_matrix <- normalized_rank_matrix(filtered_probability)

metric_rows <- list()
audit_rows <- list()
qc_rows <- list()


# Run the five family holdouts --------------------------------------------

for (heldout_family in families) {
  message("\nHeld-out family: ", heldout_family)

  heldout_technologies <- technologies[strategy_map == heldout_family]
  kept_technologies <- setdiff(technologies, heldout_technologies)
  kept_index <- match(kept_technologies, technologies)


  # Retrain ConquerBN without the held-out family ------------------------

  model <- fit_conquernet_holdout(
    dataset = dataset,
    kept_technologies = kept_technologies,
    cutoff = coverage_cutoff,
    max_iter = em_max_iteration,
    epsilon = em_convergence_epsilon,
    seed = em_seed,
    analysis_version = "Total_gt10_LOFO_v1",
    holdout_name = heldout_family
  )

  saveRDS(
    model,
    file.path(model_dir, paste0("LeaveFamily.", heldout_family, ".rds")),
    compress = TRUE
  )

  conquernet_score <- rep(NA_real_, n_sites)
  conquernet_score[model$row_index] <- model$score


  # Build the other three scores from the remaining technologies --------

  kept_probability <- filtered_probability[, kept_index, drop = FALSE]
  n_available <- rowSums(is.finite(kept_probability))
  kept_ranks <- rank_matrix[, kept_index, drop = FALSE]

  scores <- list(
    ConquerBN = conquernet_score,
    RankProduct = rank_product_score(kept_ranks),
    RRA = rra_score(kept_ranks)
  )

  idr <- read_idr_baseline(
    directory = lofo_idr_dir,
    unit = heldout_family,
    cutoff = coverage_cutoff,
    n_sites = n_sites
  )
  scores$IDR <- idr$score

  all_scores_finite <- Reduce("&", lapply(scores, is.finite))


  # Evaluate each technology inside the held-out family -----------------

  for (heldout_technology in heldout_technologies) {
    technology_index <- match(heldout_technology, technologies)

    # BBmix MAP ground truth under Total>10
    heldout_observed <- is.finite(filtered_probability[, technology_index])
    heldout_label <- rep(NA_integer_, n_sites)
    heldout_label[heldout_observed] <- as.integer(
      filtered_probability[heldout_observed, technology_index] > posterior_cutoff
    )

    # All four methods are evaluated on exactly the same sites
    common <-
      !is.na(heldout_label) &
      n_available > 0L &
      idr$coverage &
      all_scores_finite

    for (method in names(scores)) {
      current <- ranking_metrics(
        y = heldout_label[common],
        score = scores[[method]][common],
        method = method,
        unit = heldout_technology,
        cutoff = coverage_cutoff,
        family = heldout_family
      )
      current$evaluation_domain <- "four_method_cutoff_matched_IDR_common"
      metric_rows[[length(metric_rows) + 1L]] <- current
    }

    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      heldout_family = heldout_family,
      heldout_technology = heldout_technology,
      coverage_rule = paste0("Total>", coverage_cutoff),
      coverage_cutoff = coverage_cutoff,
      ground_truth = "heldout_BBmix_MAP_m6ASiteProb_gt_0.5",
      n_total_sites = n_sites,
      n_heldout_total_passing_coverage = sum(
        is.finite(total[, technology_index]) &
          total[, technology_index] > coverage_cutoff
      ),
      n_heldout_probability_observed = sum(heldout_observed),
      n_heldout_positive = sum(heldout_label == 1L, na.rm = TRUE),
      n_heldout_negative = sum(heldout_label == 0L, na.rm = TRUE),
      n_training_sites_in_family_model = length(model$row_index),
      n_remaining_probability_available = sum(n_available > 0L),
      n_idr_coverage = sum(idr$coverage),
      n_common = sum(common),
      n_positive_common = sum(heldout_label[common] == 1L),
      n_negative_common = sum(heldout_label[common] == 0L),
      valid_binary_technology =
        length(unique(heldout_label[common])) == 2L,
      identifiable_at_1000 =
        sum(heldout_label[common] == 1L) >= 1000L &&
        sum(heldout_label[common] == 0L) >= 1000L,
      idr_cutoff_matched = TRUE,
      idr_rowmin_file = idr$files[["rowmin"]],
      idr_label_file = idr$files[["label"]],
      idr_coverage_file = idr$files[["coverage"]],
      idr_score_definition =
        "-rowmin (smaller rowmin means stronger m6A rank)",
      stringsAsFactors = FALSE
    )
  }


  # Record one QC row for the family-level model -------------------------

  qc_rows[[length(qc_rows) + 1L]] <- data.frame(
    heldout_family = heldout_family,
    coverage_rule = paste0("Total>", coverage_cutoff),
    coverage_cutoff = coverage_cutoff,
    n_heldout_technologies = length(heldout_technologies),
    n_kept_technologies = length(kept_technologies),
    n_training_sites = model$n_training_sites,
    n_kept_entries_passing_coverage = model$n_kept_entries_passing_coverage,
    converged = model$converged,
    iterations = model$iterations,
    runtime_minutes = model$runtime_minutes,
    score_min = min(model$score),
    score_max = max(model$score),
    passed = isTRUE(model$converged) && all(is.finite(model$score)),
    stringsAsFactors = FALSE
  )

  rm(model, conquernet_score, idr, scores)
  gc(verbose = FALSE)
}


# Save technology-level and family-level results -------------------------

metrics <- do.call(rbind, metric_rows)
names(metrics)[names(metrics) == "fold"] <- "heldout_technology"
audit <- do.call(rbind, audit_rows)
qc <- do.call(rbind, qc_rows)

metrics <- metrics[
  order(metrics$heldout_family, metrics$heldout_technology, metrics$method),
  ,
  drop = FALSE
]
audit <- audit[
  order(audit$heldout_family, audit$heldout_technology),
  ,
  drop = FALSE
]
qc <- qc[order(qc$heldout_family), , drop = FALSE]

write.csv(
  metrics,
  file.path(output_dir, "per_technology_metrics.csv"),
  row.names = FALSE
)
write.csv(
  audit,
  file.path(output_dir, "evaluation_domain_audit.csv"),
  row.names = FALSE
)
write.csv(qc, file.path(qc_dir, "lofo_model_qc.csv"), row.names = FALSE)
summarize_lofo(metrics, output_dir)
message("\nLOFO completed.")
message("Technology-level results: ", output_dir)
message("Model QC: ", file.path(qc_dir, "lofo_model_qc.csv"))