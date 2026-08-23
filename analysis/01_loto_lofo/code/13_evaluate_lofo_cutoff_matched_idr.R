suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(matrixStats)
  library(pROC)
  library(PRROC)
})

## W6 B/C: evaluate coverage-aware LOSO models with held-out BBmix MAP labels
##
## For each Total cutoff and each held-out family:
##   - load the LOSO model trained under that same cutoff;
##   - define each held-out technology's ground truth using only sites with
##       held-out Total > cutoff and finite held-out m6ASiteProb;
##       positive: m6ASiteProb > 0.5;
##       negative: m6ASiteProb <= 0.5;
##   - recompute MeanProb, RankProduct and RRA from the same coverage-filtered
##     remaining four-family input;
##   - evaluate all methods on their common domain.
##
## Per-technology metrics are retained for transparency, but the five held-out
## families are the independent LOSO units. The primary LOSO summaries and
## method comparisons therefore use one macro result per held-out family.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot determine current script path.")
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
script_dir <- dirname(script_file)

source(file.path(script_dir, "w6_helpers.R"))
cli <- parse_cli()

data_dir <- normalizePath(
  cli_value(cli, "data-dir", "/home/xiyao/15data/ConquerNet/ConquerBN_basic"),
  winslash = "/",
  mustWork = TRUE
)
prediction_dir <- normalizePath(
  cli_value(
    cli,
    "prediction-dir",
    "/home/xiyao/15data/W6_corrected/results/02_loso_models"
  ),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- cli_value(
  cli,
  "output-dir",
  "/home/xiyao/15data/W7/02_cutoff_matched_IDR_benchmark/lofo"
)
cutoffs <- as.integer(strsplit(
  cli_value(cli, "cutoffs", "10"), ",", fixed = TRUE
)[[1]])
cutoffs <- sort(unique(cutoffs[is.finite(cutoffs)]))
if (!identical(cutoffs, 10L)) stop("This review script only permits --cutoffs=10.")

posterior_cutoff <- as.numeric(cli_value(cli, "posterior-cutoff", "0.5"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

idr_dir <- normalizePath(
  cli_value(
    cli,
    "idr-dir",
    "/home/xichen/ConquerNet/ConquerBN_basic/LOFO_IDR"
  ),
  winslash = "/",
  mustWork = TRUE
)

coverage_tag <- function(cutoff) paste0("Total_gt_", cutoff)

safe_sign_test <- function(delta) {
  delta <- delta[is.finite(delta) & delta != 0]
  if (!length(delta)) return(NA_real_)
  stats::binom.test(
    x = sum(delta > 0),
    n = length(delta),
    p = 0.5,
    alternative = "two.sided"
  )$p.value
}

safe_exact_wilcox <- function(delta) {
  delta <- delta[is.finite(delta)]
  if (length(delta) < 2L || any(delta == 0) || anyDuplicated(abs(delta))) {
    return(NA_real_)
  }
  tryCatch(
    stats::wilcox.test(
      delta,
      mu = 0,
      alternative = "two.sided",
      exact = TRUE,
      conf.int = FALSE
    )$p.value,
    error = function(e) NA_real_,
    warning = function(w) NA_real_
  )
}

family_macro_from_technology <- function(evaluable) {
  split_key <- interaction(
    evaluable$coverage_cutoff,
    evaluable$coverage_rule,
    evaluable$heldout_family,
    evaluable$method,
    drop = TRUE,
    lex.order = TRUE
  )
  groups <- split(evaluable, split_key)

  result <- do.call(rbind, lapply(groups, function(x) {
    data.frame(
      coverage_cutoff = unique(x$coverage_cutoff),
      coverage_rule = unique(x$coverage_rule),
      heldout_family = unique(x$heldout_family),
      method = unique(x$method),
      n_evaluable_technologies = nrow(x),
      technology_macro_auroc = mean(x$auroc, na.rm = TRUE),
      technology_macro_auprc = mean(x$auprc, na.rm = TRUE),
      mean_n_per_technology = mean(x$n, na.rm = TRUE),
      total_n_across_technologies = sum(x$n, na.rm = TRUE),
      total_positives_across_technologies = sum(x$positives, na.rm = TRUE),
      total_negatives_across_technologies = sum(x$negatives, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result[order(
    result$coverage_cutoff,
    result$heldout_family,
    result$method
  ), ]
}

dataset_path <- file.path(
  data_dir,
  "m6AConquer_hg38_BenchmarkDataSet_Combined_SE.rds"
)
dataset <- readRDS(dataset_path)
if (!inherits(dataset, "SummarizedExperiment")) {
  stop("Combined dataset is not a SummarizedExperiment.")
}

required_assays <- c("Total", "m6ASiteProb")
missing_assays <- setdiff(required_assays, assayNames(dataset))
if (length(missing_assays)) {
  stop("Missing required assays: ", paste(missing_assays, collapse = ", "))
}

total <- assay(dataset, "Total")
probability <- assay(dataset, "m6ASiteProb")
technologies <- colnames(dataset)
n_sites <- nrow(dataset)
strategy_map <- conquernet_strategy_map()
stopifnot(
  identical(dim(total), dim(probability)),
  identical(dimnames(total), dimnames(probability)),
  setequal(technologies, names(strategy_map))
)
families <- unique(unname(strategy_map[technologies]))

idr_paths <- function(heldout_family, coverage_cutoff) {
  result <- list(
    rowmin = file.path(
      idr_dir,
      paste0(
        "Leave.", heldout_family, ".idr.m6A.rowmin.cov",
        coverage_cutoff, ".rds"
      )
    ),
    label = file.path(
      idr_dir,
      paste0(
        "Leave.", heldout_family, ".idr.m6A.label.cov",
        coverage_cutoff, ".rds"
      )
    ),
    coverage = file.path(
      idr_dir,
      paste0(
        "Leave.", heldout_family, ".idr.idx.cov",
        coverage_cutoff, ".rds"
      )
    )
  )
  missing <- unlist(result, use.names = FALSE)
  missing <- missing[!file.exists(missing)]
  if (length(missing)) {
    stop("Missing cutoff-matched LOFO IDR files:\n", paste(missing, collapse = "\n"))
  }
  result
}

for (coverage_cutoff in cutoffs) {
  invisible(lapply(families, idr_paths, coverage_cutoff = coverage_cutoff))
}

metric_rows <- list()
audit_rows <- list()
model_qc_rows <- list()

for (coverage_cutoff in cutoffs) {
  message("=== LOSO EVALUATION: Total>", coverage_cutoff, " ===")

  threshold_probability <- probability
  threshold_probability[
    !is.finite(total) | total <= coverage_cutoff
  ] <- NA_real_
  probability_count <- rowSums(is.finite(threshold_probability))
  probability_sum <- rowSums(threshold_probability, na.rm = TRUE)
  rank_matrix <- normalized_rank_matrix(
    threshold_probability,
    universe_size = n_sites
  )

  for (heldout_family in families) {
    heldout_technologies <- names(strategy_map)[strategy_map == heldout_family]
    kept <- setdiff(technologies, heldout_technologies)
    heldout_index <- match(heldout_technologies, technologies)
    kept_index <- match(kept, technologies)

    prediction_path <- file.path(
      prediction_dir,
      coverage_tag(coverage_cutoff),
      paste0("LeaveFamily.", heldout_family, ".conquernet.rds")
    )
    if (!file.exists(prediction_path)) {
      stop("Missing coverage-matched LOSO prediction: ", prediction_path)
    }
    prediction <- readRDS(prediction_path)

    stopifnot(
      identical(prediction$analysis_version, "coverage_strict_loso_v3"),
      identical(as.character(prediction$heldout_family), heldout_family),
      setequal(as.character(prediction$heldout_technologies), heldout_technologies),
      identical(as.character(prediction$kept_technologies), kept),
      identical(as.integer(prediction$coverage_cutoff), as.integer(coverage_cutoff)),
      isTRUE(prediction$coverage_applied_to_training_entries),
      isTRUE(prediction$converged),
      length(prediction$row_index) == length(prediction$score),
      !anyDuplicated(prediction$row_index),
      all(is.finite(prediction$score)),
      all(prediction$score >= 0 & prediction$score <= 1)
    )

    model_qc_rows[[length(model_qc_rows) + 1L]] <- data.frame(
      heldout_family = heldout_family,
      coverage_rule = paste0("Total>", coverage_cutoff),
      coverage_cutoff = coverage_cutoff,
      n_heldout_technologies = length(heldout_technologies),
      n_kept_technologies = length(kept),
      n_training_sites = length(prediction$row_index),
      n_kept_entries_passing_coverage = prediction$n_kept_entries_passing_coverage,
      converged = prediction$converged,
      iterations = prediction$iterations,
      runtime_minutes = prediction$runtime_minutes,
      score_min = min(prediction$score),
      score_max = max(prediction$score),
      stringsAsFactors = FALSE
    )

    conquer_score <- rep(NA_real_, n_sites)
    conquer_score[prediction$row_index] <- prediction$score

    current_idr_paths <- idr_paths(heldout_family, coverage_cutoff)
    idr_rowmin <- readRDS(current_idr_paths$rowmin)
    idr_label <- readRDS(current_idr_paths$label)
    idr_coverage <- readRDS(current_idr_paths$coverage)
    if (
      !is.numeric(idr_rowmin) ||
        !is.logical(idr_label) ||
        !is.logical(idr_coverage) ||
        length(idr_rowmin) != n_sites ||
        length(idr_label) != n_sites ||
        length(idr_coverage) != n_sites
    ) {
      stop("Invalid LOFO IDR object structure for ", heldout_family)
    }
    if (
      !identical(is.finite(idr_rowmin), idr_coverage) ||
        !identical(idr_label, idr_coverage & idr_rowmin < 0.05)
    ) {
      stop("LOFO IDR rowmin/label/coverage relationship failed for ", heldout_family)
    }
    if (
      any(idr_rowmin[idr_coverage] < 0 | idr_rowmin[idr_coverage] > 1) ||
        any(!idr_coverage & !(is.infinite(idr_rowmin) & idr_rowmin > 0))
    ) {
      stop("Unexpected LOFO IDR rowmin range for ", heldout_family)
    }
    idr_score <- -idr_rowmin

    ## Baselines use the same coverage-filtered remaining four-family input.
    kept_probability <- threshold_probability[, kept_index, drop = FALSE]
    remaining_count <- rowSums(is.finite(kept_probability))
    remaining_sum <- rowSums(kept_probability, na.rm = TRUE)
    mean_score <- remaining_sum / remaining_count
    mean_score[remaining_count == 0L] <- NA_real_

    family_ranks <- rank_matrix[, kept_index, drop = FALSE]
    score_list <- list(
      ConquerBN = conquer_score,
      MeanProb = mean_score,
      RankProduct = rank_product_score(family_ranks),
      RRA = rra_score(family_ranks),
      IDR = idr_score
    )
    all_scores_finite <- Reduce(`&`, lapply(score_list, is.finite)) &
      idr_coverage

    for (j in heldout_index) {
      heldout <- technologies[j]

      ## Held-out BBmix MAP label under the same coverage rule.
      heldout_observed <- is.finite(threshold_probability[, j])
      label <- rep(NA_integer_, n_sites)
      label[heldout_observed] <- as.integer(
        threshold_probability[heldout_observed, j] > posterior_cutoff
      )

      common <- !is.na(label) &
        remaining_count > 0L &
        all_scores_finite

      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        heldout_family = heldout_family,
        heldout_technology = heldout,
        coverage_rule = paste0("Total>", coverage_cutoff),
        coverage_cutoff = coverage_cutoff,
        ground_truth = "heldout_BBmix_MAP_m6ASiteProb_gt_0.5",
        n_total_sites = n_sites,
        n_heldout_total_passing_coverage = sum(
          is.finite(total[, j]) & total[, j] > coverage_cutoff
        ),
        n_heldout_probability_observed = sum(heldout_observed),
        n_heldout_positive = sum(label == 1L, na.rm = TRUE),
        n_heldout_negative = sum(label == 0L, na.rm = TRUE),
        n_training_sites_in_family_model = length(prediction$row_index),
        n_remaining_probability_available = sum(remaining_count > 0L),
        n_idr_coverage = sum(idr_coverage),
        n_common = sum(common),
        n_positive_common = sum(label[common] == 1L),
        n_negative_common = sum(label[common] == 0L),
        valid_binary_technology = sum(label[common] == 1L) > 0L &&
          sum(label[common] == 0L) > 0L,
        identifiable_at_1000 = sum(label[common] == 1L) >= 1000L &&
          sum(label[common] == 0L) >= 1000L,
        idr_cutoff_matched = TRUE,
        idr_rowmin_file = basename(current_idr_paths$rowmin),
        idr_label_file = basename(current_idr_paths$label),
        idr_coverage_file = basename(current_idr_paths$coverage),
        idr_score_definition =
          "-rowmin (smaller rowmin means stronger m6A rank)",
        stringsAsFactors = FALSE
      )

      for (method in names(score_list)) {
        metric_result <- ranking_metrics(
          label[common],
          score_list[[method]][common],
          method,
          fold = heldout,
          heldout_family = heldout_family,
          cutoff = coverage_cutoff
        )
        metric_result$evaluation_domain <-
          "five_method_cutoff_matched_IDR_common"
        metric_rows[[length(metric_rows) + 1L]] <- metric_result
      }
    }

    rm(
      prediction,
      conquer_score,
      kept_probability,
      family_ranks,
      score_list,
      idr_rowmin,
      idr_label,
      idr_coverage,
      idr_score
    )
    gc()
  }

  rm(
    threshold_probability,
    probability_count,
    probability_sum,
    rank_matrix
  )
  gc()
}

metrics <- do.call(rbind, metric_rows)
audit <- do.call(rbind, audit_rows)
model_qc <- do.call(rbind, model_qc_rows)
names(metrics)[names(metrics) == "fold"] <- "heldout_technology"

metrics <- metrics[order(
  metrics$coverage_cutoff,
  metrics$heldout_family,
  metrics$heldout_technology,
  metrics$method
), ]
audit <- audit[order(
  audit$coverage_cutoff,
  audit$heldout_family,
  audit$heldout_technology
), ]
model_qc <- model_qc[order(
  model_qc$coverage_cutoff,
  model_qc$heldout_family
), ]

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
write.csv(
  model_qc,
  file.path(output_dir, "loso_model_qc.csv"),
  row.names = FALSE
)

## Primary independent units: five held-out families.
evaluable <- metrics[is.finite(metrics$auroc) & is.finite(metrics$auprc), ]
if (!nrow(evaluable)) stop("No evaluable LOSO technology results.")
family_macro <- family_macro_from_technology(evaluable)
write.csv(
  family_macro,
  file.path(output_dir, "family_macro_summary.csv"),
  row.names = FALSE
)

family_completeness <- do.call(rbind, lapply(
  split(
    family_macro,
    interaction(
      family_macro$coverage_cutoff,
      family_macro$method,
      drop = TRUE,
      lex.order = TRUE
    )
  ),
  function(x) {
    missing_families <- setdiff(families, x$heldout_family)
    data.frame(
      coverage_cutoff = unique(x$coverage_cutoff),
      coverage_rule = unique(x$coverage_rule),
      method = unique(x$method),
      n_expected_families = length(families),
      n_evaluable_families = length(unique(x$heldout_family)),
      missing_families = paste(missing_families, collapse = ";"),
      complete_five_family_summary = length(missing_families) == 0L,
      stringsAsFactors = FALSE
    )
  }
))
rownames(family_completeness) <- NULL
write.csv(
  family_completeness,
  file.path(output_dir, "five_family_completeness_audit.csv"),
  row.names = FALSE
)

five_family_macro <- do.call(rbind, lapply(
  split(
    family_macro,
    interaction(
      family_macro$coverage_cutoff,
      family_macro$method,
      drop = TRUE,
      lex.order = TRUE
    )
  ),
  function(x) {
    data.frame(
      coverage_cutoff = unique(x$coverage_cutoff),
      coverage_rule = unique(x$coverage_rule),
      method = unique(x$method),
      n_independent_families = length(unique(x$heldout_family)),
      family_macro_auroc = mean(x$technology_macro_auroc, na.rm = TRUE),
      family_macro_auprc = mean(x$technology_macro_auprc, na.rm = TRUE),
      sd_family_auroc = stats::sd(x$technology_macro_auroc, na.rm = TRUE),
      sd_family_auprc = stats::sd(x$technology_macro_auprc, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
rownames(five_family_macro) <- NULL
five_family_macro <- five_family_macro[order(
  five_family_macro$coverage_cutoff,
  five_family_macro$method
), ]
write.csv(
  five_family_macro,
  file.path(output_dir, "five_family_macro_summary.csv"),
  row.names = FALSE
)

## ConquerBN-vs-baseline comparisons: exactly one paired value per family.
comparison_rows <- list()
comparison_summary_rows <- list()
for (coverage_cutoff in cutoffs) {
  current <- family_macro[family_macro$coverage_cutoff == coverage_cutoff, ]
  conquer <- current[
    current$method == "ConquerBN",
    c("heldout_family", "technology_macro_auroc", "technology_macro_auprc")
  ]
  names(conquer)[2:3] <- c("conquernet_auroc", "conquernet_auprc")

  for (baseline in setdiff(unique(current$method), "ConquerBN")) {
    other <- current[
      current$method == baseline,
      c("heldout_family", "technology_macro_auroc", "technology_macro_auprc")
    ]
    names(other)[2:3] <- c("baseline_auroc", "baseline_auprc")
    paired <- merge(conquer, other, by = "heldout_family")
    paired$coverage_cutoff <- coverage_cutoff
    paired$coverage_rule <- paste0("Total>", coverage_cutoff)
    paired$baseline_method <- baseline
    paired$delta_auroc_conquernet_minus_baseline <-
      paired$conquernet_auroc - paired$baseline_auroc
    paired$delta_auprc_conquernet_minus_baseline <-
      paired$conquernet_auprc - paired$baseline_auprc
    comparison_rows[[length(comparison_rows) + 1L]] <- paired

    comparison_summary_rows[[length(comparison_summary_rows) + 1L]] <- data.frame(
      coverage_cutoff = coverage_cutoff,
      coverage_rule = paste0("Total>", coverage_cutoff),
      baseline_method = baseline,
      n_independent_families = nrow(paired),
      n_auroc_wins = sum(
        paired$delta_auroc_conquernet_minus_baseline > 0,
        na.rm = TRUE
      ),
      mean_delta_auroc = mean(
        paired$delta_auroc_conquernet_minus_baseline,
        na.rm = TRUE
      ),
      median_delta_auroc = median(
        paired$delta_auroc_conquernet_minus_baseline,
        na.rm = TRUE
      ),
      auroc_exact_sign_test_p = safe_sign_test(
        paired$delta_auroc_conquernet_minus_baseline
      ),
      auroc_exact_wilcox_p = safe_exact_wilcox(
        paired$delta_auroc_conquernet_minus_baseline
      ),
      n_auprc_wins = sum(
        paired$delta_auprc_conquernet_minus_baseline > 0,
        na.rm = TRUE
      ),
      mean_delta_auprc = mean(
        paired$delta_auprc_conquernet_minus_baseline,
        na.rm = TRUE
      ),
      median_delta_auprc = median(
        paired$delta_auprc_conquernet_minus_baseline,
        na.rm = TRUE
      ),
      auprc_exact_sign_test_p = safe_sign_test(
        paired$delta_auprc_conquernet_minus_baseline
      ),
      auprc_exact_wilcox_p = safe_exact_wilcox(
        paired$delta_auprc_conquernet_minus_baseline
      ),
      stringsAsFactors = FALSE
    )
  }
}

family_method_differences <- do.call(rbind, comparison_rows)
family_method_comparison <- do.call(rbind, comparison_summary_rows)
write.csv(
  family_method_differences,
  file.path(output_dir, "family_method_paired_differences.csv"),
  row.names = FALSE
)
write.csv(
  family_method_comparison,
  file.path(output_dir, "five_family_method_comparison.csv"),
  row.names = FALSE
)


message("Completed coverage-aware LOSO evaluation.")
message("Primary independent units: five held-out families.")
message("Results: ", output_dir)
