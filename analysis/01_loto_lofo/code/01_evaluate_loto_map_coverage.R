suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(bnlearn)
  library(extraDistr)
  library(matrixStats)
  library(pROC)
  library(PRROC)
})

## W6 A/B: coverage-aware strict LOTO + BBmix MAP evaluation
##
## Final formal route: Total > 10 for each held-out technology:
##   1. Remove every kept technology-site observation with Total <= cutoff.
##      Removed observations are represented as m6A=0, Total=0 and weight=NA,
##      so they contribute no emission evidence.
##   2. Retrain a new strict LOTO ConquerBN model from the remaining 13 technologies.
##   3. Define the held-out ground truth as:
##        held-out Total > cutoff AND finite m6ASiteProb
##        positive: m6ASiteProb > 0.5
##        negative: m6ASiteProb <= 0.5
##   4. Recompute MeanProb, RankProduct and RRA from the same coverage-filtered
##      remaining-technology probability matrix.
##   5. Evaluate all four methods on exactly the same sites.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot determine current script path.")
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")
script_dir <- dirname(script_file)

source(file.path(script_dir, "w6_helpers.R"))
source(file.path(script_dir, "latent_switch_fixed.R"))
cli <- parse_cli()

data_dir <- normalizePath(
  cli_value(cli, "data-dir", "/home/xiyao/15data/ConquerNet/ConquerBN_basic"),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(data_dir, "function.R"))

output_dir <- cli_value(
  cli,
  "output-dir",
  "/home/xiyao/15data/W6_corrected/results/01_loto_map_coverage"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

idr_dir <- normalizePath(
  cli_value(
    cli,
    "idr-dir",
    "/home/xiyao/15data/W6_corrected/inputs/LOTO_IDR"
  ),
  winslash = "/",
  mustWork = TRUE
)

cutoffs <- as.integer(strsplit(
  cli_value(cli, "cutoffs", "10"), ",", fixed = TRUE
)[[1]])
cutoffs <- sort(unique(cutoffs[is.finite(cutoffs)]))
if (!identical(cutoffs, 10L)) stop("This review script only permits --cutoffs=10.")

posterior_cutoff <- as.numeric(cli_value(cli, "posterior-cutoff", "0.5"))
max_iter <- as.integer(cli_value(cli, "max-iter", "500"))
epsilon <- as.numeric(cli_value(cli, "epsilon", "1e-6"))
seed <- as.integer(cli_value(cli, "seed", "20260725"))
smoke_n <- as.integer(cli_value(cli, "smoke-n", "0"))
overwrite <- cli_flag(cli, "overwrite")
train_only <- cli_flag(cli, "train-only")
evaluate_only <- cli_flag(cli, "evaluate-only")

if (train_only && evaluate_only) {
  stop("--train-only and --evaluate-only cannot be used together.")
}

fold_arg <- cli_value(cli, "fold", NULL)
run_all <- cli_flag(cli, "all")
if (is.null(fold_arg) && !run_all) {
  stop("Specify --all or exactly one --fold=TECHNOLOGY.")
}
if (!is.null(fold_arg) && run_all) {
  stop("Use either --all or --fold=TECHNOLOGY, not both.")
}

dataset_path <- file.path(
  data_dir,
  "m6AConquer_hg38_BenchmarkDataSet_Combined_SE.rds"
)
dataset <- readRDS(dataset_path)
if (!inherits(dataset, "SummarizedExperiment")) {
  stop("Combined dataset is not a SummarizedExperiment.")
}

required_assays <- c("m6A", "Total", "m6ASiteProb")
missing_assays <- setdiff(required_assays, assayNames(dataset))
if (length(missing_assays)) {
  stop("Missing required assays: ", paste(missing_assays, collapse = ", "))
}

m6a_all <- assay(dataset, "m6A")
total_all <- assay(dataset, "Total")
probability_all <- assay(dataset, "m6ASiteProb")
technologies <- colnames(dataset)
n_sites <- nrow(dataset)
strategy_map <- conquernet_strategy_map()

stopifnot(
  identical(dim(m6a_all), dim(total_all)),
  identical(dim(probability_all), dim(total_all)),
  identical(dimnames(m6a_all), dimnames(total_all)),
  identical(dimnames(probability_all), dimnames(total_all)),
  setequal(technologies, names(strategy_map))
)

idr_rowmin_paths <- setNames(
  file.path(
    idr_dir,
    paste0("Leave.", technologies, ".idr.m6A.rowmin.rds")
  ),
  technologies
)
idr_label_paths <- setNames(
  file.path(
    idr_dir,
    paste0("Leave.", technologies, ".idr.m6A.label.rds")
  ),
  technologies
)
idr_coverage_paths <- setNames(
  file.path(
    idr_dir,
    paste0("Leave.", technologies, ".idr.coverage.idx.rds")
  ),
  technologies
)
missing_idr <- c(
  idr_rowmin_paths[!file.exists(idr_rowmin_paths)],
  idr_label_paths[!file.exists(idr_label_paths)],
  idr_coverage_paths[!file.exists(idr_coverage_paths)]
)
if (length(missing_idr)) {
  stop("Missing IDR input files:\n", paste(missing_idr, collapse = "\n"))
}

if (!is.null(fold_arg) && !fold_arg %in% technologies) {
  stop("Unknown fold: ", fold_arg)
}
folds <- if (run_all) technologies else fold_arg

message("Dataset: ", n_sites, " sites x ", length(technologies), " technologies")
message("Coverage cutoffs: ", paste0("Total>", cutoffs, collapse = ", "))
message("LOTO folds requested: ", paste(folds, collapse = ", "))

coverage_tag <- function(cutoff) paste0("Total_gt_", cutoff)

make_dirs <- function(cutoff) {
  cutoff_dir <- file.path(output_dir, coverage_tag(cutoff))
  result <- list(
    cutoff = cutoff_dir,
    predictions = file.path(cutoff_dir, "predictions"),
    per_fold = file.path(cutoff_dir, "per_fold")
  )
  lapply(result[c("cutoff", "predictions", "per_fold")], dir.create,
         recursive = TRUE, showWarnings = FALSE)
  result
}

extract_root_score <- function(root_belief) {
  if (ncol(root_belief) != 2L) {
    stop("TrueM6A belief must have exactly two states.")
  }
  if (!is.null(colnames(root_belief)) && "1" %in% colnames(root_belief)) {
    score <- exp(root_belief[, "1"])
    method <- "named_state_1"
  } else {
    score <- exp(root_belief[, 2L])
    method <- "second_column_after_fixed_orientation"
  }
  if (!all(is.finite(score)) || any(score < 0 | score > 1)) {
    stop("Invalid TrueM6A posterior probability.")
  }
  list(score = score, method = method)
}

fit_loto_fold <- function(heldout, coverage_cutoff, dirs) {
  result_path <- file.path(
    dirs$predictions,
    paste0("Leave.", heldout, ".conquernet.rds")
  )

  if (file.exists(result_path) && !overwrite) {
    existing <- readRDS(result_path)
    valid_existing <- identical(existing$analysis_version, "coverage_strict_loto_v3") &&
      identical(as.character(existing$heldout), heldout) &&
      identical(as.integer(existing$coverage_cutoff), as.integer(coverage_cutoff)) &&
      identical(as.integer(existing$smoke_n), as.integer(smoke_n)) &&
      isTRUE(existing$converged) &&
      length(existing$row_index) == length(existing$score) &&
      !anyDuplicated(existing$row_index) &&
      all(is.finite(existing$score)) &&
      all(existing$score >= 0 & existing$score <= 1)
    if (!valid_existing) {
      stop("Existing prediction is incompatible; rerun with --overwrite: ", result_path)
    }
    message("Using existing model: ", result_path)
    return(invisible(result_path))
  }

  heldout_index <- match(heldout, technologies)
  kept <- setdiff(technologies, heldout)
  kept_index <- match(kept, technologies)

  ## Coverage filtering is applied before model fitting.
  kept_coverage <- is.finite(total_all[, kept_index, drop = FALSE]) &
    total_all[, kept_index, drop = FALSE] > coverage_cutoff
  training_available <- rowSums(kept_coverage) > 0L
  row_index <- which(training_available)

  if (smoke_n > 0L && length(row_index) > smoke_n) {
    set.seed(seed + 1000L * coverage_cutoff + heldout_index)
    row_index <- sort(sample(row_index, smoke_n))
  }

  message(
    "Fitting ", heldout, " at Total>", coverage_cutoff,
    ": ", length(row_index), " sites, ", length(kept), " kept technologies"
  )

  m6a <- m6a_all[row_index, kept_index, drop = FALSE]
  total <- total_all[row_index, kept_index, drop = FALSE]
  weight <- probability_all[row_index, kept_index, drop = FALSE]
  observed <- is.finite(total) & total > coverage_cutoff
  if (any(observed & (!is.finite(m6a) | !is.finite(weight)))) {
    stop(
      "Found Total>", coverage_cutoff,
      " observations with missing m6A or m6ASiteProb in ", heldout
    )
  }

  ## Treat below-threshold observations as absent from the likelihood.
  m6a[!observed] <- 0
  total[!observed] <- 0
  weight[!observed] <- NA_real_

  if (any(m6a < 0, na.rm = TRUE) || any(total < 0, na.rm = TRUE) ||
      any(m6a > total, na.rm = TRUE)) {
    stop("Invalid m6A/Total counts after coverage filtering: ", heldout)
  }

  meta <- data.frame(
    Technique = kept,
    Strategy = unname(strategy_map[kept]),
    stringsAsFactors = FALSE
  )
  strategy_nodes <- unique(meta$Strategy)
  technology_nodes <- paste0(meta$Technique, "_k")
  continuous.nodes <- paste0(meta$Technique, "_x")

  arc_set <- rbind(
    data.frame(from = "TrueM6A", to = strategy_nodes),
    data.frame(from = meta$Strategy, to = technology_nodes),
    data.frame(from = technology_nodes, to = continuous.nodes)
  )
  dag <- empty.graph(unique(c(arc_set$from, arc_set$to)))
  arcs(dag) <- as.matrix(arc_set)

  colnames(m6a) <- continuous.nodes
  colnames(total) <- continuous.nodes
  colnames(weight) <- continuous.nodes
  assign("continuous.nodes", continuous.nodes, envir = .GlobalEnv)

  checked <- format_check(m6a, total, meta)
  m6a <- checked$m6A
  total <- checked$Total
  meta <- checked$meta

  beta_param <- build_beta_params(dataset, kept)
  emission_loglik <- calculate_emission_loglik(m6a, total, beta_param)
  cpt_param <- initialization_cpt_params(m6a, total, beta_param, weight, dag)

  log_likelihood <- rep(NA_real_, max_iter)
  converged <- FALSE
  previous <- NA_real_
  started <- Sys.time()

  for (iteration in seq_len(max_iter)) {
    iteration_started <- Sys.time()
    expectation <- conquerNet_expectation(emission_loglik, cpt_param, dag)
    current <- expectation$loglik
    belief <- expectation$belief
    if (!is.finite(current)) {
      stop("Non-finite log-likelihood in ", heldout,
           " at Total>", coverage_cutoff, ", iteration ", iteration)
    }

    log_likelihood[iteration] <- current
    cpt_param <- conquerNet_maximization_safe(cpt_param, belief, dag)
    relative_change <- if (iteration > 1L) {
      abs(current - previous) / max(abs(previous), .Machine$double.eps)
    } else {
      NA_real_
    }

    message(sprintf(
      "  cutoff=%d heldout=%s iteration=%d loglik=%.6f relative_change=%s elapsed_min=%.3f",
      coverage_cutoff,
      heldout,
      iteration,
      current,
      ifelse(is.na(relative_change), "NA",
             format(relative_change, scientific = TRUE, digits = 4)),
      as.numeric(difftime(Sys.time(), iteration_started, units = "mins"))
    ))

    if (iteration > 1L && relative_change < epsilon) {
      converged <- TRUE
      break
    }
    previous <- current
    rm(expectation, belief)
    gc()
  }

  iterations <- max(which(!is.na(log_likelihood)))
  log_likelihood <- log_likelihood[seq_len(iterations)]
  if (!converged) {
    stop("LOTO model did not converge: ", heldout, " at Total>", coverage_cutoff)
  }

  path <- sapply(
    bnlearn::leaf.nodes(dag),
    function(node) get_parent_chain_bn(node, dag)
  )
  cpt_param <- latent_param_label_switch_fixed(cpt_param, path)
  final_expectation <- conquerNet_expectation(emission_loglik, cpt_param, dag)
  root_result <- extract_root_score(final_expectation$belief[["TrueM6A"]])

  result <- list(
    analysis_version = "coverage_strict_loto_v3",
    heldout = heldout,
    kept_technologies = kept,
    coverage_rule = paste0("Total>", coverage_cutoff),
    coverage_cutoff = coverage_cutoff,
    coverage_applied_to_training_entries = TRUE,
    below_threshold_representation = "m6A=0; Total=0; m6ASiteProb=NA",
    row_index = row_index,
    score = root_result$score,
    cpt_param = cpt_param,
    beta_param = beta_param,
    dag = dag,
    log_likelihood = log_likelihood,
    converged = converged,
    iterations = iterations,
    epsilon = epsilon,
    max_iter = max_iter,
    smoke_n = smoke_n,
    seed = seed,
    n_total_sites = n_sites,
    n_training_sites = length(row_index),
    n_kept_entries_total = length(observed),
    n_kept_entries_passing_coverage = sum(observed),
    started = started,
    finished = Sys.time(),
    runtime_minutes = as.numeric(difftime(Sys.time(), started, units = "mins")),
    dataset_path = dataset_path,
    orientation_method = "latent_param_label_switch_fixed_v1",
    root_state_extraction = root_result$method,
    session_info = sessionInfo()
  )
  saveRDS(result, result_path)
  message("Saved model: ", result_path)
  invisible(result_path)
}

evaluate_loto_fold <- function(heldout, coverage_cutoff, dirs,
                               threshold_probability, probability_count,
                               probability_sum, rank_matrix) {
  heldout_index <- match(heldout, technologies)
  kept_index <- setdiff(seq_along(technologies), heldout_index)

  prediction_path <- file.path(
    dirs$predictions,
    paste0("Leave.", heldout, ".conquernet.rds")
  )
  if (!file.exists(prediction_path)) stop("Missing LOTO prediction: ", prediction_path)
  prediction <- readRDS(prediction_path)

  stopifnot(
    identical(prediction$analysis_version, "coverage_strict_loto_v3"),
    identical(as.character(prediction$heldout), heldout),
    identical(as.integer(prediction$coverage_cutoff), as.integer(coverage_cutoff)),
    identical(as.character(prediction$kept_technologies), technologies[kept_index]),
    isTRUE(prediction$converged),
    length(prediction$row_index) == length(prediction$score),
    !anyDuplicated(prediction$row_index),
    all(is.finite(prediction$score)),
    all(prediction$score >= 0 & prediction$score <= 1)
  )

  conquer_score <- rep(NA_real_, n_sites)
  conquer_score[prediction$row_index] <- prediction$score

  ## Held-out BBmix MAP ground truth with the same coverage requirement.
  heldout_observed <- is.finite(threshold_probability[, heldout_index])
  label <- rep(NA_integer_, n_sites)
  label[heldout_observed] <- as.integer(
    threshold_probability[heldout_observed, heldout_index] > posterior_cutoff
  )

  ## Baselines use the same coverage-filtered remaining technologies.
  remaining_count <- probability_count - as.integer(heldout_observed)
  remaining_sum <- probability_sum - ifelse(
    heldout_observed,
    threshold_probability[, heldout_index],
    0
  )
  mean_score <- remaining_sum / remaining_count
  mean_score[remaining_count == 0L] <- NA_real_

  fold_ranks <- rank_matrix[, kept_index, drop = FALSE]
  four_method_scores <- list(
    ConquerBN = conquer_score,
    MeanProb = mean_score,
    RankProduct = rank_product_score(fold_ranks),
    RRA = rra_score(fold_ranks)
  )

  four_method_common <- !is.na(label) &
    remaining_count > 0L &
    Reduce(`&`, lapply(four_method_scores, is.finite))

  ## IDR rowmin is a continuous LOTO prediction supplied by the original
  ## IDR workflow. Smaller rowmin means a more reproducible/m6A-like site.
  ## It is a prediction baseline, never a ground-truth label.
  idr_rowmin <- readRDS(idr_rowmin_paths[[heldout]])
  idr_label <- readRDS(idr_label_paths[[heldout]])
  idr_coverage <- readRDS(idr_coverage_paths[[heldout]])
  if (
    !is.numeric(idr_rowmin) ||
      !is.logical(idr_label) ||
      !is.logical(idr_coverage) ||
      length(idr_rowmin) != n_sites ||
      length(idr_label) != n_sites ||
      length(idr_coverage) != n_sites
  ) {
    stop("Invalid IDR object structure for ", heldout)
  }
  if (
    !identical(is.finite(idr_rowmin), idr_coverage) ||
      !identical(idr_label, idr_coverage & idr_rowmin < 0.05)
  ) {
    stop("IDR rowmin/label/coverage relationship failed for ", heldout)
  }
  if (
    any(idr_rowmin[idr_coverage] < 0 | idr_rowmin[idr_coverage] > 1) ||
      any(!idr_coverage & !(is.infinite(idr_rowmin) & idr_rowmin > 0))
  ) {
    stop("Unexpected IDR rowmin value range for ", heldout)
  }
  idr_score <- -idr_rowmin
  five_method_scores <- c(four_method_scores, list(IDR = idr_score))
  five_method_common <- four_method_common & idr_coverage

  audit <- data.frame(
    heldout = heldout,
    coverage_rule = paste0("Total>", coverage_cutoff),
    coverage_cutoff = coverage_cutoff,
    ground_truth = "heldout_BBmix_MAP_m6ASiteProb_gt_0.5",
    n_total_sites = n_sites,
    n_heldout_total_passing_coverage = sum(
      is.finite(total_all[, heldout_index]) & total_all[, heldout_index] > coverage_cutoff
    ),
    n_heldout_probability_observed = sum(heldout_observed),
    n_heldout_positive = sum(label == 1L, na.rm = TRUE),
    n_heldout_negative = sum(label == 0L, na.rm = TRUE),
    n_training_sites_in_model = length(prediction$row_index),
    n_remaining_probability_available = sum(remaining_count > 0L),
    n_four_method_common = sum(four_method_common),
    n_positive_four_method_common = sum(label[four_method_common] == 1L),
    n_negative_four_method_common = sum(label[four_method_common] == 0L),
    four_method_valid_binary = sum(label[four_method_common] == 1L) > 0L &&
      sum(label[four_method_common] == 0L) > 0L,
    four_method_identifiable_at_1000 =
      sum(label[four_method_common] == 1L) >= 1000L &&
      sum(label[four_method_common] == 0L) >= 1000L,
    n_idr_coverage = sum(idr_coverage),
    n_five_method_common = sum(five_method_common),
    n_positive_five_method_common = sum(label[five_method_common] == 1L),
    n_negative_five_method_common = sum(label[five_method_common] == 0L),
    five_method_valid_binary = sum(label[five_method_common] == 1L) > 0L &&
      sum(label[five_method_common] == 0L) > 0L,
    five_method_identifiable_at_1000 =
      sum(label[five_method_common] == 1L) >= 1000L &&
      sum(label[five_method_common] == 0L) >= 1000L,
    idr_score_definition = "-rowmin (smaller rowmin means stronger m6A rank)",
    stringsAsFactors = FALSE
  )

  four_metric_rows <- lapply(names(four_method_scores), function(method) {
    result <- ranking_metrics(
      label[four_method_common],
      four_method_scores[[method]][four_method_common],
      method,
      heldout,
      cutoff = coverage_cutoff
    )
    result$evaluation_domain <- "four_method_common_primary"
    result
  })
  five_metric_rows <- lapply(names(five_method_scores), function(method) {
    result <- ranking_metrics(
      label[five_method_common],
      five_method_scores[[method]][five_method_common],
      method,
      heldout,
      cutoff = coverage_cutoff
    )
    result$evaluation_domain <- "five_method_IDR_common_sensitivity"
    result
  })
  metrics <- do.call(rbind, c(four_metric_rows, five_metric_rows))

  metric_path <- file.path(
    dirs$per_fold,
    paste0("Leave.", heldout, ".metrics.csv")
  )
  audit_path <- file.path(
    dirs$per_fold,
    paste0("Leave.", heldout, ".audit.csv")
  )
  write.csv(metrics, metric_path, row.names = FALSE)
  write.csv(audit, audit_path, row.names = FALSE)
  message("Saved evaluation: ", metric_path)

  invisible(list(metrics = metrics, audit = audit))
}

## ------------------------------------------------------------------
## Phase 1: retrain strict LOTO models under each coverage cutoff.
## Rank matrices are not held in memory during EM fitting.
## ------------------------------------------------------------------
if (!evaluate_only) {
  for (coverage_cutoff in cutoffs) {
    dirs <- make_dirs(coverage_cutoff)
    message("=== TRAINING: Total>", coverage_cutoff, " ===")
    for (heldout in folds) {
      fit_loto_fold(heldout, coverage_cutoff, dirs)
      gc()
    }
  }
}

## ------------------------------------------------------------------
## Phase 2: evaluate with held-out BBmix MAP labels under the same cutoff.
## ------------------------------------------------------------------
if (!train_only) {
  for (coverage_cutoff in cutoffs) {
    dirs <- make_dirs(coverage_cutoff)
    message("=== EVALUATION: Total>", coverage_cutoff, " ===")

    threshold_probability <- probability_all
    threshold_probability[
      !is.finite(total_all) | total_all <= coverage_cutoff
    ] <- NA_real_
    probability_count <- rowSums(is.finite(threshold_probability))
    probability_sum <- rowSums(threshold_probability, na.rm = TRUE)
    rank_matrix <- normalized_rank_matrix(
      threshold_probability,
      universe_size = n_sites
    )

    for (heldout in folds) {
      evaluate_loto_fold(
        heldout,
        coverage_cutoff,
        dirs,
        threshold_probability,
        probability_count,
        probability_sum,
        rank_matrix
      )
      gc()
    }

    rm(threshold_probability, probability_count, probability_sum, rank_matrix)
    gc()
  }
}

if (train_only) {
  message("Training-only mode completed.")
  quit(save = "no", status = 0L)
}

## ------------------------------------------------------------------
## Phase 3: summarize every completed fold currently in output_dir.
## ------------------------------------------------------------------
metric_files <- list.files(
  output_dir,
  pattern = "\\.metrics\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
audit_files <- list.files(
  output_dir,
  pattern = "\\.audit\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
if (!length(metric_files) || !length(audit_files)) {
  stop("No completed per-fold evaluation files found.")
}

metrics <- do.call(rbind, lapply(metric_files, read.csv, stringsAsFactors = FALSE))
audit <- do.call(rbind, lapply(audit_files, read.csv, stringsAsFactors = FALSE))
metrics <- metrics[order(metrics$coverage_cutoff, metrics$fold, metrics$method), ]
audit <- audit[order(audit$coverage_cutoff, audit$heldout), ]

if (run_all) {
  expected_metric_rows <- length(cutoffs) * length(technologies) * 9L
  expected_audit_rows <- length(cutoffs) * length(technologies)
  metric_key <- paste(
    metrics$coverage_cutoff,
    metrics$fold,
    metrics$evaluation_domain,
    metrics$method,
    sep = "\t"
  )
  audit_key <- paste(audit$coverage_cutoff, audit$heldout, sep = "\t")
  stopifnot(
    nrow(metrics) == expected_metric_rows,
    nrow(audit) == expected_audit_rows,
    !anyDuplicated(metric_key),
    !anyDuplicated(audit_key),
    setequal(
      unique(metrics$evaluation_domain),
      c(
        "four_method_common_primary",
        "five_method_IDR_common_sensitivity"
      )
    ),
    setequal(
      unique(metrics$method[metrics$evaluation_domain == "four_method_common_primary"]),
      c("ConquerBN", "MeanProb", "RankProduct", "RRA")
    ),
    setequal(
      unique(metrics$method[metrics$evaluation_domain == "five_method_IDR_common_sensitivity"]),
      c("ConquerBN", "MeanProb", "RankProduct", "RRA", "IDR")
    )
  )
}

write.csv(metrics, file.path(output_dir, "per_fold_metrics.csv"), row.names = FALSE)
write.csv(audit, file.path(output_dir, "evaluation_domain_audit.csv"), row.names = FALSE)

evaluable <- metrics[is.finite(metrics$auroc) & is.finite(metrics$auprc), ]
macro <- aggregate(
  cbind(auroc, auprc, n, positives, negatives) ~
    coverage_cutoff + coverage_rule + evaluation_domain + method,
  data = evaluable,
  FUN = mean
)
fold_counts <- aggregate(
  fold ~ coverage_cutoff + coverage_rule + evaluation_domain + method,
  data = evaluable,
  FUN = length
)
names(fold_counts)[names(fold_counts) == "fold"] <- "n_evaluable_folds"
macro <- merge(
  macro,
  fold_counts,
  by = c("coverage_cutoff", "coverage_rule", "evaluation_domain", "method")
)
macro <- macro[order(macro$method, macro$coverage_cutoff), ]
write.csv(macro, file.path(output_dir, "macro_average_summary.csv"), row.names = FALSE)

## ConquerBN versus every baseline, paired by held-out technology.
method_difference_rows <- list()
method_summary_rows <- list()
for (coverage_cutoff in cutoffs) {
  for (domain in unique(metrics$evaluation_domain)) {
    current <- metrics[
      metrics$coverage_cutoff == coverage_cutoff &
        metrics$evaluation_domain == domain,
      ,
      drop = FALSE
    ]
    conquer <- current[
      current$method == "ConquerBN",
      c("fold", "auroc", "auprc")
    ]
    names(conquer)[2:3] <- c("conquernet_auroc", "conquernet_auprc")
    for (baseline in setdiff(unique(current$method), "ConquerBN")) {
      other <- current[
        current$method == baseline,
        c("fold", "auroc", "auprc")
      ]
      names(other)[2:3] <- c("baseline_auroc", "baseline_auprc")
      paired <- merge(conquer, other, by = "fold")
      paired$coverage_cutoff <- coverage_cutoff
      paired$coverage_rule <- paste0("Total>", coverage_cutoff)
      paired$evaluation_domain <- domain
      paired$baseline_method <- baseline
      paired$delta_auroc_conquernet_minus_baseline <-
        paired$conquernet_auroc - paired$baseline_auroc
      paired$delta_auprc_conquernet_minus_baseline <-
        paired$conquernet_auprc - paired$baseline_auprc
      method_difference_rows[[length(method_difference_rows) + 1L]] <- paired

      finite_roc <- paired$delta_auroc_conquernet_minus_baseline
      finite_roc <- finite_roc[is.finite(finite_roc)]
      finite_pr <- paired$delta_auprc_conquernet_minus_baseline
      finite_pr <- finite_pr[is.finite(finite_pr)]
      safe_wilcox <- function(delta) {
        if (length(delta) < 2L || all(delta == 0)) return(NA_real_)
        suppressWarnings(tryCatch(
          stats::wilcox.test(
            delta,
            mu = 0,
            alternative = "two.sided",
            exact = FALSE
          )$p.value,
          error = function(e) NA_real_
        ))
      }
      method_summary_rows[[length(method_summary_rows) + 1L]] <- data.frame(
        coverage_cutoff = coverage_cutoff,
        coverage_rule = paste0("Total>", coverage_cutoff),
        evaluation_domain = domain,
        baseline_method = baseline,
        n_paired_auroc = length(finite_roc),
        mean_delta_auroc = mean(finite_roc),
        median_delta_auroc = median(finite_roc),
        n_auroc_wins = sum(finite_roc > 0),
        n_auroc_losses = sum(finite_roc < 0),
        wilcoxon_p_auroc = safe_wilcox(finite_roc),
        n_paired_auprc = length(finite_pr),
        mean_delta_auprc = mean(finite_pr),
        median_delta_auprc = median(finite_pr),
        n_auprc_wins = sum(finite_pr > 0),
        n_auprc_losses = sum(finite_pr < 0),
        wilcoxon_p_auprc = safe_wilcox(finite_pr),
        stringsAsFactors = FALSE
      )
    }
  }
}
method_differences <- do.call(rbind, method_difference_rows)
method_summary <- do.call(rbind, method_summary_rows)
adjustment_group <- interaction(
  method_summary$coverage_cutoff,
  method_summary$evaluation_domain,
  drop = TRUE,
  lex.order = TRUE
)
method_summary$wilcoxon_p_auroc_BH <- ave(
  method_summary$wilcoxon_p_auroc,
  adjustment_group,
  FUN = function(x) p.adjust(x, method = "BH")
)
method_summary$wilcoxon_p_auprc_BH <- ave(
  method_summary$wilcoxon_p_auprc,
  adjustment_group,
  FUN = function(x) p.adjust(x, method = "BH")
)
write.csv(
  method_differences,
  file.path(output_dir, "method_paired_differences.csv"),
  row.names = FALSE
)
write.csv(
  method_summary,
  file.path(output_dir, "method_paired_summary.csv"),
  row.names = FALSE
)

message("Completed coverage-aware strict LOTO training and evaluation.")
message("Results: ", output_dir)
