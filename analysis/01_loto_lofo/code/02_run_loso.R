suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(bnlearn)
  library(extraDistr)
  library(matrixStats)
})

## W6 C: coverage-aware Leave-One-Strategy/Family-Out training
##
## Final formal route: Total > 10 for each held-out technology family:
##   1. Remove the complete held-out family from model fitting.
##   2. Among the remaining four families, retain a technology-site observation
##      only when Total > cutoff.
##   3. Represent below-threshold observations as m6A=0, Total=0 and
##      m6ASiteProb=NA so that they contribute no emission evidence.
##   4. Retrain a new ConquerBN model and save P(TrueM6A=1).
##
## This produces 5 independently fitted family-out models. Evaluation is
## performed separately by
## 03_evaluate_loso_map_coverage.R.

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
  "/home/xiyao/15data/W6_corrected/results/02_loso_models"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

family_arg <- cli_value(cli, "family", NULL)
run_all <- cli_flag(cli, "all")
if (is.null(family_arg) && !run_all) {
  stop("Specify --all or exactly one --family=FAMILY.")
}
if (!is.null(family_arg) && run_all) {
  stop("Use either --all or --family=FAMILY, not both.")
}

cutoff_arg <- cli_value(cli, "cutoff", NULL)
cutoffs_arg <- cli_value(cli, "cutoffs", NULL)
if (!is.null(cutoff_arg) && !is.null(cutoffs_arg)) {
  stop("Use --cutoff or --cutoffs, not both.")
}
cutoff_text <- if (!is.null(cutoff_arg)) {
  cutoff_arg
} else if (!is.null(cutoffs_arg)) {
  cutoffs_arg
} else {
  "10"
}
cutoffs <- as.integer(strsplit(cutoff_text, ",", fixed = TRUE)[[1]])
cutoffs <- sort(unique(cutoffs[is.finite(cutoffs)]))
if (!identical(cutoffs, 10L)) stop("This review script only permits cutoff 10.")

max_iter <- as.integer(cli_value(cli, "max-iter", "500"))
epsilon <- as.numeric(cli_value(cli, "epsilon", "1e-6"))
seed <- as.integer(cli_value(cli, "seed", "20260725"))
smoke_n <- as.integer(cli_value(cli, "smoke-n", "0"))
overwrite <- cli_flag(cli, "overwrite")

coverage_tag <- function(cutoff) paste0("Total_gt_", cutoff)

extract_root_score <- function(root_belief) {
  if (ncol(root_belief) != 2L) {
    stop("TrueM6A belief must contain exactly two states.")
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

families <- unique(unname(strategy_map[technologies]))
if (!is.null(family_arg) && !family_arg %in% families) {
  stop("Unknown family: ", family_arg)
}
heldout_families <- if (run_all) families else family_arg

message("Dataset: ", n_sites, " sites x ", length(technologies), " technologies")
message("Coverage cutoffs: ", paste0("Total>", cutoffs, collapse = ", "))
message("Held-out families: ", paste(heldout_families, collapse = ", "))

run_family <- function(heldout_family, coverage_cutoff) {
  cutoff_dir <- file.path(output_dir, coverage_tag(coverage_cutoff))
  dir.create(cutoff_dir, recursive = TRUE, showWarnings = FALSE)

  result_path <- file.path(
    cutoff_dir,
    paste0("LeaveFamily.", heldout_family, ".conquernet.rds")
  )

  if (file.exists(result_path) && !overwrite) {
    existing <- readRDS(result_path)
    valid_existing <-
      identical(existing$analysis_version, "coverage_strict_loso_v3") &&
      identical(as.character(existing$heldout_family), heldout_family) &&
      identical(as.integer(existing$coverage_cutoff), as.integer(coverage_cutoff)) &&
      identical(as.integer(existing$smoke_n), as.integer(smoke_n)) &&
      isTRUE(existing$converged) &&
      length(existing$row_index) == length(existing$score) &&
      !anyDuplicated(existing$row_index) &&
      all(is.finite(existing$score)) &&
      all(existing$score >= 0 & existing$score <= 1)
    if (!valid_existing) {
      stop("Existing LOSO result is incompatible; rerun with --overwrite: ", result_path)
    }
    message("Using existing model: ", result_path)
    return(invisible(result_path))
  }

  heldout_technologies <- names(strategy_map)[strategy_map == heldout_family]
  kept <- setdiff(technologies, heldout_technologies)
  kept_index <- match(kept, technologies)
  family_index <- match(heldout_family, families)

  ## Apply the coverage rule before selecting training sites.
  kept_coverage <- is.finite(total_all[, kept_index, drop = FALSE]) &
    total_all[, kept_index, drop = FALSE] > coverage_cutoff
  training_available <- rowSums(kept_coverage) > 0L
  row_index <- which(training_available)

  if (smoke_n > 0L && length(row_index) > smoke_n) {
    set.seed(seed + 1000L * coverage_cutoff + family_index)
    row_index <- sort(sample(row_index, smoke_n))
  }

  message(
    "Fitting held-out family ", heldout_family,
    " at Total>", coverage_cutoff,
    ": holding out ", paste(heldout_technologies, collapse = ", "),
    "; fitting ", length(row_index), " sites with ", length(kept),
    " technologies from four families."
  )

  m6a <- m6a_all[row_index, kept_index, drop = FALSE]
  total <- total_all[row_index, kept_index, drop = FALSE]
  weight <- probability_all[row_index, kept_index, drop = FALSE]
  observed <- is.finite(total) & total > coverage_cutoff & is.finite(m6a)

  ## At Total>cutoff, m6ASiteProb is expected to be finite. Stop rather than
  ## silently using a different input domain for ConquerBN and the baselines.
  if (any(observed & (!is.finite(m6a) | !is.finite(weight)))) {
    stop(
      "Found Total>", coverage_cutoff,
      " observations with missing m6A or m6ASiteProb in family ",
      heldout_family
    )
  }

  ## Below-threshold observations are absent from the likelihood.
  m6a[!observed] <- 0
  total[!observed] <- 0
  weight[!observed] <- NA_real_

  if (any(m6a < 0, na.rm = TRUE) ||
      any(total < 0, na.rm = TRUE) ||
      any(m6a > total, na.rm = TRUE)) {
    stop(
      "Invalid m6A/Total counts after coverage filtering: ",
      heldout_family, " at Total>", coverage_cutoff
    )
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
      stop(
        "Non-finite log-likelihood in ", heldout_family,
        " at Total>", coverage_cutoff,
        ", iteration ", iteration
      )
    }

    log_likelihood[iteration] <- current
    cpt_param <- conquerNet_maximization_safe(cpt_param, belief, dag)
    relative_change <- if (iteration > 1L) {
      abs(current - previous) / max(abs(previous), .Machine$double.eps)
    } else {
      NA_real_
    }

    message(sprintf(
      "  cutoff=%d family=%s iteration=%d loglik=%.6f relative_change=%s elapsed_min=%.3f",
      coverage_cutoff,
      heldout_family,
      iteration,
      current,
      ifelse(
        is.na(relative_change),
        "NA",
        format(relative_change, scientific = TRUE, digits = 4)
      ),
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
    stop(
      "LOSO model did not converge: ", heldout_family,
      " at Total>", coverage_cutoff
    )
  }

  path <- sapply(
    bnlearn::leaf.nodes(dag),
    function(node) get_parent_chain_bn(node, dag)
  )
  cpt_param <- latent_param_label_switch_fixed(cpt_param, path)
  final_expectation <- conquerNet_expectation(emission_loglik, cpt_param, dag)
  root_result <- extract_root_score(final_expectation$belief[["TrueM6A"]])

  result <- list(
    analysis_version = "coverage_strict_loso_v3",
    heldout_family = heldout_family,
    heldout_technologies = heldout_technologies,
    kept_families = setdiff(families, heldout_family),
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

for (coverage_cutoff in cutoffs) {
  message("=== LOSO TRAINING: Total>", coverage_cutoff, " ===")
  for (heldout_family in heldout_families) {
    run_family(heldout_family, coverage_cutoff)
    gc()
  }
}

message("Completed coverage-aware LOSO model fitting.")
message("Results: ", output_dir)
