## Shared data, ranking and reporting utilities.

ensure_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing)) {
    stop("Missing required files:\n", paste(missing, collapse = "\n"))
  }
  invisible(paths)
}

conquernet_strategy_map <- function() {
  c(
    DART_seq = "EnzymeConversion",
    eTAM_seq = "EnzymeConversion",
    GLORI = "ChemicalConversion",
    m6A_REF_seq = "EnzymeCleavage",
    m6A_SAC_seq = "EnzymeConversion",
    m6ACE_seq = "Antibody",
    MAZTER_seq = "EnzymeCleavage",
    MeRIP_seq = "Antibody",
    Oxford_nanopore_rna002_m6Anet = "Nanopore",
    scDART_seq = "EnzymeConversion",
    CAM_seq = "ChemicalConversion",
    Oxford_nanopore_rna004_DeepRM = "Nanopore",
    Oxford_nanopore_rna004_SingleMod = "Nanopore",
    Oxford_nanopore_rna004_Dorado = "Nanopore"
  )
}

load_benchmark_dataset <- function(config) {
  path <- config$benchmark_file
  require_files(path)
  dataset <- readRDS(path)
  if (!inherits(dataset, "SummarizedExperiment")) {
    stop("Benchmark input is not a SummarizedExperiment.")
  }

  required <- c("m6A", "Total", "m6ASiteProb")
  missing <- setdiff(required, SummarizedExperiment::assayNames(dataset))
  if (length(missing)) {
    stop("Missing assays: ", paste(missing, collapse = ", "))
  }

  m6a <- SummarizedExperiment::assay(dataset, "m6A")
  total <- SummarizedExperiment::assay(dataset, "Total")
  probability <- SummarizedExperiment::assay(dataset, "m6ASiteProb")
  if (!identical(dim(m6a), dim(total)) ||
      !identical(dim(probability), dim(total)) ||
      !identical(dimnames(m6a), dimnames(total)) ||
      !identical(dimnames(probability), dimnames(total))) {
    stop("Benchmark assays are not aligned.")
  }

  technologies <- colnames(dataset)
  strategy_map <- conquernet_strategy_map()
  if (!setequal(technologies, names(strategy_map))) {
    stop("Dataset technologies do not match the fixed strategy map.")
  }
  dataset
}

normalized_rank_matrix <- function(score_matrix) {
  result <- matrix(
    1,
    nrow = nrow(score_matrix),
    ncol = ncol(score_matrix),
    dimnames = dimnames(score_matrix)
  )
  universe_size <- nrow(score_matrix)
  for (j in seq_len(ncol(score_matrix))) {
    observed <- is.finite(score_matrix[, j])
    if (any(observed)) {
      result[observed, j] <- rank(
        -score_matrix[observed, j],
        ties.method = "average"
      ) / universe_size
    }
  }
  result
}

rank_product_score <- function(rank_matrix) {
  -rowMeans(log10(pmax(rank_matrix, .Machine$double.xmin)))
}

rra_score <- function(rank_matrix) {
  k <- ncol(rank_matrix)
  minimum <- rep(1, nrow(rank_matrix))
  for (order_statistic in seq_len(k)) {
    ordered_rank <- matrixStats::rowOrderStats(
      rank_matrix,
      which = order_statistic
    )
    minimum <- pmin(
      minimum,
      stats::pbeta(
        ordered_rank,
        shape1 = order_statistic,
        shape2 = k - order_statistic + 1L
      )
    )
  }
  -log10(pmax(pmin(1, minimum * k), .Machine$double.xmin))
}

ranking_metrics <- function(y, score, method, unit, cutoff, family = NA_character_) {
  keep <- !is.na(y) & is.finite(score)
  y <- as.integer(y[keep])
  score <- as.numeric(score[keep])
  positives <- sum(y == 1L)
  negatives <- sum(y == 0L)
  auroc <- auprc <- NA_real_

  if (positives > 0L && negatives > 0L) {
    auroc <- as.numeric(pROC::auc(pROC::roc(
      response = y,
      predictor = score,
      direction = "<",
      quiet = TRUE
    )))
    auprc <- as.numeric(PRROC::pr.curve(
      scores.class0 = score[y == 1L],
      scores.class1 = score[y == 0L],
      curve = FALSE
    )$auc.integral)
  }

  data.frame(
    fold = unit,
    heldout_family = family,
    coverage_rule = paste0("Total>", cutoff),
    coverage_cutoff = cutoff,
    ground_truth = "bbmix_MAP_p_gt_0.5",
    method = method,
    n = length(y),
    positives = positives,
    negatives = negatives,
    prevalence = if (length(y)) mean(y) else NA_real_,
    auroc = auroc,
    auprc = auprc,
    stringsAsFactors = FALSE
  )
}

read_idr_baseline <- function(directory, unit, cutoff, n_sites) {
  paths <- c(
    rowmin = file.path(
      directory,
      paste0("Leave.", unit, ".idr.m6A.rowmin.cov", cutoff, ".rds")
    ),
    label = file.path(
      directory,
      paste0("Leave.", unit, ".idr.m6A.label.cov", cutoff, ".rds")
    ),
    coverage = file.path(
      directory,
      paste0("Leave.", unit, ".idr.idx.cov", cutoff, ".rds")
    )
  )
  require_files(paths)

  rowmin <- readRDS(paths[["rowmin"]])
  label <- readRDS(paths[["label"]])
  coverage <- readRDS(paths[["coverage"]])
  if (!is.numeric(rowmin) || !is.logical(label) || !is.logical(coverage) ||
      length(rowmin) != n_sites || length(label) != n_sites ||
      length(coverage) != n_sites) {
    stop("Invalid IDR objects for ", unit)
  }
  if (!identical(is.finite(rowmin), coverage) ||
      !identical(label, coverage & rowmin < 0.05)) {
    stop("IDR rowmin/label/coverage relationship failed for ", unit)
  }
  if (any(rowmin[coverage] < 0 | rowmin[coverage] > 1)) {
    stop("IDR rowmin is outside [0,1] for ", unit)
  }

  list(
    score = -rowmin,
    coverage = coverage,
    files = stats::setNames(basename(paths), names(paths))
  )
}

safe_wilcox <- function(delta, exact = FALSE) {
  delta <- delta[is.finite(delta)]
  if (length(delta) < 2L || all(delta == 0)) return(NA_real_)
  if (exact && (any(delta == 0) || anyDuplicated(abs(delta)))) return(NA_real_)
  suppressWarnings(tryCatch(
    stats::wilcox.test(
      delta,
      mu = 0,
      alternative = "two.sided",
      exact = exact
    )$p.value,
    error = function(e) NA_real_
  ))
}

safe_sign_test <- function(delta) {
  delta <- delta[is.finite(delta) & delta != 0]
  if (!length(delta)) return(NA_real_)
  stats::binom.test(sum(delta > 0), length(delta), p = 0.5)$p.value
}
