#!/usr/bin/env Rscript

# Week8 B supplement: effect sizes for posterior-derived FP versus TP
# meta-transcript and normalized full-gene distributions.
#
# Contract:
#   * class labels are the Week8 direct four-class posterior/global-FDR calls
#   * 0=TN, 1=FN, 2=FP, 3=TP, NA=unclassified
#   * FP=2 and TP=3 only; no MAP/Viterbi fallback
#   * every analysed site is independently rechecked to have Total > 10
#   * model/data source is the no-cali IVT-matched SummarizedExperiment
#
# Point estimates use every assigned site. Confidence intervals use a fixed,
# balanced, without-replacement subsample (at most 10,000 sites per class),
# followed by 500 percentile-bootstrap resamples within that fixed subsample.

options(stringsAsFactors = FALSE)

parse_cli <- function(args) {
  defaults <- list(
    class_matrix = paste0(
      "/home/xiyao/15data/W8/results/01_four_class_global_fdr/",
      "four_class_matrix_Total_gt_10.rds"
    ),
    features = paste0(
      "/home/xiyao/15data/W8/input_data/features/",
      "m6Aconquer_omicsFeatures_hg38.rds"
    ),
    benchmark = paste0(
      "/home/xichen/ConquerNet/ConquerBN_basic/",
      "m6AConquer_hg38_BenchmarkDataSet_Combined_SE_",
      "no_cali_IVT_matched.rds"
    ),
    output_dir = paste0(
      "/home/xiyao/15data/W8/results/02_fp_tp_meta_distributions/",
      "Total_gt_10_no_cali_IVT_matched/effect_sizes"
    ),
    coverage_cutoff = "10",
    reliability_min_per_class = "1000",
    bootstrap_max_per_class = "10000",
    bootstrap_repeats = "500",
    bootstrap_seed = "20260812"
  )

  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Arguments must use --name=value syntax: ", arg)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pieces[[1L]], fixed = TRUE)
    if (!key %in% names(defaults)) stop("Unknown argument: --", key)
    defaults[[key]] <- paste(pieces[-1L], collapse = "=")
  }
  defaults
}

extract_class_matrix <- function(object) {
  if (is.matrix(object) || is.data.frame(object)) return(as.matrix(object))
  candidates <- c(
    "unique_class_matrix", "class_matrix", "global_class_matrix",
    "global_fdr_class_matrix", "score"
  )
  for (name in candidates) {
    value <- object[[name]]
    if (is.matrix(value) || is.data.frame(value)) return(as.matrix(value))
  }
  stop("Could not find a four-class matrix in the class object.")
}

normalize_class_codes <- function(x) {
  original_dimnames <- dimnames(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    labels <- toupper(trimws(x))
    map <- c(TN = 0L, FN = 1L, FP = 2L, TP = 3L)
    unexpected <- setdiff(unique(labels[!is.na(labels)]), names(map))
    if (length(unexpected)) {
      stop("Unexpected class labels: ", paste(unexpected, collapse = ", "))
    }
    out <- matrix(unname(map[labels]), nrow = nrow(x), ncol = ncol(x))
  } else {
    out <- matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
    unexpected <- setdiff(unique(out[!is.na(out)]), 0:3)
    if (length(unexpected)) {
      stop("Unexpected class codes: ", paste(unexpected, collapse = ", "))
    }
  }
  dimnames(out) <- original_dimnames
  out
}

extract_feature_table <- function(object) {
  if (inherits(object, "GRanges")) return(S4Vectors::mcols(object))
  if (inherits(object, "SummarizedExperiment")) {
    return(SummarizedExperiment::rowData(object))
  }
  if (is.matrix(object) || is.data.frame(object)) return(object)
  stop("Unsupported OmicsFeatures class: ", paste(class(object), collapse = ", "))
}

extract_feature_ranges <- function(object) {
  if (inherits(object, "GRanges")) return(object)
  if (inherits(object, "SummarizedExperiment")) {
    return(SummarizedExperiment::rowRanges(object))
  }
  NULL
}

coordinate_signature_equal <- function(left, right) {
  if (is.null(left) || is.null(right)) return(NA)
  if (length(left) != length(right)) return(FALSE)
  identical(as.character(GenomicRanges::seqnames(left)),
            as.character(GenomicRanges::seqnames(right))) &&
    identical(GenomicRanges::start(left), GenomicRanges::start(right)) &&
    identical(GenomicRanges::end(left), GenomicRanges::end(right)) &&
    identical(as.character(GenomicRanges::strand(left)),
              as.character(GenomicRanges::strand(right)))
}

feature_assignment_mask <- function(feature_name, feature_table) {
  if (feature_name == "meta_tx_topology") {
    requested <- c(
      "overlap_exonicFivePrimeUTR", "overlap_exonicCDS",
      "overlap_exonicThreePrimeUTR"
    )
  } else if (feature_name == "relativePOS_fullGenes") {
    requested <- "overlap_fullGenes"
  } else {
    stop("No explicit assignment-mask rule for feature: ", feature_name)
  }

  missing <- setdiff(requested, colnames(feature_table))
  if (length(missing)) {
    stop(
      "Required annotation-mask column(s) missing for ", feature_name, ": ",
      paste(missing, collapse = ", ")
    )
  }
  annotation <- as.data.frame(feature_table[, requested, drop = FALSE])
  annotation[] <- lapply(annotation, function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[!is.finite(x)] <- 0
    x
  })
  list(
    mask = rowSums(as.matrix(annotation) > 0) > 0,
    method = paste0("explicit overlap mask: ", paste(requested, collapse = " OR "))
  )
}

safe_quantile <- function(x, probability) {
  if (!length(x)) return(NA_real_)
  as.numeric(stats::quantile(x, probability, names = FALSE, type = 7))
}

# Positive Cliff's delta means TP values tend to be greater than FP values.
cliffs_delta <- function(fp, tp) {
  n_fp <- length(fp)
  n_tp <- length(tp)
  if (!n_fp || !n_tp) return(NA_real_)
  pooled_ranks <- rank(c(fp, tp), ties.method = "average")
  rank_sum_tp <- sum(tail(pooled_ranks, n_tp))
  ## Coerce before multiplication: several formal classes have >50,000
  ## observations and integer n_tp*n_fp can exceed R's 32-bit integer limit.
  n_tp_double <- as.double(n_tp)
  n_fp_double <- as.double(n_fp)
  u_tp <- rank_sum_tp - n_tp_double * (n_tp_double + 1) / 2
  auc_tp_greater <- u_tp / (n_tp_double * n_fp_double)
  2 * auc_tp_greater - 1
}

ks_d <- function(fp, tp) {
  if (!length(fp) || !length(tp)) return(NA_real_)
  unname(suppressWarnings(stats::ks.test(fp, tp, exact = FALSE)$statistic))
}

point_statistics <- function(fp, tp) {
  fp_q25 <- safe_quantile(fp, 0.25)
  fp_q75 <- safe_quantile(fp, 0.75)
  tp_q25 <- safe_quantile(tp, 0.25)
  tp_q75 <- safe_quantile(tp, 0.75)
  fp_mean <- if (length(fp)) mean(fp) else NA_real_
  tp_mean <- if (length(tp)) mean(tp) else NA_real_
  fp_median <- if (length(fp)) stats::median(fp) else NA_real_
  tp_median <- if (length(tp)) stats::median(tp) else NA_real_
  c(
    fp_mean = fp_mean,
    fp_median = fp_median,
    fp_q25 = fp_q25,
    fp_q75 = fp_q75,
    fp_iqr = fp_q75 - fp_q25,
    tp_mean = tp_mean,
    tp_median = tp_median,
    tp_q25 = tp_q25,
    tp_q75 = tp_q75,
    tp_iqr = tp_q75 - tp_q25,
    mean_difference_tp_minus_fp = tp_mean - fp_mean,
    median_difference_tp_minus_fp = tp_median - fp_median,
    ks_d = ks_d(fp, tp),
    cliffs_delta_tp_greater = cliffs_delta(fp, tp)
  )
}

balanced_bootstrap_ci <- function(fp, tp, max_per_class, repeats) {
  balanced_n <- min(length(fp), length(tp), max_per_class)
  empty <- c(
    mean_difference_ci_low = NA_real_, mean_difference_ci_high = NA_real_,
    median_difference_ci_low = NA_real_, median_difference_ci_high = NA_real_,
    ks_d_ci_low = NA_real_, ks_d_ci_high = NA_real_,
    cliffs_delta_ci_low = NA_real_, cliffs_delta_ci_high = NA_real_
  )
  if (balanced_n < 2L || repeats < 1L) {
    return(list(
      balanced_n = balanced_n,
      ci = empty,
      status = "not_computed_min_assigned_class_lt_2"
    ))
  }

  fp_fixed <- if (length(fp) > balanced_n) sample(fp, balanced_n) else fp
  tp_fixed <- if (length(tp) > balanced_n) sample(tp, balanced_n) else tp
  estimates <- matrix(NA_real_, nrow = repeats, ncol = 4L)
  colnames(estimates) <- c(
    "mean_difference", "median_difference", "ks_d", "cliffs_delta"
  )
  for (iteration in seq_len(repeats)) {
    fp_boot <- sample(fp_fixed, balanced_n, replace = TRUE)
    tp_boot <- sample(tp_fixed, balanced_n, replace = TRUE)
    estimates[iteration, ] <- c(
      mean(tp_boot) - mean(fp_boot),
      stats::median(tp_boot) - stats::median(fp_boot),
      ks_d(fp_boot, tp_boot),
      cliffs_delta(fp_boot, tp_boot)
    )
  }
  interval <- apply(
    estimates, 2L, stats::quantile,
    probs = c(0.025, 0.975), names = FALSE, type = 7, na.rm = TRUE
  )
  list(
    balanced_n = balanced_n,
    ci = c(
      ## Drop the matrix element name. Otherwise c(new_name=named_scalar)
      ## creates names such as mean_difference_ci_low.mean_difference,
      ## which differ from the non-computable-row schema and break rbind().
      mean_difference_ci_low = as.numeric(interval[1L, "mean_difference"]),
      mean_difference_ci_high = as.numeric(interval[2L, "mean_difference"]),
      median_difference_ci_low = as.numeric(interval[1L, "median_difference"]),
      median_difference_ci_high = as.numeric(interval[2L, "median_difference"]),
      ks_d_ci_low = as.numeric(interval[1L, "ks_d"]),
      ks_d_ci_high = as.numeric(interval[2L, "ks_d"]),
      cliffs_delta_ci_low = as.numeric(interval[1L, "cliffs_delta"]),
      cliffs_delta_ci_high = as.numeric(interval[2L, "cliffs_delta"])
    ),
    status = "computed_balanced_percentile_bootstrap"
  )
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
coverage_cutoff <- as.numeric(args$coverage_cutoff)
reliability_min <- as.integer(args$reliability_min_per_class)
bootstrap_max <- as.integer(args$bootstrap_max_per_class)
bootstrap_repeats <- as.integer(args$bootstrap_repeats)
bootstrap_seed <- as.integer(args$bootstrap_seed)
if (!identical(coverage_cutoff, 10)) {
  stop("This formal Week8 B analysis requires coverage_cutoff=10 exactly.")
}
if (is.na(reliability_min) || reliability_min < 1L ||
    is.na(bootstrap_max) || bootstrap_max < 2L ||
    is.na(bootstrap_repeats) || bootstrap_repeats < 1L ||
    is.na(bootstrap_seed)) {
  stop("Invalid reliability/bootstrap settings.")
}

paths <- c(
  class_matrix = args$class_matrix,
  features = args$features,
  benchmark = args$benchmark
)
missing_paths <- names(paths)[!file.exists(paths)]
if (length(missing_paths)) {
  stop(
    "Missing required input(s):\n",
    paste(paste0("  ", missing_paths, ": ", paths[missing_paths]), collapse = "\n")
  )
}
required_packages <- c("S4Vectors", "GenomicRanges", "SummarizedExperiment")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

message("Loading formal Week8 class matrix...")
class_object <- readRDS(args$class_matrix)
class_matrix <- normalize_class_codes(extract_class_matrix(class_object))
message("Loading matched OmicsFeatures and no-cali IVT-matched SE...")
feature_object <- readRDS(args$features)
feature_table <- extract_feature_table(feature_object)
benchmark <- readRDS(args$benchmark)

if (!inherits(benchmark, "SummarizedExperiment")) {
  stop("Benchmark is not a SummarizedExperiment.")
}
if (!"Total" %in% SummarizedExperiment::assayNames(benchmark)) {
  stop("Benchmark lacks a Total assay.")
}
if (nrow(class_matrix) != nrow(feature_table) ||
    nrow(class_matrix) != nrow(benchmark)) {
  stop("Class matrix, OmicsFeatures, and benchmark row counts differ.")
}
if (is.null(colnames(class_matrix)) || anyDuplicated(colnames(class_matrix))) {
  stop("Class matrix must have unique technology column names.")
}

total <- SummarizedExperiment::assay(benchmark, "Total")
technologies <- colnames(class_matrix)
if (!all(technologies %in% colnames(total))) {
  stop(
    "Total assay lacks technologies: ",
    paste(setdiff(technologies, colnames(total)), collapse = ", ")
  )
}
coordinate_identity <- coordinate_signature_equal(
  extract_feature_ranges(feature_object),
  SummarizedExperiment::rowRanges(benchmark)
)
if (!identical(coordinate_identity, TRUE)) {
  stop("OmicsFeatures and benchmark coordinates are not exactly identical.")
}

feature_names <- c("meta_tx_topology", "relativePOS_fullGenes")
missing_features <- setdiff(feature_names, colnames(feature_table))
if (length(missing_features)) {
  stop("Missing feature(s): ", paste(missing_features, collapse = ", "))
}
feature_inputs <- setNames(lapply(feature_names, function(feature_name) {
  values <- suppressWarnings(as.numeric(feature_table[[feature_name]]))
  finite <- is.finite(values)
  in_range <- finite & values >= -1e-12 & values <= 1 + 1e-12
  values[in_range] <- pmin(1, pmax(0, values[in_range]))
  assignment <- feature_assignment_mask(feature_name, feature_table)
  list(
    values = values,
    finite = finite,
    in_range = in_range,
    assignment_mask = assignment$mask,
    assignment_method = assignment$method,
    usable = in_range & assignment$mask
  )
}), feature_names)

set.seed(bootstrap_seed)
results <- list()
result_index <- 0L
expected_result_schema <- NULL

for (technology in technologies) {
  message("Effect sizes: ", technology)
  class_vector <- class_matrix[, technology]
  coverage <- suppressWarnings(as.numeric(total[, technology]))
  coverage_mask <- is.finite(coverage) & coverage > coverage_cutoff
  fp_in_matrix <- !is.na(class_vector) & class_vector == 2L
  tp_in_matrix <- !is.na(class_vector) & class_vector == 3L
  fp_base <- fp_in_matrix & coverage_mask
  tp_base <- tp_in_matrix & coverage_mask

  # The formal Total>10 matrix must already contain no sub-threshold calls.
  if (sum(fp_in_matrix) != sum(fp_base) || sum(tp_in_matrix) != sum(tp_base)) {
    stop(
      "Class matrix contains FP/TP outside Total>10 for technology: ",
      technology
    )
  }
  raw_status <- if (sum(fp_base) == 0L || sum(tp_base) == 0L) {
    "untrainable"
  } else if (sum(fp_base) >= reliability_min && sum(tp_base) >= reliability_min) {
    "reliable"
  } else {
    "exploratory"
  }

  for (feature_name in feature_names) {
    input <- feature_inputs[[feature_name]]
    fp_mask <- fp_base & input$usable
    tp_mask <- tp_base & input$usable
    fp <- input$values[fp_mask]
    tp <- input$values[tp_mask]
    statistics <- point_statistics(fp, tp)
    bootstrap <- balanced_bootstrap_ci(
      fp, tp,
      max_per_class = bootstrap_max,
      repeats = bootstrap_repeats
    )

    result_index <- result_index + 1L
    result_row <- data.frame(
      technology = technology,
      feature = feature_name,
      analysis_status = raw_status,
      coverage_rule = "Total>10",
      label_rule = "direct posterior global-FDR; FP=2; TP=3; no fallback",
      n_fp_label_matrix = sum(fp_in_matrix),
      n_tp_label_matrix = sum(tp_in_matrix),
      n_fp_after_Total_gt_10_recheck = sum(fp_base),
      n_tp_after_Total_gt_10_recheck = sum(tp_base),
      n_fp_assigned = length(fp),
      n_tp_assigned = length(tp),
      fp_assignment_rate = if (sum(fp_base)) length(fp) / sum(fp_base) else NA_real_,
      tp_assignment_rate = if (sum(tp_base)) length(tp) / sum(tp_base) else NA_real_,
      assignment_rate_difference_tp_minus_fp =
        if (sum(fp_base) && sum(tp_base)) {
          length(tp) / sum(tp_base) - length(fp) / sum(fp_base)
        } else NA_real_,
      n_fp_missing_or_nonfinite = sum(fp_base & !input$finite),
      n_tp_missing_or_nonfinite = sum(tp_base & !input$finite),
      n_fp_outside_0_1 = sum(fp_base & input$finite & !input$in_range),
      n_tp_outside_0_1 = sum(tp_base & input$finite & !input$in_range),
      n_fp_outside_annotation_mask = sum(
        fp_base & input$in_range & !input$assignment_mask
      ),
      n_tp_outside_annotation_mask = sum(
        tp_base & input$in_range & !input$assignment_mask
      ),
      assignment_mask_method = input$assignment_method,
      as.list(statistics),
      bootstrap_seed = bootstrap_seed,
      bootstrap_max_per_class = bootstrap_max,
      bootstrap_repeats = bootstrap_repeats,
      bootstrap_balanced_n_per_class = bootstrap$balanced_n,
      bootstrap_status = bootstrap$status,
      as.list(bootstrap$ci),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    if (is.null(expected_result_schema)) {
      expected_result_schema <- names(result_row)
    } else if (!identical(names(result_row), expected_result_schema)) {
      stop(
        "Internal result-schema mismatch for ", technology, " / ",
        feature_name, ". Expected: ",
        paste(expected_result_schema, collapse = " | "),
        "; observed: ", paste(names(result_row), collapse = " | ")
      )
    }
    results[[result_index]] <- result_row
  }
  rm(class_vector, coverage, coverage_mask, fp_in_matrix, tp_in_matrix)
  invisible(gc())
}

effect_sizes <- do.call(rbind, results)
if (nrow(effect_sizes) != length(technologies) * length(feature_names)) {
  stop("Unexpected number of result rows.")
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  effect_sizes,
  file.path(args$output_dir, "fp_tp_meta_effect_sizes.csv"),
  row.names = FALSE
)
write.csv(
  effect_sizes[effect_sizes$analysis_status == "reliable", , drop = FALSE],
  file.path(args$output_dir, "fp_tp_meta_effect_sizes_reliable_only.csv"),
  row.names = FALSE
)

status_summary <- aggregate(
  technology ~ analysis_status,
  data = unique(effect_sizes[c("technology", "analysis_status")]),
  FUN = length
)
names(status_summary)[2L] <- "n_technologies"
write.csv(
  status_summary,
  file.path(args$output_dir, "fp_tp_meta_effect_size_status_summary.csv"),
  row.names = FALSE
)

qc <- data.frame(
  item = c(
    "class_matrix_path", "features_path", "benchmark_path", "output_dir",
    "n_universe_sites", "n_technologies", "n_features",
    "coordinate_identity", "coverage_cutoff", "coverage_operator",
    "class_code_contract", "fallback", "reliability_min_per_class",
    "bootstrap_seed", "bootstrap_max_per_class", "bootstrap_repeats",
    "bootstrap_design", "point_estimate_design",
    "all_label_counts_unchanged_after_coverage_recheck",
    "all_assigned_values_finite_and_in_0_1"
  ),
  value = c(
    normalizePath(args$class_matrix, winslash = "/", mustWork = TRUE),
    normalizePath(args$features, winslash = "/", mustWork = TRUE),
    normalizePath(args$benchmark, winslash = "/", mustWork = TRUE),
    normalizePath(args$output_dir, winslash = "/", mustWork = TRUE),
    nrow(class_matrix), length(technologies), length(feature_names),
    coordinate_identity, coverage_cutoff, ">",
    "0=TN;1=FN;2=FP;3=TP;NA=unclassified", "none",
    reliability_min, bootstrap_seed, bootstrap_max, bootstrap_repeats,
    paste0(
      "fixed balanced without-replacement subsample, then percentile bootstrap ",
      "with replacement; 95% CI"
    ),
    "all assigned sites",
    all(
      effect_sizes$n_fp_label_matrix ==
        effect_sizes$n_fp_after_Total_gt_10_recheck &
      effect_sizes$n_tp_label_matrix ==
        effect_sizes$n_tp_after_Total_gt_10_recheck
    ),
    all(
      effect_sizes$n_fp_missing_or_nonfinite == 0 &
      effect_sizes$n_tp_missing_or_nonfinite == 0 &
      effect_sizes$n_fp_outside_0_1 == 0 &
      effect_sizes$n_tp_outside_0_1 == 0
    )
  ),
  stringsAsFactors = FALSE
)
write.csv(
  qc,
  file.path(args$output_dir, "fp_tp_meta_effect_size_qc.csv"),
  row.names = FALSE
)

message("Saved Week8 B effect-size supplement: ", args$output_dir)
print(status_summary, row.names = FALSE)
