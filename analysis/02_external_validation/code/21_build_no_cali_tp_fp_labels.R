suppressPackageStartupMessages({
  library(SummarizedExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
cli_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) default else sub(prefix, "", hit[[1]], fixed = TRUE)
}

dataset_path <- normalizePath(
  cli_value(
    "dataset",
    paste0(
      "/home/xichen/ConquerNet/ConquerBN_basic/",
      "m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali.rds"
    )
  ),
  winslash = "/",
  mustWork = TRUE
)
model_dir <- normalizePath(
  cli_value(
    "model-dir",
    "/home/xiyao/15data/W7/03_no_cali_full_model/Total_gt_10"
  ),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- cli_value(
  "output-dir",
  "/home/xiyao/15data/W7/04_no_cali_tp_fp_labels/Total_gt_10"
)
coverage_cutoff <- as.integer(cli_value("cutoff", "10"))
alpha <- as.numeric(cli_value("alpha", "0.05"))
min_reliable <- as.integer(cli_value("min-reliable", "1000"))

if (!identical(coverage_cutoff, 10L)) {
  stop("Week7 formal analysis is fixed to --cutoff=10.")
}
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
  stop("--alpha must be between 0 and 1.")
}
if (!is.finite(min_reliable) || min_reliable < 1L) {
  stop("--min-reliable must be a positive integer.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

posterior_path <- file.path(model_dir, "node_state1_posterior.rds")
viterbi_path <- file.path(model_dir, "viterbi_class_matrix.rds")
if (!all(file.exists(c(posterior_path, viterbi_path)))) {
  stop("Full no-cali posterior or Viterbi output is missing.")
}

dataset <- readRDS(dataset_path)
posterior_object <- readRDS(posterior_path)
viterbi_object <- readRDS(viterbi_path)

if (!inherits(dataset, "SummarizedExperiment")) {
  stop("Dataset is not a SummarizedExperiment.")
}

technologies <- colnames(dataset)
n_sites <- nrow(dataset)
total <- assay(dataset, "Total")
posterior <- posterior_object$node_state1_probability
viterbi_class <- viterbi_object$class_matrix
technology_nodes <- paste0(technologies, "_k")

stopifnot(
  identical(posterior_object$dataset_path, dataset_path),
  identical(viterbi_object$dataset_path, dataset_path),
  identical(as.integer(posterior_object$coverage_cutoff), coverage_cutoff),
  identical(as.integer(viterbi_object$coverage_cutoff), coverage_cutoff),
  identical(dim(posterior), c(n_sites, 1L + length(technologies))),
  identical(colnames(posterior), c("TrueM6A", technology_nodes)),
  identical(dim(viterbi_class), c(n_sites, length(technologies))),
  identical(colnames(viterbi_class), technologies)
)

if (
  any(
    !is.na(viterbi_class) &
      !viterbi_class %in% 0:3
  )
) {
  stop("Unexpected Viterbi class code.")
}

global_select <- function(local_error, alpha) {
  selected <- rep(FALSE, length(local_error))
  valid_index <- which(is.finite(local_error))
  if (!length(valid_index)) {
    return(list(
      selected = selected,
      n_selected = 0L,
      local_error_cutoff = NA_real_,
      global_error_at_cutoff = NA_real_
    ))
  }
  ordered_index <- valid_index[order(local_error[valid_index])]
  ordered_error <- local_error[ordered_index]
  global_error <- cumsum(ordered_error) / seq_along(ordered_error)
  valid_k <- which(global_error < alpha)
  if (!length(valid_k)) {
    return(list(
      selected = selected,
      n_selected = 0L,
      local_error_cutoff = NA_real_,
      global_error_at_cutoff = NA_real_
    ))
  }
  k <- max(valid_k)
  selected[ordered_index[seq_len(k)]] <- TRUE
  list(
    selected = selected,
    n_selected = k,
    local_error_cutoff = ordered_error[[k]],
    global_error_at_cutoff = global_error[[k]]
  )
}

count_classes <- function(class_vector, mask) {
  values <- class_vector[mask]
  counts <- table(factor(values, levels = 0:3))
  c(
    n_domain = sum(mask),
    n_labeled = sum(!is.na(values)),
    n_unclassified = sum(is.na(values)),
    n_TN = as.integer(counts[[1]]),
    n_FN = as.integer(counts[[2]]),
    n_FP = as.integer(counts[[3]]),
    n_TP = as.integer(counts[[4]])
  )
}

root_p1 <- posterior[, "TrueM6A"]
root_positive_result <- global_select(1 - root_p1, alpha)
root_negative_result <- global_select(root_p1, alpha)
root_positive <- root_positive_result$selected
root_negative <- root_negative_result$selected
root_conflict <- root_positive & root_negative

threshold_rows <- list(
  data.frame(
    node = "TrueM6A",
    state = "positive",
    local_error = "1-p",
    alpha = alpha,
    n_selected = root_positive_result$n_selected,
    local_error_cutoff = root_positive_result$local_error_cutoff,
    global_error_at_cutoff = root_positive_result$global_error_at_cutoff,
    n_conflict_with_opposite = sum(root_conflict),
    stringsAsFactors = FALSE
  ),
  data.frame(
    node = "TrueM6A",
    state = "negative",
    local_error = "p",
    alpha = alpha,
    n_selected = root_negative_result$n_selected,
    local_error_cutoff = root_negative_result$local_error_cutoff,
    global_error_at_cutoff = root_negative_result$global_error_at_cutoff,
    n_conflict_with_opposite = sum(root_conflict),
    stringsAsFactors = FALSE
  )
)

global_class_matrix <- matrix(
  NA_integer_,
  nrow = n_sites,
  ncol = length(technologies),
  dimnames = list(NULL, technologies)
)
selected_class_matrix <- global_class_matrix
decision_rows <- vector("list", length(technologies))

for (j in seq_along(technologies)) {
  technology <- technologies[[j]]
  technology_node <- technology_nodes[[j]]
  technology_p1 <- posterior[, technology_node]

  technology_positive_result <- global_select(1 - technology_p1, alpha)
  technology_negative_result <- global_select(technology_p1, alpha)
  technology_positive <- technology_positive_result$selected
  technology_negative <- technology_negative_result$selected
  technology_conflict <- technology_positive & technology_negative

  threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
    node = technology_node,
    state = "positive",
    local_error = "1-p",
    alpha = alpha,
    n_selected = technology_positive_result$n_selected,
    local_error_cutoff = technology_positive_result$local_error_cutoff,
    global_error_at_cutoff = technology_positive_result$global_error_at_cutoff,
    n_conflict_with_opposite = sum(technology_conflict),
    stringsAsFactors = FALSE
  )
  threshold_rows[[length(threshold_rows) + 1L]] <- data.frame(
    node = technology_node,
    state = "negative",
    local_error = "p",
    alpha = alpha,
    n_selected = technology_negative_result$n_selected,
    local_error_cutoff = technology_negative_result$local_error_cutoff,
    global_error_at_cutoff = technology_negative_result$global_error_at_cutoff,
    n_conflict_with_opposite = sum(technology_conflict),
    stringsAsFactors = FALSE
  )

  invalid <- root_conflict | technology_conflict
  global_class <- rep(NA_integer_, n_sites)
  global_class[root_negative & technology_negative & !invalid] <- 0L
  global_class[root_positive & technology_negative & !invalid] <- 1L
  global_class[root_negative & technology_positive & !invalid] <- 2L
  global_class[root_positive & technology_positive & !invalid] <- 3L
  global_class_matrix[, j] <- global_class

  observed <- is.finite(total[, j]) & total[, j] > coverage_cutoff
  global_counts <- count_classes(global_class, observed)
  map_counts <- count_classes(viterbi_class[, j], observed)

  use_global <- (
    global_counts[["n_TP"]] >= min_reliable &&
      global_counts[["n_FP"]] >= min_reliable
  )
  selected_source <- if (use_global) {
    "global_FDR_FOR"
  } else {
    "MAP_Viterbi_fallback"
  }
  selected_class <- if (use_global) global_class else viterbi_class[, j]
  selected_class[!observed] <- NA_integer_
  selected_class_matrix[, j] <- selected_class
  selected_counts <- count_classes(selected_class, observed)

  status <- if (
    selected_counts[["n_TP"]] == 0L ||
      selected_counts[["n_FP"]] == 0L
  ) {
    "untrainable_single_class"
  } else if (
    selected_counts[["n_TP"]] >= min_reliable &&
      selected_counts[["n_FP"]] >= min_reliable
  ) {
    "reliably_identifiable"
  } else {
    "exploratory_low_minority"
  }

  decision_rows[[j]] <- data.frame(
    technology = technology,
    coverage_rule = "Total>10",
    n_observed = sum(observed),
    global_n_TN = global_counts[["n_TN"]],
    global_n_FN = global_counts[["n_FN"]],
    global_n_FP = global_counts[["n_FP"]],
    global_n_TP = global_counts[["n_TP"]],
    global_n_unclassified = global_counts[["n_unclassified"]],
    map_n_TN = map_counts[["n_TN"]],
    map_n_FN = map_counts[["n_FN"]],
    map_n_FP = map_counts[["n_FP"]],
    map_n_TP = map_counts[["n_TP"]],
    selected_label_source = selected_source,
    selected_n_TN = selected_counts[["n_TN"]],
    selected_n_FN = selected_counts[["n_FN"]],
    selected_n_FP = selected_counts[["n_FP"]],
    selected_n_TP = selected_counts[["n_TP"]],
    selected_n_unclassified = selected_counts[["n_unclassified"]],
    tp_fp_status = status,
    stringsAsFactors = FALSE
  )
}

threshold_table <- do.call(rbind, threshold_rows)
decisions <- do.call(rbind, decision_rows)

saveRDS(
  list(
    analysis_version = "week7_no_cali_global_labels_Total_gt_10_v1",
    dataset_path = dataset_path,
    coverage_cutoff = coverage_cutoff,
    alpha = alpha,
    class_code_definition = c(TN = 0L, FN = 1L, FP = 2L, TP = 3L),
    class_matrix = global_class_matrix
  ),
  file.path(output_dir, "global_fdr_for_class_matrix.rds"),
  compress = TRUE
)
saveRDS(
  list(
    analysis_version = "week7_no_cali_selected_labels_Total_gt_10_v1",
    dataset_path = dataset_path,
    coverage_cutoff = coverage_cutoff,
    alpha = alpha,
    min_reliable = min_reliable,
    class_code_definition = c(TN = 0L, FN = 1L, FP = 2L, TP = 3L),
    class_matrix = selected_class_matrix,
    decisions = decisions
  ),
  file.path(output_dir, "selected_tp_fp_class_matrix.rds"),
  compress = TRUE
)
write.csv(
  threshold_table,
  file.path(output_dir, "global_fdr_for_thresholds.csv"),
  row.names = FALSE
)
write.csv(
  decisions,
  file.path(output_dir, "tp_fp_label_source_decisions.csv"),
  row.names = FALSE
)

me_rip <- decisions[decisions$technology == "MeRIP_seq", , drop = FALSE]
glori <- decisions[decisions$technology == "GLORI", , drop = FALSE]
if (nrow(me_rip) != 1L || nrow(glori) != 1L) {
  stop("Could not identify exactly one MeRIP and one GLORI decision row.")
}
write.csv(
  rbind(me_rip, glori),
  file.path(output_dir, "week7_validation_label_summary.csv"),
  row.names = FALSE
)

message("Saved no-cali TP/FP labels and source decisions to: ", output_dir)
print(rbind(me_rip, glori), row.names = FALSE)
