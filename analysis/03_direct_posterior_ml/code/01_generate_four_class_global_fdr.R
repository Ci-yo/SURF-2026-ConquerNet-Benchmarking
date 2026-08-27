options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(bnlearn)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Please run this file with Rscript.")
script_dir <- dirname(normalizePath(
  sub("^--file=", "", script_arg[[1L]]), winslash = "/"
))
model_path <- if (length(args) >= 1L) args[[1L]] else {
  paste0(
    "/home/xiyao/15data/W7/06_no_cali_IVT_matched_full_model/",
    "Global_no_additional_coverage/full_no_cali_global.conquernet.rds"
  )
}
output_dir <- if (length(args) >= 2L) args[[2L]] else {
  "/home/xiyao/15data/W8/results/01_four_class_global_fdr"
}
alpha <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 0.05
downstream_coverage_cutoff <- if (length(args) >= 4L) {
  as.integer(args[[4L]])
} else {
  10L
}
smoke_n <- if (length(args) >= 5L) as.integer(args[[5L]]) else 0L

helper_path <- if (length(args) >= 6L) args[[6L]] else {
  normalizePath(
    file.path(script_dir, "EM_function_2.R"),
    winslash = "/", mustWork = FALSE
  )
}
if (!file.exists(model_path)) stop("Missing model: ", model_path)
if (!file.exists(helper_path)) stop("Missing inference helper: ", helper_path)
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("Invalid alpha.")
if (!is.finite(downstream_coverage_cutoff) || downstream_coverage_cutoff < 0) {
  stop("Invalid downstream coverage cutoff.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
posterior_dir <- file.path(output_dir, "per_technology")
dir.create(posterior_dir, recursive = TRUE, showWarnings = FALSE)

source(helper_path)

started <- Sys.time()
model <- readRDS(model_path)
if (!isTRUE(model$converged)) stop("Saved ConquerBN model did not converge.")
if (is.null(model$dataset_path) || !file.exists(model$dataset_path)) {
  stop("The dataset recorded in the model is missing: ", model$dataset_path)
}
data <- readRDS(model$dataset_path)
if (nrow(data) != model$n_total_sites) stop("Model/data universe mismatch.")

technologies <- as.character(model$technologies)
if (!all(technologies %in% colnames(data))) {
  stop("At least one model technology is absent from the dataset.")
}
row_index <- as.integer(model$row_index)
if (smoke_n > 0L) row_index <- head(row_index, smoke_n)
if (!length(row_index) || anyNA(row_index) || anyDuplicated(row_index) ||
    min(row_index) < 1L || max(row_index) > nrow(data)) {
  stop("Invalid model row_index.")
}

dag <- model$dag
cpt_param <- model$cpt_param
beta_param <- model$beta_param
continuous.nodes <- bnlearn::leaf.nodes(dag)
assign("continuous.nodes", continuous.nodes, envir = .GlobalEnv)
expected_leaves <- paste0(technologies, "_x")
if (!setequal(expected_leaves, continuous.nodes)) {
  stop("Saved technologies do not match DAG leaves.")
}

m6a <- assays(data)$m6A[row_index, technologies, drop = FALSE]
total_raw <- assays(data)$Total[row_index, technologies, drop = FALSE]
total <- total_raw
if (isTRUE(model$coverage_applied_to_training_entries)) {
  model_missing <- !is.finite(total) | total <= model$coverage_cutoff
} else {
  model_missing <- !is.finite(total) | total <= 0
}
m6a[model_missing] <- 0
total[model_missing] <- 0
colnames(m6a) <- colnames(total) <- expected_leaves
m6a <- m6a[, continuous.nodes, drop = FALSE]
total <- total[, continuous.nodes, drop = FALSE]

message("Calculating emission likelihood for ", length(row_index), " sites...")
emission <- calculate_emission_loglik(m6a, total, beta_param)
message("Running final sum-product inference...")
inference <- conquerNet_inference(emission, cpt_param, dag, "sum-product")
rm(emission, m6a, total)
gc(verbose = FALSE)

root_node <- bnlearn::root.nodes(dag)
if (length(root_node) != 1L) stop("Expected one root node.")
root_state1 <- exp(inference$belief[[root_node]][, 2L])

calculate_four_class <- function(inference, technology_node, dag) {
  strategy <- bnlearn::parents(dag, technology_node)
  if (length(strategy) != 1L) {
    stop("Expected exactly one strategy parent for ", technology_node)
  }
  root_strategy <- inference$belief[[strategy]]
  technology_strategy <- inference$belief[[technology_node]]
  if (!identical(dim(root_strategy)[2:3], c(2L, 2L)) ||
      !identical(dim(technology_strategy)[2:3], c(2L, 2L))) {
    stop("Unexpected pairwise-belief dimensions for ", technology_node)
  }

  log_strategy0 <- matrixStats::rowLogSumExps(cbind(
    root_strategy[, 1L, 1L], root_strategy[, 1L, 2L]
  ))
  log_strategy1 <- matrixStats::rowLogSumExps(cbind(
    root_strategy[, 2L, 1L], root_strategy[, 2L, 2L]
  ))
  joint <- function(root_index, technology_index) {
    exp(matrixStats::rowLogSumExps(cbind(
      root_strategy[, 1L, root_index] +
        technology_strategy[, technology_index, 1L] - log_strategy0,
      root_strategy[, 2L, root_index] +
        technology_strategy[, technology_index, 2L] - log_strategy1
    )))
  }
  cbind(
    TN = joint(1L, 1L),
    FN = joint(2L, 1L),
    FP = joint(1L, 2L),
    TP = joint(2L, 2L)
  )
}

global_call <- function(posterior, alpha) {
  if (any(!is.finite(posterior)) || any(posterior < 0 | posterior > 1)) {
    stop("Invalid class posterior values.")
  }
  local_error <- 1 - posterior
  ordering <- base::order(local_error, seq_along(local_error))
  global_error_ordered <- cumsum(local_error[ordering]) / seq_along(ordering)
  selected_ordered <- global_error_ordered < alpha
  selected <- rep(FALSE, length(posterior))
  selected[ordering[selected_ordered]] <- TRUE
  last <- if (any(selected_ordered)) max(which(selected_ordered)) else NA_integer_
  list(
    selected = selected,
    n_selected = sum(selected),
    boundary_global_fdr = if (is.na(last)) NA_real_ else {
      global_error_ordered[[last]]
    },
    minimum_selected_posterior = if (!any(selected)) NA_real_ else {
      min(posterior[selected])
    },
    maximum_unselected_posterior = if (all(selected)) NA_real_ else {
      max(posterior[!selected])
    }
  )
}

class_names <- c("TN", "FN", "FP", "TP")
class_codes <- setNames(0:3, class_names)
n_total <- nrow(data)
class_matrix_global <- matrix(
  NA_integer_, nrow = n_total, ncol = length(technologies),
  dimnames = list(NULL, technologies)
)
class_matrix_total_gt_cutoff <- class_matrix_global
count_rows <- list()
fdr_rows <- list()
identity_rows <- list()

for (technology in technologies) {
  message("Four-class posterior: ", technology)
  technology_node <- paste0(technology, "_k")
  posterior <- calculate_four_class(inference, technology_node, dag)

  technology_state1 <- exp(matrixStats::rowLogSumExps(cbind(
    inference$belief[[technology_node]][, 2L, 1L],
    inference$belief[[technology_node]][, 2L, 2L]
  )))
  identity <- data.frame(
    technology = technology,
    n_inference_sites = nrow(posterior),
    min_posterior = min(posterior),
    max_posterior = max(posterior),
    max_abs_sum_to_one_error = max(abs(rowSums(posterior) - 1)),
    max_abs_root_state1_error = max(abs(
      posterior[, "FN"] + posterior[, "TP"] - root_state1
    )),
    max_abs_technology_state1_error = max(abs(
      posterior[, "FP"] + posterior[, "TP"] - technology_state1
    )),
    stringsAsFactors = FALSE
  )
  identity_rows[[technology]] <- identity
  if (max(unlist(identity[, 5:7])) > 1e-10) {
    stop("Four-class posterior identity failed for ", technology)
  }

  calls <- lapply(class_names, function(class_name) {
    global_call(posterior[, class_name], alpha)
  })
  names(calls) <- class_names
  selected <- do.call(cbind, lapply(calls, `[[`, "selected"))
  colnames(selected) <- class_names
  n_selected_classes <- rowSums(selected)
  label <- rep(NA_integer_, nrow(selected))
  unique_call <- n_selected_classes == 1L
  if (any(unique_call)) {
    selected_name <- class_names[max.col(selected[unique_call, , drop = FALSE])]
    label[unique_call] <- unname(class_codes[selected_name])
  }

  observed_gt_cutoff <- is.finite(total_raw[, technology]) &
    total_raw[, technology] > downstream_coverage_cutoff
  label_gt_cutoff <- label
  label_gt_cutoff[!observed_gt_cutoff] <- NA_integer_
  class_matrix_global[row_index, technology] <- label
  class_matrix_total_gt_cutoff[row_index, technology] <- label_gt_cutoff

  for (class_name in class_names) {
    fdr_rows[[paste(technology, class_name, sep = "::")]] <- data.frame(
      technology = technology,
      class = class_name,
      alpha = alpha,
      n_selected = calls[[class_name]]$n_selected,
      boundary_global_fdr = calls[[class_name]]$boundary_global_fdr,
      minimum_selected_posterior =
        calls[[class_name]]$minimum_selected_posterior,
      maximum_unselected_posterior =
        calls[[class_name]]$maximum_unselected_posterior,
      stringsAsFactors = FALSE
    )
  }
  count_rows[[technology]] <- data.frame(
    technology = technology,
    n_total_universe = n_total,
    n_inference_sites = length(label),
    n_observed_total_gt_cutoff = sum(observed_gt_cutoff),
    coverage_cutoff = downstream_coverage_cutoff,
    n_unclassified_global = sum(is.na(label)),
    n_overlap_global = sum(n_selected_classes > 1L),
    n_TN_global = sum(label == class_codes[["TN"]], na.rm = TRUE),
    n_FN_global = sum(label == class_codes[["FN"]], na.rm = TRUE),
    n_FP_global = sum(label == class_codes[["FP"]], na.rm = TRUE),
    n_TP_global = sum(label == class_codes[["TP"]], na.rm = TRUE),
    n_unclassified_Total_gt_cutoff = sum(is.na(label_gt_cutoff)),
    n_TN_Total_gt_cutoff = sum(
      label_gt_cutoff == class_codes[["TN"]], na.rm = TRUE
    ),
    n_FN_Total_gt_cutoff = sum(
      label_gt_cutoff == class_codes[["FN"]], na.rm = TRUE
    ),
    n_FP_Total_gt_cutoff = sum(
      label_gt_cutoff == class_codes[["FP"]], na.rm = TRUE
    ),
    n_TP_Total_gt_cutoff = sum(
      label_gt_cutoff == class_codes[["TP"]], na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )

  saveRDS(
    list(
      analysis_version = "week8_direct_four_class_posterior_global_fdr_v1",
      technology = technology,
      row_index = row_index,
      posterior = posterior,
      selected_global_fdr = selected,
      n_selected_classes = n_selected_classes,
      label_global = label,
      observed_total_gt_cutoff = observed_gt_cutoff,
      label_total_gt_cutoff = label_gt_cutoff,
      class_code = class_codes,
      alpha = alpha,
      downstream_coverage_cutoff = downstream_coverage_cutoff,
      model_path = model_path,
      dataset_path = model$dataset_path
    ),
    file.path(posterior_dir, paste0(technology, ".four_class.rds")),
    compress = "gzip"
  )
  rm(
    posterior, technology_state1, calls, selected, n_selected_classes,
    label, observed_gt_cutoff, label_gt_cutoff
  )
  gc(verbose = FALSE)
}

counts <- do.call(rbind, count_rows)
fdr_audit <- do.call(rbind, fdr_rows)
identity_audit <- do.call(rbind, identity_rows)

saveRDS(
  class_matrix_global,
  file.path(output_dir, "four_class_matrix_global_unique.rds"),
  compress = "gzip"
)
saveRDS(
  class_matrix_total_gt_cutoff,
  file.path(
    output_dir,
    paste0("four_class_matrix_Total_gt_", downstream_coverage_cutoff, ".rds")
  ),
  compress = "gzip"
)
write.csv(counts, file.path(output_dir, "four_class_counts.csv"), row.names = FALSE)
write.csv(
  fdr_audit,
  file.path(output_dir, "global_fdr_call_audit.csv"),
  row.names = FALSE
)
write.csv(
  identity_audit,
  file.path(output_dir, "posterior_identity_audit.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    code = unname(class_codes),
    class = names(class_codes),
    TrueM6A_state = c(0L, 1L, 0L, 1L),
    technology_state = c(0L, 0L, 1L, 1L),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "class_codebook.csv"),
  row.names = FALSE
)
finished <- Sys.time()
write.csv(
  data.frame(
    analysis_version = "week8_direct_four_class_posterior_global_fdr_v1",
    model_path = model_path,
    dataset_path = model$dataset_path,
    helper_path = helper_path,
    n_total_sites = n_total,
    n_inference_sites = length(row_index),
    smoke_n = smoke_n,
    posterior_definition =
      "direct_joint_P(TrueM6A_state,technology_state_given_all_evidence)",
    global_call_definition = "local_error=1-posterior;cumulative_mean<alpha",
    alpha = alpha,
    downstream_coverage_rule = paste0("technology Total>", downstream_coverage_cutoff),
    overlap_resolution = "only_exactly_one_selected_class_receives_label",
    fallback_rule = "none",
    started = format(started, "%Y-%m-%d %H:%M:%S %Z"),
    finished = format(finished, "%Y-%m-%d %H:%M:%S %Z"),
    runtime_minutes = as.numeric(difftime(finished, started, units = "mins")),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "protocol.csv"),
  row.names = FALSE
)

cat("Saved Week8 four-class posterior results to: ", output_dir, "\n", sep = "")
