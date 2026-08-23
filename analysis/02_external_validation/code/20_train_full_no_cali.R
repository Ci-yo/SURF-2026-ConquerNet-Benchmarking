suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(bnlearn)
  library(matrixStats)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Cannot determine current script path.")
script_dir <- dirname(normalizePath(
  sub("^--file=", "", script_arg[[1]]),
  winslash = "/"
))

source(file.path(script_dir, "w6_helpers.R"))
source(file.path(script_dir, "latent_switch_fixed.R"))
cli <- parse_cli()

function_file <- normalizePath(
  cli_value(
    cli,
    "function-file",
    "/home/xiyao/15data/ConquerNet/ConquerBN_basic/function.R"
  ),
  winslash = "/",
  mustWork = TRUE
)
source(function_file)
inference_file <- normalizePath(
  cli_value(
    cli,
    "inference-file",
    "/home/xiyao/15data/ConquerNet/ConquerBN_basic/EM_function_2.R"
  ),
  winslash = "/",
  mustWork = TRUE
)

dataset_path <- normalizePath(
  cli_value(
    cli,
    "dataset",
    paste0(
      "/home/xichen/ConquerNet/ConquerBN_basic/",
      "m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali.rds"
    )
  ),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- cli_value(
  cli,
  "output-dir",
  "/home/xiyao/15data/W7/03_no_cali_full_model/Total_gt_10"
)
coverage_cutoff <- as.integer(cli_value(cli, "cutoff", "10"))
max_iter <- as.integer(cli_value(cli, "max-iter", "500"))
epsilon <- as.numeric(cli_value(cli, "epsilon", "1e-6"))
seed <- as.integer(cli_value(cli, "seed", "20260730"))
smoke_n <- as.integer(cli_value(cli, "smoke-n", "0"))
overwrite <- cli_flag(cli, "overwrite")

if (!identical(coverage_cutoff, 10L)) {
  stop("Week7 formal analysis is fixed to --cutoff=10.")
}
if (!is.finite(max_iter) || max_iter < 2L) stop("Invalid --max-iter.")
if (!is.finite(epsilon) || epsilon <= 0) stop("Invalid --epsilon.")
if (!is.finite(smoke_n) || smoke_n < 0L) stop("Invalid --smoke-n.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
model_path <- file.path(output_dir, "full_no_cali_Total_gt_10.conquernet.rds")
posterior_path <- file.path(output_dir, "node_state1_posterior.rds")
viterbi_class_path <- file.path(output_dir, "viterbi_class_matrix.rds")
training_qc_path <- file.path(output_dir, "training_qc.csv")

if (
  all(file.exists(c(model_path, posterior_path, viterbi_class_path, training_qc_path))) &&
    !overwrite
) {
  message("All full-model outputs already exist. Use --overwrite to replace them.")
  quit(status = 0)
}

message("Reading no-calibration benchmark: ", dataset_path)
dataset <- readRDS(dataset_path)
if (!inherits(dataset, "SummarizedExperiment")) {
  stop("Dataset is not a SummarizedExperiment.")
}

required_assays <- c("m6A", "Total", "m6ASiteProb")
missing_assays <- setdiff(required_assays, assayNames(dataset))
if (length(missing_assays)) {
  stop("Missing assays: ", paste(missing_assays, collapse = ", "))
}

m6a_all <- assay(dataset, "m6A")
total_all <- assay(dataset, "Total")
weight_all <- assay(dataset, "m6ASiteProb")
technologies <- colnames(dataset)
n_sites <- nrow(dataset)
strategy_map <- conquernet_strategy_map()

stopifnot(
  identical(dim(m6a_all), dim(total_all)),
  identical(dim(weight_all), dim(total_all)),
  identical(dimnames(m6a_all), dimnames(total_all)),
  identical(dimnames(weight_all), dimnames(total_all)),
  setequal(technologies, names(strategy_map))
)
strategy_map <- strategy_map[technologies]

observed_all <- is.finite(total_all) & total_all > coverage_cutoff
bad_observed <- observed_all & (!is.finite(m6a_all) | !is.finite(weight_all))
if (any(bad_observed)) {
  stop("Total>10 entries contain missing m6A or m6ASiteProb.")
}

training_available <- rowSums(observed_all) > 0L
row_index <- which(training_available)
if (smoke_n > 0L && length(row_index) > smoke_n) {
  set.seed(seed)
  row_index <- sort(sample(row_index, smoke_n))
}
message(
  "Training full no-cali model on ", length(row_index), " / ", n_sites,
  " sites with at least one Total>10 technology."
)

m6a <- m6a_all[row_index, , drop = FALSE]
total <- total_all[row_index, , drop = FALSE]
weight <- weight_all[row_index, , drop = FALSE]
observed <- observed_all[row_index, , drop = FALSE]

m6a[!observed] <- 0
total[!observed] <- 0
weight[!observed] <- NA_real_

if (
  any(m6a < 0, na.rm = TRUE) ||
    any(total < 0, na.rm = TRUE) ||
    any(m6a > total, na.rm = TRUE)
) {
  stop("Invalid m6A/Total counts after coverage filtering.")
}

meta <- data.frame(
  Technique = technologies,
  Strategy = unname(strategy_map),
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

beta_param <- build_beta_params(dataset, technologies)
message("Computing beta-binomial emission log-likelihood.")
emission_loglik <- calculate_emission_loglik(m6a, total, beta_param)
message("Initializing full-model CPTs.")
cpt_param <- initialization_cpt_params(
  m6a,
  total,
  beta_param,
  weight,
  dag
)

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
    stop("Non-finite log-likelihood at iteration ", iteration)
  }

  log_likelihood[[iteration]] <- current
  cpt_param <- conquerNet_maximization_safe(cpt_param, belief, dag)
  relative_change <- if (iteration > 1L) {
    abs(current - previous) / max(abs(previous), .Machine$double.eps)
  } else {
    NA_real_
  }

  message(sprintf(
    "iteration=%d loglik=%.6f relative_change=%s elapsed_min=%.3f",
    iteration,
    current,
    ifelse(
      is.na(relative_change),
      "NA",
      format(relative_change, scientific = TRUE, digits = 4)
    ),
    as.numeric(difftime(Sys.time(), iteration_started, units = "mins"))
  ))

  rm(expectation, belief)
  invisible(gc())

  if (iteration > 1L && relative_change < epsilon) {
    converged <- TRUE
    break
  }
  previous <- current
}

iterations <- max(which(!is.na(log_likelihood)))
log_likelihood <- log_likelihood[seq_len(iterations)]
if (!converged) {
  stop("Full no-calibration model did not converge.")
}

path <- sapply(
  bnlearn::leaf.nodes(dag),
  function(node) get_parent_chain_bn(node, dag)
)
cpt_param <- latent_param_label_switch_fixed(cpt_param, path)

## The verified W6 training functions above do not expose the generic
## sum-product/max-sum wrapper. Load the original ConquerNet inference
## implementation only after EM convergence.
source(inference_file)

extract_state1_probability <- function(node_name, belief_node, expected_n) {
  if (identical(node_name, "TrueM6A")) {
    marginal_log <- belief_node
  } else {
    marginal_log <- apply(belief_node, 2L, matrixStats::rowLogSumExps)
  }
  if (nrow(marginal_log) != expected_n && ncol(marginal_log) == expected_n) {
    marginal_log <- t(marginal_log)
  }
  if (nrow(marginal_log) != expected_n || ncol(marginal_log) != 2L) {
    stop("Unexpected marginal dimensions for ", node_name)
  }
  marginal_log <- marginal_log - matrixStats::rowLogSumExps(marginal_log)
  probability <- exp(marginal_log)
  state_names <- colnames(marginal_log)
  state1_column <- if (!is.null(state_names) && "1" %in% state_names) {
    match("1", state_names)
  } else {
    2L
  }
  result <- probability[, state1_column]
  if (
    any(!is.finite(result)) ||
      any(result < -1e-10 | result > 1 + 1e-10)
  ) {
    stop("Invalid posterior probability for ", node_name)
  }
  pmin(1, pmax(0, result))
}

message("Running final sum-product inference.")
sum_product <- conquerNet_inference(
  emission_loglik,
  cpt_param,
  dag,
  mode = "sum-product"
)

posterior_nodes <- c("TrueM6A", technology_nodes)
posterior_training <- matrix(
  NA_real_,
  nrow = length(row_index),
  ncol = length(posterior_nodes),
  dimnames = list(NULL, posterior_nodes)
)
for (node in posterior_nodes) {
  posterior_training[, node] <- extract_state1_probability(
    node,
    sum_product$belief[[node]],
    length(row_index)
  )
}
inference_loglik <- sum_product$loglik
rm(sum_product)
invisible(gc())

posterior_full <- matrix(
  NA_real_,
  nrow = n_sites,
  ncol = length(posterior_nodes),
  dimnames = list(NULL, posterior_nodes)
)
posterior_full[row_index, ] <- posterior_training

message("Running final max-sum inference for Viterbi fallback labels.")
max_sum <- conquerNet_inference(
  emission_loglik,
  cpt_param,
  dag,
  mode = "max-sum"
)
optimal_states <- max_sum$optimal_states
if (is.null(names(optimal_states))) {
  stop("Max-sum optimal states do not have node names.")
}
missing_state_nodes <- setdiff(posterior_nodes, names(optimal_states))
if (length(missing_state_nodes)) {
  stop(
    "Max-sum result is missing nodes: ",
    paste(missing_state_nodes, collapse = ", ")
  )
}

root_state <- as.integer(optimal_states[["TrueM6A"]])
if (length(root_state) != length(row_index) || any(!root_state %in% 0:1)) {
  stop("Invalid TrueM6A Viterbi state.")
}

viterbi_class_training <- matrix(
  NA_integer_,
  nrow = length(row_index),
  ncol = length(technologies),
  dimnames = list(NULL, technologies)
)
for (j in seq_along(technologies)) {
  technology_state <- as.integer(optimal_states[[technology_nodes[[j]]]])
  if (
    length(technology_state) != length(row_index) ||
      any(!technology_state %in% 0:1)
  ) {
    stop("Invalid Viterbi state for ", technologies[[j]])
  }
  viterbi_class_training[, j] <- 2L * technology_state + root_state
}
rm(max_sum, optimal_states, root_state)
invisible(gc())

viterbi_class_full <- matrix(
  NA_integer_,
  nrow = n_sites,
  ncol = length(technologies),
  dimnames = list(NULL, technologies)
)
viterbi_class_full[row_index, ] <- viterbi_class_training

finished <- Sys.time()
model <- list(
  analysis_version = "week7_no_cali_full_Total_gt_10_v1",
  dataset_path = dataset_path,
  coverage_rule = "Total>10",
  coverage_cutoff = coverage_cutoff,
  coverage_applied_to_training_entries = TRUE,
  below_threshold_representation = "m6A=0; Total=0; m6ASiteProb=NA",
  technologies = technologies,
  strategy_map = strategy_map,
  row_index = row_index,
  cpt_param = cpt_param,
  beta_param = beta_param,
  dag = dag,
  log_likelihood = log_likelihood,
  converged = converged,
  iterations = iterations,
  epsilon = epsilon,
  max_iter = max_iter,
  seed = seed,
  smoke_n = smoke_n,
  n_total_sites = n_sites,
  n_training_sites = length(row_index),
  n_entries_passing_coverage = sum(observed_all),
  inference_loglik = inference_loglik,
  orientation_method = "latent_param_label_switch_fixed_v1",
  started = started,
  finished = finished,
  runtime_minutes = as.numeric(difftime(finished, started, units = "mins")),
  session_info = sessionInfo()
)

saveRDS(model, model_path, compress = TRUE)
saveRDS(
  list(
    analysis_version = model$analysis_version,
    dataset_path = dataset_path,
    coverage_cutoff = coverage_cutoff,
    row_index = row_index,
    node_state1_probability = posterior_full,
    inference_loglik = inference_loglik
  ),
  posterior_path,
  compress = TRUE
)
saveRDS(
  list(
    analysis_version = model$analysis_version,
    dataset_path = dataset_path,
    coverage_cutoff = coverage_cutoff,
    row_index = row_index,
    class_code_definition = c(TN = 0L, FN = 1L, FP = 2L, TP = 3L),
    class_matrix = viterbi_class_full
  ),
  viterbi_class_path,
  compress = TRUE
)

qc <- data.frame(
  analysis_version = model$analysis_version,
  dataset_path = dataset_path,
  coverage_rule = model$coverage_rule,
  n_total_sites = n_sites,
  n_training_sites = length(row_index),
  n_technologies = length(technologies),
  n_entries_passing_coverage = sum(observed_all),
  converged = converged,
  iterations = iterations,
  final_log_likelihood = tail(log_likelihood, 1L),
  inference_loglik = inference_loglik,
  runtime_minutes = model$runtime_minutes,
  posterior_min = min(posterior_training),
  posterior_max = max(posterior_training),
  orientation_method = model$orientation_method,
  stringsAsFactors = FALSE
)
write.csv(qc, training_qc_path, row.names = FALSE)

message("Saved full model: ", model_path)
message("Saved posterior: ", posterior_path)
message("Saved Viterbi classes: ", viterbi_class_path)
message("Saved QC: ", training_qc_path)
