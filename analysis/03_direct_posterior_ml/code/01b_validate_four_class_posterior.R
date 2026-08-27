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
n_smoke <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5000L
helper_path <- if (length(args) >= 3L) args[[3L]] else {
  normalizePath(
    file.path(script_dir, "EM_function_2.R"),
    winslash = "/", mustWork = FALSE
  )
}

if (!file.exists(helper_path)) stop("Missing inference helper: ", helper_path)
source(helper_path)
model <- readRDS(model_path)
data <- readRDS(model$dataset_path)
dag <- model$dag
cpt_param <- model$cpt_param
beta_param <- model$beta_param
continuous.nodes <- bnlearn::leaf.nodes(dag)
assign("continuous.nodes", continuous.nodes, envir = .GlobalEnv)

row_index <- head(model$row_index, n_smoke)
m6a <- assays(data)$m6A[row_index, model$technologies, drop = FALSE]
total <- assays(data)$Total[row_index, model$technologies, drop = FALSE]

if (isTRUE(model$coverage_applied_to_training_entries)) {
  missing_entry <- !is.finite(total) | total <= model$coverage_cutoff
} else {
  missing_entry <- !is.finite(total) | total <= 0
}
m6a[missing_entry] <- 0
total[missing_entry] <- 0

expected_leaves <- paste0(model$technologies, "_x")
if (!setequal(expected_leaves, continuous.nodes)) {
  stop("Technology leaves do not match the saved model.")
}
colnames(m6a) <- colnames(total) <- expected_leaves
m6a <- m6a[, continuous.nodes, drop = FALSE]
total <- total[, continuous.nodes, drop = FALSE]

emission <- calculate_emission_loglik(m6a, total, beta_param)
res <- conquerNet_inference(emission, cpt_param, dag, "sum-product")

calculate_eval_prob_fixed <- function(res, tech, dag) {
  if (!tech %in% unlist(lapply(bnlearn::leaf.nodes(dag), function(x) {
    bnlearn::parents(dag, x)
  }))) {
    stop("Technology latent node is absent from the network: ", tech)
  }
  strategy <- bnlearn::parents(dag, tech)
  if (length(strategy) != 1L) stop("Expected exactly one strategy parent.")

  root_strategy <- res$belief[[strategy]]
  tech_strategy <- res$belief[[tech]]
  if (length(dim(root_strategy)) != 3L || length(dim(tech_strategy)) != 3L) {
    stop("Unexpected belief dimensions.")
  }

  log_p_strategy0 <- matrixStats::rowLogSumExps(
    cbind(root_strategy[, 1L, 1L], root_strategy[, 1L, 2L])
  )
  log_p_strategy1 <- matrixStats::rowLogSumExps(
    cbind(root_strategy[, 2L, 1L], root_strategy[, 2L, 2L])
  )

  joint <- function(root_state, tech_state) {
    root_index <- as.integer(root_state) + 1L
    tech_index <- as.integer(tech_state) + 1L
    matrixStats::rowLogSumExps(cbind(
      root_strategy[, 1L, root_index] +
        tech_strategy[, tech_index, 1L] - log_p_strategy0,
      root_strategy[, 2L, root_index] +
        tech_strategy[, tech_index, 2L] - log_p_strategy1
    ))
  }

  out <- cbind(
    TN = exp(joint("0", "0")),
    FN = exp(joint("1", "0")),
    FP = exp(joint("0", "1")),
    TP = exp(joint("1", "1"))
  )
  out
}

root_probability <- exp(res$belief[[bnlearn::root.nodes(dag)]][, 2L])
technology_nodes <- paste0(model$technologies, "_k")

audit <- do.call(rbind, lapply(technology_nodes, function(tech) {
  p <- calculate_eval_prob_fixed(res, tech, dag)
  tech_probability <- exp(matrixStats::rowLogSumExps(cbind(
    res$belief[[tech]][, 2L, 1L],
    res$belief[[tech]][, 2L, 2L]
  )))
  data.frame(
    technology = sub("_k$", "", tech),
    n_sites = nrow(p),
    min_probability = min(p),
    max_probability = max(p),
    max_abs_four_class_sum_error = max(abs(rowSums(p) - 1)),
    max_abs_root_state1_marginal_error = max(abs(
      p[, "FN"] + p[, "TP"] - root_probability
    )),
    max_abs_technology_state1_marginal_error = max(abs(
      p[, "FP"] + p[, "TP"] - tech_probability
    )),
    stringsAsFactors = FALSE
  )
}))

print(audit, row.names = FALSE)
if (any(!is.finite(as.matrix(audit[, -(1:2)])))) {
  stop("Non-finite audit result.")
}
if (max(audit$max_abs_four_class_sum_error) > 1e-10 ||
    max(audit$max_abs_root_state1_marginal_error) > 1e-10 ||
    max(audit$max_abs_technology_state1_marginal_error) > 1e-10) {
  stop("Four-class posterior identities failed.")
}
cat("\nFOUR_CLASS_POSTERIOR_IDENTITIES_PASSED\n")
