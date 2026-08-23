parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  values <- list()
  flags <- character()
  for (arg in args) {
    if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      values[[key]] <- sub("^--[^=]+=", "", arg)
    } else if (grepl("^--", arg)) {
      flags <- c(flags, sub("^--", "", arg))
    }
  }
  list(values = values, flags = unique(flags))
}

cli_value <- function(cli, key, default = NULL) {
  value <- cli$values[[key]]
  if (is.null(value)) default else value
}

cli_flag <- function(cli, key) key %in% cli$flags

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

normalized_rank_matrix <- function(score_matrix, universe_size = nrow(score_matrix)) {
  result <- matrix(
    1,
    nrow = nrow(score_matrix),
    ncol = ncol(score_matrix),
    dimnames = dimnames(score_matrix)
  )
  for (column in seq_len(ncol(score_matrix))) {
    observed <- is.finite(score_matrix[, column])
    if (any(observed)) {
      result[observed, column] <- rank(
        -score_matrix[observed, column],
        ties.method = "average"
      ) / universe_size
    }
  }
  result
}

rank_product_score <- function(rank_matrix) {
  -rowMeans(log10(pmax(rank_matrix, .Machine$double.xmin)))
}

rra_pvalue_vectorized <- function(rank_matrix) {
  list_count <- ncol(rank_matrix)
  minimum_beta <- rep(1, nrow(rank_matrix))
  for (order_statistic in seq_len(list_count)) {
    ordered_rank <- matrixStats::rowOrderStats(
      rank_matrix,
      which = order_statistic
    )
    beta_probability <- pbeta(
      ordered_rank,
      shape1 = order_statistic,
      shape2 = list_count - order_statistic + 1L
    )
    minimum_beta <- pmin(minimum_beta, beta_probability)
  }
  pmin(1, minimum_beta * list_count)
}

rra_score <- function(rank_matrix) {
  -log10(pmax(rra_pvalue_vectorized(rank_matrix), .Machine$double.xmin))
}

ranking_metrics <- function(y, score, method, fold, heldout_family = NA_character_, cutoff) {
  ok <- !is.na(y) & is.finite(score)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])
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
    auprc <- PRROC::pr.curve(
      scores.class0 = score[y == 1L],
      scores.class1 = score[y == 0L],
      curve = FALSE
    )$auc.integral
  }
  data.frame(
    fold = fold,
    heldout_family = heldout_family,
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

build_beta_params <- function(dataset, kept_technologies) {
  parameter_table <- S4Vectors::metadata(dataset)[["OmixM6A_para"]]
  result <- list()
  for (technology in kept_technologies) {
    node <- paste0(technology, "_x")
    result[[node]] <- matrix(
      c(
        parameter_table[technology, "alpha_m6A_bg"],
        parameter_table[technology, "alpha_m6A_fg"],
        parameter_table[technology, "beta_m6A_bg"],
        parameter_table[technology, "beta_m6A_fg"]
      ),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(c("alpha", "beta"), c("0", "1"))
    )
  }
  result
}

conquerNet_maximization_safe <- function(cpt_param, belief_list, dag) {
  root_prob <- colMeans(exp(belief_list[["TrueM6A"]]))
  cpt_param[["TrueM6A"]][, 1] <- root_prob / sum(root_prob)
  for (node in setdiff(names(cpt_param), "TrueM6A")) {
    parent <- bnlearn::parents(dag, node)
    if (identical(parent, "TrueM6A")) {
      cpt_param[[node]] <- param_update_strategy(
        belief_list[[node]], belief_list[[parent]]
      )
    } else {
      cpt_param[[node]] <- param_update_technique(
        belief_list[[node]], belief_list[[parent]]
      )
    }
    dimnames(cpt_param[[node]]) <- list(c("0", "1"), c("0", "1"))
  }
  cpt_param
}
