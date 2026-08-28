## Minimal ConquerBN implementation used by the LOTO and LOFO benchmarks.
## This file contains the shared model mathematics; the holdout and evaluation
## steps remain visible in the two benchmark main scripts.

expand_array <- function(x, target, child = TRUE) {
  if (child) {
    return(array(rep(x, dim(target)[3]), dim = dim(target)))
  }
  result <- array(0, dim = dim(target))
  for (parent_state in seq_len(dim(target)[3])) {
    result[, , parent_state] <- x[, parent_state]
  }
  result
}
normalize_columns <- function(x, input_is_log = TRUE) {
  if (!input_is_log) x <- log(x)
  maximum <- apply(x, 2, max)
  scaled <- sweep(x, 2, maximum, FUN = "-")
  probability <- exp(scaled)
  sweep(probability, 2, colSums(probability), FUN = "/")
}
normalize_log_belief <- function(x) {
  two_dimensional <- matrix(
    x,
    nrow = dim(x)[1],
    ncol = prod(dim(x)[-1])
  )
  normalizer <- matrixStats::rowLogSumExps(two_dimensional)
  array(two_dimensional - normalizer, dim = dim(x), dimnames = dimnames(x))
}

apply_dirichlet <- function(x, alpha = 0.5) {
  smoothed <- x + alpha
  sweep(smoothed, 2, colSums(smoothed), FUN = "/")
}

build_conquernet_dag <- function(technologies, strategy_map) {
  meta <- data.frame(
    Technique = technologies,
    Strategy = unname(strategy_map[technologies]),
    stringsAsFactors = FALSE
  )
  strategies <- unique(meta$Strategy)
  technology_nodes <- paste0(meta$Technique, "_k")
  leaf_nodes <- paste0(meta$Technique, "_x")
  arcs <- rbind(
    data.frame(from = "TrueM6A", to = strategies),
    data.frame(from = meta$Strategy, to = technology_nodes),
    data.frame(from = technology_nodes, to = leaf_nodes)
  )
  dag <- bnlearn::empty.graph(unique(c(arcs$from, arcs$to)))
  bnlearn::arcs(dag) <- as.matrix(arcs)
  list(dag = dag, meta = meta, leaf_nodes = leaf_nodes)
}

build_beta_params <- function(dataset, technologies) {
  table <- S4Vectors::metadata(dataset)[["OmixM6A_para"]]
  if (is.null(table)) stop("Dataset metadata does not contain OmixM6A_para.")

  result <- vector("list", length(technologies))
  names(result) <- paste0(technologies, "_x")
  for (technology in technologies) {
    result[[paste0(technology, "_x")]] <- matrix(
      c(
        table[technology, "alpha_m6A_bg"],
        table[technology, "alpha_m6A_fg"],
        table[technology, "beta_m6A_bg"],
        table[technology, "beta_m6A_fg"]
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("alpha", "beta"), c("0", "1"))
    )
  }
  result
}

calculate_emission_loglik <- function(m6a, total, beta_param) {
  result <- vector("list", length(beta_param))
  names(result) <- names(beta_param)
  for (node in names(beta_param)) {
    parameter <- beta_param[[node]]
    result[[node]] <- cbind(
      "0" = extraDistr::dbbinom(
        m6a[, node], total[, node],
        parameter["alpha", "0"], parameter["beta", "0"],
        log = TRUE
      ),
      "1" = extraDistr::dbbinom(
        m6a[, node], total[, node],
        parameter["alpha", "1"], parameter["beta", "1"],
        log = TRUE
      )
    )
  }
  result
}

initialize_cpt <- function(m6a, weight, beta_param, dag) {
  root <- "TrueM6A"
  cpt <- list()
  cpt[[root]] <- matrix(
    c(0.9, 0.1),
    ncol = 1,
    dimnames = list(c("0", "1"), NULL)
  )

  strategy_nodes <- bnlearn::children(dag, root)
  technology_nodes <- vapply(
    bnlearn::leaf.nodes(dag),
    function(node) bnlearn::parents(dag, node)[[1]],
    character(1)
  )

  gamma_technology <- weight
  gamma_technology[is.na(gamma_technology)] <- 0
  gamma_strategy_background <- matrix(
    NA_real_,
    nrow = nrow(m6a),
    ncol = length(strategy_nodes),
    dimnames = list(NULL, strategy_nodes)
  )

  for (strategy in strategy_nodes) {
    cpt[[strategy]] <- matrix(
      c(0.9, 0.1, 0.1, 0.9),
      nrow = 2,
      dimnames = list(c("0", "1"), c("0", "1"))
    )
    children <- bnlearn::children(dag, strategy)
    background <- rep(0, nrow(m6a))
    for (technology_node in children) {
      leaf <- bnlearn::children(dag, technology_node)[[1]]
      background <- background + (1 - gamma_technology[, leaf])
    }
    gamma_strategy_background[, strategy] <- background / length(children)
  }

  rough_strategy_label <- ifelse(gamma_strategy_background <= 0.5, "1", "0")
  for (technology_node in technology_nodes) {
    strategy <- bnlearn::parents(dag, technology_node)[[1]]
    leaf <- bnlearn::children(dag, technology_node)[[1]]
    states <- colnames(beta_param[[leaf]])
    counts <- matrix(
      0,
      nrow = length(states),
      ncol = 2,
      dimnames = list(states, c("0", "1"))
    )

    for (strategy_state in c("0", "1")) {
      rows <- rough_strategy_label[, strategy] == strategy_state
      counts[, strategy_state] <- c(
        sum(1 - gamma_technology[rows, leaf]),
        sum(gamma_technology[rows, leaf])
      )
    }
    cpt[[technology_node]] <- apply_dirichlet(counts)
  }
  cpt
}

param_update_strategy <- function(current, parent) {
  unnormalized <- exp(sweep(
    log(colSums(exp(current))),
    2,
    log(colSums(exp(parent))),
    FUN = "-"
  ))
  normalize_columns(unnormalized, input_is_log = FALSE)
}

param_update_technology <- function(current, parent) {
  unnormalized <- exp(sweep(
    log(colSums(exp(current))),
    2,
    log(rowSums(colSums(exp(parent)))),
    FUN = "-"
  ))
  normalize_columns(unnormalized, input_is_log = FALSE)
}

conquernet_expectation <- function(emission_loglik, cpt, dag) {
  root <- "TrueM6A"
  n <- nrow(emission_loglik[[1]])
  leaf_nodes <- bnlearn::leaf.nodes(dag)
  technology_nodes <- vapply(
    leaf_nodes,
    function(node) bnlearn::parents(dag, node)[[1]],
    character(1)
  )
  strategy_nodes <- bnlearn::children(dag, root)

  theta <- emission_loglik
  for (node in names(cpt)) {
    if (identical(node, root)) {
      theta[[node]] <- matrix(
        rep(log(cpt[[node]]), each = n),
        nrow = n,
        dimnames = list(NULL, rownames(cpt[[node]]))
      )
    } else {
      theta[[node]] <- array(
        rep(log(cpt[[node]]), each = n),
        dim = c(n, dim(cpt[[node]]))
      )
    }
  }

  psi <- list()
  for (node in names(cpt)) {
    child <- bnlearn::children(dag, node)
    if (child[[1]] %in% names(cpt)) {
      psi[[node]] <- theta[[node]]
    } else {
      if (length(child) != 1L) {
        stop("Each technology node must have exactly one leaf child.")
      }
      psi[[node]] <- expand_array(theta[[child]], theta[[node]]) + theta[[node]]
    }
  }

  message <- list()

  ## Collect: technology -> strategy.
  for (node in technology_nodes) {
    key <- paste0(node, "|", bnlearn::parents(dag, node))
    message[[key]] <- apply(
      psi[[node]],
      3,
      matrixStats::rowLogSumExps
    )
  }

  ## Collect: strategy -> root.
  for (node in strategy_nodes) {
    key <- paste0(node, "|", root)
    child_keys <- paste0(bnlearn::children(dag, node), "|", node)
    child_sum <- Reduce("+", message[child_keys])
    temp <- expand_array(child_sum, psi[[node]]) + psi[[node]]
    message[[key]] <- apply(temp, 3, matrixStats::rowLogSumExps)
  }

  ## Distribute: root -> strategy.
  for (node in strategy_nodes) {
    key <- paste0(root, "|", node)
    other_keys <- paste0(setdiff(strategy_nodes, node), "|", root)
    other_sum <- if (length(other_keys)) {
      Reduce("+", message[other_keys])
    } else {
      matrix(0, nrow = n, ncol = 2)
    }
    message[[key]] <- psi[[root]] + other_sum
  }

  ## Distribute: strategy -> technology.
  for (strategy in strategy_nodes) {
    children <- bnlearn::children(dag, strategy)
    for (child in children) {
      key <- paste0(strategy, "|", child)
      up <- message[[paste0(root, "|", strategy)]]
      temp <- expand_array(up, psi[[strategy]], child = FALSE) + psi[[strategy]]
      siblings <- setdiff(children, child)
      if (length(siblings)) {
        down <- Reduce("+", message[paste0(siblings, "|", strategy)])
        temp <- temp + expand_array(down, psi[[strategy]])
      }
      message[[key]] <- apply(temp, 2, matrixStats::rowLogSumExps)
    }
  }

  belief <- list()
  for (node in names(cpt)) {
    child_keys <- paste0(bnlearn::children(dag, node), "|", node)
    if (identical(node, root)) {
      temp <- Reduce("+", message[child_keys]) + psi[[node]]
      loglik <- sum(matrixStats::rowLogSumExps(temp))
      belief[[node]] <- normalize_log_belief(temp)
    } else {
      parent_key <- paste0(bnlearn::parents(dag, node), "|", node)
      temp <- expand_array(message[[parent_key]], psi[[node]], child = FALSE) +
        psi[[node]]
      if (child_keys[[1]] %in% names(message)) {
        child_sum <- Reduce("+", message[child_keys])
        temp <- temp + expand_array(child_sum, psi[[node]])
      }
      belief[[node]] <- normalize_log_belief(temp)
    }
  }
  list(loglik = loglik, belief = belief)
}

conquernet_maximization <- function(cpt, belief, dag) {
  root <- "TrueM6A"
  root_probability <- colMeans(exp(belief[[root]]))
  cpt[[root]][, 1] <- root_probability / sum(root_probability)

  for (node in setdiff(names(cpt), root)) {
    parent <- bnlearn::parents(dag, node)[[1]]
    cpt[[node]] <- if (identical(parent, root)) {
      param_update_strategy(belief[[node]], belief[[parent]])
    } else {
      param_update_technology(belief[[node]], belief[[parent]])
    }
    dimnames(cpt[[node]]) <- list(c("0", "1"), c("0", "1"))
  }
  cpt
}

get_parent_chain <- function(node, dag) {
  path <- character()
  current <- node
  repeat {
    parent <- bnlearn::parents(dag, current)
    if (!length(parent)) break
    if (length(parent) > 1L) stop("Model nodes must have at most one parent.")
    path <- c(parent, path)
    current <- parent
  }
  path
}

switch_cpt_label <- function(cpt, columns = TRUE) {
  original <- cpt
  if (columns) {
    cpt[, 1] <- original[, 2]
    cpt[, 2] <- original[, 1]
  } else {
    cpt[1, ] <- original[2, ]
    cpt[2, ] <- original[1, ]
  }
  cpt
}

cpt_odds_ratio <- function(root, strategy, technology = NULL) {
  root_log <- log(root)
  root_log <- rbind(t(root_log), t(root_log))
  strategy_joint <- exp(log(strategy) + root_log)
  if (is.null(technology)) {
    return(
      (strategy_joint[1, 1] * strategy_joint[2, 2]) /
        (strategy_joint[1, 2] * strategy_joint[2, 1])
    )
  }
  strategy_marginal <- log(rowSums(strategy_joint))
  technology_joint <- exp(log(technology) + strategy_marginal)
  (technology_joint[1, 1] * technology_joint[2, 2]) /
    (technology_joint[1, 2] * technology_joint[2, 1])
}

orient_latent_labels <- function(cpt, path) {
  if (is.null(dim(path)) || nrow(path) < 3L) {
    stop("Path must contain root, strategy and technology rows.")
  }
  result <- cpt
  root <- "TrueM6A"
  strategy_by_technology <- as.character(path[2, ])
  technology_nodes <- as.character(path[3, ])
  strategy_nodes <- unique(strategy_by_technology)

  child_or <- setNames(vapply(seq_along(technology_nodes), function(i) {
    cpt_odds_ratio(
      result[[root]],
      result[[strategy_by_technology[[i]]]],
      result[[technology_nodes[[i]]]]
    )
  }, numeric(1)), technology_nodes)

  for (strategy in strategy_nodes) {
    children <- technology_nodes[strategy_by_technology == strategy]
    agreed <- sum(child_or[children] < 1)
    if (agreed != length(children)) {
      before <- sum(child_or[children])
      switched_strategy <- switch_cpt_label(result[[strategy]], columns = FALSE)
      after <- sum(vapply(children, function(child) {
        cpt_odds_ratio(
          result[[root]],
          switched_strategy,
          switch_cpt_label(result[[child]], columns = TRUE)
        )
      }, numeric(1)))
      if (after <= before) next
    }
    for (child in children) {
      result[[child]] <- switch_cpt_label(result[[child]], columns = TRUE)
    }
    result[[strategy]] <- switch_cpt_label(result[[strategy]], columns = FALSE)
  }

  strategy_or <- setNames(vapply(strategy_nodes, function(strategy) {
    cpt_odds_ratio(result[[root]], result[[strategy]])
  }, numeric(1)), strategy_nodes)

  if (!all(strategy_or >= 1)) {
    if (!all(strategy_or < 1)) {
      switched_root <- switch_cpt_label(result[[root]], columns = FALSE)
      after <- sum(vapply(strategy_nodes, function(strategy) {
        cpt_odds_ratio(
          switched_root,
          switch_cpt_label(result[[strategy]], columns = TRUE)
        )
      }, numeric(1)))
      if (after <= sum(strategy_or)) return(result)
    }
    for (strategy in strategy_nodes) {
      result[[strategy]] <- switch_cpt_label(result[[strategy]], columns = TRUE)
    }
    result[[root]] <- switch_cpt_label(result[[root]], columns = FALSE)
  }
  result
}

extract_root_probability <- function(root_belief) {
  if (ncol(root_belief) != 2L) stop("TrueM6A must have two states.")
  score <- if (!is.null(colnames(root_belief)) && "1" %in% colnames(root_belief)) {
    exp(root_belief[, "1"])
  } else {
    exp(root_belief[, 2L])
  }
  if (!all(is.finite(score)) || any(score < 0 | score > 1)) {
    stop("Invalid TrueM6A posterior probability.")
  }
  score
}

fit_conquernet_holdout <- function(
    dataset,
    kept_technologies,
    cutoff,
    max_iter,
    epsilon,
    seed,
    analysis_version,
    holdout_name) {
  technologies <- colnames(dataset)
  kept_index <- match(kept_technologies, technologies)
  m6a_all <- SummarizedExperiment::assay(dataset, "m6A")
  total_all <- SummarizedExperiment::assay(dataset, "Total")
  probability_all <- SummarizedExperiment::assay(dataset, "m6ASiteProb")

  observed_all <- is.finite(total_all[, kept_index, drop = FALSE]) &
    total_all[, kept_index, drop = FALSE] > cutoff
  row_index <- which(rowSums(observed_all) > 0L)
  m6a <- m6a_all[row_index, kept_index, drop = FALSE]
  total <- total_all[row_index, kept_index, drop = FALSE]
  weight <- probability_all[row_index, kept_index, drop = FALSE]
  observed <- is.finite(total) & total > cutoff

  if (any(observed & (!is.finite(m6a) | !is.finite(weight)))) {
    stop("Observed entries contain missing m6A or m6ASiteProb: ", holdout_name)
  }
  m6a[!observed] <- 0
  total[!observed] <- 0
  weight[!observed] <- NA_real_
  if (any(m6a < 0, na.rm = TRUE) || any(total < 0, na.rm = TRUE) ||
      any(m6a > total, na.rm = TRUE)) {
    stop("Invalid counts after coverage filtering: ", holdout_name)
  }

  model <- build_conquernet_dag(kept_technologies, conquernet_strategy_map())
  colnames(m6a) <- model$leaf_nodes
  colnames(total) <- model$leaf_nodes
  colnames(weight) <- model$leaf_nodes

  beta <- build_beta_params(dataset, kept_technologies)
  emission <- calculate_emission_loglik(m6a, total, beta)
  cpt <- initialize_cpt(m6a, weight, beta, model$dag)

  set.seed(seed)
  log_likelihood <- numeric(max_iter)
  previous <- NA_real_
  converged <- FALSE
  started <- Sys.time()

  for (iteration in seq_len(max_iter)) {
    expectation <- conquernet_expectation(emission, cpt, model$dag)
    current <- expectation$loglik
    if (!is.finite(current)) {
      stop("Non-finite log-likelihood: ", holdout_name, ", iteration ", iteration)
    }
    log_likelihood[[iteration]] <- current
    cpt <- conquernet_maximization(cpt, expectation$belief, model$dag)

    relative_change <- if (iteration > 1L) {
      abs(current - previous) / max(abs(previous), .Machine$double.eps)
    } else {
      Inf
    }
    if (iteration == 1L || iteration %% 10L == 0L || relative_change < epsilon) {
      message(sprintf(
        "%s: iteration=%d loglik=%.6f relative_change=%.3e",
        holdout_name, iteration, current, relative_change
      ))
    }
    if (iteration > 1L && relative_change < epsilon) {
      converged <- TRUE
      break
    }
    previous <- current
  }

  iterations <- iteration
  log_likelihood <- log_likelihood[seq_len(iterations)]
  if (!converged) stop("EM did not converge: ", holdout_name)

  path <- sapply(
    bnlearn::leaf.nodes(model$dag),
    get_parent_chain,
    dag = model$dag
  )
  cpt <- orient_latent_labels(cpt, path)
  final <- conquernet_expectation(emission, cpt, model$dag)
  score <- extract_root_probability(final$belief[["TrueM6A"]])

  list(
    analysis_version = analysis_version,
    holdout_name = holdout_name,
    kept_technologies = kept_technologies,
    coverage_rule = paste0("Total>", cutoff),
    coverage_cutoff = cutoff,
    row_index = row_index,
    score = score,
    cpt_param = cpt,
    beta_param = beta,
    dag = model$dag,
    log_likelihood = log_likelihood,
    converged = converged,
    iterations = iterations,
    epsilon = epsilon,
    max_iter = max_iter,
    n_total_sites = nrow(dataset),
    n_training_sites = length(row_index),
    n_kept_entries_total = length(observed),
    n_kept_entries_passing_coverage = sum(observed),
    runtime_minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))
  )
}
