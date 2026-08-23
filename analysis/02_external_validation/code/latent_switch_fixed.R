## Corrected implementation of the original latent_param_label_switch().
## It preserves the original odds-ratio/voting rule and only fixes indexing,
## stale-variable reuse, and the invalid break outside a loop.

latent_param_label_switch_fixed <- function(cpt.param, path) {
  if (is.null(dim(path)) || nrow(path) < 3L) {
    stop("path must contain root, strategy, and technology rows.")
  }
  result <- cpt.param
  root_node <- "TrueM6A"
  strategy_by_technology <- as.character(path[2, ])
  technology_nodes <- as.character(path[3, ])
  strategy_nodes <- unique(strategy_by_technology)
  required <- unique(c(root_node, strategy_nodes, technology_nodes))
  missing <- setdiff(required, names(result))
  if (length(missing)) stop("Missing CPT nodes: ", paste(missing, collapse = ", "))
  if (anyDuplicated(technology_nodes)) stop("Technology nodes in path must be unique.")

  child_or <- setNames(vapply(seq_along(technology_nodes), function(i) {
    odds_ratio(
      result[[root_node]],
      result[[strategy_by_technology[i]]],
      result[[technology_nodes[i]]]
    )
  }, numeric(1)), technology_nodes)
  vote_total <- table(factor(strategy_by_technology, levels = strategy_nodes))
  vote_agreed <- table(factor(strategy_by_technology[child_or < 1], levels = strategy_nodes))

  for (strategy in strategy_nodes) {
    children <- technology_nodes[strategy_by_technology == strategy]
    total <- as.integer(vote_total[strategy])
    agreed <- as.integer(vote_agreed[strategy])
    if (total == 0L) next
    if (total != agreed) {
      before <- sum(child_or[children])
      strategy_switched <- label_switch(result[[strategy]], column = FALSE)
      after <- sum(vapply(children, function(child) {
        child_switched <- label_switch(result[[child]], column = TRUE)
        odds_ratio(result[[root_node]], strategy_switched, child_switched)
      }, numeric(1)))
      if (after <= before) next
    }
    for (child in children) {
      result[[child]] <- label_switch(result[[child]], column = TRUE)
    }
    result[[strategy]] <- label_switch(result[[strategy]], column = FALSE)
  }

  strategy_or <- setNames(vapply(strategy_nodes, function(strategy) {
    odds_ratio(result[[root_node]], result[[strategy]])
  }, numeric(1)), strategy_nodes)
  if (!all(strategy_or >= 1)) {
    if (!all(strategy_or < 1)) {
      root_switched <- label_switch(result[[root_node]], column = FALSE)
      after <- sum(vapply(strategy_nodes, function(strategy) {
        strategy_switched <- label_switch(result[[strategy]], column = TRUE)
        odds_ratio(root_switched, strategy_switched)
      }, numeric(1)))
      if (after <= sum(strategy_or)) return(result)
    }
    for (strategy in strategy_nodes) {
      result[[strategy]] <- label_switch(result[[strategy]], column = TRUE)
    }
    result[[root_node]] <- label_switch(result[[root_node]], column = FALSE)
  }
  result
}
