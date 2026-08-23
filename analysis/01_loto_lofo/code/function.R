library(bnlearn)
library(extraDistr)
library(ggplot2)
library(tidyr)  
library(dplyr)
library(SummarizedExperiment)

has_technique_strategy_issue <- function(meta) {
  
  any_issue <- meta %>%
    group_by(Technique) %>%
    summarise(count = n_distinct(Strategy)) %>%
    filter(count > 1) %>%
    nrow() > 0
  
  return(any_issue)
}

merge_duplicate_entries <- function(meta, data_matrix, merge_function = "sum") {
  cat("\nMerging duplicates...\n")
  
  unique_combos <- unique(meta)
  unique_combos$Combo_ID <- 1:nrow(unique_combos)
  
  meta$Combo_ID <- apply(meta, 1, function(row) {
    which(unique_combos$Technique == row[1] & unique_combos$Strategy == row[2])[1]
  })
  
  merged_matrix_list <- list()
  new_meta_list <- list()
  
  for (combo_id in unique(meta$Combo_ID)) {
    row_indices <- which(meta$Combo_ID == combo_id)
    
    if (length(row_indices) == 1) {
      merged_matrix_list[[as.character(combo_id)]] <- data_matrix[, row_indices, drop = FALSE]
    } else {
      cat(sprintf("Merging %d duplicated samples from '%s' data\n", 
                  length(row_indices), meta$Technique[row_indices[1]]))
      
      cols_to_merge <- data_matrix[, row_indices, drop = FALSE]
      
      if (merge_function == "sum") {
        merged_col <- rowSums(cols_to_merge)
      } else if (merge_function == "mean") {
        merged_col <- rowMeans(cols_to_merge)
      } else if (merge_function == "median") {
        merged_col <- apply(cols_to_merge, 1, median)
      } else {
        stop("No valid merging function。Please use 'sum', 'mean', or 'median'")
      }
      
      merged_matrix_list[[as.character(combo_id)]] <- matrix(merged_col, ncol = 1)
    }
    
    new_meta_list[[as.character(combo_id)]] <- meta[row_indices[1], c("Technique", "Strategy")]
  }
  
  new_meta <- do.call(rbind, new_meta_list)
  rownames(new_meta) <- NULL
  
  new_data_matrix <- do.call(cbind, merged_matrix_list)
  colnames(new_data_matrix) <- new_meta$Technique
  
  return(list(
    meta = new_meta,
    data_matrix = new_data_matrix
  ))
}

format_check <- function(m6A, Total, meta){
  if (!all(c("Technique", "Strategy") %in% colnames(meta))) {
    stop("Metadata should contain 'Technique' and 'Strategy' column")
  }
  
  if (any(colnames(m6A) != colnames(Total))){
    stop("Column names of m6A and Total count matrices are not matched")
  }
  
  if(any(colnames(m6A) != paste0(meta$Technique,"_x"))) {
    stop("Column names of m6A and Total count matrices are not matched to metadata")
  }
  
  if(has_technique_strategy_issue(meta)) {
    stop("At least one technique is assigned with more than one strategy.")
  }
  
  if(any(duplicated(meta) | duplicated(meta, fromLast = TRUE))){
    res <- merge_duplicate_entries(meta, m6A)
    m6A <- res$data_matrix
    res <- merge_duplicate_entries(meta, Total)
    Total <- res$data_matrix
    meta <- res$meta
  }
  
  return(list(m6A = m6A, Total = Total, meta = meta))
}

generate_random_cpt <- function(is_root = FALSE) {
  if (is_root) {
    # 根节点：生成一个2维概率向量
    probs <- runif(2) # 生成两个随机数
    cpt <- matrix(probs / sum(probs),
                  ncol = 1, byrow = TRUE,
                  dimnames = list(c(0, 1), NULL)) # 归一化确保和为1
    # cpt <- c(0.5,0.5)
    return(cpt)
  } else {
    # 非根节点（单父节点，二值）：生成一个2x2概率矩阵
    # 每行对应父节点的一个状态，需要独立归一化
    cpt <- matrix(runif(4), ncol = 2, byrow = TRUE, dimnames = list(c(0, 1), c(0, 1))) # 生成4个随机数填充2x2矩阵
    cpt <- t(t(cpt) / colSums(cpt)) # 按列归一化
    return(cpt)
  }
}

initialization_discrete <- function(seed=NULL, discrete.nodes){
  is.root <- c(TRUE, rep(FALSE, length(discrete.nodes)-1))
  if(is.null(seed)){
    cpt.list <- lapply(is.root, generate_uniform_cpt)
  } else {
    set.seed(seed)
    cpt.list <- lapply(is.root, generate_random_cpt)
  }
  names(cpt.list) <- discrete.nodes
  return(cpt.list)
}

apply_dirichlet <- function(mat, alpha = 0.1) {
  mat_smoothed <- mat + alpha
  sweep(mat_smoothed, MARGIN = 2, STATS = colSums(mat_smoothed), FUN = "/")
}

initialization_cpt_params <- function(m6A, Total, beta.param, weight, dag){
  cpt <- list()
  cpt[[root.nodes(dag)]] <- matrix(c(0.9, 0.1), nrow = 2, byrow = TRUE,
                                   dimnames = list(c("0", "1"), NULL))
  k.nodes <- sapply(leaf.nodes(dag), function(i) parents(dag,i))
  c.nodes <- children(dag, root.nodes(dag))
  
  gamma_k <- weight
  gamma_k[is.na(gamma_k)] <- 0

  # background gamma mean
  gamma_c <- matrix(NA, ncol = length(c.nodes), 
                    nrow = nrow(m6A),
                    dimnames = list(NULL, c.nodes))
  
  for (c in c.nodes){
    cpt[[c]] <- matrix(c(0.9,0.1,
                         0.1,0.9),
                       nrow = 2, byrow = TRUE,
                       dimnames = list(c("0","1"),c("0","1")))
    c.k <- children(dag,c)
    
    gamma_c_0 <- rep(0, nrow(m6A))
    
    for (k in c.k) {
      x <- children(dag, k)
      gamma_c_0 <- gamma_c_0 + (1-gamma_k[,x])
    }
    
    gamma_c_0 <- gamma_c_0/length(c.k)
    
    gamma_c[,c] <- gamma_c_0
  }
  
  rough_labels <- ifelse(gamma_c <= 0.5, "1", "0")
  
  for (k in k.nodes) {
    k.c <- parents(dag, k)
    k.x <- children(dag, k)
    K <- ncol(beta.param[[k.x]])
    
    states <- switch (as.character(K),
                      "2" = c("0","1"),
                      "3" = c("-1","0","1"),
                      as.character(seq(0,K-1,1))
    )
    
    count.k.given.c <- matrix(0, nrow = K, ncol = 2, dimnames = list(states, c("0","1")))
    
    count.k.given.c[, "0"] <- c(sum(1-gamma_k[rough_labels[,k.c] == "0", k.x]), sum(gamma_k[rough_labels[,k.c] == "0", k.x]))
    count.k.given.c[, "1"] <- c(sum(1-gamma_k[rough_labels[,k.c] == "1", k.x]), sum(gamma_k[rough_labels[,k.c] == "1", k.x]))
    
    cpt[[k]] <- apply_dirichlet(count.k.given.c, alpha = 0.5)
  }
  
  return(cpt)
}


get_parent_chain_bn <- function(node_name, bn_object) {
  path <- c()
  current <- node_name
  
  while(TRUE) {
    parents <- bnlearn::parents(bn_object, current)
    if(length(parents) == 0) break
    if(length(parents) > 1) warning(paste("node", current, "has more than one parent node; only use the first parent node."))
    path <- c(parents[1], path)
    current <- parents[1]
  }
  
  return(path)
}

expand_array <- function(x, target, child = TRUE){
  # child indicates whether x is the child of target
  if (child) {
    res <- array(rep(x, dim(target)[3]), dim = c(dim(target)))
  } else {
    # res <- aperm(array(rep(x, dim(target)[2]), dim = c(dim(target))), c(1,3,2))
    res <- array(0, dim = dim(target))
    for(p in 1:dim(target)[3]) {
      # R 会自动把 N 长度的向量，重复填充为一个 N x Child_states 的矩阵
      res[,,p] <- x[,p]
    }
  }
  return(res)
}

param_update_strategy <- function(belief_current, belief_parent){
  unnorm <- exp(sweep(log(colSums(exp(belief_current))),2,log(colSums(exp(belief_parent))),FUN = "-"))
  norm <- normalization(unnorm, log = FALSE)
  return(norm)
}

param_update_technique <- function(belief_current, belief_parent){
  unnorm <- exp(sweep(log(colSums(exp(belief_current))),2,log(rowSums(colSums(exp(belief_parent)))),FUN = "-"))
  norm <- normalization(unnorm, log = FALSE)
  return(norm)
}

normalization <- function(p, log=TRUE){
  if (log == FALSE) {
    p <- log(p)
  }
  p.max <- apply(p, 2, max)
  p.scaled <- sweep(p, 2, p.max, FUN = "-")
  p.exp <- exp(p.scaled)
  p.norm <- sweep(p.exp, 2, colSums(p.exp), FUN = "/")
  return(p.norm)
}

log_belief_normalization <- function(beta){
  beta_2d <- matrix(beta, nrow = dim(beta)[1], ncol = prod(dim(beta)[2:length(dim(beta))]))
  norm_const <- matrixStats::rowLogSumExps(beta_2d)
  beta_2d <- beta_2d - norm_const
  log_belief <- array(beta_2d, dim = dim(beta))
}

calculate_emission_loglik <- function(m6A, Total, beta.param){
  theta <- list()
  for (nodes in names(beta.param)) {
    theta[[nodes]] <- cbind(dbbinom(m6A[,nodes], Total[,nodes],
                                    beta.param[[nodes]]["alpha","0"],
                                    beta.param[[nodes]]["beta","0"], log = TRUE),
                            dbbinom(m6A[,nodes], Total[,nodes],
                                    beta.param[[nodes]]["alpha","1"],
                                    beta.param[[nodes]]["beta","1"], log = TRUE))
  }
  return(theta)
}

conquerNet_expectation <- function(emission_loglik, cpt.param, dag){
  trail <- sapply(continuous.nodes, function(i) ancestors(dag, i))
  N <- nrow(emission_loglik[[1]])
  theta <- emission_loglik
  for (nodes in names(cpt.param)) {
    if (nodes == "TrueM6A") {
      theta[[nodes]] <- matrix(rep(log(cpt.param[[nodes]]), each = N), nrow = N, byrow = FALSE)
    } else {
      theta[[nodes]] <- array(rep(log(cpt.param[[nodes]]), each = N), dim = c(N, dim(cpt.param[[nodes]])))
    }
  }
  
  psi <- list()
  for (nodes in names(cpt.param)) {
    child <- children(dag, nodes)
    if (child[1] %in% names(cpt.param)) {
      psi[[nodes]] <- theta[[nodes]]
    } else {
      if (length(child) != 1) {
        stop("More than one replicates for one techniques. Do replicates aggregation first.")
      }
      theta_child <- theta[[child]]
      theta_parent <- theta[[nodes]]
      psi[[nodes]] <- expand_array(theta_child, theta_parent) + theta_parent
    }
  }
  
  message <- list()
  
  for (nodes in trail[1,]){
    passage <- paste0(nodes,"|",parents(dag, nodes))
    message[[passage]] <- apply(psi[[nodes]], 3, matrixStats::rowLogSumExps)
  }
  
  for (nodes in unique(c(trail[2,]))) {
    passage <- paste0(nodes,"|",parents(dag, nodes))
    message_sum <- Reduce("+",message[paste0(children(dag, nodes),"|",nodes)])
    temp <- expand_array(message_sum, psi[[nodes]]) + psi[[nodes]]
    message[[passage]] <- apply(temp, 3, matrixStats::rowLogSumExps)
  }
  
  
  for (nodes in children(dag, "TrueM6A")){
    passage <- paste0("TrueM6A|",nodes)
    passage_down <- paste0(children(dag, "TrueM6A")[children(dag, "TrueM6A") != nodes], "|TrueM6A")
    message[[passage]] <- psi[["TrueM6A"]] + Reduce("+", message[passage_down])
  }
  
  for (nodes in children(dag, "TrueM6A")){
    for (child in children(dag, nodes)) {
      passage <- paste0(nodes, "|", child)
      passage_up <- paste0(parents(dag, nodes),"|",nodes)
      message_up <- message[[passage_up]]
      temp <- expand_array(message_up, psi[[nodes]], FALSE) + psi[[nodes]]
      if (length(children(dag, nodes)) != 1) {
        passage_down <- paste0(children(dag, nodes)[children(dag, nodes) != child], "|", nodes)
        message_down <- Reduce("+", message[passage_down])
        temp <- temp + expand_array(message_down, psi[[nodes]])
      }
      message[[passage]] <- apply(temp, 2, matrixStats::rowLogSumExps)
    }
  }
  
  belief <- list()
  
  for (nodes in names(cpt.param)){
    passage_down <- paste0(children(dag, nodes), "|", nodes)
    if (nodes == "TrueM6A") {
      temp <- Reduce("+", message[passage_down]) + psi[[nodes]]
      ll <- sum(matrixStats::rowLogSumExps(temp))
      belief[[nodes]] <- log_belief_normalization(temp)
    } else {
      passage_up <- paste0(parents(dag, nodes),"|",nodes)
      message_up <- expand_array(message[[passage_up]], psi[[nodes]], FALSE)
      temp <- message_up + psi[[nodes]]
      if (passage_down[1] %in% names(message)) {
        message_down <- expand_array(Reduce("+", message[passage_down]), psi[[nodes]])
        temp <- temp + message_down
      }
      belief[[nodes]] <- log_belief_normalization(temp)
    }
  }
  
  res <- list(loglik = ll, belief = belief)
  
  return(res)
}

conquerNet_maximization <- function(cpt.param, belief.list, dag){
  # Parameter update
  
  temp <- colMeans(exp(belief.list[[1]]))
  cpt.param[["TrueM6A"]][,1] <- temp/sum(temp)
  
  
  for (nodes in names(cpt.param)[2:length(cpt.param)]) {
    parent <- parents(dag, nodes)
    if (parent == "TrueM6A") {
      cpt.param[[nodes]] <- param_update_strategy(belief[[nodes]], belief[[parent]])
    } else {
      cpt.param[[nodes]] <- param_update_technique(belief[[nodes]], belief[[parent]])
    }
    dimnames(cpt.param[[nodes]]) <- list(c("0","1"), c("0","1"))
  }
  
  return(cpt.param)
}

plot_likelihood <- function(path_file){
  log.likelihood.track <- readRDS(paste0(path_file,"log.likelihood.track.rds"))
  p_ll <- ggplot(data.frame(
    log_likelihood = unlist(log.likelihood.track),
    iteration = seq(1:length(log.likelihood.track))),
    aes(x = iteration, y = log_likelihood)) +
    geom_line(color = "black", linewidth = 1) +  # 折线
    labs(
      title = "Log likelihood change",
      x = "Iteration",
      y = "Log-likihood"
    ) +
    theme_classic() +  # 简洁主题
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12))
  ggsave(paste0(path_file, "LogLikelihood.pdf"), p_ll, height = 3.5, width = 4, dpi = 320)
  ggsave(paste0(path_file, "LogLikelihood.png"), p_ll, height = 3.5, width = 4, dpi = 320)
}

label_switch <- function(cpt, column = TRUE){
  temp <- cpt
  if (column) {
    cpt[,1] <- temp[,2]
    cpt[,2] <- temp[,1]
  } else {
    cpt[1,] <- temp[2,]
    cpt[2,] <- temp[1,]
  }
  return(cpt)
}

odds_ratio <- function(z3.param, z2.param, z1.param=NULL){
  z3.marginal.log <- log(z3.param)
  z3.marginal.log <- rbind(t(z3.marginal.log),t(z3.marginal.log))
  z2.given.z3.log <- log(z2.param)
  z2.joint.z3 <- exp(z2.given.z3.log+z3.marginal.log)
  if (is.null(z1.param)) {
    or <- (z2.joint.z3[1,1]*z2.joint.z3[2,2])/(z2.joint.z3[1,2]*z2.joint.z3[2,1]) 
  } else {
    z2.marginal.log <- log(rowSums(z2.joint.z3))
    z1.given.z2.log <- log(z1.param)
    z1.joint.z2 <- exp(z1.given.z2.log+z2.marginal.log)
    or <- (z1.joint.z2[1,1]*z1.joint.z2[2,2])/(z1.joint.z2[1,2]*z1.joint.z2[2,1])
  }
  return(or)
}

latent_param_label_switch <- function(cpt.param, path){
  cpt.param.estimated <- cpt.param
  or.z1 <- matrix(NA, ncol = ncol(path))
  colnames(or.z1) <- path[3,]
  vote.z2 <- matrix(0, nrow = 2, ncol = length(unique(path[2,])))
  colnames(vote.z2) <- unique(path[2,])
  rownames(vote.z2) <- c("Total", "Agreed")
  vote.z2["Total", ] <- table(path[2,])[colnames(vote.z2)]
  # Calculate or
  for (i in 1:length(unique(path[3,]))) {
    z1 <- path[3,i]
    z2 <- path[2,i]
    or.z1[i] <- odds_ratio(cpt.param.estimated[["TrueM6A"]], cpt.param.estimated[[z2]], cpt.param.estimated[[z1]])
  }
  
  z1.switch.col <- path[3, (or.z1 < 1)]
  z2.switch.row <- path[2, (or.z1 < 1)]
  vote.z2["Agreed", ] <- table(z2.switch.row)[colnames(vote.z2)]
  vote.z2[is.na(vote.z2)] <- 0
  
  for (z2 in unique(path[2,])){
    if(vote.z2["Total", z2] == 0){
      next
    }
    if (vote.z2["Total", z2] != vote.z2["Agreed", z2]) {
      or.sum.before <- sum(or.z1[1, path[3, (path[2, ] == z2)]])
      or.sum.after <- 0
      z2.temp <- label_switch(cpt.param.estimated[[z2]], column = FALSE)
      for (z1 in path[3, (path[2, ] == z2)]){
        z1.temp <- label_switch(cpt.param.estimated[[z1]], column = TRUE)
        or.sum.after <- or.sum.after + odds_ratio(cpt.param.estimated[["TrueM6A"]], z2.temp, z1.temp)
      }
      if (or.sum.after <= or.sum.before) {
        next
      }
    }
    for (z1 in path[3, (path[2, ] == z2)]){
      cpt.param.estimated[[z1]] <- label_switch(cpt.param.estimated[[z1]], column = TRUE)
    }
    cpt.param.estimated[[z2]] <- label_switch(cpt.param.estimated[[z2]], column = FALSE)
  }
  
  or.z2 <- matrix(NA, ncol = length(unique(path[2,])))
  colnames(or.z2) <- unique(path[2,])
  
  for (i in 1:length(unique(path[2,]))) {
    z2 <- path[2,i]
    or.z2[i] <- odds_ratio(cpt.param.estimated[["TrueM6A"]], cpt.param.estimated[[z2]])
  }
  
  if(!(all(or.z2 >= 1))){
    if(!(all(or.z2 < 1))) {
      or.sum.after <- 0
      z3.temp <- label_switch(cpt.param.estimated[["TrueM6A"]], column = FALSE)
      for (i in 1:length(unique(path[2,]))) {
        z2.temp <- label_switch(cpt.param.estimated[[z2]], column = TRUE)
        or.sum.after <- or.sum.after + odds_ratio(z3.temp, z2.temp)
      }
      if(or.sum.after <= sum(or.z2)) {
        break
      }
    }
    
    for(z2 in unique(path[2,])) {
      cpt.param.estimated[[z2]] <- label_switch(cpt.param.estimated[[z2]], column = TRUE)
    }
    cpt.param.estimated[["TrueM6A"]] <- label_switch(cpt.param.estimated[["TrueM6A"]], column = FALSE)
  }
  return(cpt.param.estimated)
  # saveRDS(cpt.param.estimated, paste0(path_file, "cpt.param.estimated.switched.rds"))
}