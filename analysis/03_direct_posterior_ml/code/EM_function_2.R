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

meta_format_check <- function(meta){
  if (!all(c("Technique", "Strategy") %in% colnames(meta))) {
    stop("Metadata should contain 'Technique' and 'Strategy' column")
  }
  
  if(has_technique_strategy_issue(meta)) {
    stop("At least one technique is assigned with more than one strategy.")
  }
  return(meta)
}

input_format_check <- function(m6A, Total, meta){
  if ((any(colnames(m6A) != colnames(Total)))){
    stop("Column names of m6A and Total count matrices are not matched")
  }
  
  if(ncol(m6A) != nrow(meta)){
    stop("Column numbers of m6A and Total count matrices are not matched to metadata row numbers")
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
    theta[[nodes]] <- cbind("0" = dbbinom(m6A[,nodes], Total[,nodes],
                                          beta.param[[nodes]]["alpha","0"],
                                          beta.param[[nodes]]["beta","0"], log = TRUE),
                            "1" = dbbinom(m6A[,nodes], Total[,nodes],
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
    if (nodes == root.nodes(dag)) {
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
  
  
  for (nodes in children(dag, root.nodes(dag))){
    passage <- paste0("TrueM6A|",nodes)
    passage_down <- paste0(children(dag, root.nodes(dag))[children(dag, root.nodes(dag)) != nodes], "|TrueM6A")
    message[[passage]] <- psi[[root.nodes(dag)]] + Reduce("+", message[passage_down])
  }
  
  for (nodes in children(dag, root.nodes(dag))){
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
    if (nodes == root.nodes(dag)) {
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
  cpt.param[[root.nodes(dag)]][,1] <- temp/sum(temp)
  
  
  for (nodes in names(cpt.param)[2:length(cpt.param)]) {
    parent <- parents(dag, nodes)
    if (parent == root.nodes(dag)) {
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
    or.z1[i] <- odds_ratio(cpt.param.estimated[[root.nodes(dag)]], cpt.param.estimated[[z2]], cpt.param.estimated[[z1]])
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
        or.sum.after <- or.sum.after + odds_ratio(cpt.param.estimated[[root.nodes(dag)]], z2.temp, z1.temp)
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
    or.z2[i] <- odds_ratio(cpt.param.estimated[[root.nodes(dag)]], cpt.param.estimated[[z2]])
  }
  
  if(!(all(or.z2 >= 1))){
    if(!(all(or.z2 < 1))) {
      or.sum.after <- 0
      z3.temp <- label_switch(cpt.param.estimated[[root.nodes(dag)]], column = FALSE)
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
    cpt.param.estimated[[root.nodes(dag)]] <- label_switch(cpt.param.estimated[[root.nodes(dag)]], column = FALSE)
  }
  return(cpt.param.estimated)
  # saveRDS(cpt.param.estimated, paste0(path_file, "cpt.param.estimated.switched.rds"))
}








conquerNet_inference <- function(emission_loglik, cpt.param, dag, mode = c("sum-product", "max-sum")){
  
  # 1. 验证参数并匹配模式
  mode <- match.arg(mode)
  
  # 2. 动态绑定核心算子
  # Sum-product 使用 rowLogSumExps; Max-sum 使用 rowMaxs
  agg_fun <- if (mode == "sum-product") matrixStats::rowLogSumExps else matrixStats::rowMaxs
  
  # ==========================================
  # 阶段 1：数据与网络拓扑初始化 (两者完全共享)
  # ==========================================
  trail <- sapply(continuous.nodes, function(i) ancestors(dag, i))
  N <- nrow(emission_loglik[[1]])
  theta <- emission_loglik
  
  for (nodes in names(cpt.param)) {
    if (nodes == root.nodes(dag)) {
      theta[[nodes]] <- matrix(rep(log(cpt.param[[nodes]]), each = N), 
                               nrow = N, byrow = FALSE, 
                               dimnames = list(NULL, dimnames(cpt.param[[root.nodes(dag)]])[[1]]))
    } else {
      theta[[nodes]] <- array(rep(log(cpt.param[[nodes]]), each = N), 
                              dim = c(N, dim(cpt.param[[nodes]])), 
                              dimnames = list(NULL, 
                                              dimnames(cpt.param[[root.nodes(dag)]])[[1]],
                                              dimnames(cpt.param[[root.nodes(dag)]])[[2]]))
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
  
  # ==========================================
  # 阶段 2：双向消息传递 (仅算子不同，逻辑共享)
  # ==========================================
  
  # --- 自下而上 (Collect) ---
  k.nodes <- sapply(leaf.nodes(dag), function(i) parents(dag, i))
  for (nodes in k.nodes){
    passage <- paste0(nodes,"|",parents(dag, nodes))
    # 使用动态绑定的 agg_fun 替代写死的算子
    message[[passage]] <- apply(psi[[nodes]], 3, agg_fun)
  }
  
  c.nodes <- children(dag, root.nodes(dag))
  for (nodes in c.nodes) {
    passage <- paste0(nodes,"|",parents(dag, nodes))
    message_sum <- Reduce("+", message[paste0(children(dag, nodes),"|",nodes)])
    temp <- expand_array(message_sum, psi[[nodes]]) + psi[[nodes]]
    message[[passage]] <- apply(temp, 3, agg_fun)
  }
  
  # --- 自上而下 (Distribute) ---
  for (nodes in children(dag, root.nodes(dag))){
    passage <- paste0(root.nodes(dag),"|",nodes)
    passage_down <- paste0(children(dag, root.nodes(dag))[children(dag, root.nodes(dag)) != nodes], "|", root.nodes(dag))
    message[[passage]] <- psi[[root.nodes(dag)]] + Reduce("+", message[passage_down])
  }
  
  for (nodes in children(dag, root.nodes(dag))){
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
      # 使用动态绑定的 agg_fun
      message[[passage]] <- apply(temp, 2, agg_fun)
    }
  }
  
  # ==========================================
  # 阶段 3：解码输出 (终点逻辑分化)
  # ==========================================
  
  res <- list()
  
  if (mode == "sum-product") {
    # ---------------------------------------------
    # 模式 A: 计算后验边缘概率 (用于 EM 算法 E-step)
    # ---------------------------------------------
    belief <- list()
    for (nodes in names(cpt.param)){
      passage_down <- paste0(children(dag, nodes), "|", nodes)
      if (nodes == root.nodes(dag)) {
        temp <- Reduce("+", message[passage_down]) + psi[[nodes]]
        ll <- sum(agg_fun(temp)) # 这里的 agg_fun 自动是 rowLogSumExps
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
      rm(temp)
    }
    res$loglik <- ll
    res$belief <- belief
    rm(theta, psi, message)
    
  } else if (mode == "max-sum") {
    # ---------------------------------------------
    # 模式 B: 提取全局最优组合 (维特比硬解码)
    # ---------------------------------------------
    max_marginals <- list()
    optimal_states <- list()
    
    for (nodes in names(cpt.param)){
      passage_down <- paste0(children(dag, nodes), "|", nodes)
      if (nodes == root.nodes(dag)) {
        temp <- Reduce("+", message[passage_down]) + psi[[nodes]]
        ll <- sum(agg_fun(temp)) # 这里的 agg_fun 自动是 rowMaxs (即 MAP 对数似然)
        max_marginals[[nodes]] <- temp
        optimal_states[[nodes]] <- colnames(temp)[max.col(temp, ties.method = "first")]
      } else {
        passage_up <- paste0(parents(dag, nodes),"|",nodes)
        message_up <- expand_array(message[[passage_up]], psi[[nodes]], FALSE)
        temp <- message_up + psi[[nodes]]
        if (passage_down[1] %in% names(message)) {
          message_down <- expand_array(Reduce("+", message[passage_down]), psi[[nodes]])
          temp <- temp + message_down
        }
        max_marginals[[nodes]] <- temp
        node_temp_2d <- apply(temp, 2, matrixStats::rowMaxs)
        optimal_states[[nodes]] <- colnames(node_temp_2d)[max.col(node_temp_2d, ties.method = "first")]
        rm(node_temp_2d)
      }
      rm(temp)
    }
    res$loglik <- ll     # MAP Score
    res$max_marginals <- max_marginals # 最大势函数矩阵
    res$optimal_states <- optimal_states # 最优路径 (Viterbi Path)
    rm(theta, psi, message)
  }
  gc(verbose = FALSE)
  return(res)
}

global_fdr <- function(local.fdr, order = FALSE){
  idx <- order(local.fdr)
  lfdr.ordered <- local.fdr[idx]
  gfdr <- cumsum(lfdr.ordered) / seq_along(lfdr.ordered)
  df <- data.frame(
    idx = idx,
    local = lfdr.ordered,
    global = gfdr
  )
  if (!order){
    df <- df[order(df$idx),]
  }
  return(df)
}

integrate_tri_decoding <- function(viterbi_mat, fdr_mat, for_mat) {
  
  # 1. 维度安全检查
  if (!all(dim(viterbi_mat) == dim(fdr_mat)) || !all(dim(viterbi_mat) == dim(for_mat))) {
    stop("错误：三个矩阵的维度必须完全一致！")
  }
  
  # 2. 核心魔法：三维二进制加权编码 (向量化极速运算)
  # Viterbi 占 4分，FDR 占 2分，FOR 占 1分
  score_mat <- (viterbi_mat * 4) + (fdr_mat * 2) + (for_mat * 1)
  
  # 3. 初始化文本矩阵
  class_mat <- matrix("Unknown", nrow = nrow(score_mat), ncol = ncol(score_mat))
  dimnames(class_mat) <- dimnames(viterbi_mat)
  
  # 4. 状态映射 (完美穷举)
  
  anno <- data.frame(
    Class = as.character(c(seq(0,5), -1)),
    Annotation = c("High-Conf 0 (True Background)",
                   "High-Conf 1 (Gold Standard)",
                   "Low-Conf 1 (Viterbi Guess)",
                   "Low-Conf 0 (Viterbi Guess)",
                   "Rescue 1 (Viterbi Suppressed)",
                   "Reject 1 (Viterbi Hallucination)",
                   "Error: FDR/FOR Paradox"),
    Viterbi = c("0","1","1","0","0","1","0/1"),
    FDR = as.character(c(0,1,0,0,1,0,1)),
    FOR = as.character(c(1,0,0,0,0,1,1))
  )
  
  # --- 完美一致 (双重认证) ---
  # 4+2+0 = 6: Viterbi=1, FDR=1, FOR=0 
  # class_mat[score_mat == 6] <- "High-Conf 1 (Gold Standard)"  
  class_mat[score_mat == 6] <- 1 
  
  # 0+0+1 = 1: Viterbi=0, FDR=0, FOR=1
  # class_mat[score_mat == 1] <- "High-Conf 0 (True Background)" 
  class_mat[score_mat == 1] <- 0
  
  # --- 缺乏强证 (Viterbi 孤立推演) ---
  # 4+0+0 = 4: Viterbi=1, 但 FDR 和 FOR 都不置信
  # class_mat[score_mat == 4] <- "Low-Conf 1 (Viterbi Guess)" 
  class_mat[score_mat == 4] <- 2 
  
  # 0+0+0 = 0: Viterbi=0, 但 FDR 和 FOR 都不置信
  # class_mat[score_mat == 0] <- "Low-Conf 0 (Viterbi Guess)"
  class_mat[score_mat == 0] <- 3
  
  # --- 价值连城的冲突点 (网络拓扑 vs 局部信号) ---
  # 0+2+0 = 2: Viterbi=0, 但 FDR 说它是 1！(假阴性漏网之鱼)
  # class_mat[score_mat == 2] <- "Rescue 1 (Viterbi Suppressed)" 
  class_mat[score_mat == 2] <- 4
  
  # 4+0+1 = 5: Viterbi=1, 但 FOR 说它是 0！(假阳性系统幻觉)
  # class_mat[score_mat == 5] <- "Reject 1 (Viterbi Hallucination)"
  class_mat[score_mat == 5] <- 5
  
  # --- 悖论预警 ---
  # FDR=1 且 FOR=1 (概率上不可能同时是高置信1和高置信0)
  # class_mat[score_mat == 3 | score_mat == 7] <- "Error: FDR/FOR Paradox"
  class_mat[score_mat == 3 | score_mat == 7] <- -1
  
  
  # 5. 汇总统计
  summary_stats <- table(class_mat)
  stats.class <- names(summary_stats)
  stat.anno <- sapply(stats.class, function(i) anno$Annotation[anno$Class == i])
  
  names(summary_stats) <- stat.anno
  cat("=== 三路联合解码状态分布 ===\n")
  print(summary_stats)
  
  mode(class_mat) <- "numeric"
  
  return(list(
    class_matrix = class_mat,     # 类别鉴定矩阵
    score_matrix = score_mat,    # 数字分类矩阵 (用于 Heatmap)
    stats = summary_stats,
    anno = anno)
  )
}

technique_evaluation <- function(viterbi.res, dag) {
  root.nodes <- root.nodes(dag)
  leaf.nodes <- leaf.nodes(dag)
  tech.nodes <- sapply(leaf.nodes, function(i) parents(dag, i))
  
  eval.mat <- matrix(NA, nrow = nrow(viterbi.res), ncol = length(tech.nodes))
  root.score <- viterbi.res[,root.nodes]
  tech.score <- viterbi.res[,tech.nodes] * 2
  eval.score <- tech.score + root.score
  
  anno <- data.frame(
    Class = as.character(c(seq(0,3))),
    Annotation = c("TN",
                   "FN",
                   "FP",
                   "TP"),
    m6A = as.character(c(0,1,0,1)),
    Tech = as.character(c(0,0,1,1))
  )
  
  return(list(score = eval.score, anno = anno))
}
