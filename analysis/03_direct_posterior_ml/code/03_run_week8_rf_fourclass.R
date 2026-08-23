options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(S4Vectors)
  library(ranger)
  library(pROC)
  library(PRROC)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript 03_run_week8_rf_fourclass.R TECHNOLOGY TASK ",
    "[--include-exploratory=false] [--ntree=300] ",
    "[--max-per-class=100000] [--nthread=8] [--smoke-n=0]"
  )
}

technology <- args[[1L]]
task <- args[[2L]]

parse_options <- function(x) {
  x <- x[grepl("^--[^=]+=", x)]
  if (!length(x)) return(list())
  key <- sub("^--([^=]+)=.*$", "\\1", x)
  value <- sub("^--[^=]+=", "", x)
  stats::setNames(as.list(value), key)
}

cli <- parse_options(args[-seq_len(2L)])
value_or <- function(name, default) {
  value <- cli[[name]]
  if (is.null(value) || !nzchar(value)) default else value
}
as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

num_trees <- as.integer(value_or("ntree", "300"))
max_per_class <- as.integer(value_or("max-per-class", "100000"))
n_thread <- as.integer(value_or("nthread", "8"))
min_reliable <- as.integer(value_or("min-reliable", "1000"))
smoke_n <- as.integer(value_or("smoke-n", "0"))
include_exploratory <- as_flag(value_or("include-exploratory", "false"))
seed <- as.integer(value_or("seed", "20260812"))

if (!is.finite(num_trees) || num_trees < 1L) stop("Invalid --ntree.")
if (!is.finite(max_per_class) || max_per_class < 1L) {
  stop("Invalid --max-per-class.")
}
if (!is.finite(n_thread) || n_thread < 1L) stop("Invalid --nthread.")
if (!is.finite(min_reliable) || min_reliable < 1L) {
  stop("Invalid --min-reliable.")
}
if (!is.finite(smoke_n) || smoke_n < 0L) stop("Invalid --smoke-n.")

task_definitions <- list(
  calibration_tp_fp = c(negative = 2L, positive = 3L),
  error_fp_tn = c(negative = 0L, positive = 2L),
  artifact_tp_tn = c(negative = 0L, positive = 3L)
)
if (!task %in% names(task_definitions)) {
  stop("TASK must be calibration_tp_fp, error_fp_tn, or artifact_tp_tn.")
}

class_path <- value_or(
  "class-matrix",
  paste0(
    "/home/xiyao/15data/W8/results/01_four_class_global_fdr/",
    "four_class_matrix_Total_gt_10.rds"
  )
)
dataset_path <- value_or(
  "dataset",
  paste0(
    "/home/xichen/ConquerNet/ConquerBN_basic/",
    "m6AConquer_hg38_BenchmarkDataSet_Combined_SE_",
    "no_cali_IVT_matched.rds"
  )
)
feature_path <- value_or(
  "features",
  "/home/xiyao/15data/W8/input_data/features/m6Aconquer_omicsFeatures_hg38.rds"
)
design_path <- value_or(
  "rf-design",
  "/home/xiyao/15data/W8/input_data/rf_design/tp_fp_tp_tn_metadata.rds"
)
output_root <- value_or(
  "output-root",
  "/home/xiyao/15data/W8/results/03_rf_four_class_global_fdr"
)

for (path in c(class_path, dataset_path, feature_path, design_path)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

class_matrix <- readRDS(class_path)
dataset <- readRDS(dataset_path)
features <- readRDS(feature_path)
design <- readRDS(design_path)

if (!inherits(dataset, "SummarizedExperiment")) {
  stop("Dataset is not a SummarizedExperiment.")
}
if (!"Total" %in% assayNames(dataset)) stop("Dataset lacks Total assay.")
if (is.null(dim(class_matrix)) || is.null(colnames(class_matrix))) {
  stop("Class matrix must be a named matrix.")
}
technologies <- colnames(dataset)
if (!identical(colnames(class_matrix), technologies)) {
  stop("Class-matrix technology order differs from dataset.")
}
if (!identical(dim(class_matrix), c(nrow(dataset), ncol(dataset)))) {
  stop("Class-matrix dimensions differ from dataset.")
}
if (length(features) != nrow(dataset)) {
  stop("OmicsFeatures length differs from dataset rows.")
}
if (!technology %in% technologies) stop("Unknown technology: ", technology)
if (any(!is.na(class_matrix) & !class_matrix %in% 0:3)) {
  stop("Class matrix contains values outside 0:3/NA.")
}

## Full coordinate audit: RF rows must refer to the same genomic sites as labels.
dataset_ranges <- rowRanges(dataset)
coordinates_identical <- (
  identical(as.character(seqnames(features)), as.character(seqnames(dataset_ranges))) &&
    identical(start(features), start(dataset_ranges)) &&
    identical(end(features), end(dataset_ranges)) &&
    identical(as.character(strand(features)), as.character(strand(dataset_ranges)))
)
if (!coordinates_identical) stop("Feature/dataset genomic coordinates are misaligned.")

screening_feature_names <- as.character(design[["screening_features"]])
chromosome_fold <- design[["chromosome_fold"]]
if (length(screening_feature_names) != 185L) {
  stop("Expected the frozen 185-feature screening set; found ",
       length(screening_feature_names), ".")
}
if (!all(screening_feature_names %in% names(mcols(features)))) {
  stop("At least one frozen screening feature is absent from OmicsFeatures.")
}
if (is.null(names(chromosome_fold)) || !all(chromosome_fold %in% 1:5)) {
  stop("Invalid chromosome-fold mapping in frozen RF design.")
}

technology_column <- match(technology, technologies)
total <- assay(dataset, "Total")[, technology_column]
observed <- is.finite(total) & total > 10
class_vector <- as.integer(class_matrix[, technology_column])

## The upstream file is already Total>10. Recheck rather than silently masking a
## mismatched label file.
if (any(!is.na(class_vector) & !observed)) {
  stop("The supplied Total>10 class matrix has labels outside Total>10.")
}

codes <- task_definitions[[task]]
eligible_full <- which(observed & class_vector %in% unname(codes))
y_full <- as.integer(class_vector[eligible_full] == codes[["positive"]])
class_counts <- table(factor(y_full, levels = 0:1))
n_negative_full <- as.integer(class_counts[[1L]])
n_positive_full <- as.integer(class_counts[[2L]])

identifiability_status <- if (
  n_negative_full == 0L || n_positive_full == 0L
) {
  "untrainable_single_class"
} else if (
  n_negative_full >= min_reliable && n_positive_full >= min_reliable
) {
  "reliable"
} else {
  "exploratory_low_class_count"
}

strategy <- NA_character_
if (!is.null(design[["label_source_decisions"]])) {
  old_decision <- design[["label_source_decisions"]]
  hit <- old_decision[old_decision$technology == technology, , drop = FALSE]
  if (nrow(hit) == 1L && "strategy" %in% names(hit)) {
    strategy <- as.character(hit$strategy)
  }
}

result_dir <- file.path(output_root, task, technology)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
profile <- paste0("week8_direct_joint_globalFDR_Total_gt_10_185_", task)
status_path <- file.path(result_dir, "status.csv")

class_names <- c("TN", "FN", "FP", "TP")
status <- data.frame(
  technology = technology,
  strategy = strategy,
  task = task,
  profile = profile,
  class_definition = paste0(
    class_names[codes[["positive"]] + 1L], "=1; ",
    class_names[codes[["negative"]] + 1L], "=0"
  ),
  negative_code = unname(codes[["negative"]]),
  positive_code = unname(codes[["positive"]]),
  label_source = "direct_joint_four_class_posterior_global_FDR_lt_0.05",
  map_viterbi_fallback = FALSE,
  coverage_rule = "Total>10",
  coverage_cutoff = 10L,
  min_reliable = min_reliable,
  identifiability_status = identifiability_status,
  n_negative_full = n_negative_full,
  n_positive_full = n_positive_full,
  n_total_full = length(y_full),
  n_screening_features = length(screening_feature_names),
  num_trees = num_trees,
  max_per_class = max_per_class,
  n_thread = n_thread,
  include_exploratory = include_exploratory,
  smoke_n = smoke_n,
  run_status = "prepared",
  stringsAsFactors = FALSE
)

write.csv(status, status_path, row.names = FALSE)

if (identifiability_status == "untrainable_single_class") {
  status$run_status <- "not_trainable_single_class"
  status$completed <- as.character(Sys.time())
  write.csv(status, status_path, row.names = FALSE)
  cat("Skipped single-class task: ", technology, " / ", task, "\n", sep = "")
  quit(save = "no", status = 0L)
}
if (
  identifiability_status == "exploratory_low_class_count" &&
    !include_exploratory
) {
  status$run_status <- "skipped_exploratory_by_policy"
  status$completed <- as.character(Sys.time())
  write.csv(status, status_path, row.names = FALSE)
  cat("Skipped exploratory task by default: ", technology, " / ", task, "\n", sep = "")
  quit(save = "no", status = 0L)
}

eligible <- eligible_full
if (smoke_n > 0L && length(eligible) > smoke_n) {
  set.seed(seed)
  model_pick <- unlist(lapply(0:1, function(label) {
    candidates <- which(y_full == label)
    target <- min(length(candidates), ceiling(smoke_n / 2))
    if (length(candidates) > target) sample(candidates, target) else candidates
  }), use.names = FALSE)
  model_pick <- sort(unique(model_pick))
  eligible <- eligible_full[model_pick]
}
y_all <- as.integer(class_vector[eligible] == codes[["positive"]])
chromosome <- as.character(seqnames(features))[eligible]
fold_id <- unname(chromosome_fold[chromosome])
if (anyNA(fold_id) || !all(fold_id %in% 1:5)) {
  stop("Chromosome fold assignment failed.")
}

status$n_model_domain <- length(eligible)
status$n_model_negative <- sum(y_all == 0L)
status$n_model_positive <- sum(y_all == 1L)
write.csv(status, status_path, row.names = FALSE)

fold_precheck <- do.call(rbind, lapply(1:5, function(fold) {
  train_pool <- which(fold_id != fold)
  validation_pool <- which(fold_id == fold)
  data.frame(
    technology = technology,
    strategy = strategy,
    task = task,
    fold = fold,
    n_train_pool_negative = sum(y_all[train_pool] == 0L),
    n_train_pool_positive = sum(y_all[train_pool] == 1L),
    n_validation_negative = sum(y_all[validation_pool] == 0L),
    n_validation_positive = sum(y_all[validation_pool] == 1L),
    stringsAsFactors = FALSE
  )
}))
write.csv(fold_precheck, file.path(result_dir, "fold_precheck.csv"), row.names = FALSE)
if (any(
  fold_precheck$n_train_pool_negative == 0L |
    fold_precheck$n_train_pool_positive == 0L
)) {
  status$run_status <- "not_trainable_training_fold_single_class"
  status$completed <- as.character(Sys.time())
  write.csv(status, status_path, row.names = FALSE)
  quit(save = "no", status = 0L)
}

safe_logit <- function(p) qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))
prior_correct <- function(p, target_prevalence, sample_prevalence) {
  shift <- qlogis(target_prevalence) - qlogis(sample_prevalence)
  plogis(safe_logit(p) + shift)
}
metric_row <- function(y, p, fold) {
  valid <- is.finite(p) & !is.na(y)
  y <- y[valid]
  p <- p[valid]
  prevalence <- if (length(y)) mean(y) else NA_real_
  if (!length(y) || length(unique(y)) < 2L) {
    return(data.frame(
      technology = technology, strategy = strategy, task = task,
      fold = fold, n = length(y), n_negative = sum(y == 0L),
      n_positive = sum(y == 1L), prevalence = prevalence,
      auroc = NA_real_, auprc = NA_real_, stringsAsFactors = FALSE
    ))
  }
  roc_result <- pROC::roc(y, p, quiet = TRUE, direction = "<")
  pr_result <- PRROC::pr.curve(
    scores.class0 = p[y == 1L],
    scores.class1 = p[y == 0L],
    curve = FALSE
  )
  data.frame(
    technology = technology, strategy = strategy, task = task,
    fold = fold, n = length(y), n_negative = sum(y == 0L),
    n_positive = sum(y == 1L), prevalence = prevalence,
    auroc = as.numeric(pROC::auc(roc_result)),
    auprc = as.numeric(pr_result$auc.integral), stringsAsFactors = FALSE
  )
}

set.seed(seed)
metric_rows <- vector("list", 5L)
importance_rows <- vector("list", 5L)
fold_audit_rows <- vector("list", 5L)
oof <- data.frame(
  technology = technology,
  strategy = strategy,
  task = task,
  row_index = eligible,
  chromosome = chromosome,
  fold = fold_id,
  class_code = class_vector[eligible],
  y = y_all,
  raw_probability = NA_real_,
  prior_corrected_probability = NA_real_,
  stringsAsFactors = FALSE
)

for (fold in 1:5) {
  cat("Running ", technology, " / ", task, " / fold ", fold, "\n", sep = "")
  train_pool <- which(fold_id != fold)
  validation_pool <- which(fold_id == fold)
  train0 <- train_pool[y_all[train_pool] == 0L]
  train1 <- train_pool[y_all[train_pool] == 1L]
  sampled0 <- if (length(train0) > max_per_class) sample(train0, max_per_class) else train0
  sampled1 <- if (length(train1) > max_per_class) sample(train1, max_per_class) else train1
  train_index <- sort(c(sampled0, sampled1))

  train_prevalence <- mean(y_all[train_pool])
  sample_prevalence <- mean(y_all[train_index])
  train_x <- as.data.frame(
    mcols(features)[eligible[train_index], screening_feature_names, drop = FALSE]
  )
  validation_x <- as.data.frame(
    mcols(features)[eligible[validation_pool], screening_feature_names, drop = FALSE]
  )
  for (name in screening_feature_names) {
    train_x[[name]] <- as.numeric(train_x[[name]])
    validation_x[[name]] <- as.numeric(validation_x[[name]])
  }
  keep_feature <- vapply(train_x, function(x) {
    finite <- x[is.finite(x)]
    length(finite) > 1L && length(unique(finite)) > 1L
  }, logical(1))
  if (!any(keep_feature)) stop("All features are constant in fold ", fold)
  train_x <- train_x[, keep_feature, drop = FALSE]
  validation_x <- validation_x[, keep_feature, drop = FALSE]
  if (anyNA(train_x) || anyNA(validation_x)) {
    stop("Missing predictor values detected after feature conversion.")
  }

  train_data <- data.frame(
    y = factor(y_all[train_index], levels = 0:1),
    train_x,
    check.names = FALSE
  )
  model <- ranger(
    y ~ .,
    data = train_data,
    probability = TRUE,
    num.trees = num_trees,
    mtry = max(1L, floor(sqrt(ncol(train_x)))),
    min.node.size = 20L,
    importance = "permutation",
    num.threads = n_thread,
    seed = seed + fold,
    write.forest = TRUE
  )
  raw_probability <- predict(
    model, data = validation_x, num.threads = n_thread
  )$predictions[, "1"]
  corrected_probability <- prior_correct(
    raw_probability, train_prevalence, sample_prevalence
  )
  oof$raw_probability[validation_pool] <- raw_probability
  oof$prior_corrected_probability[validation_pool] <- corrected_probability
  metric_rows[[fold]] <- metric_row(
    y_all[validation_pool], corrected_probability, fold
  )

  importance_rows[[fold]] <- data.frame(
    technology = technology,
    strategy = strategy,
    task = task,
    fold = fold,
    feature = screening_feature_names,
    entered_model = screening_feature_names %in% names(model$variable.importance),
    importance = 0,
    stringsAsFactors = FALSE
  )
  importance_match <- match(
    names(model$variable.importance), importance_rows[[fold]]$feature
  )
  importance_rows[[fold]]$importance[importance_match] <- unname(
    model$variable.importance
  )

  fold_audit_rows[[fold]] <- data.frame(
    technology = technology, strategy = strategy, task = task, fold = fold,
    n_train_pool_negative = length(train0),
    n_train_pool_positive = length(train1),
    n_train_sample_negative = length(sampled0),
    n_train_sample_positive = length(sampled1),
    train_pool_prevalence = train_prevalence,
    train_sample_prevalence = sample_prevalence,
    n_validation_negative = sum(y_all[validation_pool] == 0L),
    n_validation_positive = sum(y_all[validation_pool] == 1L),
    n_features_candidate = length(screening_feature_names),
    n_features_entered = ncol(train_x), stringsAsFactors = FALSE
  )
  rm(train_x, validation_x, train_data, model)
  gc()
}

if (
  anyNA(oof$raw_probability) ||
    anyNA(oof$prior_corrected_probability) ||
    any(!is.finite(oof$raw_probability)) ||
    any(!is.finite(oof$prior_corrected_probability)) ||
    any(oof$raw_probability < 0 | oof$raw_probability > 1) ||
    any(oof$prior_corrected_probability < 0 |
          oof$prior_corrected_probability > 1) ||
    anyDuplicated(oof$row_index)
) {
  stop("OOF integrity check failed.")
}

metrics <- do.call(rbind, metric_rows)
importance <- do.call(rbind, importance_rows)
fold_audit <- do.call(rbind, fold_audit_rows)
write.csv(metrics, file.path(result_dir, "fold_metrics.csv"), row.names = FALSE)
write.csv(
  importance, file.path(result_dir, "permutation_importance.csv"), row.names = FALSE
)
write.csv(fold_audit, file.path(result_dir, "fold_audit.csv"), row.names = FALSE)
saveRDS(oof, file.path(result_dir, "oof_predictions.rds"), compress = TRUE)

finite_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
status$run_status <- "complete"
status$completed <- as.character(Sys.time())
status$n_evaluable_folds <- sum(is.finite(metrics$auroc) & is.finite(metrics$auprc))
status$mean_auroc <- finite_mean(metrics$auroc)
status$mean_auprc <- finite_mean(metrics$auprc)
write.csv(status, status_path, row.names = FALSE)
cat(
  "Completed: ", technology, " / ", task,
  "; AUROC=", sprintf("%.4f", status$mean_auroc),
  "; AUPRC=", sprintf("%.4f", status$mean_auprc), "\n", sep = ""
)
