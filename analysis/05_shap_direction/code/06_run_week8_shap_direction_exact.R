options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(S4Vectors)
  library(ranger)
  library(fastshap)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript 06_run_week8_shap_direction_exact.R TECHNOLOGY TASK ",
    "[--output-root=PATH] [--nsim=16] [--per-class=100] ",
    "[--background-per-class=100] [--nthread=8]"
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

output_root <- value_or(
  "output-root",
  "/home/xiyao/15data/W8/results/06_shap_direction_exact"
)
nsim <- as.integer(value_or("nsim", "16"))
per_class <- as.integer(value_or("per-class", "100"))
background_per_class <- as.integer(value_or("background-per-class", "100"))
n_thread <- as.integer(value_or("nthread", "8"))

if (!is.finite(nsim) || nsim < 2L) stop("Invalid --nsim.")
if (!is.finite(per_class) || per_class < 1L) stop("Invalid --per-class.")
if (!is.finite(background_per_class) || background_per_class < 1L) {
  stop("Invalid --background-per-class.")
}
if (!is.finite(n_thread) || n_thread < 1L) stop("Invalid --nthread.")

allowed_tasks <- c("calibration_tp_fp", "error_fp_tn")
if (!task %in% allowed_tasks) {
  stop("This exact SHAP analysis is restricted to calibration_tp_fp/error_fp_tn.")
}

task_definitions <- list(
  calibration_tp_fp = c(negative = 2L, positive = 3L),
  error_fp_tn = c(negative = 0L, positive = 2L)
)
positive_class <- unname(c(
  calibration_tp_fp = "TP",
  error_fp_tn = "FP"
)[[task]])
negative_class <- unname(c(
  calibration_tp_fp = "FP",
  error_fp_tn = "TN"
)[[task]])

rf_seed <- 20260812L
shap_seed_base <- 20260813L
num_trees_expected <- 300L
max_per_class_expected <- 100000L
prediction_tolerance <- 1e-12

class_path <- paste0(
  "/home/xiyao/15data/W8/results/01_four_class_global_fdr/",
  "four_class_matrix_Total_gt_10.rds"
)
dataset_path <- paste0(
  "/home/xichen/ConquerNet/ConquerBN_basic/",
  "m6AConquer_hg38_BenchmarkDataSet_Combined_SE_",
  "no_cali_IVT_matched.rds"
)
feature_path <- paste0(
  "/home/xiyao/15data/W8/input_data/features/",
  "m6Aconquer_omicsFeatures_hg38.rds"
)
design_path <- paste0(
  "/home/xiyao/15data/W8/input_data/rf_design/",
  "tp_fp_tp_tn_metadata.rds"
)
rf_root <- "/home/xiyao/15data/W8/results/03_rf_four_class_global_fdr"
top20_root <- paste0(
  "/home/xiyao/15data/W8/results/05_rf_feature_exports/",
  "per_model_top20"
)

original_dir <- file.path(rf_root, task, technology)
top20_path <- file.path(
  top20_root,
  paste0(task, "__", technology, "__top20.csv")
)
required_paths <- c(
  class_path, dataset_path, feature_path, design_path,
  file.path(original_dir, "status.csv"),
  file.path(original_dir, "oof_predictions.rds"),
  file.path(original_dir, "fold_audit.csv"),
  top20_path
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths)) {
  stop("Missing required inputs: ", paste(missing_paths, collapse = "; "))
}

class_matrix <- readRDS(class_path)
dataset <- readRDS(dataset_path)
features <- readRDS(feature_path)
design <- readRDS(design_path)
status <- read.csv(file.path(original_dir, "status.csv"), check.names = FALSE)
original_oof <- readRDS(file.path(original_dir, "oof_predictions.rds"))
original_fold_audit <- read.csv(
  file.path(original_dir, "fold_audit.csv"), check.names = FALSE
)
top20 <- read.csv(top20_path, check.names = FALSE)

if (nrow(status) != 1L || status$run_status[[1L]] != "complete") {
  stop("Original RF status is not complete.")
}
if (status$identifiability_status[[1L]] != "reliable") {
  stop("SHAP is restricted to reliable models.")
}
if (status$label_source[[1L]] !=
    "direct_joint_four_class_posterior_global_FDR_lt_0.05") {
  stop("Unexpected label source.")
}
if (isTRUE(status$map_viterbi_fallback[[1L]])) stop("Fallback is not allowed.")
if (status$coverage_cutoff[[1L]] != 10L) stop("Expected Total>10.")
if (status$num_trees[[1L]] != num_trees_expected) stop("Unexpected tree count.")
if (status$max_per_class[[1L]] != max_per_class_expected) {
  stop("Unexpected training class cap.")
}
if (status$n_thread[[1L]] != n_thread) {
  stop("nthread differs from the original RF run.")
}
if (!inherits(dataset, "SummarizedExperiment")) stop("Invalid dataset.")
if (!"Total" %in% assayNames(dataset)) stop("Dataset lacks Total.")
if (!identical(dim(class_matrix), c(nrow(dataset), ncol(dataset)))) {
  stop("Class-matrix dimensions differ from dataset.")
}
if (!identical(colnames(class_matrix), colnames(dataset))) {
  stop("Technology order differs between class matrix and dataset.")
}
if (length(features) != nrow(dataset)) stop("Feature row-count mismatch.")
if (!technology %in% colnames(dataset)) stop("Unknown technology.")

dataset_ranges <- rowRanges(dataset)
coordinates_identical <- (
  identical(as.character(seqnames(features)), as.character(seqnames(dataset_ranges))) &&
    identical(start(features), start(dataset_ranges)) &&
    identical(end(features), end(dataset_ranges)) &&
    identical(as.character(strand(features)), as.character(strand(dataset_ranges)))
)
if (!coordinates_identical) stop("Feature/dataset coordinates are misaligned.")

screening_features <- as.character(design[["screening_features"]])
chromosome_fold <- design[["chromosome_fold"]]
if (length(screening_features) != 185L) stop("Expected exactly 185 features.")
if (!all(screening_features %in% names(mcols(features)))) {
  stop("At least one frozen feature is missing.")
}
if (is.null(names(chromosome_fold)) || !all(chromosome_fold %in% 1:5)) {
  stop("Invalid frozen chromosome folds.")
}

top20 <- top20[
  top20$technology == technology & top20$task == task,
  , drop = FALSE
]
top20 <- top20[order(top20$top20_rank, top20$feature), , drop = FALSE]
if (nrow(top20) < 10L || !identical(top20$top20_rank[1:10], 1:10)) {
  stop("Top20 export does not contain ranks 1:10.")
}
top_features <- as.character(top20$feature[1:10])
if (anyDuplicated(top_features) || !all(top_features %in% screening_features)) {
  stop("Invalid model Top10 feature list.")
}

technology_column <- match(technology, colnames(dataset))
total <- assay(dataset, "Total")[, technology_column]
observed <- is.finite(total) & total > 10
class_vector <- as.integer(class_matrix[, technology_column])
if (any(!is.na(class_vector) & !observed)) {
  stop("Total>10 class matrix contains labels outside Total>10.")
}
codes <- task_definitions[[task]]
eligible <- which(observed & class_vector %in% unname(codes))
y_all <- as.integer(class_vector[eligible] == codes[["positive"]])
chromosome <- as.character(seqnames(features))[eligible]
fold_id <- unname(chromosome_fold[chromosome])
if (anyNA(fold_id) || !all(fold_id %in% 1:5)) stop("Fold assignment failed.")

if (!identical(as.integer(original_oof$row_index), as.integer(eligible))) {
  stop("Original OOF row order differs from reconstructed eligible rows.")
}
if (!identical(as.integer(original_oof$fold), as.integer(fold_id))) {
  stop("Original OOF folds differ from frozen design.")
}
if (!identical(as.integer(original_oof$y), as.integer(y_all))) {
  stop("Original OOF labels differ from reconstructed labels.")
}
if (!identical(as.integer(original_oof$class_code), class_vector[eligible])) {
  stop("Original OOF class codes differ from current labels.")
}
if (nrow(original_fold_audit) != 5L ||
    !identical(as.integer(original_fold_audit$fold), 1:5)) {
  stop("Original fold audit is incomplete.")
}

safe_logit <- function(p) qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))
prior_correct <- function(p, target_prevalence, sample_prevalence) {
  shift <- qlogis(target_prevalence) - qlogis(sample_prevalence)
  plogis(safe_logit(p) + shift)
}

as_numeric_frame <- function(row_positions, columns) {
  answer <- as.data.frame(
    mcols(features)[eligible[row_positions], columns, drop = FALSE],
    check.names = FALSE
  )
  for (name in columns) answer[[name]] <- as.numeric(answer[[name]])
  answer
}

sample_at_most <- function(index, n) {
  if (length(index) > n) sample(index, n) else index
}

direction_row <- function(feature_value, shap_value, feature_name, fold) {
  valid <- is.finite(feature_value) & is.finite(shap_value)
  x <- feature_value[valid]
  s <- shap_value[valid]
  unique_x <- sort(unique(x))
  is_binary <- length(unique_x) <= 2L && all(unique_x %in% c(0, 1))
  if (!length(x) || length(unique_x) < 2L) {
    direction_score <- NA_real_
    score_definition <- if (is_binary) {
      "mean_SHAP(x=1)-mean_SHAP(x=0)"
    } else {
      "Spearman(feature_value,SHAP)"
    }
  } else if (is_binary) {
    direction_score <- mean(s[x == 1]) - mean(s[x == 0])
    score_definition <- "mean_SHAP(x=1)-mean_SHAP(x=0)"
  } else {
    direction_score <- suppressWarnings(cor(x, s, method = "spearman"))
    score_definition <- "Spearman(feature_value,SHAP)"
  }
  direction_sign <- if (!is.finite(direction_score) || direction_score == 0) {
    0L
  } else {
    as.integer(sign(direction_score))
  }
  data.frame(
    technology = technology,
    strategy = as.character(status$strategy[[1L]]),
    task = task,
    positive_class = positive_class,
    negative_class = negative_class,
    fold = fold,
    feature = feature_name,
    top10_rank = match(feature_name, top_features),
    feature_type = if (is_binary) "binary_0_1" else "continuous_or_ordinal",
    n_explained = length(x),
    n_unique_values = length(unique_x),
    mean_shap = if (length(s)) mean(s) else NA_real_,
    mean_abs_shap = if (length(s)) mean(abs(s)) else NA_real_,
    median_abs_shap = if (length(s)) median(abs(s)) else NA_real_,
    direction_score = direction_score,
    direction_score_definition = score_definition,
    direction_sign = direction_sign,
    direction_label = if (direction_sign > 0L) {
      paste0("higher_value_pushes_", positive_class)
    } else if (direction_sign < 0L) {
      paste0("higher_value_pushes_", negative_class)
    } else {
      "neutral_or_not_estimable"
    },
    stringsAsFactors = FALSE
  )
}

## SHAP is deliberately run inside a function that restores the exact RF RNG
## stream on exit. Consequently the next fold's class-capped training sample is
## identical to the original sequential five-fold run.
run_fold_shap <- function(
    model, model_features, train_pool, validation_pool,
    train_prevalence, sample_prevalence, fold) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  technology_hash <- sum(utf8ToInt(technology))
  task_index <- match(task, allowed_tasks)
  set.seed(shap_seed_base + 10000L * task_index + 100L * fold + technology_hash)

  validation0 <- validation_pool[y_all[validation_pool] == 0L]
  validation1 <- validation_pool[y_all[validation_pool] == 1L]
  explain_positions <- sort(c(
    sample_at_most(validation0, per_class),
    sample_at_most(validation1, per_class)
  ))
  train0 <- train_pool[y_all[train_pool] == 0L]
  train1 <- train_pool[y_all[train_pool] == 1L]
  background_positions <- sort(c(
    sample_at_most(train0, background_per_class),
    sample_at_most(train1, background_per_class)
  ))

  explain_x <- as_numeric_frame(explain_positions, screening_features)
  background_x <- as_numeric_frame(background_positions, screening_features)
  if (anyNA(explain_x) || anyNA(background_x)) {
    stop("Missing SHAP predictor value in fold ", fold)
  }
  explain_features <- intersect(top_features, model_features)
  if (!length(explain_features)) stop("No Top10 feature entered fold ", fold)

  pred_wrapper <- function(object, newdata) {
    model_data <- newdata[, model_features, drop = FALSE]
    raw <- predict(
      object, data = model_data, num.threads = n_thread
    )$predictions[, "1"]
    prior_correct(raw, train_prevalence, sample_prevalence)
  }

  shap_matrix <- fastshap::explain(
    object = model,
    feature_names = explain_features,
    X = background_x,
    nsim = nsim,
    pred_wrapper = pred_wrapper,
    newdata = explain_x,
    adjust = FALSE,
    shap_only = TRUE,
    parallel = FALSE
  )
  shap_matrix <- as.matrix(shap_matrix)
  if (!identical(colnames(shap_matrix), explain_features)) {
    stop("fastshap returned an unexpected feature set/order in fold ", fold)
  }
  if (any(!is.finite(shap_matrix))) stop("Non-finite SHAP value in fold ", fold)

  sampled_metadata <- data.frame(
    technology = technology,
    strategy = as.character(status$strategy[[1L]]),
    task = task,
    positive_class = positive_class,
    negative_class = negative_class,
    fold = fold,
    row_index = eligible[explain_positions],
    y = y_all[explain_positions],
    class_code = class_vector[eligible[explain_positions]],
    prior_corrected_probability = pred_wrapper(model, explain_x),
    stringsAsFactors = FALSE
  )
  fold_direction <- do.call(rbind, lapply(colnames(shap_matrix), function(name) {
    direction_row(explain_x[[name]], shap_matrix[, name], name, fold)
  }))
  list(
    sample = list(
      metadata = sampled_metadata,
      feature_values = as.matrix(explain_x[, colnames(shap_matrix), drop = FALSE]),
      shap_values = shap_matrix,
      background_row_index = eligible[background_positions]
    ),
    fold_direction = fold_direction,
    audit = data.frame(
      technology = technology,
      strategy = as.character(status$strategy[[1L]]),
      task = task,
      fold = fold,
      n_explained = length(explain_positions),
      n_explained_negative = sum(y_all[explain_positions] == 0L),
      n_explained_positive = sum(y_all[explain_positions] == 1L),
      n_background = length(background_positions),
      n_background_negative = sum(y_all[background_positions] == 0L),
      n_background_positive = sum(y_all[background_positions] == 1L),
      n_top10_requested = length(top_features),
      n_top10_entered_and_explained = ncol(shap_matrix),
      omitted_constant_top10 = paste(setdiff(top_features, colnames(shap_matrix)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

result_dir <- file.path(output_root, task, technology)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(result_dir, "COMPLETE.ok"))

set.seed(rf_seed)
sample_rows <- vector("list", 5L)
direction_rows <- vector("list", 5L)
qc_rows <- vector("list", 5L)

for (fold in 1:5) {
  message("Exact RF + SHAP: ", technology, " / ", task, " / fold ", fold)
  train_pool <- which(fold_id != fold)
  validation_pool <- which(fold_id == fold)
  train0 <- train_pool[y_all[train_pool] == 0L]
  train1 <- train_pool[y_all[train_pool] == 1L]
  sampled0 <- if (length(train0) > max_per_class_expected) {
    sample(train0, max_per_class_expected)
  } else train0
  sampled1 <- if (length(train1) > max_per_class_expected) {
    sample(train1, max_per_class_expected)
  } else train1
  train_index <- sort(c(sampled0, sampled1))

  train_prevalence <- mean(y_all[train_pool])
  sample_prevalence <- mean(y_all[train_index])
  train_x <- as_numeric_frame(train_index, screening_features)
  validation_x <- as_numeric_frame(validation_pool, screening_features)
  keep_feature <- vapply(train_x, function(x) {
    finite <- x[is.finite(x)]
    length(finite) > 1L && length(unique(finite)) > 1L
  }, logical(1))
  if (!any(keep_feature)) stop("All features constant in fold ", fold)
  train_x <- train_x[, keep_feature, drop = FALSE]
  validation_x <- validation_x[, keep_feature, drop = FALSE]
  if (anyNA(train_x) || anyNA(validation_x)) {
    stop("Missing predictor after conversion in fold ", fold)
  }

  original_audit_row <- original_fold_audit[original_fold_audit$fold == fold, ]
  expected_counts <- c(
    length(train0), length(train1), length(sampled0), length(sampled1),
    sum(y_all[validation_pool] == 0L), sum(y_all[validation_pool] == 1L),
    ncol(train_x)
  )
  recorded_counts <- as.integer(c(
    original_audit_row$n_train_pool_negative,
    original_audit_row$n_train_pool_positive,
    original_audit_row$n_train_sample_negative,
    original_audit_row$n_train_sample_positive,
    original_audit_row$n_validation_negative,
    original_audit_row$n_validation_positive,
    original_audit_row$n_features_entered
  ))
  if (!identical(as.integer(expected_counts), recorded_counts)) {
    stop("Fold-audit reconstruction mismatch in fold ", fold)
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
    num.trees = num_trees_expected,
    mtry = max(1L, floor(sqrt(ncol(train_x)))),
    min.node.size = 20L,
    importance = "permutation",
    num.threads = n_thread,
    seed = rf_seed + fold,
    write.forest = TRUE
  )

  raw_probability <- predict(
    model, data = validation_x, num.threads = n_thread
  )$predictions[, "1"]
  corrected_probability <- prior_correct(
    raw_probability, train_prevalence, sample_prevalence
  )
  original_position <- match(
    eligible[validation_pool], original_oof$row_index
  )
  if (anyNA(original_position)) stop("OOF row match failed in fold ", fold)
  raw_error <- max(abs(
    raw_probability - original_oof$raw_probability[original_position]
  ))
  corrected_error <- max(abs(
    corrected_probability -
      original_oof$prior_corrected_probability[original_position]
  ))
  if (!is.finite(raw_error) || !is.finite(corrected_error) ||
      raw_error > prediction_tolerance ||
      corrected_error > prediction_tolerance) {
    stop(
      "Exact OOF reproduction failed in fold ", fold,
      ": raw=", format(raw_error, digits = 17),
      ", corrected=", format(corrected_error, digits = 17)
    )
  }

  shap_result <- run_fold_shap(
    model = model,
    model_features = colnames(train_x),
    train_pool = train_pool,
    validation_pool = validation_pool,
    train_prevalence = train_prevalence,
    sample_prevalence = sample_prevalence,
    fold = fold
  )
  shap_result$audit$n_validation_full <- length(validation_pool)
  shap_result$audit$n_train_pool <- length(train_pool)
  shap_result$audit$n_train_sample <- length(train_index)
  shap_result$audit$n_model_features <- ncol(train_x)
  shap_result$audit$raw_probability_max_abs_error <- raw_error
  shap_result$audit$prior_corrected_probability_max_abs_error <- corrected_error
  shap_result$audit$prediction_tolerance <- prediction_tolerance
  shap_result$audit$exact_reproduction_pass <- TRUE
  shap_result$audit$importance_mode <- "permutation"
  shap_result$audit$rf_seed <- rf_seed + fold
  shap_result$audit$shap_rng_restored_before_next_fold <- TRUE
  sample_rows[[fold]] <- shap_result$sample
  direction_rows[[fold]] <- shap_result$fold_direction
  qc_rows[[fold]] <- shap_result$audit

  rm(
    train_x, validation_x, train_data, model, raw_probability,
    corrected_probability, shap_result
  )
  gc()
}

fold_direction <- do.call(rbind, direction_rows)
qc <- do.call(rbind, qc_rows)
if (nrow(qc) != 5L || !all(qc$exact_reproduction_pass)) {
  stop("Five-fold exact reproduction audit is incomplete.")
}

summary_rows <- lapply(top_features, function(name) {
  x <- fold_direction[fold_direction$feature == name, , drop = FALSE]
  n_positive <- sum(x$direction_sign > 0L, na.rm = TRUE)
  n_negative <- sum(x$direction_sign < 0L, na.rm = TRUE)
  stable_sign <- if (n_positive >= 4L) 1L else if (n_negative >= 4L) -1L else 0L
  data.frame(
    technology = technology,
    strategy = as.character(status$strategy[[1L]]),
    task = task,
    positive_class = positive_class,
    negative_class = negative_class,
    feature = name,
    top10_rank = match(name, top_features),
    n_folds_requested = 5L,
    n_folds_estimable = sum(x$direction_sign != 0L & is.finite(x$direction_score)),
    n_positive_direction_folds = n_positive,
    n_negative_direction_folds = n_negative,
    n_neutral_or_missing_folds = 5L - n_positive - n_negative,
    stable_direction_4of5 = stable_sign != 0L,
    stable_direction_sign = stable_sign,
    stable_direction_label = if (stable_sign > 0L) {
      paste0("higher_value_pushes_", positive_class)
    } else if (stable_sign < 0L) {
      paste0("higher_value_pushes_", negative_class)
    } else {
      "not_stable_4of5"
    },
    median_direction_score = if (any(is.finite(x$direction_score))) {
      median(x$direction_score[is.finite(x$direction_score)])
    } else NA_real_,
    mean_abs_shap = if (any(is.finite(x$mean_abs_shap))) {
      mean(x$mean_abs_shap[is.finite(x$mean_abs_shap)])
    } else NA_real_,
    stringsAsFactors = FALSE
  )
})
direction_summary <- do.call(rbind, summary_rows)

write.csv(
  fold_direction,
  file.path(result_dir, "shap_fold_directions.csv"),
  row.names = FALSE
)
write.csv(
  direction_summary,
  file.path(result_dir, "shap_direction_summary.csv"),
  row.names = FALSE
)
write.csv(qc, file.path(result_dir, "retraining_qc.csv"), row.names = FALSE)
saveRDS(
  sample_rows,
  file.path(result_dir, "sampled_shap.rds"),
  compress = TRUE
)
write.csv(
  data.frame(
    technology = technology,
    task = task,
    positive_class = positive_class,
    negative_class = negative_class,
    label_source = status$label_source[[1L]],
    coverage_rule = "Total>10",
    map_viterbi_fallback = FALSE,
    scope = "complete_and_reliable_only",
    model_features = 185L,
    explained_top_n = 10L,
    top10_source = top20_path,
    nsim = nsim,
    validation_max_per_class_per_fold = per_class,
    background_source = "corresponding_training_pool",
    background_max_per_class_per_fold = background_per_class,
    background_contains_validation_rows = FALSE,
    shap_adjust = FALSE,
    original_rf_seed = rf_seed,
    independent_shap_seed_base = shap_seed_base,
    shap_rng_restored_before_next_fold = TRUE,
    exact_full_validation_qc_tolerance = prediction_tolerance,
    importance_mode = "permutation",
    forest_saved = FALSE,
    stringsAsFactors = FALSE
  ),
  file.path(result_dir, "protocol.csv"),
  row.names = FALSE
)
writeLines(
  c(
    paste("completed", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("technology", technology),
    paste("task", task),
    paste("max_raw_reproduction_error", max(qc$raw_probability_max_abs_error)),
    paste(
      "max_corrected_reproduction_error",
      max(qc$prior_corrected_probability_max_abs_error)
    )
  ),
  file.path(result_dir, "COMPLETE.ok")
)
message("Completed exact SHAP: ", technology, " / ", task)
