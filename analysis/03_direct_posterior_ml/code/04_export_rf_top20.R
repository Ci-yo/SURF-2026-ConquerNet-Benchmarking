#!/usr/bin/env Rscript

# Export clean Week8 RF feature-ranking products from the completed summary.
#
# Scope:
#   * only identifiability_status == "reliable"
#   * only run_status == "complete"
#   * no exploratory or untrainable task contributes to any export/shared rank
#
# Main products:
#   1. one Top-20 CSV per reliable completed technology x task
#   2. one combined Top-20 CSV
#   3. cross-technology stable-feature summary by task
#   4. Top-20 shared-feature summary for each task and a combined table
#   5. manifests and completeness audits

options(stringsAsFactors = FALSE)

parse_args <- function(x) {
  out <- list(
    summary_dir = "/home/xiyao/15data/W8/results/04_rf_four_class_summary",
    output_dir = "/home/xiyao/15data/W8/results/05_rf_feature_exports",
    overwrite = FALSE
  )
  for (arg in x) {
    if (grepl("^--summary-dir=", arg)) {
      out$summary_dir <- sub("^--summary-dir=", "", arg)
    } else if (grepl("^--output-dir=", arg)) {
      out$output_dir <- sub("^--output-dir=", "", arg)
    } else if (grepl("^--overwrite=", arg)) {
      value <- tolower(sub("^--overwrite=", "", arg))
      if (!value %in% c("true", "false", "1", "0", "yes", "no")) {
        stop("--overwrite must be true or false.")
      }
      out$overwrite <- value %in% c("true", "1", "yes")
    } else if (arg %in% c("--help", "-h")) {
      cat(
        "Usage: Rscript 04_export_rf_top20.R ",
        "[--summary-dir=DIR] [--output-dir=DIR] [--overwrite=true|false]\n",
        sep = ""
      )
      quit(status = 0L)
    } else {
      stop("Unknown argument: ", arg)
    }
  }
  out
}

assert_columns <- function(x, required, object_name) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(object_name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

as_logical_strict <- function(x, column_name) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %in% c("true", "t", "1")] <- TRUE
  out[y %in% c("false", "f", "0")] <- FALSE
  if (anyNA(out)) {
    stop("Unexpected value in ", column_name, ": ", paste(unique(x[is.na(out)]), collapse = ", "))
  }
  out
}

safe_median <- function(x) {
  x <- as.numeric(x)
  if (!any(is.finite(x))) return(NA_real_)
  median(x[is.finite(x)])
}

sanitize_file_part <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("_+", "_", x)
}

order_model_features <- function(x) {
  x[order(
    -x$top20_fold_fraction,
    x$median_rank,
    -x$median_normalized_importance,
    -x$positive_fold_fraction,
    -x$entered_fold_fraction,
    x$feature,
    na.last = TRUE
  ), , drop = FALSE]
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
summary_dir <- normalizePath(args$summary_dir, winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args$output_dir, winslash = "/", mustWork = FALSE)

status_path <- file.path(summary_dir, "all_task_status.csv")
importance_path <- file.path(summary_dir, "feature_importance_stability.csv")
if (!file.exists(status_path)) stop("Missing input: ", status_path)
if (!file.exists(importance_path)) stop("Missing input: ", importance_path)

status <- read.csv(status_path, check.names = FALSE)
importance <- read.csv(importance_path, check.names = FALSE)

assert_columns(
  status,
  c("technology", "strategy", "task", "identifiability_status", "run_status",
    "n_evaluable_folds", "n_model_negative", "n_model_positive"),
  "all_task_status.csv"
)
assert_columns(
  importance,
  c("technology", "strategy", "task", "feature", "entered_fold_fraction",
    "positive_fold_fraction", "top10_fold_fraction", "top20_fold_fraction",
    "median_rank", "median_importance", "median_normalized_importance",
    "stable_4of5_top20"),
  "feature_importance_stability.csv"
)

importance$stable_4of5_top20 <- as_logical_strict(
  importance$stable_4of5_top20,
  "stable_4of5_top20"
)

key <- function(technology, task) paste(technology, task, sep = "\r")

eligible <- status[
  status$identifiability_status == "reliable" & status$run_status == "complete",
  , drop = FALSE
]
if (nrow(eligible) == 0L) stop("No reliable completed models were found.")
if (anyDuplicated(key(eligible$technology, eligible$task))) {
  stop("all_task_status.csv has duplicate reliable completed technology-task rows.")
}
if (any(!is.finite(eligible$n_evaluable_folds) | eligible$n_evaluable_folds != 5L)) {
  stop("Every reliable completed model must have exactly five evaluable folds.")
}

eligible_keys <- key(eligible$technology, eligible$task)
importance_keys <- unique(key(importance$technology, importance$task))
if (!setequal(eligible_keys, importance_keys)) {
  stop(
    "Importance technology-task keys do not exactly match reliable completed status keys. ",
    "Missing importance: ", paste(setdiff(eligible_keys, importance_keys), collapse = "; "),
    ". Unexpected importance: ", paste(setdiff(importance_keys, eligible_keys), collapse = "; ")
  )
}
if (anyDuplicated(key(key(importance$technology, importance$task), importance$feature))) {
  stop("feature_importance_stability.csv contains duplicate technology-task-feature rows.")
}

task_order <- c("calibration_tp_fp", "error_fp_tn", "artifact_tp_tn")
unknown_tasks <- setdiff(unique(eligible$task), task_order)
task_order <- c(task_order[task_order %in% unique(eligible$task)], sort(unknown_tasks))

if (dir.exists(output_dir)) {
  existing <- list.files(output_dir, all.files = TRUE, no.. = TRUE)
  if (length(existing) > 0L && !args$overwrite) {
    stop("Output directory is not empty. Re-run with --overwrite=true: ", output_dir)
  }
  if (length(existing) > 0L && args$overwrite) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
per_model_dir <- file.path(output_dir, "per_model_top20")
per_task_dir <- file.path(output_dir, "per_task_shared_top20")
dir.create(per_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(per_task_dir, recursive = TRUE, showWarnings = FALSE)

# Export each reliable completed technology-task model's Top 20.
eligible <- eligible[order(match(eligible$task, task_order), eligible$technology), , drop = FALSE]
top20_list <- vector("list", nrow(eligible))
manifest_list <- vector("list", nrow(eligible))

for (j in seq_len(nrow(eligible))) {
  technology <- eligible$technology[j]
  task <- eligible$task[j]
  z <- importance[importance$technology == technology & importance$task == task, , drop = FALSE]
  if (nrow(z) < 20L) stop("Fewer than 20 feature rows for ", technology, " / ", task)
  z <- order_model_features(z)
  z <- z[seq_len(20L), , drop = FALSE]
  z$top20_rank <- seq_len(nrow(z))
  z$identifiability_status <- "reliable"
  z$run_status <- "complete"
  z$n_evaluable_folds <- eligible$n_evaluable_folds[j]
  z$n_model_negative <- eligible$n_model_negative[j]
  z$n_model_positive <- eligible$n_model_positive[j]
  z <- z[, c(
    "technology", "strategy", "task", "identifiability_status", "run_status",
    "n_evaluable_folds", "n_model_negative", "n_model_positive", "top20_rank",
    "feature", "entered_fold_fraction", "positive_fold_fraction",
    "top10_fold_fraction", "top20_fold_fraction", "median_rank",
    "median_importance", "median_normalized_importance", "stable_4of5_top20"
  )]

  filename <- paste0(
    sanitize_file_part(task), "__",
    sanitize_file_part(technology), "__top20.csv"
  )
  write.csv(z, file.path(per_model_dir, filename), row.names = FALSE, na = "")
  top20_list[[j]] <- z
  manifest_list[[j]] <- data.frame(
    technology = technology,
    strategy = eligible$strategy[j],
    task = task,
    identifiability_status = "reliable",
    run_status = "complete",
    n_available_features = nrow(importance[
      importance$technology == technology & importance$task == task,
      , drop = FALSE
    ]),
    n_exported_top20 = nrow(z),
    n_stable_in_exported_top20 = sum(z$stable_4of5_top20),
    output_file = file.path("per_model_top20", filename),
    stringsAsFactors = FALSE
  )
}

all_top20 <- do.call(rbind, top20_list)
model_manifest <- do.call(rbind, manifest_list)
write.csv(
  all_top20,
  file.path(output_dir, "all_reliable_models_top20.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  model_manifest,
  file.path(output_dir, "reliable_complete_model_manifest.csv"),
  row.names = FALSE,
  na = ""
)

# Build cross-technology summaries using only features that are stable within a
# reliable completed technology-task model (4/5 positive and 4/5 Top-20).
stable <- importance[
  key(importance$technology, importance$task) %in% eligible_keys &
    importance$stable_4of5_top20,
  , drop = FALSE
]
if (nrow(stable) == 0L) stop("No stable_4of5_top20 features were found.")

shared_rows <- list()
row_counter <- 0L
for (task in task_order) {
  task_technologies <- sort(unique(eligible$technology[eligible$task == task]))
  task_stable <- stable[stable$task == task, , drop = FALSE]
  for (feature in sort(unique(task_stable$feature))) {
    z <- task_stable[task_stable$feature == feature, , drop = FALSE]
    row_counter <- row_counter + 1L
    shared_rows[[row_counter]] <- data.frame(
      task = task,
      feature = feature,
      n_reliable_technologies = length(unique(z$technology)),
      n_reliable_technologies_total = length(task_technologies),
      reliable_technology_fraction = length(unique(z$technology)) / length(task_technologies),
      median_rank = safe_median(z$median_rank),
      median_normalized_importance = safe_median(z$median_normalized_importance),
      technologies = paste(sort(unique(z$technology)), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}

shared <- do.call(rbind, shared_rows)
shared <- shared[order(
  match(shared$task, task_order),
  -shared$n_reliable_technologies,
  -shared$reliable_technology_fraction,
  shared$median_rank,
  -shared$median_normalized_importance,
  shared$feature
), , drop = FALSE]
rownames(shared) <- NULL
write.csv(
  shared,
  file.path(output_dir, "cross_technology_stable_feature_summary.csv"),
  row.names = FALSE,
  na = ""
)

# Export the first 20 cross-technology shared stable features for every task.
shared_top20_list <- vector("list", length(task_order))
names(shared_top20_list) <- task_order
for (task in task_order) {
  z <- shared[shared$task == task, , drop = FALSE]
  z <- z[seq_len(min(20L, nrow(z))), , drop = FALSE]
  z$shared_rank <- seq_len(nrow(z))
  z <- z[, c(
    "task", "shared_rank", "feature", "n_reliable_technologies",
    "n_reliable_technologies_total", "reliable_technology_fraction",
    "median_rank", "median_normalized_importance", "technologies"
  )]
  filename <- paste0(sanitize_file_part(task), "__shared_top20.csv")
  write.csv(z, file.path(per_task_dir, filename), row.names = FALSE, na = "")
  shared_top20_list[[task]] <- z
}
all_shared_top20 <- do.call(rbind, shared_top20_list)
rownames(all_shared_top20) <- NULL
write.csv(
  all_shared_top20,
  file.path(output_dir, "all_tasks_shared_top20.csv"),
  row.names = FALSE,
  na = ""
)

# Completeness and exclusion audits.
expected_models_by_task <- as.data.frame(table(eligible$task), stringsAsFactors = FALSE)
names(expected_models_by_task) <- c("task", "n_reliable_complete_models")
expected_models_by_task <- expected_models_by_task[
  order(match(expected_models_by_task$task, task_order)), , drop = FALSE
]

audit <- data.frame(
  check = c(
    "reliable_complete_models",
    "per_model_top20_files",
    "all_top20_rows",
    "all_top20_unique_models",
    "all_top20_each_model_has_20_rows",
    "all_top20_ranks_are_1_to_20",
    "all_top20_only_reliable_complete",
    "importance_keys_equal_reliable_complete_keys",
    "shared_summary_only_uses_reliable_complete",
    "shared_top20_tasks",
    "shared_top20_rows"
  ),
  expected = c(
    nrow(eligible),
    nrow(eligible),
    20L * nrow(eligible),
    nrow(eligible),
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    length(task_order),
    sum(pmin(20L, as.integer(table(shared$task)[task_order])))
  ),
  observed = c(
    nrow(eligible),
    length(list.files(per_model_dir, pattern = "\\.csv$", full.names = TRUE)),
    nrow(all_top20),
    length(unique(key(all_top20$technology, all_top20$task))),
    all(table(key(all_top20$technology, all_top20$task)) == 20L),
    all(vapply(
      split(all_top20$top20_rank, key(all_top20$technology, all_top20$task)),
      function(x) identical(sort(as.integer(x)), seq_len(20L)),
      logical(1)
    )),
    all(all_top20$identifiability_status == "reliable" & all_top20$run_status == "complete"),
    setequal(eligible_keys, importance_keys),
    all(key(stable$technology, stable$task) %in% eligible_keys),
    length(unique(all_shared_top20$task)),
    nrow(all_shared_top20)
  ),
  stringsAsFactors = FALSE
)
audit$pass <- as.character(audit$expected) == as.character(audit$observed)
write.csv(audit, file.path(output_dir, "export_completeness_audit.csv"), row.names = FALSE)

task_audit <- merge(
  expected_models_by_task,
  aggregate(
    technology ~ task,
    data = stable,
    FUN = function(x) length(unique(x))
  ),
  by = "task",
  all.x = TRUE,
  sort = FALSE
)
names(task_audit)[names(task_audit) == "technology"] <- "n_technologies_with_any_stable_feature"
task_audit$n_cross_technology_stable_features <- vapply(
  task_audit$task,
  function(task) sum(shared$task == task),
  integer(1)
)
task_audit$n_shared_top20_exported <- vapply(
  task_audit$task,
  function(task) sum(all_shared_top20$task == task),
  integer(1)
)
task_audit$all_reliable_models_represented <-
  task_audit$n_reliable_complete_models == task_audit$n_technologies_with_any_stable_feature
task_audit <- task_audit[order(match(task_audit$task, task_order)), , drop = FALSE]
write.csv(task_audit, file.path(output_dir, "task_shared_summary_audit.csv"), row.names = FALSE)

excluded <- status[
  !(status$identifiability_status == "reliable" & status$run_status == "complete"),
  c("technology", "strategy", "task", "identifiability_status", "run_status",
    "n_model_negative", "n_model_positive"),
  drop = FALSE
]
excluded$included_in_top20_or_shared <- FALSE
excluded <- excluded[order(match(excluded$task, task_order), excluded$technology), , drop = FALSE]
write.csv(excluded, file.path(output_dir, "excluded_task_audit.csv"), row.names = FALSE, na = "")

readme <- c(
  "Week8 RF feature exports",
  "==========================",
  "",
  paste0("Source summary: ", summary_dir),
  "Eligibility: identifiability_status == reliable AND run_status == complete.",
  "Exploratory and untrainable tasks are excluded from every Top20/shared summary.",
  "",
  "Per-model ranking order:",
  "  1) top20_fold_fraction descending",
  "  2) median_rank ascending",
  "  3) median_normalized_importance descending",
  "  4) positive_fold_fraction descending",
  "  5) entered_fold_fraction descending",
  "  6) feature name",
  "",
  "A stable shared occurrence requires stable_4of5_top20 == TRUE within a model:",
  "at least 4/5 folds with positive importance and at least 4/5 folds in Top20.",
  "Shared denominators are the reliable completed technologies within each task.",
  "These screening Top20 lists are not yet the final de-correlated reduced feature sets.",
  "",
  paste0("Reliable completed models exported: ", nrow(eligible)),
  paste0("Combined per-model Top20 rows: ", nrow(all_top20)),
  paste0("Tasks: ", paste(task_order, collapse = ", "))
)
writeLines(readme, file.path(output_dir, "README.txt"))

if (!all(audit$pass)) {
  print(audit)
  stop("At least one export completeness check failed.")
}
if (!all(task_audit$all_reliable_models_represented)) {
  print(task_audit)
  stop("At least one reliable completed model has no stable feature representation.")
}

cat("Week8 feature export completed successfully.\n")
cat("Output:", output_dir, "\n")
cat("Reliable completed models:", nrow(eligible), "\n")
print(expected_models_by_task, row.names = FALSE)
cat("Per-model Top20 rows:", nrow(all_top20), "\n")
cat("Cross-technology stable rows:", nrow(shared), "\n")
cat("Shared Top20 rows:", nrow(all_shared_top20), "\n")
cat("All completeness checks passed.\n")
