options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
result_root <- if (length(args) >= 1L) args[[1L]] else {
  "/home/xiyao/15data/W8/results/03_rf_four_class_global_fdr"
}
summary_dir <- if (length(args) >= 2L) args[[2L]] else {
  "/home/xiyao/15data/W8/results/04_rf_four_class_summary"
}
if (!dir.exists(result_root)) stop("Missing result root: ", result_root)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

files_named <- function(name) {
  list.files(result_root, pattern = paste0("^", name, "$"), recursive = TRUE,
             full.names = TRUE)
}
read_csv_rows <- function(paths) {
  if (!length(paths)) return(NULL)
  tables <- lapply(paths, read.csv, check.names = FALSE)
  ## Completed, skipped-exploratory, and single-class status files have
  ## intentionally different optional fields.  Bind by the union of names
  ## instead of requiring identical column counts/order.
  all_names <- unique(unlist(lapply(tables, names), use.names = FALSE))
  tables <- lapply(tables, function(x) {
    missing_names <- setdiff(all_names, names(x))
    for (name in missing_names) x[[name]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, tables)
}

status <- read_csv_rows(files_named("status.csv"))
if (is.null(status) || !nrow(status)) stop("No status.csv files found.")
write.csv(status, file.path(summary_dir, "all_task_status.csv"), row.names = FALSE)

metrics <- read_csv_rows(files_named("fold_metrics.csv"))
fold_audit <- read_csv_rows(files_named("fold_audit.csv"))
importance <- read_csv_rows(files_named("permutation_importance.csv"))
if (!is.null(metrics)) {
  write.csv(metrics, file.path(summary_dir, "all_fold_metrics.csv"), row.names = FALSE)
}
if (!is.null(fold_audit)) {
  write.csv(
    fold_audit, file.path(summary_dir, "all_fold_training_audit.csv"),
    row.names = FALSE
  )
}

complete_status <- status[status$run_status == "complete", , drop = FALSE]
oof_rows <- vector("list", nrow(complete_status))
for (i in seq_len(nrow(complete_status))) {
  row <- complete_status[i, , drop = FALSE]
  path <- file.path(result_root, row$task, row$technology, "oof_predictions.rds")
  ok <- file.exists(path)
  n <- n_unique <- n_missing <- n_bad_probability <- NA_integer_
  if (ok) {
    x <- readRDS(path)
    n <- nrow(x)
    n_unique <- length(unique(x$row_index))
    n_missing <- sum(
      !is.finite(x$raw_probability) |
        !is.finite(x$prior_corrected_probability)
    )
    n_bad_probability <- sum(
      x$raw_probability < 0 | x$raw_probability > 1 |
        x$prior_corrected_probability < 0 |
        x$prior_corrected_probability > 1,
      na.rm = TRUE
    )
    ok <- (
      n == row$n_model_domain && n_unique == n && n_missing == 0L &&
        n_bad_probability == 0L && all(x$fold %in% 1:5)
    )
  }
  oof_rows[[i]] <- data.frame(
    technology = row$technology, task = row$task, n_oof = n,
    n_unique_row_index = n_unique, n_missing_probability = n_missing,
    n_bad_probability = n_bad_probability, passed = ok,
    stringsAsFactors = FALSE
  )
}
oof_audit <- if (length(oof_rows)) do.call(rbind, oof_rows) else data.frame()
write.csv(oof_audit, file.path(summary_dir, "oof_integrity_audit.csv"), row.names = FALSE)
if (nrow(oof_audit) && !all(oof_audit$passed)) stop("OOF audit failed.")

if (!is.null(metrics) && nrow(metrics)) {
  keys <- unique(metrics[c("technology", "strategy", "task")])
  performance_rows <- lapply(seq_len(nrow(keys)), function(i) {
    z <- metrics[
      metrics$technology == keys$technology[i] &
        metrics$task == keys$task[i], , drop = FALSE
    ]
    decision <- status[
      status$technology == keys$technology[i] & status$task == keys$task[i],
      , drop = FALSE
    ]
    data.frame(
      technology = keys$technology[i], strategy = keys$strategy[i],
      task = keys$task[i],
      identifiability_status = decision$identifiability_status[1L],
      n_negative = decision$n_negative_full[1L],
      n_positive = decision$n_positive_full[1L],
      n_evaluable_folds = sum(is.finite(z$auroc) & is.finite(z$auprc)),
      mean_auroc = mean(z$auroc[is.finite(z$auroc)]),
      mean_auprc = mean(z$auprc[is.finite(z$auprc)]),
      stringsAsFactors = FALSE
    )
  })
  performance <- do.call(rbind, performance_rows)
  write.csv(
    performance, file.path(summary_dir, "per_technology_performance.csv"),
    row.names = FALSE
  )
  reliable <- performance[
    performance$identifiability_status == "reliable", , drop = FALSE
  ]
  if (nrow(reliable)) {
    macro <- do.call(rbind, lapply(split(reliable, reliable$task), function(z) {
      data.frame(
        task = z$task[1L], n_technologies = nrow(z),
        macro_auroc = mean(z$mean_auroc, na.rm = TRUE),
        macro_auprc = mean(z$mean_auprc, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    write.csv(
      macro, file.path(summary_dir, "reliable_macro_performance.csv"),
      row.names = FALSE
    )
  }
}

if (!is.null(importance) && nrow(importance)) {
  importance$rank <- ave(
    -importance$importance,
    interaction(importance$technology, importance$task, importance$fold),
    FUN = function(x) rank(x, ties.method = "min")
  )
  group_id <- interaction(
    importance$technology, importance$task, importance$fold, drop = TRUE
  )
  positive_importance <- pmax(importance$importance, 0)
  denom <- ave(positive_importance, group_id, FUN = sum)
  importance$normalized_importance <- ifelse(
    denom > 0, positive_importance / denom, 0
  )
  keys <- unique(importance[c("technology", "strategy", "task", "feature")])
  stability_rows <- lapply(seq_len(nrow(keys)), function(i) {
    z <- importance[
      importance$technology == keys$technology[i] &
        importance$task == keys$task[i] &
        importance$feature == keys$feature[i], , drop = FALSE
    ]
    data.frame(
      technology = keys$technology[i], strategy = keys$strategy[i],
      task = keys$task[i], feature = keys$feature[i],
      entered_fold_fraction = mean(z$entered_model),
      positive_fold_fraction = mean(z$importance > 0),
      top10_fold_fraction = mean(z$rank <= 10),
      top20_fold_fraction = mean(z$rank <= 20),
      median_rank = median(z$rank),
      median_importance = median(z$importance),
      median_normalized_importance = median(z$normalized_importance),
      stable_4of5_top20 = (
        mean(z$importance > 0) >= 0.8 && mean(z$rank <= 20) >= 0.8
      ),
      stringsAsFactors = FALSE
    )
  })
  stability <- do.call(rbind, stability_rows)
  stability <- stability[order(
    stability$technology, stability$task, -stability$top20_fold_fraction,
    stability$median_rank, -stability$median_normalized_importance
  ), ]
  write.csv(
    stability, file.path(summary_dir, "feature_importance_stability.csv"),
    row.names = FALSE
  )
}

cat("Week8 RF summary saved to: ", summary_dir, "\n", sep = "")
