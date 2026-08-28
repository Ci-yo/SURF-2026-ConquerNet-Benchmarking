## Summary-table functions used by the expanded LOTO and LOFO scripts.
## Fold construction, model fitting and evaluation-domain selection are kept
## in the main scripts so that the code can be reviewed from top to bottom.

summarize_loto <- function(metrics, output_dir) {
  evaluable <- metrics[is.finite(metrics$auroc) & is.finite(metrics$auprc), ]

  ## All-evaluable macro summary for the IDR-common primary domain and the
  ## no-IDR sensitivity domain.
  macro_groups <- split(
    evaluable,
    interaction(
      evaluable$coverage_cutoff,
      evaluable$evaluation_domain,
      evaluable$method,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  macro <- do.call(rbind, lapply(macro_groups, function(group) {
    data.frame(
      evaluation_domain = unique(group$evaluation_domain),
      coverage_cutoff = unique(group$coverage_cutoff),
      coverage_rule = unique(group$coverage_rule),
      method = unique(group$method),
      n_evaluable_folds = nrow(group),
      mean_auroc = mean(group$auroc),
      mean_auprc = mean(group$auprc),
      mean_n = mean(group$n),
      mean_positives = mean(group$positives),
      mean_negatives = mean(group$negatives),
      stringsAsFactors = FALSE
    )
  }))
  rownames(macro) <- NULL
  write.csv(
    macro,
    file.path(output_dir, "macro_average_summary.csv"),
    row.names = FALSE
  )

  analysis_sets <- list(
    all_14_technologies = evaluable,
    both_classes_ge_1000 = evaluable[
      evaluable$positives >= 1000L & evaluable$negatives >= 1000L,
      ,
      drop = FALSE
    ]
  )
  summary_rows <- list()
  for (set_name in names(analysis_sets)) {
    current_set <- analysis_sets[[set_name]]
    groups <- split(
      current_set,
      interaction(
        current_set$coverage_cutoff,
        current_set$evaluation_domain,
        current_set$method,
        drop = TRUE,
        lex.order = TRUE
      )
    )
    for (group in groups) {
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        analysis_set = set_name,
        evaluation_domain = unique(group$evaluation_domain),
        coverage_cutoff = unique(group$coverage_cutoff),
        coverage_rule = unique(group$coverage_rule),
        method = unique(group$method),
        n_folds = nrow(group),
        mean_auroc = mean(group$auroc),
        mean_auprc = mean(group$auprc),
        mean_n = mean(group$n),
        mean_positives = mean(group$positives),
        mean_negatives = mean(group$negatives),
        stringsAsFactors = FALSE
      )
    }
  }
  summary <- do.call(rbind, summary_rows)
  rownames(summary) <- NULL
  write.csv(
    summary[
      summary$evaluation_domain == "four_method_IDR_common_primary",
      ,
      drop = FALSE
    ],
    file.path(output_dir, "four_method_IDR_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    summary[
      summary$evaluation_domain == "three_method_no_IDR_common_sensitivity",
      ,
      drop = FALSE
    ],
    file.path(output_dir, "three_method_no_IDR_summary.csv"),
    row.names = FALSE
  )

  comparisons <- list()
  for (cutoff in unique(evaluable$coverage_cutoff)) {
    cutoff_rows <- evaluable[evaluable$coverage_cutoff == cutoff, ]
    for (domain in unique(cutoff_rows$evaluation_domain)) {
      current <- cutoff_rows[cutoff_rows$evaluation_domain == domain, ]
      conquer <- current[
        current$method == "ConquerBN",
        c("fold", "auroc", "auprc")
      ]
      names(conquer)[2:3] <- c("conquernet_auroc", "conquernet_auprc")
      for (baseline in setdiff(unique(current$method), "ConquerBN")) {
        other <- current[
          current$method == baseline,
          c("fold", "auroc", "auprc")
        ]
        names(other)[2:3] <- c("baseline_auroc", "baseline_auprc")
        paired <- merge(conquer, other, by = "fold")
        delta_roc <- paired$conquernet_auroc - paired$baseline_auroc
        delta_pr <- paired$conquernet_auprc - paired$baseline_auprc
        comparisons[[length(comparisons) + 1L]] <- data.frame(
          coverage_cutoff = cutoff,
          coverage_rule = unique(current$coverage_rule),
          evaluation_domain = domain,
          baseline_method = baseline,
          n_paired_auroc = sum(is.finite(delta_roc)),
          mean_delta_auroc = mean(delta_roc, na.rm = TRUE),
          median_delta_auroc = median(delta_roc, na.rm = TRUE),
          n_auroc_wins = sum(delta_roc > 0, na.rm = TRUE),
          n_auroc_losses = sum(delta_roc < 0, na.rm = TRUE),
          wilcoxon_p_auroc = safe_wilcox(delta_roc),
          n_paired_auprc = sum(is.finite(delta_pr)),
          mean_delta_auprc = mean(delta_pr, na.rm = TRUE),
          median_delta_auprc = median(delta_pr, na.rm = TRUE),
          n_auprc_wins = sum(delta_pr > 0, na.rm = TRUE),
          n_auprc_losses = sum(delta_pr < 0, na.rm = TRUE),
          wilcoxon_p_auprc = safe_wilcox(delta_pr),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  comparisons <- do.call(rbind, comparisons)
  adjustment_group <- interaction(
    comparisons$coverage_cutoff,
    comparisons$evaluation_domain,
    drop = TRUE
  )
  comparisons$wilcoxon_p_auroc_BH <- ave(
    comparisons$wilcoxon_p_auroc,
    adjustment_group,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  comparisons$wilcoxon_p_auprc_BH <- ave(
    comparisons$wilcoxon_p_auprc,
    adjustment_group,
    FUN = function(x) p.adjust(x, method = "BH")
  )
  write.csv(comparisons, file.path(output_dir, "paired_method_tests.csv"), row.names = FALSE)
}

family_macro_summary <- function(metrics) {
  evaluable <- metrics[is.finite(metrics$auroc) & is.finite(metrics$auprc), ]
  groups <- split(
    evaluable,
    interaction(evaluable$heldout_family, evaluable$method, drop = TRUE)
  )
  result <- do.call(rbind, lapply(groups, function(x) {
    data.frame(
      coverage_cutoff = unique(x$coverage_cutoff),
      coverage_rule = unique(x$coverage_rule),
      heldout_family = unique(x$heldout_family),
      method = unique(x$method),
      n_evaluable_technologies = nrow(x),
      technology_macro_auroc = mean(x$auroc),
      technology_macro_auprc = mean(x$auprc),
      mean_n_per_technology = mean(x$n),
      total_n_across_technologies = sum(x$n),
      total_positives_across_technologies = sum(x$positives),
      total_negatives_across_technologies = sum(x$negatives),
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result[order(result$heldout_family, result$method), ]
}

summarize_lofo <- function(metrics, output_dir) {
  family_macro <- family_macro_summary(metrics)
  write.csv(
    family_macro,
    file.path(output_dir, "per_family_method_summary.csv"),
    row.names = FALSE
  )

  groups <- split(family_macro, family_macro$method)
  overall <- do.call(rbind, lapply(groups, function(x) {
    data.frame(
      coverage_cutoff = unique(x$coverage_cutoff),
      coverage_rule = unique(x$coverage_rule),
      method = unique(x$method),
      n_independent_families = length(unique(x$heldout_family)),
      family_macro_auroc = mean(x$technology_macro_auroc),
      family_macro_auprc = mean(x$technology_macro_auprc),
      sd_family_auroc = stats::sd(x$technology_macro_auroc),
      sd_family_auprc = stats::sd(x$technology_macro_auprc),
      stringsAsFactors = FALSE
    )
  }))
  rownames(overall) <- NULL
  write.csv(
    overall,
    file.path(output_dir, "four_method_IDR_overall_summary.csv"),
    row.names = FALSE
  )

  conquer <- family_macro[
    family_macro$method == "ConquerBN",
    c("heldout_family", "technology_macro_auroc", "technology_macro_auprc")
  ]
  names(conquer)[2:3] <- c("conquernet_auroc", "conquernet_auprc")
  comparisons <- list()
  for (baseline in setdiff(unique(family_macro$method), "ConquerBN")) {
    other <- family_macro[
      family_macro$method == baseline,
      c("heldout_family", "technology_macro_auroc", "technology_macro_auprc")
    ]
    names(other)[2:3] <- c("baseline_auroc", "baseline_auprc")
    paired <- merge(conquer, other, by = "heldout_family")
    delta_roc <- paired$conquernet_auroc - paired$baseline_auroc
    delta_pr <- paired$conquernet_auprc - paired$baseline_auprc
    comparisons[[length(comparisons) + 1L]] <- data.frame(
      coverage_cutoff = unique(family_macro$coverage_cutoff),
      coverage_rule = unique(family_macro$coverage_rule),
      baseline_method = baseline,
      n_independent_families = nrow(paired),
      n_auroc_wins = sum(delta_roc > 0),
      mean_delta_auroc = mean(delta_roc),
      median_delta_auroc = median(delta_roc),
      auroc_exact_sign_test_p = safe_sign_test(delta_roc),
      auroc_exact_wilcox_p = safe_wilcox(delta_roc, exact = TRUE),
      n_auprc_wins = sum(delta_pr > 0),
      mean_delta_auprc = mean(delta_pr),
      median_delta_auprc = median(delta_pr),
      auprc_exact_sign_test_p = safe_sign_test(delta_pr),
      auprc_exact_wilcox_p = safe_wilcox(delta_pr, exact = TRUE),
      stringsAsFactors = FALSE
    )
  }
  write.csv(
    do.call(rbind, comparisons),
    file.path(output_dir, "paired_method_tests.csv"),
    row.names = FALSE
  )
}
