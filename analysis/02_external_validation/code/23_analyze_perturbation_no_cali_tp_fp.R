suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
cli_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) default else sub(prefix, "", hit[[1]], fixed = TRUE)
}

benchmark_path <- normalizePath(
  cli_value(
    "benchmark",
    paste0(
      "/home/xichen/ConquerNet/ConquerBN_basic/",
      "m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali.rds"
    )
  ),
  winslash = "/",
  mustWork = TRUE
)
label_path <- normalizePath(
  cli_value(
    "labels",
    paste0(
      "/home/xiyao/15data/W7/04_no_cali_tp_fp_labels/Total_gt_10/",
      "selected_tp_fp_class_matrix.rds"
    )
  ),
  winslash = "/",
  mustWork = TRUE
)
perturbation_dir <- normalizePath(
  cli_value(
    "perturbation-dir",
    "/home/xiyao/15data/W7/purturbation"
  ),
  winslash = "/",
  mustWork = TRUE
)
output_dir <- cli_value(
  "output-dir",
  "/home/xiyao/15data/W7/06_perturbation_no_cali_TP_FP"
)
coverage_cutoff <- as.integer(cli_value("cutoff", "10"))
posterior_cutoff <- as.numeric(cli_value("posterior-cutoff", "0.5"))
minimum_replicates <- as.integer(cli_value("min-replicates", "2"))

if (!identical(coverage_cutoff, 10L)) {
  stop("Week7 formal analysis is fixed to Total>10.")
}
if (!is.finite(minimum_replicates) || minimum_replicates < 1L) {
  stop("--min-replicates must be positive.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

benchmark <- readRDS(benchmark_path)
label_object <- readRDS(label_path)
if (!inherits(benchmark, "SummarizedExperiment")) {
  stop("Benchmark is not a SummarizedExperiment.")
}
if (
  !identical(
    label_object$analysis_version,
    "week7_no_cali_selected_labels_Total_gt_10_v1"
  ) ||
    !identical(label_object$dataset_path, benchmark_path) ||
    !identical(as.integer(label_object$coverage_cutoff), coverage_cutoff)
) {
  stop("Label object is not the formal no-cali Total>10 version.")
}

coordinate_key <- function(ranges) {
  paste(
    as.character(seqnames(ranges)),
    start(ranges),
    end(ranges),
    as.character(strand(ranges)),
    sep = ":"
  )
}
benchmark_key <- coordinate_key(rowRanges(benchmark))
n_sites <- nrow(benchmark)
rm(benchmark)
invisible(gc())

## Only contrasts whose treatment and control can be identified from one
## accession series/study are included. Ambiguous or control-free samples
## are intentionally excluded.
contrast_manifest <- data.frame(
  contrast = c(
    "GLORI_siMETTL3_vs_siNC_writer",
    "GLORI_siMETTL14_vs_siNC_writer",
    "GLORI_siWTAP_vs_siNC_writer",
    "m6ACE_Mettl3_KO100_vs_WT100",
    "m6ACE_FTO_OE_vs_WT100",
    "MeRIP_A549_shMETTL3_vs_shGFP",
    "MeRIP_A549_shMETTL14_vs_shGFP",
    "MeRIP_A549_shWTAP_vs_shGFP"
  ),
  technology = c(
    rep("GLORI", 3L),
    rep("m6ACE_seq", 2L),
    rep("MeRIP_seq", 3L)
  ),
  perturbation_file = c(
    rep("GLORI_hg38_Perturb_m6A_SE.rds", 3L),
    rep("m6ACE_hg38_Perturb_m6A_SE.rds", 2L),
    rep("MeRIP_hg38_Perturb_m6A_SE.rds", 3L)
  ),
  perturbation_type = c(
    rep("writer_loss", 4L),
    "eraser_overexpression",
    rep("writer_loss", 3L)
  ),
  expected_probability_change = rep("decrease", 8L),
  control_sample_ids = c(
    rep("SRR21356221;SRR21356222", 3L),
    rep("SRR9060612;SRR9060622;SRR9060632", 2L),
    rep("SRR1182619;SRR1182621;SRR1182623", 3L)
  ),
  treatment_sample_ids = c(
    "SRR21356233;SRR21356234",
    "SRR21356229;SRR21356230",
    "SRR21356225;SRR21356226",
    "SRR9060620;SRR9060630;SRR9060640",
    "SRR9060606;SRR9060608;SRR9060610",
    "SRR1182615;SRR1182617",
    "SRR1182607;SRR1182609;SRR1182611;SRR1182613",
    "SRR1182603;SRR1182605"
  ),
  stringsAsFactors = FALSE
)

split_ids <- function(value) strsplit(value, ";", fixed = TRUE)[[1]]

group_probability <- function(probability, total, sample_indices, minimum) {
  group_probability <- probability[, sample_indices, drop = FALSE]
  group_total <- total[, sample_indices, drop = FALSE]
  valid <- (
    is.finite(group_total) &
      group_total > coverage_cutoff &
      is.finite(group_probability)
  )
  group_probability[!valid] <- NA_real_
  n_valid <- rowSums(valid)
  mean_probability <- rowMeans(group_probability, na.rm = TRUE)
  mean_probability[n_valid < minimum] <- NA_real_
  list(mean_probability = mean_probability, n_valid = n_valid)
}

safe_fisher <- function(tp_event, tp_nonevent, fp_event, fp_nonevent) {
  table <- matrix(
    c(tp_event, tp_nonevent, fp_event, fp_nonevent),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("TP", "FP"), c("event", "no_event"))
  )
  fit <- fisher.test(table)
  c(
    odds_ratio_tp_vs_fp = unname(fit$estimate),
    ci_low = unname(fit$conf.int[[1]]),
    ci_high = unname(fit$conf.int[[2]]),
    p_value = fit$p.value
  )
}

summarize_binary <- function(
  contrast,
  endpoint,
  label,
  eligible,
  event
) {
  use <- eligible & label %in% c("TP", "FP") & !is.na(event)
  tp <- use & label == "TP"
  fp <- use & label == "FP"
  n_tp <- sum(tp)
  n_fp <- sum(fp)
  tp_event <- sum(event[tp])
  fp_event <- sum(event[fp])
  if (n_tp == 0L || n_fp == 0L) {
    return(data.frame(
      contrast = contrast,
      endpoint = endpoint,
      n_tp = n_tp,
      n_fp = n_fp,
      tp_event = tp_event,
      fp_event = fp_event,
      tp_rate = NA_real_,
      fp_rate = NA_real_,
      risk_difference_tp_minus_fp = NA_real_,
      odds_ratio_tp_vs_fp = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  statistics <- safe_fisher(
    tp_event,
    n_tp - tp_event,
    fp_event,
    n_fp - fp_event
  )
  data.frame(
    contrast = contrast,
    endpoint = endpoint,
    n_tp = n_tp,
    n_fp = n_fp,
    tp_event = tp_event,
    fp_event = fp_event,
    tp_rate = tp_event / n_tp,
    fp_rate = fp_event / n_fp,
    risk_difference_tp_minus_fp = tp_event / n_tp - fp_event / n_fp,
    odds_ratio_tp_vs_fp = statistics[["odds_ratio_tp_vs_fp"]],
    ci_low = statistics[["ci_low"]],
    ci_high = statistics[["ci_high"]],
    p_value = statistics[["p_value"]],
    stringsAsFactors = FALSE
  )
}

continuous_rows <- list()
binary_rows <- list()
audit_rows <- list()
object_cache <- list()

for (i in seq_len(nrow(contrast_manifest))) {
  definition <- contrast_manifest[i, , drop = FALSE]
  technology <- definition$technology[[1]]
  object_name <- definition$perturbation_file[[1]]
  object_path <- file.path(perturbation_dir, object_name)

  if (is.null(object_cache[[object_name]])) {
    object <- readRDS(object_path)
    if (!inherits(object, "SummarizedExperiment")) {
      stop("Perturbation object is not a SummarizedExperiment: ", object_path)
    }
    if (
      nrow(object) != n_sites ||
        !identical(coordinate_key(rowRanges(object)), benchmark_key)
    ) {
      stop("Perturbation coordinates/order differ: ", object_path)
    }
    object_cache[[object_name]] <- object
  }
  object <- object_cache[[object_name]]

  metadata <- as.data.frame(colData(object))
  if (!"SampleID" %in% names(metadata)) {
    stop("SampleID is missing from ", object_name)
  }
  control_ids <- split_ids(definition$control_sample_ids[[1]])
  treatment_ids <- split_ids(definition$treatment_sample_ids[[1]])
  control_indices <- match(control_ids, metadata$SampleID)
  treatment_indices <- match(treatment_ids, metadata$SampleID)
  if (anyNA(control_indices) || anyNA(treatment_indices)) {
    stop("Contrast sample IDs were not found: ", definition$contrast[[1]])
  }

  probability <- assay(object, "m6ASiteProb")
  total <- assay(object, "Total")
  control <- group_probability(
    probability,
    total,
    control_indices,
    minimum_replicates
  )
  treatment <- group_probability(
    probability,
    total,
    treatment_indices,
    minimum_replicates
  )

  ## All currently included contrasts are expected to reduce m6A signal:
  ## writer loss or eraser overexpression. Positive response is therefore
  ## control probability minus treatment probability.
  response <- control$mean_probability - treatment$mean_probability
  eligible <- is.finite(response)

  class_code <- as.integer(label_object$class_matrix[, technology])
  label <- rep(NA_character_, n_sites)
  label[class_code == 2L] <- "FP"
  label[class_code == 3L] <- "TP"
  use <- eligible & label %in% c("TP", "FP")
  tp_response <- response[use & label == "TP"]
  fp_response <- response[use & label == "FP"]
  n_tp <- length(tp_response)
  n_fp <- length(fp_response)
  if (n_tp > 0L && n_fp > 0L) {
    wilcoxon_two_sided <- wilcox.test(
      tp_response,
      fp_response,
      alternative = "two.sided",
      exact = FALSE
    )
    wilcoxon_greater <- wilcox.test(
      tp_response,
      fp_response,
      alternative = "greater",
      exact = FALSE
    )
    ## R reports the rank-sum statistic after subtracting its minimum,
    ## i.e. the Mann-Whitney U statistic on the 0..(n_tp*n_fp) scale.
    u_statistic <- unname(wilcoxon_two_sided$statistic)
    probability_tp_greater_fp <- (
      as.numeric(u_statistic) /
        (as.numeric(n_tp) * as.numeric(n_fp))
    )
    continuous_row <- data.frame(
      contrast = definition$contrast[[1]],
      technology = technology,
      perturbation_type = definition$perturbation_type[[1]],
      expected_probability_change = "decrease",
      response_definition =
        "mean_control_probability_minus_mean_treatment_probability",
      n_tp = n_tp,
      n_fp = n_fp,
      tp_mean_response = mean(tp_response),
      fp_mean_response = mean(fp_response),
      tp_median_response = median(tp_response),
      fp_median_response = median(fp_response),
      median_difference_tp_minus_fp =
        median(tp_response) - median(fp_response),
      probability_tp_response_greater_than_fp =
        probability_tp_greater_fp,
      rank_biserial_effect = 2 * probability_tp_greater_fp - 1,
      wilcoxon_p_two_sided = wilcoxon_two_sided$p.value,
      wilcoxon_p_tp_greater = wilcoxon_greater$p.value,
      stringsAsFactors = FALSE
    )
  } else {
    continuous_row <- data.frame(
      contrast = definition$contrast[[1]],
      technology = technology,
      perturbation_type = definition$perturbation_type[[1]],
      expected_probability_change = "decrease",
      response_definition =
        "mean_control_probability_minus_mean_treatment_probability",
      n_tp = n_tp,
      n_fp = n_fp,
      tp_mean_response = NA_real_,
      fp_mean_response = NA_real_,
      tp_median_response = NA_real_,
      fp_median_response = NA_real_,
      median_difference_tp_minus_fp = NA_real_,
      probability_tp_response_greater_than_fp = NA_real_,
      rank_biserial_effect = NA_real_,
      wilcoxon_p_two_sided = NA_real_,
      wilcoxon_p_tp_greater = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  continuous_rows[[length(continuous_rows) + 1L]] <- continuous_row

  binary_rows[[length(binary_rows) + 1L]] <- summarize_binary(
    definition$contrast[[1]],
    "response_in_expected_direction",
    label,
    eligible,
    response > 0
  )
  binary_rows[[length(binary_rows) + 1L]] <- summarize_binary(
    definition$contrast[[1]],
    "MAP_positive_to_negative",
    label,
    eligible,
    control$mean_probability > posterior_cutoff &
      treatment$mean_probability <= posterior_cutoff
  )

  audit_rows[[length(audit_rows) + 1L]] <- data.frame(
    contrast = definition$contrast[[1]],
    technology = technology,
    perturbation_file = object_path,
    n_control_replicates = length(control_indices),
    n_treatment_replicates = length(treatment_indices),
    minimum_valid_replicates_per_arm = minimum_replicates,
    n_sites_valid_control = sum(is.finite(control$mean_probability)),
    n_sites_valid_treatment = sum(is.finite(treatment$mean_probability)),
    n_sites_valid_both = sum(eligible),
    n_tp_valid_both = n_tp,
    n_fp_valid_both = n_fp,
    label_source = label_object$decisions$selected_label_source[
      match(technology, label_object$decisions$technology)
    ],
    stringsAsFactors = FALSE
  )

  rm(
    probability,
    total,
    control,
    treatment,
    response,
    eligible,
    class_code,
    label,
    use,
    tp_response,
    fp_response
  )
  invisible(gc())
}

continuous_results <- do.call(rbind, continuous_rows)
binary_results <- do.call(rbind, binary_rows)
audit <- do.call(rbind, audit_rows)

continuous_results$wilcoxon_p_two_sided_bh <- p.adjust(
  continuous_results$wilcoxon_p_two_sided,
  method = "BH"
)
continuous_results$wilcoxon_p_tp_greater_bh <- p.adjust(
  continuous_results$wilcoxon_p_tp_greater,
  method = "BH"
)
binary_results$p_adj_bh <- p.adjust(binary_results$p_value, method = "BH")

write.csv(
  contrast_manifest,
  file.path(output_dir, "contrast_manifest.csv"),
  row.names = FALSE
)
write.csv(
  audit,
  file.path(output_dir, "contrast_eligibility_audit.csv"),
  row.names = FALSE
)
write.csv(
  continuous_results,
  file.path(output_dir, "TP_vs_FP_continuous_response.csv"),
  row.names = FALSE
)
write.csv(
  binary_results,
  file.path(output_dir, "TP_vs_FP_binary_response.csv"),
  row.names = FALSE
)

protocol <- c(
  "Week7 C: perturbation response of no-calibration ConquerBN TP versus FP labels",
  paste0("Benchmark: ", benchmark_path),
  paste0("Labels: ", label_path),
  paste0("Perturbation directory: ", perturbation_dir),
  "Only contrasts with identifiable same-series controls are included.",
  "Each arm requires at least two replicates with Total>10 and finite m6ASiteProb at a site.",
  "Primary continuous response is mean(control m6ASiteProb) - mean(treatment m6ASiteProb).",
  "This sign is positive in the expected direction for writer loss and eraser overexpression.",
  "Primary hypothesis: TP response is greater than FP response.",
  "Secondary endpoints: positive signed response and MAP positive-to-negative transition.",
  "Control-free eTAM, GLORI FTO+, and unmatched perturbation samples are excluded."
)
writeLines(protocol, file.path(output_dir, "ANALYSIS_PROTOCOL.txt"))

message("Saved no-cali perturbation analysis to: ", output_dir)
print(continuous_results, row.names = FALSE)
