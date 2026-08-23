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

ivt_path <- normalizePath(
  cli_value(
    "ivt",
    paste0(
      "/home/xiyao/15data/W7/IVT_MeRIP-seq/",
      "MeRIP_hg38_WT_m6A_SE_raw_sample_aggregated_m6A_SE.rds"
    )
  ),
  winslash = "/",
  mustWork = TRUE
)
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
output_dir <- cli_value(
  "output-dir",
  "/home/xiyao/15data/W7/05_MeRIP_IVT_no_cali_TP_FP"
)
coverage_cutoff <- as.integer(cli_value("cutoff", "10"))
posterior_cutoff <- as.numeric(cli_value("posterior-cutoff", "0.5"))

if (!identical(coverage_cutoff, 10L)) {
  stop("Week7 formal analysis is fixed to Total>10.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ivt_se <- readRDS(ivt_path)
benchmark <- readRDS(benchmark_path)
label_object <- readRDS(label_path)

if (
  !inherits(ivt_se, "SummarizedExperiment") ||
    !inherits(benchmark, "SummarizedExperiment")
) {
  stop("IVT or benchmark input is not a SummarizedExperiment.")
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
if (
  !identical(
    coordinate_key(rowRanges(ivt_se)),
    coordinate_key(rowRanges(benchmark))
  )
) {
  stop("IVT and no-cali benchmark coordinates/order differ.")
}
rm(benchmark)
invisible(gc())

class_matrix <- label_object$class_matrix
if (
  !is.matrix(class_matrix) ||
    nrow(class_matrix) != nrow(ivt_se) ||
    !"MeRIP_seq" %in% colnames(class_matrix)
) {
  stop("No-cali class matrix is missing aligned MeRIP labels.")
}

class_code <- as.integer(class_matrix[, "MeRIP_seq"])
label <- rep(NA_character_, length(class_code))
label[class_code == 2L] <- "FP"
label[class_code == 3L] <- "TP"
eligible_label <- label %in% c("TP", "FP")

expected_samples <- c(
  "h293_MeRIP_Abcam_IVT_IP",
  "h293_MeRIP_Abcam_mRNA_IP",
  "h293_MeRIP_NEB_IVT_IP",
  "h293_MeRIP_NEB_mRNA_IP",
  "h293_MeRIP_SYSY_IVT_IP",
  "h293_MeRIP_SYSY_mRNA_IP"
)
if (!identical(colnames(ivt_se), expected_samples)) {
  stop("Unexpected IVT/mRNA sample order.")
}

total <- assay(ivt_se, "Total")
probability <- assay(ivt_se, "m6ASiteProb")
observed <- (
  is.finite(total) &
    total > coverage_cutoff &
    is.finite(probability)
)
positive <- observed & probability > posterior_cutoff

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

summarize_endpoint <- function(scope, endpoint, eligible, event) {
  use <- eligible_label & eligible & !is.na(event)
  tp <- use & label == "TP"
  fp <- use & label == "FP"
  tp_event <- sum(event[tp])
  fp_event <- sum(event[fp])
  tp_total <- sum(tp)
  fp_total <- sum(fp)

  if (tp_total == 0L || fp_total == 0L) {
    return(data.frame(
      scope = scope,
      endpoint = endpoint,
      n_tp = tp_total,
      n_fp = fp_total,
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

  stats <- safe_fisher(
    tp_event,
    tp_total - tp_event,
    fp_event,
    fp_total - fp_event
  )
  data.frame(
    scope = scope,
    endpoint = endpoint,
    n_tp = tp_total,
    n_fp = fp_total,
    tp_event = tp_event,
    fp_event = fp_event,
    tp_rate = tp_event / tp_total,
    fp_rate = fp_event / fp_total,
    risk_difference_tp_minus_fp =
      tp_event / tp_total - fp_event / fp_total,
    odds_ratio_tp_vs_fp = stats[["odds_ratio_tp_vs_fp"]],
    ci_low = stats[["ci_low"]],
    ci_high = stats[["ci_high"]],
    p_value = stats[["p_value"]],
    stringsAsFactors = FALSE
  )
}

pairs <- list(
  Abcam = c(
    IVT = "h293_MeRIP_Abcam_IVT_IP",
    mRNA = "h293_MeRIP_Abcam_mRNA_IP"
  ),
  NEB = c(
    IVT = "h293_MeRIP_NEB_IVT_IP",
    mRNA = "h293_MeRIP_NEB_mRNA_IP"
  ),
  SYSY = c(
    IVT = "h293_MeRIP_SYSY_IVT_IP",
    mRNA = "h293_MeRIP_SYSY_mRNA_IP"
  )
)

rows <- list()
pair_calls <- list()

for (brand in names(pairs)) {
  ivt_index <- match(pairs[[brand]][["IVT"]], expected_samples)
  mrna_index <- match(pairs[[brand]][["mRNA"]], expected_samples)
  both_observed <- observed[, ivt_index] & observed[, mrna_index]
  ivt_positive <- positive[, ivt_index]
  mrna_positive <- positive[, mrna_index]

  endpoints <- list(
    IVT_positive = ivt_positive,
    mRNA_positive = mrna_positive,
    mRNA_specific = mrna_positive & !ivt_positive,
    IVT_specific = ivt_positive & !mrna_positive,
    both_positive = ivt_positive & mrna_positive
  )
  for (endpoint in names(endpoints)) {
    rows[[length(rows) + 1L]] <- summarize_endpoint(
      scope = paste0(brand, "_pair"),
      endpoint = endpoint,
      eligible = both_observed,
      event = endpoints[[endpoint]]
    )
  }

  pair_calls[[brand]] <- list(
    eligible = both_observed,
    IVT_positive = ivt_positive,
    mRNA_positive = mrna_positive,
    mRNA_specific = mrna_positive & !ivt_positive,
    IVT_specific = ivt_positive & !mrna_positive,
    both_positive = ivt_positive & mrna_positive
  )
}

all_pairs_observed <- Reduce(
  `&`,
  lapply(pair_calls, function(value) value$eligible)
)
for (
  endpoint in c(
    "IVT_positive",
    "mRNA_positive",
    "mRNA_specific",
    "IVT_specific",
    "both_positive"
  )
) {
  endpoint_matrix <- do.call(
    cbind,
    lapply(pair_calls, function(value) value[[endpoint]])
  )
  count <- rowSums(endpoint_matrix)
  for (threshold in c(1L, 2L, 3L)) {
    rows[[length(rows) + 1L]] <- summarize_endpoint(
      scope = "three_antibody_consensus",
      endpoint = paste0(endpoint, "_in_at_least_", threshold, "_of_3"),
      eligible = all_pairs_observed,
      event = count >= threshold
    )
  }
}

results <- do.call(rbind, rows)
results$coverage_rule <- "Total>10 in both IVT and mRNA"
results$call_rule <- "m6ASiteProb>0.5"
results$label_source <- label_object$decisions$selected_label_source[
  match("MeRIP_seq", label_object$decisions$technology)
]
results$p_adj_bh <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$scope, results$endpoint), ]

label_counts <- data.frame(
  label = c("TP", "FP"),
  n = c(sum(label == "TP", na.rm = TRUE), sum(label == "FP", na.rm = TRUE)),
  label_source = results$label_source[[1]],
  stringsAsFactors = FALSE
)

sample_counts <- data.frame(
  sample = expected_samples,
  n_total_gt_10 = colSums(observed),
  n_map_positive = colSums(positive),
  stringsAsFactors = FALSE
)

write.csv(
  results,
  file.path(output_dir, "MeRIP_no_cali_TP_FP_IVT_endpoints.csv"),
  row.names = FALSE
)
write.csv(
  label_counts,
  file.path(output_dir, "MeRIP_no_cali_TP_FP_label_counts.csv"),
  row.names = FALSE
)
write.csv(
  sample_counts,
  file.path(output_dir, "MeRIP_IVT_sample_call_counts.csv"),
  row.names = FALSE
)

primary <- results[
  results$scope == "three_antibody_consensus" &
    results$endpoint %in% c(
      "IVT_positive_in_at_least_2_of_3",
      "mRNA_specific_in_at_least_2_of_3"
    ),
  ,
  drop = FALSE
]
write.csv(
  primary,
  file.path(output_dir, "MeRIP_IVT_primary_2of3_summary.csv"),
  row.names = FALSE
)

protocol <- c(
  "Week7 B: MeRIP IVT validation against no-calibration ConquerBN TP/FP labels",
  paste0("Benchmark: ", benchmark_path),
  paste0("Labels: ", label_path),
  paste0("IVT data: ", ivt_path),
  "Coverage/call rule: each IVT and mRNA aggregate must have Total>10; MAP positive is m6ASiteProb>0.5.",
  "Primary summary: at least 2 of 3 antibodies.",
  "Supportive direction: TP enriched over FP among mRNA-specific calls.",
  "Artifact direction: FP enriched over TP among IVT-positive calls.",
  "Per-antibody endpoints are retained because antibody effects may differ.",
  "The inherited high_confidence_union/specific rowData annotations are not used as independent validation endpoints."
)
writeLines(protocol, file.path(output_dir, "ANALYSIS_PROTOCOL.txt"))

message("Saved formal no-cali MeRIP IVT analysis to: ", output_dir)
print(primary, row.names = FALSE)
