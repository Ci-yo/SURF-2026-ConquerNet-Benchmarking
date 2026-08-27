#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[1L]])
}

class_path <- get_arg(
  "class-matrix",
  "/home/xiyao/15data/W8/results/01_four_class_global_fdr/four_class_matrix_Total_gt_10.rds"
)
benchmark_path <- get_arg(
  "benchmark",
  "/home/xichen/ConquerNet/ConquerBN_basic/m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali_IVT_matched.rds"
)
features_path <- get_arg(
  "features",
  "/home/xiyao/15data/W8/input_data/features/m6Aconquer_omicsFeatures_hg38.rds"
)
output_dir <- get_arg(
  "output-dir",
  "/home/xiyao/15data/W8/results/02_fp_tp_meta_distributions/Total_gt_10_no_cali_IVT_matched/standard_metagene"
)
coverage_cutoff <- as.numeric(get_arg("coverage-cutoff", "10"))
min_class <- as.integer(get_arg("min-class", "51"))
bw <- as.numeric(get_arg("bandwidth", "0.02"))
n_grid <- as.integer(get_arg("n-grid", "512"))

if (!is.finite(min_class) || min_class < 1L) stop("Invalid --min-class.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading class matrix ...")
class_matrix <- readRDS(class_path)
message("Loading benchmark ...")
benchmark <- readRDS(benchmark_path)
message("Loading OmicsFeatures ...")
features <- readRDS(features_path)

if (!is.matrix(class_matrix)) stop("Class object must be a matrix.")
if (nrow(class_matrix) != length(features) || nrow(class_matrix) != nrow(benchmark)) {
  stop("Universe-size mismatch among class matrix, features, and benchmark.")
}

total <- SummarizedExperiment::assay(benchmark, "Total")
if (!identical(colnames(class_matrix), colnames(total))) {
  stop("Technology columns differ between class matrix and Total assay.")
}

rr <- SummarizedExperiment::rowRanges(benchmark)
coordinate_identity <- identical(as.character(seqnames(rr)), as.character(seqnames(features))) &&
  identical(start(rr), start(features)) &&
  identical(end(rr), end(features)) &&
  identical(as.character(strand(rr)), as.character(strand(features)))
if (!coordinate_identity) stop("Benchmark and feature coordinates are not identical.")

feature_df <- as.data.frame(S4Vectors::mcols(features))
required <- c(
  "meta_tx_topology",
  "overlap_exonicFivePrimeUTR",
  "overlap_exonicCDS",
  "overlap_exonicThreePrimeUTR"
)
missing_features <- setdiff(required, names(feature_df))
if (length(missing_features)) {
  stop("Missing feature columns: ", paste(missing_features, collapse = ", "))
}

meta_tx <- as.numeric(feature_df$meta_tx_topology)
mask_5utr <- as.numeric(feature_df$overlap_exonicFivePrimeUTR) > 0
mask_cds <- as.numeric(feature_df$overlap_exonicCDS) > 0
mask_3utr <- as.numeric(feature_df$overlap_exonicThreePrimeUTR) > 0
annotation_mask <- mask_5utr | mask_cds | mask_3utr
valid_meta <- is.finite(meta_tx) & meta_tx >= 0 & meta_tx <= 1 & annotation_mask

region_specs <- data.frame(
  region = c("5UTR", "CDS", "3UTR"),
  expected_min = c(0, 1 / 3, 2 / 3),
  expected_max = c(1 / 3, 2 / 3, 1),
  stringsAsFactors = FALSE
)
region_masks <- list(
  mask_5utr & !mask_cds & !mask_3utr,
  mask_cds & !mask_5utr & !mask_3utr,
  mask_3utr & !mask_5utr & !mask_cds
)

encoding_audit <- do.call(rbind, lapply(seq_len(nrow(region_specs)), function(i) {
  values <- meta_tx[region_masks[[i]] & is.finite(meta_tx)]
  in_expected <- values >= region_specs$expected_min[i] & values <= region_specs$expected_max[i]
  data.frame(
    region = region_specs$region[i],
    audit_mask = "exclusive annotation only",
    n_annotated = length(values),
    value_min = if (length(values)) min(values) else NA_real_,
    value_q01 = if (length(values)) unname(quantile(values, 0.01)) else NA_real_,
    value_median = if (length(values)) median(values) else NA_real_,
    value_q99 = if (length(values)) unname(quantile(values, 0.99)) else NA_real_,
    value_max = if (length(values)) max(values) else NA_real_,
    expected_min = region_specs$expected_min[i],
    expected_max = region_specs$expected_max[i],
    fraction_in_expected_third = if (length(values)) mean(in_expected) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
write.csv(encoding_audit, file.path(output_dir, "meta_tx_encoding_audit.csv"), row.names = FALSE)

if (any(encoding_audit$fraction_in_expected_third < 0.95, na.rm = TRUE)) {
  warning("At least one exclusively annotated region is not almost entirely within its expected meta_tx third; inspect meta_tx_encoding_audit.csv.")
}

technology_status <- lapply(seq_len(ncol(class_matrix)), function(j) {
  technology <- colnames(class_matrix)[j]
  covered <- is.finite(total[, j]) & total[, j] > coverage_cutoff
  n_fp <- sum(class_matrix[, j] == 2L & covered, na.rm = TRUE)
  n_tp <- sum(class_matrix[, j] == 3L & covered, na.rm = TRUE)
  status <- if (n_fp >= min_class && n_tp >= min_class) {
    "reliable"
  } else if (n_fp > 0L && n_tp > 0L) {
    "exploratory"
  } else {
    "untrainable"
  }
  data.frame(technology, n_fp, n_tp, status, stringsAsFactors = FALSE)
})
technology_status <- do.call(rbind, technology_status)
write.csv(technology_status, file.path(output_dir, "technology_status.csv"), row.names = FALSE)

safe_density <- function(x, bw, n_grid) {
  if (length(x) < 2L || length(unique(x)) < 2L) return(NULL)
  density(x, bw = bw, from = 0, to = 1, n = n_grid, cut = 0)
}

density_rows <- list()
sample_rows <- list()

for (j in seq_len(ncol(class_matrix))) {
  technology <- colnames(class_matrix)[j]
  covered <- is.finite(total[, j]) & total[, j] > coverage_cutoff
  fp_idx <- which(class_matrix[, j] == 2L & covered & valid_meta)
  tp_idx <- which(class_matrix[, j] == 3L & covered & valid_meta)

  sample_rows[[length(sample_rows) + 1L]] <- data.frame(
    technology = technology,
    status = technology_status$status[technology_status$technology == technology],
    n_fp_labels = technology_status$n_fp[technology_status$technology == technology],
    n_tp_labels = technology_status$n_tp[technology_status$technology == technology],
    n_fp_meta_tx = length(fp_idx),
    n_tp_meta_tx = length(tp_idx),
    fp_assignment_rate = if (technology_status$n_fp[technology_status$technology == technology] > 0) length(fp_idx) / technology_status$n_fp[technology_status$technology == technology] else NA_real_,
    tp_assignment_rate = if (technology_status$n_tp[technology_status$technology == technology] > 0) length(tp_idx) / technology_status$n_tp[technology_status$technology == technology] else NA_real_,
    stringsAsFactors = FALSE
  )

  for (class_name in c("FP", "TP")) {
    idx <- if (class_name == "FP") fp_idx else tp_idx
    d <- safe_density(meta_tx[idx], bw, n_grid)
    if (!is.null(d)) {
      density_rows[[length(density_rows) + 1L]] <- data.frame(
        technology = technology,
        status = technology_status$status[technology_status$technology == technology],
        class = class_name,
        x = d$x,
        density = d$y,
        n = length(idx),
        bandwidth = bw,
        stringsAsFactors = FALSE
      )
    }
  }
}

sample_audit <- do.call(rbind, sample_rows)
density_data <- do.call(rbind, density_rows)
write.csv(sample_audit, file.path(output_dir, "metagene_sample_audit.csv"), row.names = FALSE)
write.csv(density_data, file.path(output_dir, "metagene_density_data.csv"), row.names = FALSE)

tp_col <- "#C83B35"
fp_col <- "#2B6CB0"

draw_panel <- function(technology, show_legend = TRUE) {
  dsub <- density_data[density_data$technology == technology, , drop = FALSE]
  nsub <- sample_audit[sample_audit$technology == technology, , drop = FALSE]
  if (!nrow(dsub) || !nrow(nsub)) {
    plot.new()
    title(main = technology)
    text(0.5, 0.5, "FP or TP unavailable")
    return(invisible(NULL))
  }

  ymax <- max(dsub$density, na.rm = TRUE) * 1.08
  plot(
    NA_real_, NA_real_,
    xlim = c(0, 1), ylim = c(0, ymax),
    xaxt = "n", xlab = "", ylab = "Density",
    main = technology, bty = "o", las = 1
  )
  axis(1, at = c(1 / 6, 1 / 2, 5 / 6), labels = c("5'UTR", "CDS", "3'UTR"), tick = FALSE, line = 0.5)
  abline(v = c(1 / 3, 2 / 3), lty = 2, col = "#777777", lwd = 1.2)

  for (class_name in c("FP", "TP")) {
    one <- dsub[dsub$class == class_name, , drop = FALSE]
    if (nrow(one)) {
      lines(one$x, one$density, col = if (class_name == "TP") tp_col else fp_col, lwd = 2.2)
    }
  }

  usr <- par("usr")
  yrange <- usr[4] - usr[3]
  y0 <- usr[3] - 0.18 * yrange
  y1 <- usr[3] - 0.145 * yrange
  rect(0, y0, 1 / 3, y1, col = "#202020", border = NA, xpd = NA)
  rect(1 / 3, y0, 2 / 3, y1, col = "#8C8C8C", border = "#202020", xpd = NA)
  rect(2 / 3, y0, 1, y1, col = "#202020", border = NA, xpd = NA)

  subtitle <- sprintf(
    "FP n=%s; TP n=%s",
    format(nsub$n_fp_meta_tx, big.mark = ","),
    format(nsub$n_tp_meta_tx, big.mark = ",")
  )
  mtext(subtitle, side = 3, line = 0.2, cex = 0.72)
  if (show_legend) {
    legend("topleft", legend = c("TP", "FP"), col = c(tp_col, fp_col), lwd = 2.2, bty = "n", cex = 0.85)
  }
}

reliable <- technology_status$technology[technology_status$status == "reliable"]

if (length(reliable)) {
  panel_ncol <- if (length(reliable) > 6L) 4L else 3L
  panel_nrow <- ceiling(length(reliable) / panel_ncol)
  combined_png <- file.path(output_dir, "reliable_technologies_FP_vs_TP_standard_metagene.png")
  png(combined_png, width = 800 * panel_ncol, height = 800 * panel_nrow, res = 220)
  par(mfrow = c(panel_nrow, panel_ncol), mar = c(4.6, 4.2, 3.2, 1.0), oma = c(0.5, 0.5, 2.5, 0.5))
  for (technology in reliable) draw_panel(technology, show_legend = TRUE)
  mtext("Posterior-derived FP vs TP: standardized 5'UTR-CDS-3'UTR topology", outer = TRUE, side = 3, line = 0.7, cex = 1.25, font = 2)
  dev.off()

  combined_pdf <- file.path(output_dir, "reliable_technologies_FP_vs_TP_standard_metagene.pdf")
  pdf(combined_pdf, width = 4 * panel_ncol, height = 4 * panel_nrow)
  par(mfrow = c(panel_nrow, panel_ncol), mar = c(4.6, 4.2, 3.2, 1.0), oma = c(0.5, 0.5, 2.5, 0.5))
  for (technology in reliable) draw_panel(technology, show_legend = TRUE)
  mtext("Posterior-derived FP vs TP: standardized 5'UTR-CDS-3'UTR topology", outer = TRUE, side = 3, line = 0.7, cex = 1.25, font = 2)
  dev.off()
}

individual_dir <- file.path(output_dir, "per_technology")
dir.create(individual_dir, recursive = TRUE, showWarnings = FALSE)
for (technology in technology_status$technology[technology_status$status != "untrainable"]) {
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", technology)
  png(file.path(individual_dir, paste0(safe_name, "_FP_vs_TP_standard_metagene.png")), width = 1500, height = 1100, res = 220)
  par(mar = c(5.2, 4.6, 3.6, 1.0))
  draw_panel(technology, show_legend = TRUE)
  dev.off()
  pdf(file.path(individual_dir, paste0(safe_name, "_FP_vs_TP_standard_metagene.pdf")), width = 6.5, height = 4.8)
  par(mar = c(5.2, 4.6, 3.6, 1.0))
  draw_panel(technology, show_legend = TRUE)
  dev.off()
}

qc <- data.frame(
  item = c(
    "class_matrix_path", "benchmark_path", "features_path", "n_universe_sites",
    "coordinate_identity", "coverage_rule", "class_codes", "fallback",
    "meta_tx_annotation_mask", "region_boundaries", "density_bandwidth",
    "density_grid_points", "minimum_class_count", "n_reliable", "n_exploratory", "n_untrainable"
  ),
  value = c(
    class_path, benchmark_path, features_path, nrow(class_matrix),
    coordinate_identity, paste0("Total > ", coverage_cutoff), "FP=2;TP=3", "none",
    "exonicFivePrimeUTR OR exonicCDS OR exonicThreePrimeUTR",
    "5UTR:[0,1/3];CDS:[1/3,2/3];3UTR:[2/3,1]",
    bw, n_grid, min_class,
    sum(technology_status$status == "reliable"),
    sum(technology_status$status == "exploratory"),
    sum(technology_status$status == "untrainable")
  ),
  stringsAsFactors = FALSE
)
write.csv(qc, file.path(output_dir, "standard_metagene_qc.csv"), row.names = FALSE)

message("Completed standard metagene plotting: ", output_dir)
