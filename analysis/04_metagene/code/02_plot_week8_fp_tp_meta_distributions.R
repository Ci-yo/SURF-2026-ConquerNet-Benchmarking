#!/usr/bin/env Rscript

# Week8 B: compare posterior-derived FP and TP site distributions.
#
# This script is intentionally downstream-only. It does not call sites and it
# does not alter the Week8 A labels. It expects the unique four-class matrix
# produced by Week8 A (0=TN, 1=FN, 2=FP, 3=TP), the model-to-universe row
# mapping, and the original OmicsFeatures object.
#
# Primary outputs, for every technology:
#   plots/<technology>_FP_vs_TP_meta_distributions.{png,pdf}
#   binned/<technology>_<feature>_FP_vs_TP_bins.csv
#   all_technologies_FP_vs_TP_binned_distributions.csv
#   fp_tp_meta_distribution_qc.csv
#   input_alignment_qc.csv
#
# The plotted y-axis is the within-class fraction of assigned sites in each
# bin. Consequently, the binned fractions sum to one separately for FP and TP.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
})

parse_cli <- function(args) {
  defaults <- list(
    class_matrix = paste0(
      "/home/xiyao/15data/W8/results/",
      "01_four_class_global_fdr/",
      "four_class_matrix_Total_gt_10.rds"
    ),
    # The current A output is full-universe, so no mapping file is needed.
    # For a future model-row-only matrix, pass --row-mapping=<RDS>.
    row_mapping = "none",
    features = paste0(
      "/home/xiyao/15data/W8/input_data/features/",
      "m6Aconquer_omicsFeatures_hg38.rds"
    ),
    # Must be the exact dataset recorded by the Week7 no-cali IVT-matched
    # full model used to generate the Week8 four-class posterior labels.
    # The calibrated/original Combined_SE has identical coordinates but
    # different assay values and therefore cannot be used for Total>10 QC.
    benchmark = paste0(
      "/home/xichen/ConquerNet/ConquerBN_basic/",
      "m6AConquer_hg38_BenchmarkDataSet_Combined_SE_",
      "no_cali_IVT_matched.rds"
    ),
    output_dir = paste0(
      "/home/xiyao/15data/W8/results/",
      "02_fp_tp_meta_distributions/Total_gt_10_no_cali_IVT_matched"
    ),
    coverage_cutoff = "10",
    bins = "100",
    technologies = "all"
  )

  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Arguments must use --name=value syntax: ", arg)
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    key <- gsub("-", "_", pieces[[1L]], fixed = TRUE)
    if (!key %in% names(defaults)) stop("Unknown argument: --", key)
    defaults[[key]] <- paste(pieces[-1L], collapse = "=")
  }
  defaults
}

stop_missing <- function(paths) {
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing)) {
    stop(
      "Missing required input(s):\n",
      paste(paste0("  ", missing, ": ", paths[missing]), collapse = "\n")
    )
  }
}

extract_class_matrix <- function(object) {
  if (is.matrix(object) || is.data.frame(object)) {
    return(as.matrix(object))
  }
  candidates <- c(
    "unique_class_matrix", "class_matrix", "global_class_matrix",
    "global_fdr_class_matrix", "score"
  )
  for (name in candidates) {
    value <- object[[name]]
    if (is.matrix(value) || is.data.frame(value)) return(as.matrix(value))
  }
  stop(
    "Could not find a class matrix. Expected a matrix or a list element named ",
    paste(candidates, collapse = ", "), "."
  )
}

extract_row_index <- function(object) {
  if (is.atomic(object) && is.null(dim(object))) return(as.integer(object))
  if (!is.list(object)) return(NULL)
  candidates <- c(
    "row_index", "model_row_index", "universe_row_index",
    "training_row_index", "site_index"
  )
  for (name in candidates) {
    value <- object[[name]]
    if (is.atomic(value) && is.null(dim(value))) return(as.integer(value))
  }
  NULL
}

normalize_class_codes <- function(x) {
  original_dimnames <- dimnames(x)
  if (is.factor(x)) x <- as.character(x)

  if (is.character(x)) {
    labels <- toupper(trimws(x))
    map <- c(TN = 0L, FN = 1L, FP = 2L, TP = 3L)
    unexpected <- setdiff(unique(labels[!is.na(labels)]), names(map))
    if (length(unexpected)) {
      stop("Unexpected character class label(s): ", paste(unexpected, collapse = ", "))
    }
    out <- base::matrix(unname(map[labels]), nrow = nrow(x), ncol = ncol(x))
  } else {
    out <- base::matrix(as.integer(x), nrow = nrow(x), ncol = ncol(x))
    unexpected <- setdiff(unique(out[!is.na(out)]), 0:3)
    if (length(unexpected)) {
      stop("Unexpected numeric class code(s): ", paste(unexpected, collapse = ", "))
    }
  }
  dimnames(out) <- original_dimnames
  out
}

extract_feature_table <- function(object) {
  # Keep the S4 DataFrame lazy here. Converting all 224 columns to a base
  # data.frame would make an unnecessary second copy of the ~2.8M-row object.
  if (inherits(object, "GRanges")) return(S4Vectors::mcols(object))
  if (inherits(object, "SummarizedExperiment")) {
    return(SummarizedExperiment::rowData(object))
  }
  if (is.matrix(object) || is.data.frame(object)) return(object)
  stop("Unsupported OmicsFeatures class: ", paste(class(object), collapse = ", "))
}

extract_feature_ranges <- function(object) {
  if (inherits(object, "GRanges")) return(object)
  if (inherits(object, "SummarizedExperiment")) {
    return(SummarizedExperiment::rowRanges(object))
  }
  NULL
}

coordinate_signature_equal <- function(left, right) {
  if (is.null(left) || is.null(right)) return(NA)
  if (length(left) != length(right)) return(FALSE)
  identical(as.character(GenomicRanges::seqnames(left)),
            as.character(GenomicRanges::seqnames(right))) &&
    identical(GenomicRanges::start(left), GenomicRanges::start(right)) &&
    identical(GenomicRanges::end(left), GenomicRanges::end(right)) &&
    identical(as.character(GenomicRanges::strand(left)),
              as.character(GenomicRanges::strand(right)))
}

feature_assignment_mask <- function(feature_name, feature_table, row_index) {
  n <- length(row_index)
  if (feature_name == "meta_tx_topology") {
    annotation_columns <- c(
      "overlap_exonicFivePrimeUTR",
      "overlap_exonicCDS",
      "overlap_exonicThreePrimeUTR"
    )
  } else if (feature_name == "relativePOS_fullGenes") {
    annotation_columns <- "overlap_fullGenes"
  } else {
    annotation_columns <- character()
  }

  available <- intersect(annotation_columns, colnames(feature_table))
  if (!length(available)) {
    return(list(
      mask = rep(TRUE, n),
      method = "finite_and_in_[0,1]_only; explicit_overlap_column_unavailable"
    ))
  }

  annotation <- as.data.frame(feature_table[row_index, available, drop = FALSE])
  annotation[] <- lapply(annotation, function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[!is.finite(x)] <- 0
    x
  })
  list(
    mask = rowSums(as.matrix(annotation) > 0) > 0,
    method = paste0("explicit overlap mask: ", paste(available, collapse = " OR "))
  )
}

bin_distribution <- function(values, class_name, feature_name, technology,
                             n_bins) {
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  index <- findInterval(
    values,
    breaks,
    all.inside = TRUE,
    rightmost.closed = TRUE
  )
  counts <- tabulate(index, nbins = n_bins)
  denominator <- sum(counts)
  data.frame(
    technology = technology,
    feature = feature_name,
    class = class_name,
    bin = seq_len(n_bins),
    bin_start = head(breaks, -1L),
    bin_end = tail(breaks, -1L),
    bin_midpoint = (head(breaks, -1L) + tail(breaks, -1L)) / 2,
    count = counts,
    within_class_fraction = if (denominator > 0L) counts / denominator else NA_real_,
    n_assigned_and_binned = denominator,
    stringsAsFactors = FALSE
  )
}

safe_file_name <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

write_plot <- function(data, technology, output_base) {
  display_names <- c(
    meta_tx_topology = "Meta-transcript topology",
    relativePOS_fullGenes = "Relative position across full gene"
  )
  data$feature_display <- unname(display_names[data$feature])
  data$class <- factor(data$class, levels = c("TP", "FP"))

  plot <- ggplot(
    data,
    aes(x = bin_midpoint, y = within_class_fraction, colour = class)
  ) +
    geom_line(linewidth = 0.75, na.rm = TRUE) +
    facet_wrap(~feature_display, ncol = 1L, scales = "free_y") +
    scale_colour_manual(values = c(TP = "#C43C39", FP = "#2C6FB7")) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      title = paste0(technology, ": posterior-derived FP vs TP"),
      subtitle = "Total > 10; unique global-FDR calls; normalized within class",
      x = "Normalized position",
      y = "Fraction of assigned class sites per bin",
      colour = "Class"
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )

  ggsave(paste0(output_base, ".png"), plot, width = 7.2, height = 7.6,
         units = "in", dpi = 300, bg = "white")
  ggsave(paste0(output_base, ".pdf"), plot, width = 7.2, height = 7.6,
         units = "in", device = grDevices::cairo_pdf, bg = "white")
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
coverage_cutoff <- as.numeric(args$coverage_cutoff)
n_bins <- as.integer(args$bins)
if (!is.finite(coverage_cutoff) || coverage_cutoff < 0) {
  stop("coverage_cutoff must be a non-negative finite number.")
}
if (is.na(n_bins) || n_bins < 10L) stop("bins must be an integer >= 10.")

required_paths <- c(class_matrix = args$class_matrix, features = args$features)
stop_missing(required_paths)

if (!requireNamespace("S4Vectors", quietly = TRUE) ||
    !requireNamespace("GenomicRanges", quietly = TRUE)) {
  stop("S4Vectors and GenomicRanges are required.")
}

class_object <- readRDS(args$class_matrix)
feature_object <- readRDS(args$features)

class_matrix <- normalize_class_codes(extract_class_matrix(class_object))
row_index_from_class <- extract_row_index(class_object)
mapping_is_supplied <- !args$row_mapping %in% c("", "none", "NA", "auto")
if (mapping_is_supplied && !file.exists(args$row_mapping)) {
  stop("The requested row mapping does not exist: ", args$row_mapping)
}
row_index_from_mapping <- if (mapping_is_supplied) {
  extract_row_index(readRDS(args$row_mapping))
} else {
  NULL
}

if (!is.null(row_index_from_class) && !is.null(row_index_from_mapping) &&
    !identical(row_index_from_class, row_index_from_mapping)) {
  stop("row_index in the class object differs from the independent row mapping.")
}
row_index <- if (!is.null(row_index_from_class)) {
  row_index_from_class
} else {
  row_index_from_mapping
}

feature_table <- extract_feature_table(feature_object)
n_universe <- nrow(feature_table)

if (nrow(class_matrix) == n_universe) {
  if (!is.null(row_index) && length(row_index) == n_universe &&
      !identical(row_index, seq_len(n_universe))) {
    stop(
      "The class matrix has full-universe rows, but the supplied row mapping is ",
      "not the full ordered sequence. Refusing ambiguous alignment."
    )
  }
  row_index <- seq_len(n_universe)
  row_mapping_mode <- "full_universe_matrix"
} else if (!is.null(row_index) && nrow(class_matrix) == length(row_index)) {
  row_mapping_mode <- "model_rows_mapped_to_feature_universe"
} else {
  stop(
    "Class matrix rows (", nrow(class_matrix), ") equal neither feature rows (",
    n_universe, ") nor a supplied row_index length (",
    if (is.null(row_index)) "none" else length(row_index), ")."
  )
}

if (anyNA(row_index) || any(row_index < 1L) || any(row_index > n_universe) ||
    anyDuplicated(row_index)) {
  stop("row_index must be unique, non-missing, and within the feature universe.")
}
if (is.null(colnames(class_matrix)) || any(!nzchar(colnames(class_matrix))) ||
    anyDuplicated(colnames(class_matrix))) {
  stop("The class matrix must have unique, non-empty technology column names.")
}

technologies <- colnames(class_matrix)
if (!identical(args$technologies, "all")) {
  requested <- trimws(strsplit(args$technologies, ",", fixed = TRUE)[[1L]])
  missing <- setdiff(requested, technologies)
  if (length(missing)) stop("Unknown technology/technologies: ", paste(missing, collapse = ", "))
  technologies <- requested
}

feature_names <- c("meta_tx_topology", "relativePOS_fullGenes")
missing_features <- setdiff(feature_names, colnames(feature_table))
if (length(missing_features)) {
  stop("Missing required OmicsFeature(s): ", paste(missing_features, collapse = ", "))
}

# Recheck Total > cutoff whenever the benchmark is available. The Week8 A
# matrix should already be cutoff-matched; this independent check prevents an
# incorrectly named or unfiltered matrix from silently entering the plots.
coverage_matrix <- NULL
coverage_verification <- "not_rechecked; benchmark path missing"
coordinate_match <- NA
if (file.exists(args$benchmark)) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment is required to verify the coverage mask.")
  }
  benchmark <- readRDS(args$benchmark)
  if (!inherits(benchmark, "SummarizedExperiment")) {
    stop("The benchmark input is not a SummarizedExperiment.")
  }
  if (!"Total" %in% SummarizedExperiment::assayNames(benchmark)) {
    stop("The benchmark input has no Total assay.")
  }
  total <- SummarizedExperiment::assay(benchmark, "Total")
  if (nrow(total) != n_universe) {
    stop("Benchmark rows and OmicsFeatures rows differ.")
  }
  if (!all(technologies %in% colnames(total))) {
    stop(
      "Benchmark Total assay lacks technology/technologies: ",
      paste(setdiff(technologies, colnames(total)), collapse = ", ")
    )
  }
  coverage_matrix <- total[row_index, technologies, drop = FALSE]
  rm(total)

  coordinate_match <- coordinate_signature_equal(
    extract_feature_ranges(feature_object),
    SummarizedExperiment::rowRanges(benchmark)
  )
  if (identical(coordinate_match, FALSE)) {
    stop("Benchmark and OmicsFeatures genomic coordinates are not identical.")
  }
  coverage_verification <- paste0("independently_rechecked_Total>", coverage_cutoff)
  rm(benchmark)
  invisible(gc())
}

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir <- file.path(args$output_dir, "plots")
bin_dir <- file.path(args$output_dir, "binned")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bin_dir, recursive = TRUE, showWarnings = FALSE)

# Feature values and annotation masks are technology-independent. Compute them
# once and reuse them for all 14 technologies.
feature_inputs <- setNames(lapply(feature_names, function(feature_name) {
  values <- suppressWarnings(as.numeric(feature_table[[feature_name]][row_index]))
  finite <- is.finite(values)
  in_range <- finite & values >= -1e-12 & values <= 1 + 1e-12
  values[in_range] <- pmin(1, pmax(0, values[in_range]))
  assignment <- feature_assignment_mask(feature_name, feature_table, row_index)
  list(
    values = values,
    finite = finite,
    in_range = in_range,
    assignment_mask = assignment$mask,
    assignment_method = assignment$method,
    usable = in_range & assignment$mask
  )
}), feature_names)
rm(feature_table, feature_object)
invisible(gc())

all_binned <- list()
all_qc <- list()
result_index <- 0L
qc_index <- 0L

for (technology in technologies) {
  message("Processing ", technology)
  tech_column <- match(technology, colnames(class_matrix))
  class_vector <- class_matrix[, tech_column]
  class_label_mask <- !is.na(class_vector) & class_vector %in% c(2L, 3L)

  if (is.null(coverage_matrix)) {
    coverage_mask <- rep(TRUE, length(class_vector))
  } else {
    coverage <- suppressWarnings(as.numeric(coverage_matrix[, technology]))
    coverage_mask <- is.finite(coverage) & coverage > coverage_cutoff
  }

  fp_base <- class_label_mask & class_vector == 2L & coverage_mask
  tp_base <- class_label_mask & class_vector == 3L & coverage_mask
  n_fp_base <- sum(fp_base)
  n_tp_base <- sum(tp_base)

  technology_binned <- list()
  technology_bin_index <- 0L

  for (feature_name in feature_names) {
    feature_input <- feature_inputs[[feature_name]]
    values <- feature_input$values
    finite <- feature_input$finite
    in_range <- feature_input$in_range
    usable_feature <- feature_input$usable

    fp_mask <- fp_base & usable_feature
    tp_mask <- tp_base & usable_feature

    fp_bins <- bin_distribution(
      values[fp_mask], "FP", feature_name, technology, n_bins
    )
    tp_bins <- bin_distribution(
      values[tp_mask], "TP", feature_name, technology, n_bins
    )
    binned <- rbind(fp_bins, tp_bins)

    result_index <- result_index + 1L
    all_binned[[result_index]] <- binned
    technology_bin_index <- technology_bin_index + 1L
    technology_binned[[technology_bin_index]] <- binned

    file_stem <- paste0(
      safe_file_name(technology), "_", feature_name, "_FP_vs_TP_bins.csv"
    )
    write.csv(binned, file.path(bin_dir, file_stem), row.names = FALSE)

    for (class_name in c("FP", "TP")) {
      base_mask <- if (class_name == "FP") fp_base else tp_base
      final_mask <- if (class_name == "FP") fp_mask else tp_mask
      denominator <- sum(base_mask)
      qc_index <- qc_index + 1L
      all_qc[[qc_index]] <- data.frame(
        technology = technology,
        class = class_name,
        feature = feature_name,
        coverage_cutoff = coverage_cutoff,
        class_code = if (class_name == "FP") 2L else 3L,
        n_class_before_coverage_recheck = sum(
          class_label_mask & class_vector == if (class_name == "FP") 2L else 3L
        ),
        n_class_after_coverage_recheck = denominator,
        n_feature_assigned_and_in_range = sum(final_mask),
        feature_assignment_rate = if (denominator > 0L) {
          sum(final_mask) / denominator
        } else {
          NA_real_
        },
        n_feature_missing_or_nonfinite = sum(base_mask & !finite),
        n_feature_outside_0_1 = sum(base_mask & finite & !in_range),
        n_not_in_explicit_annotation_mask = sum(
          base_mask & in_range & !feature_input$assignment_mask
        ),
        assignment_mask_method = feature_input$assignment_method,
        stringsAsFactors = FALSE
      )
    }
  }

  technology_plot_data <- do.call(rbind, technology_binned)
  output_base <- file.path(
    plot_dir,
    paste0(safe_file_name(technology), "_FP_vs_TP_meta_distributions")
  )
  write_plot(technology_plot_data, technology, output_base)
}

binned_table <- do.call(rbind, all_binned)
qc_table <- do.call(rbind, all_qc)

write.csv(
  binned_table,
  file.path(args$output_dir, "all_technologies_FP_vs_TP_binned_distributions.csv"),
  row.names = FALSE
)
write.csv(
  qc_table,
  file.path(args$output_dir, "fp_tp_meta_distribution_qc.csv"),
  row.names = FALSE
)

fraction_check <- aggregate(
  within_class_fraction ~ technology + feature + class,
  data = binned_table,
  FUN = function(x) sum(x, na.rm = TRUE)
)
names(fraction_check)[names(fraction_check) == "within_class_fraction"] <-
  "sum_within_class_fraction"
fraction_check$expected_sum <- ifelse(
  is.finite(fraction_check$sum_within_class_fraction) &
    fraction_check$sum_within_class_fraction > 0,
  1,
  NA_real_
)
fraction_check$absolute_sum_error <- abs(
  fraction_check$sum_within_class_fraction - fraction_check$expected_sum
)
write.csv(
  fraction_check,
  file.path(args$output_dir, "binned_fraction_sum_audit.csv"),
  row.names = FALSE
)

class_metadata <- if (is.list(class_object)) {
  c(
    coverage_cutoff = if (!is.null(class_object$coverage_cutoff)) {
      as.character(class_object$coverage_cutoff)
    } else NA_character_,
    criterion = if (!is.null(class_object$criterion)) {
      as.character(class_object$criterion)
    } else NA_character_
  )
} else {
  c(coverage_cutoff = NA_character_, criterion = NA_character_)
}

alignment_qc <- data.frame(
  item = c(
    "class_matrix_path", "row_mapping_path", "features_path", "benchmark_path",
    "class_matrix_rows", "class_matrix_columns", "feature_universe_rows",
    "row_index_length", "row_index_min", "row_index_max",
    "row_index_unique", "row_mapping_mode", "coordinate_identity",
    "coverage_verification", "requested_coverage_cutoff",
    "class_object_coverage_cutoff", "class_object_criterion",
    "class_code_contract", "normalization_contract", "n_bins"
  ),
  value = c(
    normalizePath(args$class_matrix, winslash = "/", mustWork = TRUE),
    if (mapping_is_supplied) {
      normalizePath(args$row_mapping, winslash = "/", mustWork = TRUE)
    } else {
      "not_required_for_full_universe_matrix"
    },
    normalizePath(args$features, winslash = "/", mustWork = TRUE),
    if (file.exists(args$benchmark)) {
      normalizePath(args$benchmark, winslash = "/", mustWork = TRUE)
    } else args$benchmark,
    nrow(class_matrix), ncol(class_matrix), n_universe, length(row_index),
    min(row_index), max(row_index), !anyDuplicated(row_index), row_mapping_mode,
    if (is.na(coordinate_match)) "not_available" else coordinate_match,
    coverage_verification, coverage_cutoff,
    class_metadata[["coverage_cutoff"]], class_metadata[["criterion"]],
    "0=TN;1=FN;2=FP;3=TP;NA=unclassified",
    "within each technology x feature x class, bin fractions sum to one",
    n_bins
  ),
  stringsAsFactors = FALSE
)
write.csv(
  alignment_qc,
  file.path(args$output_dir, "input_alignment_qc.csv"),
  row.names = FALSE
)

if (any(qc_table$n_class_after_coverage_recheck == 0L)) {
  warning("At least one technology/class has zero sites after coverage filtering.")
}
nonempty_fraction_error <- fraction_check$absolute_sum_error[
  is.finite(fraction_check$absolute_sum_error)
]
if (length(nonempty_fraction_error) && max(nonempty_fraction_error) > 1e-12) {
  stop("Binned within-class fractions failed their sum-to-one audit.")
}

cat("\nWeek8 B outputs saved to:\n  ", args$output_dir, "\n", sep = "")
cat("Technologies: ", length(technologies), "\n", sep = "")
cat("Bins per feature: ", n_bins, "\n", sep = "")
cat("Coverage verification: ", coverage_verification, "\n", sep = "")
