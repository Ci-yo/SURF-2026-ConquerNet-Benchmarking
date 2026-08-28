## Frozen configuration for Benchmark_CodeReview2_20260827.
##
## The LOTO/LOFO scripts use the original combined benchmark dataset and are
## unchanged from the finalized Total>10 workflow. The newer four-class and
## figure scripts use the no-cali IVT-matched dataset through their explicit
## command-line arguments; those paths are recorded here for provenance.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Please run this file through an Rscript entry point.")
code_dir <- dirname(normalizePath(
  sub("^--file=", "", script_arg[[1L]]), winslash = "/"
))
project_root <- dirname(code_dir)
result_dir <- file.path(project_root, "results", "reproduced_run")

## Final LOTO/LOFO inputs -------------------------------------------------

xichen_data_dir <- "/home/xichen/ConquerNet/ConquerBN_basic"
benchmark_file <- file.path(
  xichen_data_dir,
  "m6AConquer_hg38_BenchmarkDataSet_Combined_SE.rds"
)
loto_idr_dir <- file.path(xichen_data_dir, "LOTO_IDR")
lofo_idr_dir <- file.path(xichen_data_dir, "LOFO_IDR")

coverage_cutoff <- 10L
posterior_cutoff <- 0.5
em_max_iteration <- 500L
em_convergence_epsilon <- 1e-6
em_seed <- 20260725L

## Current no-cali four-class analysis inputs ----------------------------

no_cali_ivt_matched_benchmark_file <- file.path(
  xichen_data_dir,
  "m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali_IVT_matched.rds"
)
no_cali_full_model_file <- paste0(
  "/home/xiyao/15data/W7/06_no_cali_IVT_matched_full_model/",
  "Global_no_additional_coverage/full_no_cali_global.conquernet.rds"
)
omics_feature_file <- "/home/xiyao/15data/W8/input_data/features/m6Aconquer_omicsFeatures_hg38.rds"
rf_design_file <- "/home/xiyao/15data/W8/input_data/rf_design/tp_fp_tp_tn_metadata.rds"

minimum_class_size <- 51L       # implements the approved rule: class size > 50
number_of_folds <- 5L
number_of_trees <- 300L
maximum_training_samples_per_class <- 100000L
number_of_threads <- 8L
number_of_screening_features <- 185L
rf_seed <- 20260812L

class_code <- c(TN = 0L, FN = 1L, FP = 2L, TP = 3L)

CONFIG <- list(
  project_root = project_root,
  result_dir = result_dir,
  benchmark_file = benchmark_file,
  coverage_cutoff = coverage_cutoff,
  posterior_cutoff = posterior_cutoff,
  em_max_iteration = em_max_iteration,
  em_convergence_epsilon = em_convergence_epsilon,
  em_seed = em_seed,
  no_cali_ivt_matched_benchmark_file = no_cali_ivt_matched_benchmark_file,
  no_cali_full_model_file = no_cali_full_model_file,
  omics_feature_file = omics_feature_file,
  rf_design_file = rf_design_file,
  minimum_class_size = minimum_class_size,
  number_of_folds = number_of_folds,
  number_of_trees = number_of_trees,
  maximum_training_samples_per_class = maximum_training_samples_per_class,
  number_of_threads = number_of_threads,
  number_of_screening_features = number_of_screening_features,
  rf_seed = rf_seed,
  class_code = class_code
)
