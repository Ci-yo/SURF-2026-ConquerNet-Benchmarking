# Verified R object contracts

This page records the concrete R representations used by the frozen Week10
workflow. The metadata below were checked directly with `readRDS()` on
BioXplore on 2026-08-27. They are object contracts, not public download links.

## 1. Original combined benchmark

```text
/home/xichen/ConquerNet/ConquerBN_basic/
  m6AConquer_hg38_BenchmarkDataSet_Combined_SE.rds
```

| Property | Verified value |
|---|---|
| R class | `RangedSummarizedExperiment` |
| Dimensions | 2,788,705 genomic sites × 14 technologies |
| Row coordinates | `GRanges` |
| Assays | `m6A`, `Total`, `m6ASiteProb`, `AdjPvalue` |
| Assay representation | Base R `matrix`/`array`, storage mode `double` |

The assays have different meanings even though all four are stored as numeric
matrices:

| Assay | Meaning |
|---|---|
| `m6A` | Technology-specific modified/read-signal count used by the emission model |
| `Total` | Technology-specific total coverage count |
| `m6ASiteProb` | Single-technology BBmix posterior probability of the foreground state |
| `AdjPvalue` | Binomial-test Benjamini-Hochberg-adjusted P value |

`m6ASiteProb` is **not** the ConquerNet root posterior. It supplies the
held-out BBmix MAP label in LOTO/LOFO: positive when `m6ASiteProb > 0.5` and
negative otherwise. The frozen coverage condition is `Total > 10`—strictly
greater than 10, not `Total >= 10` and not “remove values greater than 10.”

The verified object has no missing entries in `m6A` or `Total`.
`m6ASiteProb` and `AdjPvalue` each contain 23,730,999 `NA` entries, reflecting
technology-site combinations without a defined single-technology result.

## 2. LOTO and LOFO model outputs

The score saved for ConquerBN is a numeric vector in `[0,1]`:

```text
P_holdout(TrueM6A = 1 | retained technology evidence)
```

It is a model posterior probability. In the evaluation tables it is called a
`score` only because AUROC/AUPRC functions accept a common continuous-score
column for all methods.

The other final comparators are not probabilities:

| Comparator | Stored representation | Evaluation score |
|---|---|---|
| Continuous IDR | `rowmin`: numeric/double vector; `idx` and `label`: logical vectors | `-rowmin` |
| RRA | Numeric rank-derived vector | Higher values indicate stronger positive ranking evidence |
| Rank Product | Numeric rank-derived vector | Higher values indicate stronger positive ranking evidence |

The LOTO IDR directory contains 84 RDS files: 14 held-out technologies × two
historical coverage cutoffs × `rowmin`/`idx`/`label`. The LOFO IDR directory
contains 30 files: five held-out families × two cutoffs × three file types.
Only the matching `Total > 10` continuous `rowmin` score enters the frozen
four-method comparison. `idx` defines IDR coverage; the binary `label` is not
used as benchmark truth.

## 3. No-calibration IVT-matched branch

```text
/home/xichen/ConquerNet/ConquerBN_basic/
  m6AConquer_hg38_BenchmarkDataSet_Combined_SE_no_cali_IVT_matched.rds
```

This is also a `RangedSummarizedExperiment` with 2,788,705 rows, 14 columns,
`GRanges` row coordinates and the same four double-matrix assays. It is paired
with:

```text
/home/xiyao/15data/W7/06_no_cali_IVT_matched_full_model/
  Global_no_additional_coverage/full_no_cali_global.conquernet.rds
```

The fitted model is a named R `list` with 28 components, including the fitted
CPT parameters, beta-binomial parameters, DAG, inference row indices and input
dataset provenance. This branch is used for direct TN/FN/FP/TP inference and
the downstream IVT comparison. It is intentionally separate from the original
combined benchmark so IVT-derived calibration information is not reused in an
analysis that later compares labels with IVT evidence.

## 4. Direct four-class posterior and calls

For each technology-site pair, the updated inference returns four posterior
probabilities:

```text
P(TN | evidence)
P(FN | evidence)
P(FP | evidence)
P(TP | evidence)
```

They sum to one, subject only to floating-point error. Classes are called
separately by sorting local error `1 - P(class | evidence)` and retaining the
largest prefix whose cumulative mean error is strictly below `0.05`. A final
label is assigned only when exactly one class passes; the technology-specific
`Total > 10` requirement is then applied.

The frozen class matrix is:

```text
/home/xiyao/15data/W8/results/01_four_class_global_fdr/
  four_class_matrix_Total_gt_10.rds
```

| Property | Verified value |
|---|---|
| R type | Integer matrix |
| Dimensions | 2,788,705 × 14 |
| Codes | `0=TN`, `1=FN`, `2=FP`, `3=TP` |
| Uncalled sites | Missing (`NA`) |

This is a direct joint-posterior classification. It does not use the previous
marginal root/technology threshold combination, a Viterbi/MAP fallback, or a
1,000-per-class fallback rule.

## 5. Feature and fold objects

### OmicsFeatures

```text
/home/xiyao/15data/W8/input_data/features/
  m6Aconquer_omicsFeatures_hg38.rds
```

| Property | Verified value |
|---|---|
| R class | `GRanges` |
| Length | 2,788,705 ranges |
| Metadata columns | 224 |
| Screening predictors | 185 |

The 185-predictor set excludes RBP and histone-modification feature groups.
Raw chromosome, absolute coordinate and strand fields are not predictors.
Transcript-relative quantities such as exon length, distance to a junction and
meta-transcript topology are biological annotation features, not raw genomic
identifiers.

### RF design metadata

```text
/home/xiyao/15data/W8/input_data/rf_design/
  tp_fp_tp_tn_metadata.rds
```

This is a named R `list` with 27 components. It includes the 185 predictor
names, a 224-row feature manifest, class coding (`TN=0`, `FN=1`, `FP=2`,
`TP=3`) and a named integer chromosome-to-fold map for the 25 chromosome
labels. The map assigns chromosomes to five blocked folds; chromosome is not
fed to the random forest as a predictor.

## 6. RF output semantics

Each eligible technology-task model uses chromosome-blocked five-fold
out-of-fold prediction. Eligibility requires both compared classes to contain
more than 50 sites (`>=51`). The saved OOF object contains raw and
prevalence-corrected positive-class predictions. The reported AUROC/AUPRC uses
the prevalence-corrected OOF score.

These values are classifier scores used for discrimination assessment. They
are not a completed formal probability-calibration analysis: no independent
Platt/isotonic stage has yet established that, for example, a score of `0.8`
corresponds to an 80% empirical event rate.

## 7. Alignment invariants

Before combining any objects, the workflow requires:

1. identical site-universe length;
2. identical chromosome, coordinate and strand ordering where available;
3. identical technology names and intended column order;
4. held-out labels and all compared scores restricted to one common finite
   evaluation domain;
5. `Total > 10` applied to the correct technology-site observation.

The exact source and generated-object paths are listed in
[`FINAL_DATA_PATHS.csv`](FINAL_DATA_PATHS.csv).
