# Data availability and object contracts

Raw and intermediate data are not stored in this repository.

## Core benchmark object

Expected format: an R `SummarizedExperiment` with aligned rows and technology columns. Required assays include:

```text
m6A
Total
m6ASiteProb
AdjPvalue
```

The formal LOTO/LOFO benchmark uses the agreed benchmark dataset with `Total > 10`. The Week7/Week9 external-validation and direct-posterior branches use the no-IVT-calibration, IVT-matched benchmark and its matching fitted full model.

## OmicsFeatures

Expected format: an R `GRanges` object with exactly the same coordinates and row order as the benchmark universe. The frozen screening set contains 185 predictors after removing RBP and histone-modification features. Raw chromosome, coordinate and strand fields are not predictors; chromosome is used only to define blocked folds.

## IDR

Continuous IDR outputs must be aligned to the same site universe and generated from the matching LOTO or LOFO exclusion design. Binary IDR labels are not used as ground truth.

## External validation

MeRIP IVT and writer/eraser perturbation inputs require an explicit sample/contrast manifest and coordinate alignment to the benchmark universe. The repository contains only compact aggregate results and audit tables.

## Not included

```text
*.rds / *.RDS raw inputs
fitted ConquerBN model objects
per-site four-class posterior objects
OOF prediction objects
normalized rank matrices
full IDR matrices
simulation objects without sampled site-level truth
```

`data_catalogue/FINAL_DATA_PATHS.csv` records the frozen inputs and outputs used
by the current public snapshot. `INITIAL_DATA_PATHS.csv` is retained as a
historical inventory and may include superseded branches. Internal paths are
provenance notes, not public download links.
