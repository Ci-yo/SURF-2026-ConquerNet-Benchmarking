# Data availability and object contracts

Raw and large intermediate objects are not stored in this repository. Their
verified representations and internal provenance paths are documented so that
an authorized collaborator can reproduce the analyses without guessing what a
number in an RDS file means.

## Core benchmark objects

Both the original combined benchmark and the no-calibration IVT-matched
benchmark are R `RangedSummarizedExperiment` objects with 2,788,705 aligned
genomic rows and 14 technology columns. Their row coordinates are `GRanges` and
their four assays are base R double matrices:

| Assay | Statistical meaning |
|---|---|
| `m6A` | Modified/read-signal count |
| `Total` | Total coverage count |
| `m6ASiteProb` | Single-technology BBmix foreground posterior |
| `AdjPvalue` | Binomial-test BH-adjusted P value |

The last two assays are not interchangeable. `m6ASiteProb > 0.5` supplies the
held-out BBmix MAP truth in the final LOTO/LOFO benchmark; `AdjPvalue` is not
used as that truth. The frozen coverage rule is strictly `Total > 10`.

LOTO/LOFO use the original combined benchmark. Direct four-class inference,
feature modelling and IVT comparison use the no-calibration IVT-matched object
and its matching fitted full ConquerNet model. This separation avoids reusing
IVT-derived calibration information when later checking IVT consistency.

## OmicsFeatures and fold design

OmicsFeatures is a `GRanges` object with 2,788,705 ranges and 224 metadata
columns in the same row order as the benchmark universe. The frozen screening
uses 185 predictors after excluding RBP and histone-modification groups. Raw
chromosome, absolute coordinate and strand are not predictors; a separate
named integer chromosome map defines the five blocked folds.

## Direct four-class labels

The frozen `Total > 10` label object is an integer matrix with 2,788,705 rows
and 14 technology columns: `0=TN`, `1=FN`, `2=FP`, `3=TP`, with `NA` for an
uncalled technology-site pair. Labels come from the direct four-state joint
posterior and class-specific cumulative posterior global FDR `<0.05`; there is
no MAP/Viterbi fallback.

## IDR

Each matching IDR result has a numeric/double continuous `rowmin` vector and
logical `idx`/`label` vectors. The final benchmark uses `-rowmin` as the
continuous score and `idx` as its coverage domain. The binary IDR label is not
used as benchmark truth.

## External consistency

MeRIP IVT and writer/eraser perturbation inputs require an explicit
sample/contrast manifest and coordinate alignment to the benchmark universe.
They are independent consistency checks, not the LOTO/LOFO label source and
not assumed to be perfect biological truth. The repository contains only
compact aggregate results and audit tables.

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

See [`data_catalogue/OBJECT_CONTRACTS.md`](data_catalogue/OBJECT_CONTRACTS.md)
for the complete type and meaning contract and
[`data_catalogue/FINAL_DATA_PATHS.csv`](data_catalogue/FINAL_DATA_PATHS.csv) for
the frozen paths. `INITIAL_DATA_PATHS.csv` is a historical inventory and may
contain superseded branches. Internal paths are provenance notes, not public
download links.
