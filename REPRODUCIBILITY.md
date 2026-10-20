# Reproducibility notes

Start with [the step-by-step guide](docs/GETTING_STARTED.md) for environment
setup, the public archive check, expected outputs and troubleshooting.
The execution order below assumes the external research inputs are already
available and all server paths have been configured.

## Execution order

1. Run `Rscript analysis/01_loto_lofo/code/run_loto_benchmark.R`.
2. Run `Rscript analysis/01_loto_lofo/code/run_lofo_benchmark.R`; evaluate at
   technology level and summarize at the five-family level.
3. Fit or load the matching no-cali full model for external validation.
4. Run `analysis/03_direct_posterior_ml/code/run_calibration_min51.sh` to
   generate direct four-class posteriors, unique global-FDR calls,
   chromosome-blocked RF screening, OOF summaries and Top-20 exports.
5. Run `analysis/04_metagene/code/run_B_min51.sh` to generate FP/TP effect sizes
   and the final eight-technology standard-metagene figures.
6. Reconstruct the selected RF models exactly before computing SHAP directions.

## Frozen constants

```text
coverage: Total > 10
BBmix MAP: m6ASiteProb > 0.5 vs <= 0.5
four-class global FDR: < 0.05
ML eligibility: both task classes > 50 (integer minimum 51)
RF folds: chromosome-blocked 5-fold
RF trees: 300
training cap: 100,000 per class
Vote: excluded
MeanProb: excluded from the frozen Week10 benchmark comparison
```

## Path configuration

Several scripts preserve the exact server defaults used to generate the archived results. Override their CLI paths when running elsewhere. The direct-posterior scripts now default to the included `EM_function_2.R` helper in their own code directory.

## Integrity checks

- The primary LOTO/LOFO comparison contains ConquerBN, continuous IDR, RRA and
  Rank Product on common site domains.
- LOFO inferential units are five families, not 14 technologies.
- Four-class posteriors sum to one and recover the relevant marginal posteriors.
- Unique global-FDR calls contain no class overlaps in the formal run.
- Every completed RF model has one OOF prediction per eligible site, with no duplicate row index.
- SHAP analysis proceeds only when exact full-fold OOF probability reproduction passes.
