# Reproducibility notes

## Execution order

1. Run the LOTO models and evaluation in `analysis/01_loto_lofo/code/`.
2. Run the five LOFO models, then evaluate at technology level and summarize at family level.
3. Fit or load the matching no-cali full model for external validation.
4. Generate direct four-class posteriors and unique global-FDR calls.
5. Run chromosome-blocked RF screening, summarize OOF results and export Top20 features.
6. Generate metagene distributions and effect-size tables.
7. Reconstruct the selected RF models exactly before computing SHAP directions.

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
```

## Path configuration

Several scripts preserve the exact server defaults used to generate the archived results. Override their CLI paths when running elsewhere. The direct-posterior scripts now default to the included `EM_function_2.R` helper in their own code directory.

## Integrity checks

- LOTO/LOFO methods are compared on common site domains.
- LOFO inferential units are five families, not 14 technologies.
- Four-class posteriors sum to one and recover the relevant marginal posteriors.
- Unique global-FDR calls contain no class overlaps in the formal run.
- Every completed RF model has one OOF prediction per eligible site, with no duplicate row index.
- SHAP analysis proceeds only when exact full-fold OOF probability reproduction passes.
