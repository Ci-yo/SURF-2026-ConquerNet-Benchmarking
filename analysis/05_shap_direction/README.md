# Exact RF reconstruction and SHAP direction

Permutation importance ranks features but does not indicate direction. This module reconstructs each selected random forest using the frozen fold order, random-number sequence, sampling and hyperparameters, then verifies the reconstructed full validation-fold probabilities against the archived OOF probabilities.

SHAP is computed only after exact reproduction passes.

## Scope

- Six core TP-vs-FP models.
- Five core FP-vs-TN models.
- Five chromosome folds per model.
- Top10 features interpreted per model; the forest itself still uses the full 185-feature set.

## Direction semantics

- TP-vs-FP: positive SHAP increases predicted TP probability; negative SHAP supports FP.
- FP-vs-TN: positive SHAP increases predicted FP probability; negative SHAP supports TN.

These are model-association directions, not causal effects. Correlated GC windows, distance variables and transcript features can share or substitute contribution.

## QC

All 11 models completed. Maximum raw and prior-corrected OOF reproduction error was zero. Of 110 model-feature summaries, 105 had a stable sign in at least four of five folds.
