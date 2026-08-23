# Direct four-class posterior and ML screening

For root state `R`, strategy state `S`, technology state `T` and evidence `e`:

```text
P(R=r,T=t|e) = sum_s P(R=r,S=s|e) * P(T=t,S=s|e) / P(S=s|e)
```

This yields `TN=(0,0)`, `FN=(1,0)`, `FP=(0,1)` and `TP=(1,1)` without thresholding root and technology marginal posteriors separately.

Each class is called with local error `1-P(class)` and cumulative global FDR `<0.05`. Only unique calls are retained, followed by technology-specific `Total>10`. No MAP/Viterbi fallback is used.

## ML tasks

| Task | Positive | Negative |
|---|---|---|
| calibration_tp_fp | TP | FP |
| error_fp_tn | FP | TN |
| artifact_tp_tn | TP | TN |

A model runs when both classes contain more than 50 sites. Validation uses chromosome-blocked five-fold OOF predictions. Training is capped at 100,000 sites per class; the held-out fold keeps natural prevalence.

## Outputs

- `results/posterior_audit/`: class counts and probability/calling identities.
- `results/rf_summary/`: task status, fold metrics, performance, OOF QC and importance stability.
- `results/feature_exports/`: per-model Top20 and cross-technology stable features.

AUPRC must be interpreted against prevalence. In particular, a TP-vs-FP model with only tens of FP can have AUPRC near one even with modest minority-class ranking.
