# Direct four-class posterior and ML screening

## Label construction

For root state `R`, strategy state `S`, technology state `T` and evidence `e`:

```text
P(R=r,T=t|e) = sum_s P(R=r,S=s|e) * P(T=t,S=s|e) / P(S=s|e)
```

This yields `TN=(0,0)`, `FN=(1,0)`, `FP=(0,1)` and `TP=(1,1)` without
thresholding root and technology marginal posteriors separately. The saved audit
confirms sum-to-one and marginal reconstruction to numerical precision.

Each class is called with local error `1-P(class)` and cumulative global FDR
`<0.05`. Only unique calls are retained, followed by technology-specific
`Total>10`. No MAP/Viterbi fallback is used.

## ML tasks

| Task | Positive | Negative | Meaning |
|---|---|---|---|
| `calibration_tp_fp` | TP | FP | reliability of a technology-positive call |
| `error_fp_tn` | FP | TN | feature-dependent false-positive mechanism |
| `artifact_tp_tn` | TP | TN | broader positive-versus-negative contrast |

A model runs when both classes contain more than 50 sites (`>=51`). Validation
uses chromosome-blocked five-fold OOF predictions. Training is capped at 100,000
sites per class; every validation fold retains natural prevalence.

## Clean execution route

From the repository root:

```bash
bash analysis/03_direct_posterior_ml/code/run_calibration_min51.sh
```

Set `OUTPUT_ROOT=/path/to/new/output` to override the safe reproduction output
directory. The runner does not overwrite the archived W8 result folders.

The runner executes, in order:

1. `01_generate_four_class_global_fdr.R`
2. `02_run_rf_fourclass.R` for 14 technologies × three tasks
3. `03_summarize_rf_fourclass.R`
4. `04_export_rf_top20.R`

`01b_validate_four_class_posterior.R` is the independent posterior algebra audit.

## Compact outputs

- `results/posterior_audit/`: class counts and probability/calling identities.
- `results/rf_summary/`: task status, fold metrics, performance, OOF QC and
  importance stability.
- `results/feature_exports/`: 26 per-model Top-20 tables and cross-technology
  stable-feature summaries.

Current completed scope: 8 TP-vs-FP, 7 FP-vs-TN and 11 TP-vs-TN models. AUPRC
must be interpreted against prevalence, especially when FP contains only tens or
hundreds of sites.
