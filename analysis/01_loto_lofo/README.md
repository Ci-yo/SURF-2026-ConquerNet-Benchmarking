# Final LOTO and LOFO benchmark

## LOTO

Each fold removes one technology completely, masks retained observations with `Total <= 10`, rebuilds the 13-technology DAG/CPT initialization and reruns EM. The score is the root posterior `P_LOTO(TrueM6A=1)`.

Held-out truth is the held-out technology BBmix MAP label:

```text
positive: m6ASiteProb > 0.5 and Total > 10
negative: m6ASiteProb <= 0.5 and Total > 10
```

ConquerBN, MeanProb, RankProduct, RRA and matched continuous IDR are evaluated on the same site domain.

## LOFO

Each model removes one complete technology family and trains on the remaining four families. Technology metrics are first averaged within the held-out family. Statistical comparison uses five independent family-level values.

## Key files

- `code/01_evaluate_loto_map_coverage.R`: retraining and MAP evaluation.
- `code/02_run_loso.R`: historical filename for LOFO model training.
- `code/03_evaluate_loso_map_coverage.R`: historical filename for LOFO evaluation.
- `code/12_evaluate_loto_cutoff_matched_idr.R`: five-method LOTO common-domain evaluation.
- `code/13_evaluate_lofo_cutoff_matched_idr.R`: five-method LOFO common-domain evaluation.
- `results/LOTO/loto_five_method_cutoff_matched_summary.csv`: main LOTO ranking.
- `results/LOFO/lofo_five_method_overall_summary.csv`: main five-family ranking.
- `results/QC/`: convergence and model audits.

The `loso` text in legacy filenames means the finalized leave-one-family-out analysis and is described as LOFO in reports.
