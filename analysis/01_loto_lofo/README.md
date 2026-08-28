# Final LOTO and LOFO benchmark

This module contains the frozen `Total>10` held-out benchmark described in the
Week10 summary. The input is a `RangedSummarizedExperiment` with 2,788,705
aligned sites and 14 technology columns. Its `m6A`, `Total`, `m6ASiteProb` and
`AdjPvalue` assays are base R double matrices; the final benchmark label uses
only the held-out technology's `m6ASiteProb` and `Total` assays.

## LOTO

Each fold removes one technology, applies `Total>10` to the retained
technology-site observations, retrains ConquerBN on the remaining 13
technologies and outputs

```text
P_LOTO(TrueM6A = 1 | evidence from the remaining technologies).
```

The held-out BBmix MAP label is:

```text
positive: held-out m6ASiteProb > 0.5 and held-out Total > 10
negative: held-out m6ASiteProb <= 0.5 and held-out Total > 10
```

`m6ASiteProb` is the held-out **single-technology BBmix posterior**, not the
ConquerBN `TrueM6A` posterior. The final primary comparison contains four
continuous scores on one common evaluation domain:

- ConquerBN root posterior probability;
- continuous IDR score `-rowmin`;
- Rank Product score;
- RRA score.

IDR, Rank Product and RRA are ranking scores, not calibrated probabilities.
MeanProb and Vote are not part of the frozen Week10 comparison.

## LOFO

Each fold removes one complete technology family and retrains ConquerBN on the
other four families. Every technology inside the held-out family is evaluated
against its own BBmix MAP label under `Total>10`. Technology metrics are then
averaged within family, and the five held-out families are the independent
benchmark units.

The IDR input must be the matching family-out result. Technology-out IDR is not
substituted into LOFO.

## Clean execution

```bash
Rscript analysis/01_loto_lofo/code/run_loto_benchmark.R
Rscript analysis/01_loto_lofo/code/run_lofo_benchmark.R
```

Shared validation, ConquerBN fitting and summary functions are under
`code/common_functions/`. New runs write to
`analysis/01_loto_lofo/results/reproduced_run/` by default, leaving the compact
reference results untouched.

## Main result tables

- `results/LOTO/four_method_IDR_summary.csv`
- `results/LOTO/per_fold_metrics.csv`
- `results/LOTO/evaluation_domain_audit.csv`
- `results/LOFO/four_method_IDR_overall_summary.csv`
- `results/LOFO/per_family_method_summary.csv`
- `results/LOFO/per_technology_metrics.csv`
- `results/QC/`

The LOTO output also includes a three-method no-IDR sensitivity domain. The
reported main LOTO numbers use the four-method IDR-common domain.
