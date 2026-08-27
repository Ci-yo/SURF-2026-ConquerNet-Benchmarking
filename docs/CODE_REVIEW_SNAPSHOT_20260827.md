# Code-review snapshot — 2026-08-27

This snapshot freezes the current public route after the direct four-class
posterior and `>50` class-size updates.

## Preserved benchmark branch

- LOTO: 14 held-out technologies.
- LOFO: five independent held-out technology families.
- `Total>10`, held-out BBmix MAP truth, matching continuous IDR and common-site
  method comparisons.
- Vote and Total>20 sensitivity outputs remain excluded.

## Updated downstream branch

- Direct joint posterior code: TN=0, FN=1, FP=2, TP=3.
- Class-specific cumulative global FDR `<0.05` and exactly one selected class.
- Technology-specific `Total>10`; no MAP/Viterbi fallback.
- ML eligibility: both compared classes `>50` (integer minimum 51).
- 26 completed technology-task RF models; all 26 OOF audits passed.
- 26 Top-20 tables, totalling 520 model-feature rows.
- Eight eligible FP/TP standard-metagene panels; stale SingleMod panel removed.

## Numerical checks

| Audit | Result |
|---|---:|
| Maximum four-class sum-to-one error | 7.32747196252603e-15 |
| Maximum reconstructed TrueM6A marginal error | 1.69864122767649e-14 |
| Maximum reconstructed Technology marginal error | 3.33066907387547e-16 |
| Final overlapping class calls | 0 |
| LOTO technologies | 14 |
| LOFO independent families | 5 |
| Completed RF models / OOF passes | 26 / 26 |
| Metagene technologies | 8 |

## Release contents

Only small source code, audit tables, aggregate results, publication-quality
figures and the SURF poster are tracked. Large raw objects, fitted models,
per-site posterior matrices and OOF predictions are omitted; exact internal
paths are recorded in `data_catalogue/FINAL_DATA_PATHS.csv`.
