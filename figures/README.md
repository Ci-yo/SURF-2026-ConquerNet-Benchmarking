# README result figures

`readme_results_summary` is generated entirely in R from `source_data/readme_summary_metrics.csv` by `code/plot_readme_summary.R`.

The primary PDF uses R's native Helvetica device so that text size can be audited directly. A Cairo PDF is also exported for platform compatibility.

The three panels report descriptive macro averages:

- LOTO: 14 held-out technologies on the common ConquerBN/MeanProb/IDR/RRA/RankProduct domain;
- LOFO: five independent held-out technology families;
- direct-posterior RF: eight TP-vs-FP, seven FP-vs-TN and eleven TP-vs-TN completed technology models.

The figure does not show calibrated probabilities or uncertainty intervals. Per-unit metrics and QC tables remain in the corresponding `analysis/` modules.
