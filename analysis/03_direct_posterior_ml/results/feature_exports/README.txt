Week8 RF feature exports
==========================

Source summary: /home/xiyao/15data/W8/results/04_rf_four_class_summary_min51
Eligibility: identifiability_status == reliable AND run_status == complete.
Exploratory and untrainable tasks are excluded from every Top20/shared summary.

Per-model ranking order:
  1) top20_fold_fraction descending
  2) median_rank ascending
  3) median_normalized_importance descending
  4) positive_fold_fraction descending
  5) entered_fold_fraction descending
  6) feature name

A stable shared occurrence requires stable_4of5_top20 == TRUE within a model:
at least 4/5 folds with positive importance and at least 4/5 folds in Top20.
Shared denominators are the reliable completed technologies within each task.
These screening Top20 lists are not yet the final de-correlated reduced feature sets.

Reliable completed models exported: 26
Combined per-model Top20 rows: 520
Tasks: calibration_tp_fp, error_fp_tn, artifact_tp_tn
