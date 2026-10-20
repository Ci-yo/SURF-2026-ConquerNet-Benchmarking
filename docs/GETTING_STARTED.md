# Getting started

This repository publishes compact results and analysis source. Choose the
public walkthrough if you do not have the original research inputs.

## 1. Download the repository

Install Git and Python 3.10 or later. In a terminal:

```bash
git clone https://github.com/Ci-yo/SURF-2026-ConquerNet-Benchmarking.git
cd SURF-2026-ConquerNet-Benchmarking
python --version
```

Alternatively, use Code → Download ZIP and extract the entire archive. Run
subsequent commands from the folder containing `README.md`.

## 2. Verify the public snapshot

```bash
python scripts/check_public_snapshot.py
```

There are no Python package dependencies. Success reports the checked manifest
file count and 11 summary rows, then prints the LOTO, LOFO and RF summary
metrics. This checks file integrity and checks the eight benchmark summary
rows against their retained source tables. It does not refit any model or
reconstruct per-site predictions.

`MANIFEST_SHA256.tsv` records SHA-256 and byte size. Text files are compared
with LF line endings so Windows checkout conversion does not create false
failures; binary files are compared byte-for-byte. A deliberate source edit or
plot rerun can change a checksum. Preserve the published manifest for comparison;
do not dismiss an unexplained mismatch.

## 3. Read the results in the correct order

| Question | File to open |
|---|---|
| What are the headline metrics? | [`figures/source_data/readme_summary_metrics.csv`](../figures/source_data/readme_summary_metrics.csv) |
| How did each held-out technology perform? | [`analysis/01_loto_lofo/results/LOTO/per_fold_metrics.csv`](../analysis/01_loto_lofo/results/LOTO/per_fold_metrics.csv) |
| What are the five independent LOFO units? | [`analysis/01_loto_lofo/results/LOFO/per_family_method_summary.csv`](../analysis/01_loto_lofo/results/LOFO/per_family_method_summary.csv) |
| Which RF tasks completed? | [`analysis/03_direct_posterior_ml/results/rf_summary/`](../analysis/03_direct_posterior_ml/results/rf_summary/) |
| Which external checks disagreed with the internal labels? | [`analysis/02_external_validation/README.md`](../analysis/02_external_validation/README.md) |

The archived ConquerBN summary is approximately `0.7873 / 0.6785`
(LOTO AUROC/AUPRC) and `0.7679 / 0.6638` (LOFO). LOFO is a macro-average
over **five families**, not 14 independent technology observations.
Internal TP/FP/TN/FN labels are posterior-derived model states, not
experimentally established biological truth.

## 4. Optional: redraw the summary figure

Install R and, in the R console, run:

```r
install.packages(c("ggplot2", "patchwork", "svglite", "ragg"),
                 repos = "https://cloud.r-project.org")
```

Then run in a terminal:

```bash
Rscript figures/code/plot_readme_summary.R
```

This reads the included summary CSV and replaces the
`figures/readme_results_summary` exports. It redraws retained metrics; it does
not reproduce the benchmark from raw data. Rendering can vary by font and
package version, so resulting figure hashes may differ from the archive.

## 5. Full analysis: prepare the missing research objects first

Full reproduction requires access to the objects described in
[`DATA_AVAILABILITY.md`](../DATA_AVAILABILITY.md) and
[`data_catalogue/OBJECT_CONTRACTS.md`](../data_catalogue/OBJECT_CONTRACTS.md).
No public download links for those objects are supplied here.

Before running the computational pipeline:

1. Obtain the original combined benchmark, technology-out and family-out IDR
   inputs for LOTO/LOFO.
2. Separately obtain the no-calibration IVT-matched benchmark and its matching
   full ConquerNet model for the four-class/external-consistency branch.
3. Obtain aligned OmicsFeatures, the chromosome-fold design, and any required
   OOF/model objects. Check their genomic row identity and dimensions against
   the object contracts.
4. Use a Linux research environment and install the packages loaded by the
   chosen module. The repository does not ship a complete version-locked
   training environment or a universal one-command installer.
5. Set the input locations in
   [`analysis/01_loto_lofo/code/config.R`](../analysis/01_loto_lofo/code/config.R).
   Those LOTO/LOFO entry points use configuration variables; they do not
   provide a general `--input` option.
6. Inspect the shell runners for the downstream modules and replace their
   input defaults with the accessible local paths. Setting `OUTPUT_ROOT`
   alone does not replace all input paths.

Only after these prerequisites are met, follow this execution order:

| Order | Entry point | Output / prerequisite |
|---|---|---|
| 1 | `Rscript analysis/01_loto_lofo/code/run_loto_benchmark.R` | Technology-held-out benchmark; default new output under `analysis/01_loto_lofo/results/reproduced_run/`. |
| 2 | `Rscript analysis/01_loto_lofo/code/run_lofo_benchmark.R` | Family-held-out benchmark using matching family-out IDR. |
| 3 | `bash analysis/03_direct_posterior_ml/code/run_calibration_min51.sh` | Direct posterior, unique global-FDR calls, RF screening, summaries and Top-20 exports. Requires the separate no-calibration input branch. |
| 4 | `bash analysis/04_metagene/code/run_B_min51.sh` | FP/TP topology after valid labels/features exist. |
| 5 | Follow [`analysis/05_shap_direction/README.md`](../analysis/05_shap_direction/README.md) and its runner | SHAP only after exact OOF reconstruction succeeds. |

These are large research jobs. The public walkthrough does not start them.
The frozen constants and required audits are in
[`REPRODUCIBILITY.md`](../REPRODUCIBILITY.md).

## Troubleshooting

| Symptom | Action |
|---|---|
| `python` is not found on Windows | Try `py --version` or add Python to PATH. |
| A manifest file is missing | Download the complete repository; do not copy just the README and scripts. |
| Hash mismatch after your edits | Compare with a fresh clone; the check validates the published snapshot. |
| `Rscript` is not found | Install R and add its `bin` directory to PATH. |
| `/home/...` path does not exist | You reached a full-analysis entry point with server defaults. Prepare the missing objects and configure paths first. |
| Public checks pass but model fitting cannot run | These are different levels of verification; the archive does not include the training objects. |
