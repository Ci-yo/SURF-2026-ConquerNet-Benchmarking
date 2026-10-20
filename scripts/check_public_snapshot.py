#!/usr/bin/env python3
"""Check the public archive and summary tables, without private RDS inputs."""
import csv
import hashlib
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".r", ".sh", ".md", ".csv", ".tsv", ".txt", ".svg", ".py", ".yml"}


def canonical_bytes(path):
    data = path.read_bytes()
    if path.suffix.lower() in TEXT_SUFFIXES or path.name in {".gitignore", ".gitattributes"}:
        data = data.replace(b"\r\n", b"\n")
    return data


def read_csv(relative):
    with (ROOT / relative).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    with (ROOT / "MANIFEST_SHA256.tsv").open(encoding="utf-8-sig", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    for row in manifest:
        path = ROOT / row["relative_path"]
        if not path.resolve().is_relative_to(ROOT):
            raise ValueError(f"Manifest path leaves repository: {row['relative_path']}")
        data = canonical_bytes(path)
        if len(data) != int(row["size_bytes"]) or hashlib.sha256(data).hexdigest() != row["sha256"]:
            raise ValueError(f"Manifest mismatch: {row['relative_path']}")

    summary = read_csv("figures/source_data/readme_summary_metrics.csv")
    if len(summary) != 11:
        raise ValueError("Expected 8 benchmark rows and 3 RF task summaries.")
    for row in summary:
        for key in ("AUROC", "AUPRC"):
            value = float(row[key])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"Invalid {key}: {row}")
    definitions = [
        ("LOTO", "analysis/01_loto_lofo/results/LOTO/four_method_IDR_summary.csv",
         "mean_auroc", "mean_auprc", "n_folds", 14),
        ("LOFO", "analysis/01_loto_lofo/results/LOFO/four_method_IDR_overall_summary.csv",
         "family_macro_auroc", "family_macro_auprc", "n_independent_families", 5),
    ]
    for evidence, relative, auc_col, pr_col, n_col, expected_n in definitions:
        archived = read_csv(relative)
        if evidence == "LOTO":
            archived = [r for r in archived if r["analysis_set"] == "all_14_technologies"]
        by_method = {r["method"]: r for r in archived}
        selected = [r for r in summary if r["evidence"] == evidence]
        expected_methods = {"ConquerBN", "IDR", "RRA", "RankProduct"}
        if len(selected) != 4 or {r["category"] for r in selected} != expected_methods:
            raise ValueError(f"Unexpected {evidence} method set.")
        for row in selected:
            source = by_method[row["category"]]
            if int(row["n_units"]) != expected_n or int(source[n_col]) != expected_n:
                raise ValueError(f"Unexpected {evidence} unit count.")
            if not (math.isclose(float(row["AUROC"]), float(source[auc_col]), abs_tol=1e-12)
                    and math.isclose(float(row["AUPRC"]), float(source[pr_col]), abs_tol=1e-12)):
                raise ValueError(f"Summary disagrees with source: {evidence}/{row['category']}")
    print(f"Public snapshot checks passed: {len(manifest)} files; 11 summary rows.")
    for row in summary:
        print(f"{row['evidence']:4} {row['category']:12} n={int(row['n_units']):2d} "
              f"AUROC={float(row['AUROC']):.4f} AUPRC={float(row['AUPRC']):.4f}")
    print("This verifies retained files/tables; it does not rerun model fitting or per-site inference.")


if __name__ == "__main__":
    main()
