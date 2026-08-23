# README summary-figure QA notes

## Figure contract

- Core conclusion: ConquerBN ranks first in both technology-level and family-level held-out reconstruction, while discrimination from posterior-derived classes is task dependent.
- Archetype: quantitative grid with two primary benchmark panels and one supporting machine-learning panel.
- Backend: R only (`ggplot2` and `patchwork`).
- Final dimensions: 183 mm × 135 mm.
- Exports: editable SVG, standard PDF, Cairo PDF and 300-dpi PNG preview.

## Data integrity

- All values are copied from the finalized Total>10 result summaries in this repository.
- LOTO uses 14 technology folds; LOFO uses five independent family units.
- RF macro values use the eight TP-vs-FP, seven FP-vs-TN and eleven TP-vs-TN completed technology models.
- No rows were sampled or excluded for plotting.
- Values are descriptive macro averages; no uncertainty interval or probability-calibration claim is implied.

## Automated and visual checks

- The R source passed server-side `parse()` and ran without warnings.
- Static source preflight: 16 PASS, 4 documented WARN, 0 FAIL.
- The PNG was inspected panel by panel at its native 183 mm layout: no title, subtitle, value label, axis label or legend was clipped or overlapping.
- The four static warnings are intentional for a compact GitHub deliverable: no TIFF, 300-dpi rather than 600-dpi PNG, width stored in a variable, and R syntax confirmed by the actual server runtime.
- The PDF text-audit utility reports 1-pt text for both R PDF devices because they encode text using a 1-pt base font plus a transformation matrix. This value is not the rendered glyph size. The plotting source has a 6.5-pt minimum, and both the R-generated PNG and SVG were visually checked at final size. The PDF glyph-size result is therefore recorded as a device/auditor compatibility limitation rather than treated as a biological or numerical failure.
