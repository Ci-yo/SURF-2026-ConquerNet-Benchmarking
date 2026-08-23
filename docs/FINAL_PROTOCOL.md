# Frozen final analysis protocol

## Benchmark branch

### LOTO

- Apply `Total>10` to retained technology-site observations before fitting.
- Retrain 14 strict leave-one-technology-out models.
- Define held-out BBmix MAP truth with `m6ASiteProb>0.5` versus `<=0.5`, restricted to held-out `Total>10`.
- Compare ConquerBN, MeanProb, RankProduct, RRA and matching continuous IDR on one common site domain.
- Exclude Vote.

### LOFO

- Hold out one complete technology family and retrain on the other four.
- Use each held-out technology's BBmix MAP plus `Total>10` as its evaluation label.
- Average technology metrics within family and treat five families as the independent units.
- Use family-out IDR, not technology-out IDR.

## Direct-posterior ML branch

- Use the no-IVT-calibration, IVT-matched full-model/data pair.
- Reconstruct direct joint TN/FN/FP/TP posteriors.
- Apply class-specific global FDR `<0.05`, retain unique class calls and require technology-specific `Total>10`.
- Do not replace small classes with MAP/Viterbi calls.
- Run a task when both classes are `>50` (`min_reliable=51`).
- Use 185 screened OmicsFeatures, chromosome-blocked five-fold OOF, 300 trees and training-only sampling to at most 100,000 per class.
- Report permutation importance and exact-reconstruction SHAP direction separately.

## External consistency branch

- Build TP/FP labels from the no-cali full model before IVT and perturbation analyses.
- Treat IVT and writer/eraser perturbation as independent consistency evidence, not perfect ground truth.
- Retain negative or mixed findings.

## Excluded from the final route

```text
Total>20 sensitivity outputs
Vote
old separately-thresholded marginal four-class labels
MAP/Viterbi fallback for the direct-posterior branch
simulation ranking without sampled site-level TrueM6A truth
large model, cache, posterior and OOF objects
```
