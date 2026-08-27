# Posterior-derived FP/TP metagene distributions

FP and TP sites are mapped to a standardized transcript axis:

```text
[0, 1/3]   5′UTR
[1/3, 2/3] CDS
[2/3, 1]   3′UTR
```

The final main figure includes exactly eight technologies with both FP and TP
counts above 50. Red is TP and blue is FP; density curves are normalized within
class. SingleMod has only one FP call and is recorded in the QC table but does
not receive a panel.

```bash
bash analysis/04_metagene/code/run_B_min51.sh
```

`OUTPUT_ROOT` and `CLASS_MATRIX` can be overridden as environment variables; the
defaults write a new reproduction directory and read the frozen W8 class matrix.

The runner uses the same direct-posterior class matrix, the no-cali IVT-matched
benchmark, technology-specific `Total>10`, and integer minimum `51` as the RF
branch. It generates effect-size tables plus the combined and per-technology
standard metagene figures.

The observed shifts are technology-specific. m6A-SAC and DeepRM have only 62 and
236 full RF-domain FP sites; the plotted topology subset contains 57 and 232,
respectively, after requiring valid metagene coordinates.

Main figure:
`results/main_figure/reliable_technologies_FP_vs_TP_standard_metagene.png`.

Effect sizes include mean/median differences, KS D and Cliff's delta. Bootstrap
intervals use the fixed balanced capped subsample recorded in the QC table.
