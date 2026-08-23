# Posterior-derived FP/TP metagene distributions

FP and TP sites are mapped to a standardized transcript axis:

```text
[0, 1/3]   5′UTR
[1/3, 2/3] CDS
[2/3, 1]   3′UTR
```

The main figure includes eight technologies with both FP and TP counts above 50. Red is TP and blue is FP. Curves are normalized within class.

The observed shifts are technology-specific. m6A-SAC and DeepRM have only 62 and 236 FP sites, respectively, so their curves and effect sizes are more sensitive to individual sites than the large-class technologies.

Main figure: `results/main_figure/reliable_technologies_FP_vs_TP_standard_metagene.png`.

Effect sizes include full-data mean/median differences, KS D and Cliff's delta. Bootstrap intervals use a fixed balanced capped subsample and should be interpreted according to the protocol recorded in the QC table.
