# No-cali external consistency analyses

The full ConquerBN model was refit without IVT calibration so that IVT could be used downstream without directly recycling calibration information.

## MeRIP IVT

The primary endpoint uses a 2-of-3 antibody consensus and requires `Total>10` in both IVT and mRNA inputs. Contrary to the initial expectation, IVT-positive calls were more frequent among model-derived FP than TP sites (16.7% vs 11.4%). The result is retained as a negative consistency finding.

## Writer/eraser perturbation

Matched contrasts were evaluated with continuous posterior change and binary expected-direction/MAP-switch endpoints. Responsiveness was contrast-dependent; the current data do not show a universal TP-over-FP response advantage.

## Boundaries

- TP/FP labels in this branch come from the no-cali full model.
- IVT and perturbation are external consistency checks, not unquestioned gold standards.
- Contrasts with zero or very small eligible classes remain descriptive or non-evaluable.

See `results/MeRIP_IVT/` and `results/perturbation/` for compact tables and protocol text.
