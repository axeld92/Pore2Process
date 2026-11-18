# Conditional depth transfer of a shared adsorption equilibrium

[Protocol](protocol.md) · [Source-equilibrium audit](../REPORT.md)

The two printed equilibrium laws passed the original finite-bath
checks, yet their 25 cm column predictions gave opposite retention
errors. This test chooses one mixture weight using only retained
phenol at 10 and 15 cm and then carries it unchanged to three
held-out conditions.

Use the existing Python–Julia setup, the Alhamed PDF and the local
records linked from the [batch reference](../reference_summary.json).
The full run qualifies the candidate weights in the bath, solves and
seals the two calibration depths, selects one weight, and verifies
the holdouts. The depth transfer uses the same saved 42.7%-void N60
flow-path and wall-contact rules.

```sh
source scripts/env.sh
uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.conditional.run
```

The command prints a timestamped directory in
`results/adsorption/breakthrough_hypotheses/batch/conditional/`.
Pass that directory and `--check` to replay the evidence and
numerical checks from the saved trajectories.

## Selection from the 10 and 15 cm beds

Every declared weight passes the four Figure 2 batch comparisons.
The same two Figure 5 curves then select one common weight from
the grid by squared last-marker retention error.

| Freundlich weight | Joint squared error ((mg/g)²) |
| ---: | ---: |
| 0 | 134.202 |
| 0.25 | 37.430 |
| **0.50** | **4.776** |
| 0.75 | 16.099 |
| 1 | 54.607 |

The selected weight of 0.50 gives an effective equilibrium loading
of 129.78 mg/g at 150 mg/L. Its batch-compatible particle diffusivity
is 1.5115058150757535 × 10⁻¹² m²/s. Both timing screens and the
curve error screen pass at both calibration depths. At 10 cm
the prediction captures 62.08 mg/g through the last marker
against 64.05 mg/g from Figure 5, outside the reading range.
At 15 cm it captures 61.27 mg/g against 60.32 mg/g, inside the
reading range. Selecting the shared alternative minimizes the
two retention errors together. [Alhamed, Figure 5](https://doi.org/10.1016/j.jhazmat.2009.05.002)

## Prediction at the held-out conditions

| Bed (cm) | Flow (mL/min) | Curve | 5% time (min) | 50% time (min) | Outlet/feed RMSE | Retained (mg/g) | Retained after published midpoint (mg/g) |
| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 20 | 23.3 | Published | 63.85 | 192.59 | — | 56.37 | 17.62 |
| 20 | 23.3 | Prediction | 97.03 | 210.32 | 0.066 | 61.55 | 19.47 |
| 25 | 23.3 | Published | 125.74 | 274.40 | — | 57.13 | 10.71 |
| 25 | 23.3 | Prediction | 151.55 | 280.85 | 0.046 | 60.28 | 12.32 |
| 25 | 62.5 | Published | 16.16 | 71.63 | — | 55.67 | 26.08 |
| 25 | 62.5 | Prediction | 22.94 | 68.38 | 0.027 | 56.89 | 26.63 |

The high-flow 25 cm case now agrees with the last-marker reading
range of 54.28–57.06 mg/g. It captures 0.55 mg/g more than the
published curve after the published midpoint. The fixed Langmuir
prediction captured 4.90 mg/g too little during that period.
The 5% crossing remains 42% late. At 20 cm and low flow the
5% crossing is 52% late. Both held-out low-flow cases capture
too much phenol by the last marker. Every held-out midpoint
and full-curve RMSE passes its screen. The complete prediction
still fails the joint timing and retention checks.
[Alhamed, Figures 4–5](https://doi.org/10.1016/j.jhazmat.2009.05.002)

All three held-out curves pass the numerical comparison with
48 shells and a 4 s step. The two higher-weight 10 cm candidates
needed 96 shells at 2 s to meet the early timing tolerance.
The second 25 cm low-flow curve in Figure 5 is an alternate reading
of the case and supports the same conclusion. The mixture remains a
conditional model for these columns.

```sh
CONDITIONAL_RUN=results/adsorption/breakthrough_hypotheses/batch/conditional/20261118T000130Z-e9054190
uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.conditional.run "$CONDITIONAL_RUN" --check
```

[Compact verification record](reference_summary.json)
