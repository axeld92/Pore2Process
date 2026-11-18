# Figure 7 equilibrium sensitivity

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the project root with the configured Python and Julia
environments and the saved runs for
[particle equilibrium](../breakthrough_hypotheses/batch/README.md),
[common carbon access](../carbon_access/README.md),
[N160 travel-time memory](../disordered_travel_memory/README.md) and
[DWSIM bed pumping](../column_transport/README.md).

```sh
source scripts/env.sh
uv run --locked python -m examples.adsorption.figure7_equilibrium.run
```

The command prints the timestamped directory. Use its path for
`$size_run` below. Replay verifies the saved forecasts, reruns the small
GAC4 batch check, and recomputes the curve and process metrics.

```sh
uv run --locked python -m examples.adsorption.figure7_equilibrium.run "$size_run" --check
uv run --locked python -m examples.adsorption.figure7_equilibrium.clean_limit "$size_run"
uv run --locked python -m examples.adsorption.figure7_equilibrium.plot "$size_run"
```

Each run saves a plan with fixed flow, tracer, carbon and equilibrium
inputs, the saved history at each tested resolution, a forecast seal
and the scored results. Run data stay under `results/`. The process
calculation reads the saved DWSIM shaft powers and multiplies them by
service times at the published flow. Both stages use the same dry mass.
The pump stage uses a 14.325 cm synthetic hydraulic bed, and the
breakthrough comparison uses Alhamed's 10 cm bed.
The clean-bed diagnostic uses the same plan and saves the wall-allocated
and flow-proportional limits in `clean_limit.json`.
