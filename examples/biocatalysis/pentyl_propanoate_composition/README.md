# Composition and catalyst-history kinetic audit

[Engineering conclusions](REPORT.md) · [Frozen protocol](protocol.json)

All sixteen published enzyme batch assays are independently
traceable in the literature source files. Python prepares the
whole-curve splits and exports; Julia fits the constrained
mass-specific rates, integrates the activity state and reuses
the checked existing film and finite-volume reactor solver.
No new Python solver or per-curve empirical activity is used.
Temperatures are in °C at the interface, absolute Kelvin in
inverse-temperature relations; concentrations and rates use
the documented mol/kg liquid and mol/(kg_cat s) bases.

From the repository root after [setup](../../../docs/guides/setup.md):

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_composition.run
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_composition/test_composition.py
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_composition/test_composition.jl
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_composition.check
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_composition.verify_primary
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_composition.plot \
  --output results/biocatalysis/pentyl_propanoate_composition/composition_audit.png
```

The full campaign writes to a fresh timestamped results directory,
or use `--output` with a fresh destination. It saves the batch fit
and all pointwise batch predictions, then the frozen reactor
predictions, activity-age scenarios and comparison. Every numerical
output is hashed. The [retention utility](retain.py) regenerates
the compact tracked reference in a new destination:

```sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_composition.retain \
  --run results/biocatalysis/pentyl_propanoate_composition/verification_3 \
  --output results/biocatalysis/pentyl_propanoate_composition/compact_replay
```

Replace the run directory above with the newly generated run.
The compact reference retains primary parameters/holdouts,
sensitivity results, selected batch progress, all reactor
predictions and uniform-activity age checks. The full detailed
pointwise and optimizer ledgers remain only in routine results.

Verified on macOS arm64 with Python 3.12.9 and Julia 1.11.6.
All 50 new supplement rows matched Table S2 exactly.
Three new Python tests and the four preceding loading-evidence
tests passed. The 47 direct Julia assertions passed, including
mass-time invariance, an independent exact constant-hazard
recovery/retention integral, synthetic rate-parameter recovery,
and existing film/axial flux consistency. Tightened batch
integration tolerances changed primary conversion by at most
4.57e-7. Refining 512 to 1,024 axial cells changed conversion
by at most 1.17e-5; normalized reaction/flux balance error
below 4.94e-10 and maximum sampled film residual below
6.68e-13 mol/L/min. The previous central film model replays
within 1e-15. All twenty original reactor observations remain
outside calibration. The numerical solves passed; some
nonlinear candidate fits have weak/boundary or unresolved
parameters (recorded), and the predictive transfer screen fails.

A cold-process [primary replay](reference/replay.json) repeated
ten fitting cases (all data, eight paired folds, and the hot-acid
challenge), thirteen models per case, and 3,915 retained batch
predictions. All repeated parameters, scores, conversion and
activity predictions were identical. Repeating all 440 frozen
reactor predictions and 280 age-screen conditions also produced
exactly the same conversions. The sensitivity campaigns were
completed in the full run; the replay specifically rechecks the
primary fits and frozen bed outputs. Nine repository layout
checks additionally passed, and the previous two catalyst-loading
audits still match their frozen source hashes.
