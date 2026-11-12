# Stationary circular-wall transport check

This check follows the [rotating-disk and stretching benchmarks](../transport_verification/REPORT.md).
The disk radius, peak velocity, molecular diffusivity, initial mode, observation
window and stock scalar TRT are retained. The prescribed angular shear goes
to zero at the stationary smooth wall. The independent Julia reference solves
the radial concentration evolution, so it also accounts for stretching of the
mode by the new velocity. The [protocol](protocol.md) records the screens and
the reason for extending the independent reference by one grid.

Native application code is separate from the pinned OpenLB distribution and
follows its advection–diffusion and cylinder example lifecycle. The common
project build helper compiles in `.local/adsorption/stationary_curved_wall`.
No vendor edits or dependency changes are needed. Tested platform: macOS
26.6.2 ARM64, Python 3.12.9, NumPy 2.5.4, Julia 1.11.6, GCC 15.2, OpenLB 1.9r0.

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/adsorption/stationary_curved_wall/test_reference.jl
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stationary_curved_wall.run \
  --rotation-summary results/adsorption/transport_verification/curved/20261111T011153Z-1a03a0b9/summary.json
```

Omit `--rotation-summary` for a standalone run without the archived constant
comparison. The runner selects the existing Julia executable and project
recorded by bootstrap; it does not import JuliaCall or change manifests. It
qualifies the reference before any native simulation, validates a short-mesh
case first, then runs the remaining cases with at most two processes, each
using two stock OpenMP threads.

Exact commands, source snapshots, native fields and reference refinement
records are stored in a timestamped folder below
`results/adsorption/stationary_curved_wall/runs/`. Native files follow the
[previous binary layout](../transport_verification/README.md#saved-data-contract).
The Julia CSVs contain dimensionless radius and complex normalized mode
components at midpoint and final. `--check <summary.json>` replays the native
analyses against the saved radial mode without any new simulations.
`python -m examples.adsorption.stationary_curved_wall.plot` uses the latest
completed summary to produce the figures and a compact scientific record.

This check does not validate hydraulic flow, a full bed, column dispersion,
breakthrough or service time.

See the [completed report](REPORT.md) and [compact record](reference_summary.json).
The tested replay and plot commands are:

```sh
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stationary_curved_wall.run \
  --check results/adsorption/stationary_curved_wall/runs/20261111T122659Z-56628d90/summary.json
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stationary_curved_wall.plot
```
