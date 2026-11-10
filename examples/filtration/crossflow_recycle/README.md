# Crossflow recycle and circulation energy

The verified OpenLB resistance ratios and Julia membrane closure now supply a
stock DWSIM recycle flowsheet. Separate fresh-feed and circulation pumps,
channel loss, water recovery, purge and temperature control are explicit.
This is a lumped, clean-water hydraulic scenario; channel resistance and
circulation targets are declared assumptions.

From the project root on the tested Mac:

```sh
source scripts/env.sh
uv sync --locked
uv run python -m examples.filtration.crossflow_recycle.run
uv run python -m examples.filtration.crossflow_recycle.check
uv run pytest -q
/Users/dula/.julia/juliaup/julia-1.11.6+0.aarch64.apple.darwin14/bin/julia --startup-file=no --project=julia -e 'using Pkg; Pkg.test()'
```

Existing Python/Julia environments, locks, DWSIM 10.2.15.0, private .NET SDK
10.0.403 / runtime 10.0.14 and the shared bridge are reused. New builds live
under `.local/crossflow_recycle/`; no original simulator files are modified.
The Julia user compilation cache must be writable. The original VTK regression
tests require normal macOS graphics access in a restrictive automation sandbox.
No new package or privileged installation is required for this milestone.

The default parent is `results/process_filtration/latest.json`. Its saved
filtration cases and original OpenLB study must remain intact. Reproduction
from a clean checkout first follows [the graded-support instructions](../graded_support/README.md)
and [the once-through filtration instructions](../once_through/README.md).
The new stage verifies both ancestors before reusing the measured data. It
runs no new OpenLB simulations and makes no new literature fit.

Optional explicit paths:

```sh
uv run python -m examples.filtration.crossflow_recycle.run --parent results/process_filtration/runs/20261108T231207Z-316d61cb/summary.json
uv run python -m examples.filtration.crossflow_recycle.check results/crossflow_recycle/runs/RUN_ID/summary.json
```

Every run uses a fresh directory. `latest.json` is updated only after passing
the checks. Failed runs and logs are retained. The matrix contains 21 design
scenarios and four seed/tolerance variants, each solved, saved, loaded and
solved in actual DWSIM. A separate one-iteration case must fail convergence
and produce no valid summary. The checker inspects saved evidence and replays
Julia; it performs **no fresh external solves**.

Exports include input/source/parent snapshots, property and runtime versions,
implementation and artifact hashes, actual `crossflow.dwxmz` cases, DWSIM
renders, material/energy states, convergence diagnostics, unit/overall balance
checks, `comparison.csv`, `circulation_energy.png` and timings. All volume-based
metrics explicitly identify their common fresh-feed reference density; actual
flashed stream volumes are retained separately. Cooling heat is exported
separately from shaft work.

See [the model, results and limitations](../../../docs/crossflow_recycle.md).
The original Panel app remains independent of this CLI workflow.
