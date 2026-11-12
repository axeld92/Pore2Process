# Dispersion transfer from one cell to two

The [protocol](protocol.md) fixes one dispersion coefficient from the
[N384 one-cell outlet](../porous_outlet_tracer/N384.md) and tests its
two-cell prediction against the N288-per-period Julia tracer solver.
The [report](REPORT.md) records the outcome. The original N96 TRT
flow supplies both repeats, and the calculation keeps the 80 µs
step, 4 ms pulse, and open boundary convention of the one-cell test.

From the configured project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.dispersion_holdout.run
uv run python -m examples.adsorption.dispersion_holdout.run \
  --check results/adsorption/dispersion_holdout/runs/20261112T033234Z-32e27f5e/summary.json
uv run python -m examples.adsorption.dispersion_holdout.plot
```

A fresh run prints its own summary path. The model and coefficient
are hashed before the held-out computation begins. The case archives
all simulated steps and may terminate when its full-window cumulative
profile error is already mathematically above the existing limit.
The read-only replay reconstructs the repeat, verifies the frozen
model and recomputes the stopping event and comparison. It launches
no simulations.
