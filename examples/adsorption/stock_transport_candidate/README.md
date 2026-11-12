# Stock TRT alternatives at a stationary circular wall

This study repeats the [stationary-wall tracer benchmark](../stationary_curved_wall/REPORT.md)
with three combinations of stock TRT collision and wall reflection. Its
[protocol](protocol.md) defines the matched inputs, comparison cases and
acceptance limits. The [report](REPORT.md) gives the scientific results and the
[compact record](reference_summary.json) gives the numerical comparisons.

## Run and replay

The commands use the project's pinned OpenLB 1.9r0 setup on macOS 26.6.2 ARM64,
with GCC 15.2, Python 3.12.9, NumPy 2.5.4, Matplotlib 3.11.2 and Julia 1.11.6
for the reference. The comparison requires the archived stationary disk
summary and its files under
`results/adsorption/stationary_curved_wall/runs/20261111T122659Z-56628d90/`.
That archive supplies the exact physical input snapshots and qualified radial
mode used here. Its summary and mode digests are checked before native runs.

From the project root:

```sh
source scripts/env.sh
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stock_transport_candidate.run
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stock_transport_candidate.run \
  --check results/adsorption/stock_transport_candidate/runs/20261111T133159Z-7f69a322/summary.json
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stock_transport_candidate.plot
julia --startup-file=no --project=julia examples/adsorption/stationary_curved_wall/test_reference.jl
```

For a fresh execution, pass the new run's `summary.json` to `--check`. The
plot command uses the latest completed run by default and accepts
`--summary <path>` for a specific run. The original
[stationary-wall run guide](../stationary_curved_wall/README.md) documents how
the comparison archive was generated.

The external build lives in `.local/adsorption/stock_transport_candidate/`.
Each completed study stores native inputs, field snapshots, full extrema and
inventory histories, link coverage, command lines, binary, vendor integrity
and a source snapshot in its timestamped results directory. Replay verifies
hashes and recalculates the screens from the saved fields and histories. Field
layout and SI units follow the
[transport data contract](../transport_verification/README.md#saved-data-contract).
Solid cells in the half-way setup contain temporary population copies. All
tracer statistics and axial comparisons use the liquid mask.
