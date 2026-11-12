# Stationary circular wall at half the time step

The [protocol](protocol.md) repeats the regularized TRT, half-way reflection
configuration from the [stock OpenLB comparison](../stock_transport_candidate/REPORT.md).
The test uses the same 128-cell disk and the same OpenLB 1.9r0 executable.
The [report](REPORT.md) and [compact record](reference_summary.json) give the
comparison to the qualified Julia reference.

Run these commands from the project root on the configured macOS ARM64 machine.

```sh
source scripts/env.sh
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stationary_wall_time_step.run
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stationary_wall_time_step.run \
  --check results/adsorption/stationary_wall_time_step/runs/20261111T144732Z-c53cd814/summary.json
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.stationary_wall_time_step.plot
julia --startup-file=no --project=julia examples/adsorption/stationary_curved_wall/test_reference.jl
```

The run uses the exact stock comparison archive and original Julia reference
listed in the protocol and verified by digest. Their generation commands are
in the [stock comparison run guide](../stock_transport_candidate/README.md)
and the [Julia reference run guide](../stationary_curved_wall/README.md).
Each new timestamped directory under `results/adsorption/stationary_wall_time_step/runs/`
contains the binary, native command and output for each case, full histories,
snapshots, inputs, checks and source hashes. Use the new summary path with
`--check` to replay a fresh run. The plot command defaults to the latest
completed summary. The [saved data contract](../transport_verification/README.md#saved-data-contract)
describes binary layout and units.
