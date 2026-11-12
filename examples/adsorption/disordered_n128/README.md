# N128 disordered packing

The [protocol](protocol.md) defines the grid limit and the matched
tracer comparison. This case requires the [N96 refinement archive](../disordered_refinement/REPORT.md)
and the existing project OpenLB, Julia and Python environments.
[Results](REPORT.md)

From the project root, build the external stock TRT application
if needed and run the N128 solve:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.disordered_packing.native build
uv run python -m examples.adsorption.disordered_n128.native
```

Use the printed cohort path as `$n128_run`. After the native field
passes verification:

```sh
uv run python -m examples.adsorption.disordered_n128.run "$n128_run" handoff
uv run python -m examples.adsorption.disordered_n128.run "$n128_run" trace
uv run python -m examples.adsorption.disordered_n128.run "$n128_run" seal
uv run python -m examples.adsorption.disordered_n128.run "$n128_run/summary.json" check
uv run python -m examples.adsorption.disordered_n128.plot "$n128_run/summary.json"
```

The native warm field, converged field, face flows and raw marker
arrays use the contracts in the [N64 guide](../disordered_packing/README.md).
The compact report and reference summary are versioned. Full local
run directories are kept under the ignored `results/` tree.
