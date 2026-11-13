# N160 disordered packing

Use the [N128 archive](../disordered_n128/REPORT.md) and the existing
stock OpenLB 1.9r0, Python and Julia environments. The [protocol](protocol.md)
defines the matched comparison and acceptance limits. [Results](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.disordered_packing.native build
uv run python -m examples.adsorption.disordered_n160.native
```

Use the cohort path printed by the native run as `$n160_run`. After
the field passes verification:

```sh
uv run python -m examples.adsorption.disordered_n160.run "$n160_run" handoff
uv run python -m examples.adsorption.disordered_n160.run "$n160_run" trace
uv run python -m examples.adsorption.disordered_n160.run "$n160_run" seal
uv run python -m examples.adsorption.disordered_n160.run "$n160_run/summary.json" check
uv run python -m examples.adsorption.disordered_n160.plot "$n160_run/summary.json"
```

The field, face and marker files use the [packing example's data contracts](../disordered_packing/README.md).
The retained Git evidence is the compact summary and plot. Full local
cohort artifacts are in the ignored `results/` tree.
