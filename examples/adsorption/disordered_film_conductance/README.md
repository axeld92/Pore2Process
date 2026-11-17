# Liquid-film conductance from the saved N160 flow

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the project root with the accepted
[N160 stock TRT field](../disordered_n160/REPORT.md) and the
[early-escape study](../disordered_early_bypass/REPORT.md) available.
The archive hashes and paths are in the compact record. The first command
runs Julia tracer probes at half the original time step and prints
their output directory. The finish step uses that directory, runs the
quarter step if the protocol calls for it and seals the film estimates.

```sh
source scripts/env.sh
uv run python -m examples.adsorption.disordered_film_conductance.run --probe
uv run python -m examples.adsorption.disordered_film_conductance.run --finish "$film_run"
uv run python -m examples.adsorption.disordered_film_conductance.run --check "$film_run/summary.json"
uv run pytest -q tests/adsorption/test_disordered_film_conductance.py
```

The check command verifies the input and output hashes and recalculates the
estimates from the saved trajectories. Film coefficients use m/s and
conductances use m³/s. The observation arrays have the same first-passage
and wall-contact layout as in the original
[disordered packing example](../disordered_packing/README.md).
