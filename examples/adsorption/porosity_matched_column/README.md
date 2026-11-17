# Rerun the 25 cm adsorption cases with porosity-matched film coefficients

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the project root with the
[42.7%-void pack archive](../porosity_matched_film/README.md),
the [old column comparison](../carbon_access_transfer/README.md)
and the private Alhamed PDF preserved at their documented paths.
The prediction step uses the real Julia solver. It prints the run
directory and seals the predictions before the separate score step
loads the published markers.

```sh
source scripts/env.sh
uv run python -m examples.adsorption.porosity_matched_column.run
uv run python -m examples.adsorption.porosity_matched_column.run "$column_run" --score
uv run python -m examples.adsorption.porosity_matched_column.run "$column_run" --check
uv run python -m examples.adsorption.porosity_matched_column.plot "$column_run"
uv run pytest -q tests/adsorption/test_porosity_matched_column.py
```

Set `column_run` to the directory printed by the prediction step.
The run stores the Julia source, fixed geometry and mass inputs,
effective film conductance, time and shell refinement, mass balances,
and full histories under `results/adsorption/porosity_matched_column/`.
Curve fractions are outlet concentration divided by feed concentration.
The score keeps the recorded time windows and expresses retained mass
as kg phenol per kg dry carbon.
