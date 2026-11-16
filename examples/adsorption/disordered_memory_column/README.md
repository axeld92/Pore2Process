# N160 travel-time memory and phenol adsorption

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the checkout root. The runner requires the saved
[N160 probes](../disordered_section_exchange/REPORT.md) and the
[qualified memory model](../disordered_travel_memory/REPORT.md).
The comparison uses the existing audited Alhamed PDF and digitization.
The solver uses the project's JuliaCall environment.

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/adsorption/disordered_memory_column/test_model.jl
uv run python -m examples.adsorption.disordered_memory_column.run
uv run python -m examples.adsorption.disordered_memory_column.run --check "$memory_column_summary"
uv run python -m examples.adsorption.disordered_memory_column.plot "$memory_column_summary"
```

Set `memory_column_summary` to the summary path printed by the run.
The read-only check replays hashes, conservation, numerical refinement
and scoring. Run the model again to regenerate the Julia histories.
The run saves the input model and contact shares, inert checks,
unscored curve files and the freeze digest before the scoring summary.
All concentrations are in kg/m³ and all times are in seconds in the
machine-readable records.
