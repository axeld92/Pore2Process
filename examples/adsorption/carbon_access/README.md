# Carbon accessibility sensitivity

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the project root with the [frozen N160 column](../disordered_memory_column/REPORT.md),
the [repeat-join control](../disordered_repeat_join/REPORT.md) and the
audited Figure 7 inputs available.

```sh
source scripts/env.sh
uv run pytest -q tests/adsorption/test_carbon_access.py
uv run python -m examples.adsorption.carbon_access.run
uv run python -m examples.adsorption.carbon_access.run --check "$carbon_access_summary"
uv run python -m examples.adsorption.carbon_access.plot
```

Set `carbon_access_summary` to the summary path printed by the run.
The run records the scenario plan, its generated Julia solver,
both cohort histories and the curve freeze digest before scoring.
Replay checks saved source and curve hashes, phenol balances, numerical
refinement and Figure 7 scores. The inherited N160 delay network stays
in its archived source file. The generated code replaces only its
particle-contact calculation with the same particle operator applied
to both carbon cohorts.
