# Figure 7 retained-phenol audit

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

From the project root, with the audited Alhamed PDF and
[frozen repeat-join results](../disordered_repeat_join/REPORT.md) available:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.figure7_retention.run
uv run python -m examples.adsorption.figure7_retention.run --check
uv run python -m examples.adsorption.figure7_retention.plot
uv run pytest -q tests/adsorption/test_figure7_retention.py
```

The default output is `results/adsorption/figure7_retention/summary.json`.
The check recomputes the analysis and compares it with the compact record.
The run uses the frozen cumulative outlet and inventory histories in kg
and converts units only for the report and plot.
