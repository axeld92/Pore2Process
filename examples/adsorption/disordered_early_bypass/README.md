# Early wall contact and film transfer

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the checkout root with the two archived
[carbon-access transfer](../carbon_access_transfer/README.md) summaries
identified by path and SHA-256 in the compact record. The runner reads the
saved tracer arrivals and wall contacts and the frozen column setup.
It uses the existing Julia travel-time classification. Python builds and
composes the pulse distributions. No additional OpenLB or adsorbing-bed
simulation is needed.

```sh
source scripts/env.sh
uv run python -m examples.adsorption.disordered_early_bypass.run
uv run python -m examples.adsorption.disordered_early_bypass.run --check "$bypass_summary"
uv run python -m examples.adsorption.disordered_early_bypass.plot "$bypass_summary"
uv run pytest -q tests/adsorption/test_disordered_early_bypass.py
```

Set `bypass_summary` to the final path printed by the runner. Replay checks
hashes, recalculates the held-out scores and compares the archived film pulses
with the analytical lower bound. The run writes weights and pulse data to
`results/adsorption/disordered_early_bypass/`. The plot command reproduces
the report figure in the example's `figures/` directory.
