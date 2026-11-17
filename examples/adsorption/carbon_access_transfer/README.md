# Carbon-access transfer across bed depth and flow

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the checkout root with the archived
[N160 TRT field](../disordered_n160/REPORT.md),
[carbon-access cases](../carbon_access/REPORT.md), and private Alhamed
source available. The probe step runs two tracer calculations on the
archived velocity shape. The column step runs three independent Julia
processes across operating cases.

```sh
source scripts/env.sh
uv run python -m examples.adsorption.carbon_access_transfer.probe
uv run python -m examples.adsorption.carbon_access_transfer.run
uv run python -m examples.adsorption.carbon_access_transfer.run --check "$transfer_summary"
uv run python -m examples.adsorption.carbon_access_transfer.plot "$transfer_summary"
uv run pytest -q tests/adsorption/test_carbon_access_transfer.py
```

Set `transfer_summary` to the final summary printed by the column
step. The runner uses the latest completed probe archive by default.
Use `--probes` to select an explicit probe summary for a reproduction.
Replay checks source hashes, archived histories, conservation,
refinement and observed scores. The predicted curve files, plans and
probe arrays are generated under `results/`.
