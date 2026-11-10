# Date-stone carbon: LDF literature comparison

[Scientific report: purpose, theory, results and conclusions](REPORT.md).
The comparison is numerically reproducible and fails its physical screens.
It is separate from the synthetic coconut-shell-carbon pumping example.

From the configured project root:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.alhamed.source
uv run python -m examples.adsorption.alhamed.study
uv run python -m examples.adsorption.alhamed.check
```

`source` audits the supplied original PDF and marker readings; `study` performs
new Julia calibration/forward calculations; `check` replays the saved comparison.
See the [detailed execution record](../../../docs/verification/adsorption/alhamed.md)
and [frozen protocol](protocol.md) for source
prerequisites, numerical settings and the separate failed physical criteria.
