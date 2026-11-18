# Matched-conversion water activity in batch and packed bed

[Bottom line and results](REPORT.md) · [Evidence](evidence.md) · [Protocol](protocol.md)

This example compares the water inventory and bead-surface water activity
along the same pentyl propanoate reaction path. Python uses measured
phase and activity data for liquid composition sensitivities. The
previously verified Julia rate and film models give the water flux
from the 6 g bed.

From the project root after [environment setup](../../../docs/guides/setup.md):

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_water_balance.run
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_water_balance/test_activity.py
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_water_balance.check
```

The runner writes the detailed grid and summaries to new directories
under `results/biocatalysis/pentyl_propanoate_water_balance/`. Conversion
and water activity use dimensionless fractions. Blank activity cells
indicate that an input composition is outside the projection's
available data. Use `--output` to select a fresh directory in the
checkout. The retained summary stores hashes of the source data,
models and output tables.

Verified on macOS arm64 with Python 3.12.9 and Julia 1.11.6.
The five activity tests and 16 shared Julia film checks passed.
Two runs reproduced the outputs. Refining the conversion grid
from 0.005 to 0.0025 reproduced all reported film-activity maxima.
