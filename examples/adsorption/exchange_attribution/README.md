# Separating carbon allocation from path exchange

The [protocol](protocol.md) fixes nine scenarios per particle size.
The study reuses the archived TRT carrier field, the complete tracer
exchange result, and the source and batch parameters from the
[unfitted column comparison](../exchange_column/REPORT.md).
[Results and verification](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.exchange_attribution.markers
uv run python -m examples.adsorption.exchange_attribution.run --prepared <prepared-directory>
uv run python -m examples.adsorption.exchange_attribution.run --check <summary-path>
uv run python -m examples.adsorption.exchange_attribution.plot --summary <summary-path>
```

Use the prepared cohort directory printed by the marker command and
the summary path printed by the simulation command. The preparation
runs four marker jobs on the same saved synthetic field. The second
command generates the allocation-enabled Julia solver from the pinned
original and checks that the original allocation reproduces the parent
solver. The parent implementation remains available for its historical
replay. The study saves the full source and prediction records under
`results/adsorption/exchange_attribution/`.

The two-period transition matrices use receiving classes as rows and
upstream classes as columns. Dry mass, liquid volumes and time histories
use kg, m³ and seconds. Read-only replay reconstructs the allocation and
transition inputs from the archives and checks the scores and numerical
refinement using the saved column runs.
