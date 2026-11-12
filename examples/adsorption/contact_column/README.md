# Figure 7 rerun with carbon contact and Péclet-matched exchange

The [protocol](protocol.md) uses the [carbon contact study](../carbon_contact/REPORT.md)
and the [Péclet-matched exchange study](../exchange_attribution/REPORT.md).
The model reconstructs both allocations from the wall events and uses the
existing conservative particle and liquid solver. It freezes the complete
curves and startup checks before scoring the two Figure 7 marker series.
[Results and verification](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.contact_column.run
uv run python -m examples.adsorption.contact_column.run --check <summary-path>
uv run python -m examples.adsorption.contact_column.plot --summary <summary-path>
```

Use the completed path printed by the run. Read-only replay verifies the
upstream archives, exact transition counts, mass allocation, recorded
physical and numerical checks and source scores. The saved model lists
carbon fractions in the same 40 flow quantile bins as the exchange
matrix. Each complete simulation conserves its own phenol balance.
Where an independently refined 12 s simulation is used for early scoring,
its curve and full mass history are saved separately.

The targeted 0.80 mm midpoint check uses only the frozen numerical model:

```sh
uv run python -m examples.adsorption.contact_column.plateau <cohort-directory>
uv run python -m examples.adsorption.contact_column.plateau <cohort-directory> --check
```

This compares two additional 600 s runs against the finest full-duration
prediction, with their mass balances and outputs saved separately.
