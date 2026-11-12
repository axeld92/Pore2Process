# Wall exposure among first-cell survivors

This study reuses the pinned [wall-uptake cohort](../wall_uptake/REPORT.md)
and its archived inert marker pairs and first-hit times. It requires the
prior shared-carbon, wall-contact and wall-uptake result folders. The
[holdout protocol](protocol.md) fixes the prediction and the gate. The
[column protocol](column_protocol.md) defines the transfer to the
batch-qualified shared-carbon model. [Results](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run pytest -q examples/adsorption/survivor_exposure/test_predict.py
uv run python -m examples.adsorption.survivor_exposure.run
uv run python -m examples.adsorption.survivor_exposure.run --check <holdout-summary>
uv run python -m examples.adsorption.survivor_exposure.column --holdout <holdout-summary>
uv run python -m examples.adsorption.survivor_exposure.column --check <column-summary>
uv run python -m examples.adsorption.survivor_exposure.plot --summary <column-summary>
```

Use the paths printed by the two runs. The column entry point checks
the saved holdout before loading the solver. The read-only checks
rebuild the model from the archives, verify the frozen files and replay
the scoring. The column checker also replays conservation and
discretization checks.

The holdout stores fractional carrier flows as donor-by-receiver joint
fluxes. The Julia input is a receiver-by-donor concentration matrix.
Survival and normalized transit times are dimensionless. Column film
conductance is in m³/s, film coefficients are in m/s, and area weights
are fractions of the total external particle area.
