# Predicting tracer residence times over two repeated porous cells

See the [protocol](protocol.md) for the parallel-path hypothesis,
one-cell calibration, frozen predictions and acceptance criteria. This
case preserves the pore grid and physics of the earlier
[dispersion-only comparison](../dispersion_holdout/REPORT.md) and
computes the complete tracer response.
[Results and verification](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run pytest -q examples/adsorption/residence_time_holdout/test_paths.py
uv run python -m examples.adsorption.residence_time_holdout.run
uv run python -m examples.adsorption.residence_time_holdout.run --check <summary-path>
uv run python -m examples.adsorption.residence_time_holdout.plot --summary <summary-path>
```

Replace `<summary-path>` with the completed path printed by the run.
The numerical replay checks the frozen prediction, archived field,
two-cell reconstruction and recorded results without a new solve.
The run depends on the saved one-cell N384 and N96 stock TRT archives
referenced by the prior example. It uses ten Julia threads and may
take several hours. Progress is saved in the new cohort's
`holdout/stdout.log`; the reused tracer driver has a legacy N384
startup label, while the recorded ratio of 3 corresponds to N288
tracer resolution per period. The source snapshot, fitted path weights
and full prediction are written before the numerical holdout.
