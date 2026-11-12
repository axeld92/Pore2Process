# Exchange between fast and slow paths

The [protocol](protocol.md) defines a Julia marker calculation on the saved
balanced OpenLB TRT carrier field and a Python outlet prediction using
the [existing one-cell calibration](../residence_time_holdout/REPORT.md).
Marker pairs, frozen curves, and the comparison against the archived
two-cell response are saved under `results/adsorption/path_exchange/`.
[Results and verification](REPORT.md)

From the configured project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=julia examples/adsorption/path_exchange/test_track.jl
uv run pytest -q examples/adsorption/path_exchange/test_predict.py
uv run python -m examples.adsorption.path_exchange.run
uv run python -m examples.adsorption.path_exchange.run --check <summary-path>
uv run python -m examples.adsorption.path_exchange.plot --summary <summary-path>
```

Replace `<summary-path>` with the completed path printed by the run.
Read-only replay reconstructs the compact carrier field, recalculates
predictions from saved marker pairs and checks the comparison. It
does not rerun the marker solver or either archived simulation.

Each marker record stores two little-endian Float64 durations in seconds.
Rows follow inlet sample order. The first duration ends at the first
crossing into the second period. The second ends at the outlet. NaN
identifies an uncompleted duration at the specified time cap. The model
uses complete pairs only and reports the recovered fraction and the
ten-by-ten matrix of counts, with first-period decile as row and second
period as column, fastest to slowest. Outlet curves report fractions
of the injected dose against seconds from the start of the pulse.
