# Eight-cell tracer holdout

The [protocol](protocol.md) defines the carrier, perfect-sink observable,
matched within-class reset and acceptance limits. The calculation requires
the archived [two-cell survivor exposure cohort](../survivor_exposure/REPORT.md)
and its pinned TRT and wall-uptake inputs. [Results](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=julia examples/adsorption/repeated_geometry/test_track.jl
uv run pytest -q examples/adsorption/repeated_geometry/test_predict.py
uv run python -m examples.adsorption.repeated_geometry.run
uv run python -m examples.adsorption.repeated_geometry.run --check <summary-path>
uv run python -m examples.adsorption.repeated_geometry.diagnostics --summary <summary-path>
uv run python -m examples.adsorption.repeated_geometry.plot --summary <summary-path>
```

Use the summary path printed by the run. The runner verifies the
archived field, freezes the reduced predictions and checks the
two-cell tracker against a same-seed subset before running the four
eight-cell cases. The read-only checker reconstructs predictions
and scores the saved trajectories again.

Each little-endian Float64 marker row contains the eight absolute
first-passage times followed by the eight first wall-hit times,
all in seconds after inlet. The tracer clock advances the cell number
on first forward seam crossing. A later backtrack retains that
clock's phase. NaN means the corresponding event was unobserved.
No wall contact is assigned to reflecting faces at the exterior
inlet or outlet.
