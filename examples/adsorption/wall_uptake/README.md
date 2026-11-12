# Film transfer from an absorbing-wall tracer test

The [protocol](protocol.md) fixes the marker calculation and conductance
rule. The first sampled solid-wall hit removes the marker from the
absorbing mass balance. The saved continuation provides unreacted transit
time and records contacts in the second period. Its mass is not counted
as surviving after the hit.
[Results and verification](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run pytest -q examples/adsorption/wall_uptake/test_estimate.py
julia --startup-file=no --project=julia examples/adsorption/wall_uptake/test_track_uptake.jl
uv run python -m examples.adsorption.wall_uptake.prepare
uv run python -m examples.adsorption.wall_uptake.run --prepared <cohort-directory>
uv run python -m examples.adsorption.wall_uptake.run --check <summary-path>
uv run python -m examples.adsorption.wall_uptake.plot --summary <summary-path>
```

Use the printed paths for the prepared cohort and the completed summary.
The preparation runs four sets of markers. The column run reconstructs
the absorbing estimates, preserves the Péclet exchange and full shared
carbon inventory, and scores Figure 7 after freezing the predicted curves.
Read-only replay verifies the existing archives, binary reduction, shared
carbon mass balances and both numerical and literature comparisons.

Marker files contain four little-endian Float64 columns, in seconds:
two inert transit durations and the first solid-wall hit time within each
transit segment. NaN in either hit column means the marker traversed that
segment without a sampled wall hit. Bed film conductances are in m³/s,
external surface areas are in m², and film coefficients are in m/s.
