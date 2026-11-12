# Carbon contact in the saved TRT packing

The [protocol](protocol.md) fixes the flow, marker runs, contact measure
and carbon mass accounting. These calculations use the archived carrier
flow and solid geometry. Each reflecting marker produces both transit
times and solid-wall contact records. They measure external grain-envelope
exposure in the synthetic geometry.
[Results and verification](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=julia examples/adsorption/carbon_contact/test_track_contacts.jl
uv run pytest -q examples/adsorption/carbon_contact/test_analyze.py
uv run python -m examples.adsorption.carbon_contact.geometry
uv run python -m examples.adsorption.carbon_contact.prepare
uv run python -m examples.adsorption.carbon_contact.analyze --prepared <cohort-directory>
uv run python -m examples.adsorption.carbon_contact.analyze --check <summary-path>
uv run python -m examples.adsorption.carbon_contact.plot --summary <summary-path>
```

Use the directory printed by preparation and the summary path printed by
analysis. The preparation runs seven marker calculations. Read-only
replay checks the independent geometry reconstruction, archived trajectory
times, source hashes and contact summaries.

The marker file has eight little-endian Float64 columns. It contains
first and second transit times in seconds, carbon hit counts before and
after the first seam crossing, reflected normal displacement in metres
for those two phases, the count of reflections at the exterior axial ends,
and carbon hits after the first crossing that occurred back in the first
geometry. The paired `.f64.faces` stream contains a signed little-endian
Int32 marker ID, little-endian UInt32 wall face ID and little-endian Float64
reflection distance in metres, with 16 bytes per row. Marker IDs are one
based, with negative sign for hits after the first seam crossing. Face IDs
fold the repeated geometry and are one based. For zero-based liquid voxel
row `i` and axis `a`, the negative and positive sides are `6*i+2*a+1`
and `6*i+2*a+2`. The analysis checks all event IDs and per-marker sums.
