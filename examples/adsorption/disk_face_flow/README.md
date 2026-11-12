# Conservative tracer transport on the stationary disk

This Julia calculation uses the same physical disk, velocity, diffusivity,
initial tracer profiles, and observation times as the
[stationary-wall OpenLB benchmark](../stationary_curved_wall/REPORT.md).
Its independent comparison is the archived, converged Julia radial mode.
The [protocol](protocol.md) and [spatial confirmation note](refinement_note.md)
record the mesh and time-step cases. The [report](REPORT.md) gives the outcome.

The Julia solver lives in `julia/adsorption/DiskFaceTransport.jl`. It solves
the full real-valued concentration field and uses Julia's standard
LinearAlgebra and SparseArrays solvers. Python coordinates the runs, prepares
the analytical initial field, verifies the archived reference and compares
the native results. No new dependencies are required.

Run from the project root on the configured macOS ARM64 machine. The
recorded Julia 1.11.6 runtime is selected by `scripts/bootstrap_julia.py`.
The reference and original sampling points come from the exact archived
stationary wall cohort `20261111T122659Z-56628d90`.
The [stationary wall run guide](../stationary_curved_wall/README.md) documents
the generation of that archive.

```sh
source scripts/env.sh
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.disk_face_flow.run
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.disk_face_flow.run \
  --check results/adsorption/disk_face_flow/runs/20261111T153035Z-2fda74bc/summary.json
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.disk_face_flow.confirm
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.disk_face_flow.confirm \
  --check results/adsorption/disk_face_flow/confirmations/20261111T153229Z-96c67149/summary.json
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.disk_face_flow.plot \
  --summary results/adsorption/disk_face_flow/runs/20261111T153035Z-2fda74bc/summary.json \
  --confirmation results/adsorption/disk_face_flow/confirmations/20261111T153229Z-96c67149/summary.json
julia --startup-file=no --project=julia examples/adsorption/stationary_curved_wall/test_reference.jl
```

The spatial confirmation uses the original campaign selected in its note.
Replay a confirmation with `python -m examples.adsorption.disk_face_flow.confirm
--check <confirmation summary.json>`. The plot command checks that its
confirmation used the same campaign supplied with `--summary`.

The timestamped directories in `results/adsorption/disk_face_flow/` hold the
source and result hashes, commands, SI input parameters, liquid inventories,
concentration extrema and field snapshots. Sector fields are little-endian
Float64 cell averages in kg/m³ with shape (radial rings, azimuthal sectors)
and radial index contiguous. Angles increase counterclockwise from +x. Ring
coordinates and coefficients are in `rings.csv`, with 1-based radial indices
and positive volume flow along increasing angle. Snapshot times match the
archived reference exactly. Carrier flow and diffusive conductance are both
zero at the wall. Output and reference hashes can be checked without rerunning
Julia.
