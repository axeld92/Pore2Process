# Ordered and disordered packing at matched Péclet numbers

The [protocol](protocol.md) defines the synthetic packings, hydraulic
checks and eight-cell tracer comparison. The case builds OpenLB
1.9r0 in a project-owned directory using the existing external
OpenMP build helper. It requires the pinned [ordered tracer archive](../repeated_geometry/REPORT.md),
Julia 1.11.6 and the project's Python environment. [Results](REPORT.md)

From the project root, prepare a result folder and build:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.disordered_packing.native prepare
uv run python -m examples.adsorption.disordered_packing.native build
```

Use the path printed by preparation as `$packing_run` in the commands
below. The high-Mach pilot fields initialize the lower-time-step
solves. Only fields accepted by the face-transfer step reach the
tracer.

```sh
for packing_case in "disordered 48" "disordered 64" "ordered 64"; do
  read -r label grid <<< "$packing_case"
  uv run python -m examples.adsorption.disordered_packing.native run "$packing_run" "$label" "$grid" --pilot
  uv run python -m examples.adsorption.disordered_packing.native run "$packing_run" "$label" "$grid"
  uv run python -m examples.adsorption.disordered_packing.mesh "$packing_run" "$label" "$grid"
done
uv run python -m examples.adsorption.disordered_packing.tracer "$packing_run" ordered 64
uv run python -m examples.adsorption.disordered_packing.tracer "$packing_run" disordered 48
uv run python -m examples.adsorption.disordered_packing.tracer "$packing_run" disordered 64
uv run python -m examples.adsorption.disordered_packing.extend "$packing_run/ordered_N64" d450
uv run python -m examples.adsorption.disordered_packing.extend "$packing_run/disordered_N64" d450
uv run python -m examples.adsorption.disordered_packing.compare "$packing_run"
uv run python -m examples.adsorption.disordered_packing.compare "$packing_run/summary.json" --check
uv run python -m examples.adsorption.disordered_packing.plot "$packing_run/summary.json"
uv run pytest -q examples/adsorption/disordered_packing/test_mesh.py
julia --startup-file=no --project=julia examples/adsorption/disordered_packing/test_track.jl
```

The OpenLB binary exports little-endian Float64 rows with zero-based
voxel IDs, three velocity components in m/s and density in kg/m³.
Face exports use one-based compact Julia IDs, signed rates in m³/s
and axis indices 1–3. Each marker row stores eight absolute
first-passage times and eight phase-specific first-wall-hit times in
seconds. NaN marks an unobserved event. The replay verifies stored
files and reconstructs the comparison without launching simulations.
The additional ensembles check the borderline 0.45 mm time-step
comparisons and leave the matched 8,192-marker comparison intact.
