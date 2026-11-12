# One-cell inlet–outlet inert-tracer benchmark

The [protocol](protocol.md) opens the axial seam of the archived
N96-derived N192 and N288 tracer meshes. It keeps the original
stock TRT flow, periodic transverse faces and solid mask. The Julia
solver reuses the bounded interior method and adds conservative
inlet and outlet terms. The [report](REPORT.md) records the comparison.

Run from the configured project root with Julia 1.11.6. The runner
uses the pinned parent hydraulic and tracer records. No rebuild
or dependency change is needed.

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.porous_outlet_tracer.run
uv run python -m examples.adsorption.porous_outlet_tracer.run \
  --check results/adsorption/porous_outlet_tracer/runs/20261111T194841Z-0489ef5c/summary.json
uv run python -m examples.adsorption.porous_outlet_tracer.plot
```

A fresh run prints its unique summary path for replay. The two-cell
boundary checks run before the porous computations. The N96 scout
sets the common physical end time for N192 and N288. Both cases
use the same 4 ms inlet pulse and 80 µs time step. Each case
archives its executable command, time history and final field.

The [face array contract](../porous_face_handoff/README.md) applies to
the binary arrays. In `boundary_rate.f64`, positive values denote
inlet volume flow divided by cell volume and negative values denote
outlet flow divided by cell volume. Both are in s⁻¹. Internal faces
retain the original orientations and flow values. Histories include
injected and recovered tracer in kg, liquid inventory in kg, interval
exit tracer mass and cumulative unnormalized first and second exit-time
moments. Replay reconstructs the cut faces and compares the saved
results without launching a tracer solve.

## N384 late-tail comparison

The [N384 report](N384.md) compares against the archived N288 outlet
case at the same 80 µs step and final time. It reads compact files
with one row per original N96 liquid voxel. Each row stores the three
signed incoming and outgoing face flows and the adjacent liquid row
IDs. Zero neighbor IDs identify walls and open faces.

```sh
uv run python -m examples.adsorption.porous_outlet_tracer.n384
uv run python -m examples.adsorption.porous_outlet_tracer.n384 \
  --check results/adsorption/porous_outlet_tracer/N384/20261111T211006Z-4ab30fe9/summary.json
uv run python -m examples.adsorption.porous_outlet_tracer.plot_n384
```

The runner checks concentration and outlet equivalence against the
archived solver before launching N384. Use a fresh run's printed
summary path for replay.
