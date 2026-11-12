# Tracer refinement at fixed OpenLB flow

The [protocol](protocol.md) specifies the half-step and doubled tracer
mesh cases for the archived N96 periodic FCC cell. The parent
[handoff report](../porous_face_handoff/REPORT.md) describes the saved
stock TRT field and bounded Julia solver. This refinement keeps the
original solid voxels, physical flow, diffusivity, and observation times.
The [report](REPORT.md) gives the numerical outcome.

Run from the project root on the configured macOS ARM64 environment.
The bounded solver uses Julia 1.11.6 through the archived project
configuration. No OpenLB build or new dependency is needed.

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.porous_tracer_refinement.time_half
```

Use the study directory printed by the half-step command for the next
two commands. The recorded execution used the following path.

```sh
study=results/adsorption/porous_tracer_refinement/runs/20261111T160638Z-156b2680
uv run python -m examples.adsorption.porous_tracer_refinement.mesh "$study"
uv run python -m examples.adsorption.porous_tracer_refinement.analyze "$study"
uv run python -m examples.adsorption.porous_tracer_refinement.analyze "$study/summary.json" --check
uv run python -m examples.adsorption.porous_tracer_refinement.plot --summary "$study/summary.json"
julia --startup-file=no --project=julia examples/adsorption/porous_face_handoff/test_bounded.jl
```

The mesh command saves refined faces and runs both grids at a common
physical time step when the fine stability limit requires it. The
[archived data contract](../porous_face_handoff/README.md) applies to the
binary face and tracer arrays. Fine liquid rows are ordered by original
liquid row, then by the lexicographic zero-based offsets of its eight
children. Positive flow follows the stored coordinate axis, including
periodic seams. Subface sums and the resulting local balances are
reconstructed during replay. Replay launches no Julia or OpenLB solve.

## N192 time-step check

The [separate N192 check](N192_TIME.md) reuses the recorded fine faces and
halves the effective step in each observation interval. Run from the
same configured project root.

```sh
uv run python -m examples.adsorption.porous_tracer_refinement.n192_time
uv run python -m examples.adsorption.porous_tracer_refinement.n192_time \
  --check results/adsorption/porous_tracer_refinement/N192_time/20261111T162208Z-42f74625/summary.json
```

Use the summary path printed by a fresh run to replay that run.

## N288 matched-step spatial check

The [N288 comparison](N288_SPACE.md) subdivides the saved N96 voxels
while using the archived N192 half-step time grid. Julia uses eight
threads after checking concentration equivalence against the serial
solver on saved N192 fields.

```sh
uv run python -m examples.adsorption.porous_tracer_refinement.n288
uv run python -m examples.adsorption.porous_tracer_refinement.n288 \
  --check results/adsorption/porous_tracer_refinement/N288_space/20261111T164622Z-a501da76/summary.json
```

The generated face arrays use the same binary data contract. Each
original liquid row contributes 27 child rows in lexicographic offset
order. The runner logs elapsed tracer steps and archives the throughput,
serial-equivalence and time-grid checks. Use the printed summary path
to replay a fresh run.
