# N96 disordered-packing refinement

This case uses the accepted [N64 disordered packing](../disordered_packing/REPORT.md),
its original sphere-centre file, and the existing stock OpenLB 1.9r0
external application. The [protocol](protocol.md) defines the N96
geometry, Péclet matching and the 2-point grid test. [Results](REPORT.md)

From the project root, build the pinned application if required and
run the independent N96 native solution:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
uv run python -m examples.adsorption.disordered_packing.native build
uv run python -m examples.adsorption.disordered_refinement.native
```

Use the new cohort path printed by the native run as `$n96_run`.
The N64 face projection parity check can run while the native
lattice is settling. The remaining stages require the verified
N96 field.

```sh
uv run python -m examples.adsorption.disordered_refinement.projection "$n96_run" --parity
uv run python -m examples.adsorption.disordered_refinement.run "$n96_run" handoff
uv run python -m examples.adsorption.disordered_refinement.run "$n96_run" trace
uv run python -m examples.adsorption.disordered_packing.extend "$n96_run/disordered_N96" d800
uv run python -m examples.adsorption.disordered_refinement.run "$n96_run" seal
uv run python -m examples.adsorption.disordered_refinement.run "$n96_run/summary.json" check
uv run python -m examples.adsorption.disordered_refinement.plot "$n96_run/summary.json"
```

The run uses the same binary field and tracer contracts as the
[ordered/disordered study](../disordered_packing/README.md).
The interpolated warm field uses the same rows as the converged
native field and is archived independently. The comparison records
the N64 and N96 profiles, native and conservative-flow checks,
inert recovery, and time-step and spatial screens.
The added 0.80 mm ensemble checks the borderline time-step result
separately from the matched 8,192-marker grid comparisons.
