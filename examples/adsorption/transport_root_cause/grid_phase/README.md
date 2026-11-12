# Half-cell grid-phase diagnostic

Translate the same periodic FCC spheres and initial concentration profile
by half a cell in x, independently converge the shifted hydraulic solution
at the same physical flow, and repeat the stock scalar calculation. This
tests sensitivity to voxel-wall alignment; it does not create a new material
or provide a generally valid correction.

The frozen [protocol](protocol.md) specifies the translation, fixed-flow
root, parent reference, tolerances and peak-population observation.

From the repository root:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.grid_phase.run
uv run python -m examples.adsorption.transport_root_cause.grid_phase.run --check
uv run python -m examples.adsorption.transport_root_cause.grid_phase.postcheck
uv run python -m examples.adsorption.transport_root_cause.grid_phase.plot
uv run pytest tests/adsorption/test_transport_grid_phase.py -q
```

The run requires the previously completed guarded unshifted N32 and hydraulic
parents. Their exact paths and hashes are recorded in each new run. Native
water and scalar processes execute sequentially with serial CPU_SISD builds.
Builds and archived binaries stay outside the original OpenLB distribution.
Full outputs and provenance are in
`results/adsorption/transport_root_cause/grid_phase/`.

Coordinates are metres with zero-based indices at (index+0.5)dx. The
translation is delta=dx/2; the initial sine is evaluated at x-delta.
Concentration is kg/m3, velocity m/s and signed sphere-surface gap metres.
Surface normals point from the nearest periodic sphere image into the fluid;
normal_x near +1 identifies a downstream pole. The analytical distance
utility audits locations and does not construct simulation geometry.

The postcheck reads actual initial native populations and verifies the
intended scalar phase. Conservation checks use actual population links and
include the explicitly identified full-way wall population inventory.
No new adsorption fit, dispersion coefficient or process prediction is made.

The completed [report](REPORT.md), [captured replay summary](reference.json)
and [initial-profile reference](initial_profile_reference.json) distinguish
the verified alignment sensitivity from an unresolved transport failure.
