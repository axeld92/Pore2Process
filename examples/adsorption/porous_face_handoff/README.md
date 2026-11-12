# Saved OpenLB TRT flow to Julia tracer transport

This study transfers three saved stock OpenLB 1.9r0 flow fields for the
same 0.8 mm sphere FCC cell to Julia. Python prepares the periodic face
graph and checks the data contract. Julia balances carrier flows with an
axial and transverse throughput constraint, then integrates transport.
The [centered protocol](protocol.md) and [bounded protocol](bounded_protocol.md)
define the screens and scope. Results are in the [report](REPORT.md).

The original hydraulic archives are those pinned in `run.py` under
`results/adsorption/`. Their hydraulic provenance is documented in the
[flow-refinement guide](../transport_root_cause/hydraulics/README.md).
No fresh OpenLB executable or Python/Julia dependency is needed.

Run from the project root on the configured macOS ARM64 environment,
with the Julia 1.11.6 runtime selected by project bootstrap.

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/adsorption/porous_face_handoff/test_transport.jl
julia --startup-file=no --project=julia examples/adsorption/porous_face_handoff/test_bounded.jl
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.porous_face_handoff.run
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.porous_face_handoff.bounded
OPENBLAS_NUM_THREADS=1 uv run python -m examples.adsorption.porous_face_handoff.plot
```

The two runners accept `--check <summary.json>` for read-only numerical
replay. The bounded runner accepts `--centered <summary.json>` and defaults
to the latest completed centered cohort. The plot command accepts the
matching `--centered` and `--bounded` summary paths.

The archived CSV rows set the one-based liquid-cell IDs used by
`left.i32` and `right.i32`. Each face points in the positive direction
of its `axis.u8` value (1=x, 2=y, 3=z), including periodic seams. Solid
faces are absent. Face flows are little-endian Float64 m³/s and may be
negative in this orientation. Concentration arrays are little-endian
Float64 kg/m³ in the original liquid-cell row order. Histories record
physical times, concentration extrema, liquid mean and integrated
molecular variance dissipation. The summaries include exact native
sources, viscosity conversion, flow audits, solver checks and hashes.
