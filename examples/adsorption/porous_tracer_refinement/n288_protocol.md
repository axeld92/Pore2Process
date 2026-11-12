# N288 tracer mesh with the archived N192 time grid

Subdivide each liquid voxel of the saved N96 TRT field into 3×3×3
tracer cells using the same rectangular normal-flow reconstruction.
The 3×3 subfaces sum to the archived face flow within floating-point
roundoff. Keep the parent solid mask, total liquid volume, physical
diffusivity of 1e-9 m²/s and the analytical initial sinusoid sampled
at tracer cell centres. Check local continuity, impermeable solids
and all periodic plane flows against the archived throughput.

Match the archived N192 half-step time grid exactly, with 1,314 steps
in each observation interval and identical physical endpoints. Check
the fine-grid stability cap before integrating. Use the established
bounded flux-correction method and SSP RK3. For computational efficiency,
compute independent face contributions and gather them at cells in the
original accumulation order on separate Julia threads. Verify the
threaded stage and integration against the serial bounded solver on
the archived N192 initial, midpoint and final fields before N288.
Require bitwise identical concentrations and dissipation agreement to
1e-10 relative for four matched steps from each field.

Report relative variance and fundamental amplitude changes and the
shortest angular phase shift at both observation times. Use the same
5% variance/amplitude and 0.05 rad phase limits. Audit physical bounds,
liquid inventory and the estimated numerical share of variance loss
against the existing 15% limit. Compare final changes to the previous
N96–N192 matched comparison for scale. This is one further spatial
check for the reconstructed N96 flow and voxel geometry. Execute
one N288 tracer solve, no new hydraulic or full-bed case.
