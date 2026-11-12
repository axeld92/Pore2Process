# Refinement of tracer transport on the saved N96 flow

Use the archived bounded N96 tracer result, the archived N96 TRT hydraulic
field and its exact balanced faces. Keep the same single periodic 0.8 mm
sphere FCC unit cell, liquid mask, flow, physical diffusivity of 1e-9 m²/s,
initial sinusoidal concentration and observation times. Use the existing
bounded Julia solver without changing its limiter or time integrator.

First halve the original N96 step cap with time factor 0.5. Check whole
history bounds and mass, the profile RMS at both times, and how final
variance and numerical dissipation change relative to the archived N96.

For the spatial check, subdivide each liquid voxel into eight equal cubes.
The N192 tracer grid keeps the N96 voxel wall and hydraulic field. Within
each coarse voxel, interpolate the signed normal velocity linearly along
the normal coordinate between its two face values (rectangular RT0).
Allocate a quarter of a coarse face flow to each of its four fine faces.
No carrier flow crosses a solid face. Audit each coarse-face sum, every
fine axial plane, local fine-cell continuity and the unchanged liquid
volume. Use a common physical time cap on the two meshes, the smaller of
the N96 half-step cap and the fine-grid stable cap. If that cap subdivides
either observation interval differently from the N96 half-step run,
repeat the N96 comparison at the common cap. Initialize refined cells
from the same analytical sine at their centres.

Use the established limits at both observation times, 0.1% of initial
amplitude for the time-step profile change, 5% relative for spatial
variance and fundamental amplitude, and 0.05 rad for spatial phase.
Use 1e-8 relative for inventory and 1e-8 kg/m³ for concentration bounds.
Report signed numerical variance loss (actual half-variance loss minus
integrated molecular dissipation) in concentration-squared units and as
a percentage of total loss. The existing numerical contribution screen
is 15%. Report whether it fell even if it still misses that screen.
Compare fine fields averaged onto the coarse voxels as a supplemental
spatial profile diagnostic. A single doubling at fixed hydraulic mesh is
a sensitivity check and cannot establish a continuum convergence order.
No new OpenLB, full-bed, or DWSIM run is part of this study.
