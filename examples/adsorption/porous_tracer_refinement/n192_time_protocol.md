# N192 time-step check at fixed tracer mesh

Reuse the exact N192 face arrays, liquid-cell order, tracer initial field,
Julia solver, physical diffusivity and observation times from the
completed fixed-flow refinement. Keep the N96 hydraulic geometry and
the eight tracer subcells per liquid voxel. Halve the effective step in
each observation interval by doubling its integer step count from 657
to 1,314. Use a time-step cap that selects those counts for both
intervals. Verify that every second new time equals an archived time.

Report absolute and relative changes in concentration variance, signed
shortest-angle phase shift of the longitudinal fundamental, amplitude
change, and full-profile RMS change relative to the initial 0.25 kg/m³
amplitude. Use the existing 0.1% temporal profile criterion. Compare
variance and phase changes with the preceding spatial changes for scale.
Audit physical bounds, inventory and the existing numerical variance
budget at both times. Preserve the earlier run and its qualification
status. Run one new Julia tracer case. The N192 mesh, saved OpenLB
hydraulics and full bed stay outside the refinement scope.
