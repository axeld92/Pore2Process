# Conditional extension to N96

10 November 2026, declared before the N96 FCC hydraulic or scalar execution.
The full wall-free shear series gives normalized profile RMS errors 38.58%,
8.96%, 2.23% at N32/48/64. The independently equivalence-checked thin N96
control gives 0.3634%, inside the declared 1% profile target, without changing
the physical diffusivity. This provides a concrete reason to test N96 in the
particle geometry; three finer-but-underresolved scalar cases alone would not
answer whether resolution resolves the porous peaks.

Reuse the frozen hydraulic runner and its serial/OpenMP cross-check. Start a
new actual native N96 FCC field at rest, using the accepted N64 pressure
gradient only as an initial pressure-root guess. No interpolated coarse field
is a final solution. Require the same strict fixed-flow and other hydraulic
screens. Use four stock OpenMP threads, with external source/build/output
folders. Stop only the redundant serial N64 execution once a qualified field
is safely preserved, retain its partial outputs and record the interruption.

Run the same frozen instrumented scalar application at N96, D=1e-9 m2/s,
dt=1.9017455796408345e-5 s, and the same common physical observation window.
Report changed voxel porosity, exact elapsed time, maximum lattice speed,
extrema, physical-liquid moments and peak location. Do not call this formal
asymptotic convergence if the four-grid trend does not support it.

If N96 passes or strongly improves the particle range screen, test a smaller
physical time step before claiming a transport correction. A passed periodic
case still cannot waive the failed open-boundary fronts, supply dispersion,
or authorize a calibrated adsorption/process prediction. Do not change the
physical diffusivity, clip concentrations or fit a closure to rescue a run.

This is one additional FCC grid, not an unlimited collision-method development
or a commitment to larger grids. Interpret any remaining failures with the
separate wall-free and flat-wall controls, exact population balances and
stock-equilibrium candidate. Preserve the finite study's unresolved questions.
