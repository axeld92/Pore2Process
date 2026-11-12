# Hydraulic initial-state and domain checks

Declared before these additional native executions, 10 November 2026. Verify
the higher-Reynolds 0.80 mm operating point at N32 beyond its already executed
grid/time/flow screens. Use the same native gradient, viscosity, particle size,
TRT, forcing and discretization. A cold one-cell case must agree with the
completed warm-start field's superficial velocity within 1e-5 relative. A
two-cell x-periodic domain initialized with a 5% axial sinusoidal velocity
perturbation must return to the one-cell bulk state within 1e-5 in both mean
velocity and L2 inter-cell velocity difference. Preserve the existing 1e-7
stationarity, 1e-6 mass, .05 Mach and 1e-4 final density-departure screens.

All mesh/geometry/dynamics/forcing/functors remain stock OpenLB. The separate
application extends the periodic domain and initial field, not the solver.
No scalar, adsorption, temperature or experimental calibration is introduced.
If these screens fail, the final Julia-to-DWSIM hydraulic handoff is unqualified
until that failure is resolved. This does not resolve finite column walls or
distributors, or prove stability for arbitrary domains/disturbances.
