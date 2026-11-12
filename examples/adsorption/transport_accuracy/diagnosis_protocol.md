# Bounded investigation of the periodic porous-profile failure

10 November 2026, declared before these trials. The user requests a limited
investigation, then documentation and progress elsewhere if no useful correction
is found. Preserve all earlier examples, sources, thresholds and failed data.

Focus on the .80 mm application case and its smooth periodic initial condition
C=.5+.25*sin(2*pi*x/L). Reuse the EXACT stock combined regularized TRT D3Q27
executable and source guarded by the completed profile study. Make no vendor,
dynamics, boundary, descriptor, geometry-class or numerical-update change.
Keep Dm=1e-9 m2/s, geometry, physical flow and elapsed time fixed except for the
explicit no-advection control. All new fields/runs are separate files.

Execute at most four initial native diagnoses:

1. N32, zero velocity at every fluid node: same actual geometry, Dm, initial
   profile, dt and elapsed time. The provided hydraulic CSV is copied with only
   velocity columns set to zero. This is a diffusion-only diagnostic, not the
   application flow or a new measured property.
2. N32, half dt, same native velocity and physical elapsed time.
3. N32, one-tenth dt, same native velocity and physical elapsed time.
4. N48, using its separately converged matched-flow native field, at the SAME
   physical dt as case 2 and same elapsed time. Compare case 4 with case 2 to
   assess space/voxel geometry; compare cases 2/3 with the original N32 to
   assess time. Do not call a two-grid result a formal error estimate.

Keep the existing profile screen C>=.23, C<=.77 (a .02 kg/m3 engineering
allowance around the continuum .25..75 bounds), finite values and kinetic mass
residual <=1e-8. Report actual timing and all extrema. Pure diffusion must also
be compared with its initial range; there is no analytic exact sphere-bed field.

Read-only audit: confirm native/scalar node coordinates and fluid masks, units,
velocity scales, physical Q, weak-compressibility magnitude, and stock TRT's odd
diffusive/even kinetic relaxation assignment. Calculate an explicitly labelled
Cartesian finite-volume velocity-divergence PROXY, with zero face flux across
voxel solid faces. This is not OpenLB's population mass flux, a replacement flow
solver, a correction or proof of local fluid continuity. Inspect where bounds
fail in saved fields and whether they are next to solid voxels. Simple diffusion
length/grid-spacing comparisons are scale estimates, not mesh requirements.

Stopping rule: if no clear application setup defect emerges and boundedness
still fails under the matched refinement, record a limitation of this tested
high-Peclet configuration and stop this investigation. Do not launch a new
collision-method development project. A single stock-configuration correction
may be tested if a concrete defect is found; acceptance would require new
verification for both sizes. Do not fit dispersion or adsorption to failed data.
No new experiments, unavailable raw data, reactive extension or system install
is needed. The accepted hydraulic/process and batch components remain usable.
