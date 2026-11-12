# Prospective spatial and local-transport diagnosis

10 November 2026 local time. Authorized continuation of the bounded diagnosis;
the earlier stopping rule is superseded for this separate study. Preserve the
earlier sources, evidence and failed screens. No fitted dispersion or adsorption
prediction may be extracted from an unqualified transport calculation.

Study the d800 synthetic FCC case: particle diameter 0.8 mm, nominal porosity
0.6, side 0.0013891725543141453 m, actual imposed superficial flow corresponding
to 43.8 mL/min through the 11 mm reference column. Use independently converged
stock TRT hydraulic fields at each resolution, scaled to the same carrier
viscosity as before. Report the actual voxel porosity and actual flow; these
can differ slightly with mesh. No geometry interpolation is a final flow field.

Use stock CombinedAdvectionDiffusionRLBdynamics wrapping stock TRT, D3Q27,
stock full-way BounceBack, all-periodic boundaries, D=1e-9 m2/s, initial
C=0.5+0.25*sin(2*pi*x/L) kg/m3. Keep physical dt=1.9017455796408345e-5 s
and requested end=0.2147547125969688 s fixed; rounding to a full step is explicit.
Run N32, N48 and N64; consider N96 only if practical and useful after N64.
Added diagnostics must not alter populations or the numerical update. Check
N32/N48 reproducibility against the existing matched-dt cases to roundoff.

Record extrema at every time step and preserve the entire field at the largest
excursion, including fluid/solid labels and SI velocities. Classify peak-error
locations against neighboring solid links, with background wall proximity also
reported. Save physical-liquid concentration mean, standard deviation, and
first longitudinal Fourier coefficients: these are periodic-wave diagnostics,
not a column residence-time distribution. Compare these engineering quantities
as well as extrema across meshes. Do not infer a formal convergence order from
nonmonotone voxel geometries or unqualified asymptotic behavior.

Audit exact local post-collision population streaming at steps 1 and 2 and in
a separate deterministic prefix replay at the peak step. Recover unshifted
populations before positivity tests. Include numerical wall storage in kinetic
mass closure while reporting physical-liquid inventory separately. Solute
population fluxes and water population fluxes must be identified independently;
cell-centered velocity differences are not a water-population continuity test.
Additional constant-C controls may isolate uniform-preservation mechanisms.

Separate controls use a prescribed smooth divergence-free periodic shear,
without particles, with a converged Fourier continuum reference. They test
velocity variation in isolation, not physical curved-wall geometry. Water-flow
and scalar-population audits retain unchanged stock dynamics and operators.
No custom collision, descriptor, transport solver or vendor patch is authorized
by this investigation. A candidate stock-configuration correction requires a
mechanistic justification and fresh verification; no clipping or extra hidden
diffusion is an acceptable fix.

Screens: finite, normalized total kinetic mass residual <=1e-8, concentration
in [0.23,0.77] kg/m3 for the sine wave. A short periodic pass does not qualify
the previously failed open-boundary front, long-time dispersion or adsorption
service prediction. Seek the root cause; report unresolved limitations honestly.
