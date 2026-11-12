# Half-cell FCC translation: frozen prospective diagnostic

This tests sensitivity to the alignment of a fixed continuum geometry with
the voxel lattice. It is not a new material design, a fitted geometry, or a
proposed universal fix. Keep N=32, the FCC side and 0.8 mm sphere diameter,
nominal porosity 0.6, periodic domain and physical superficial flow unchanged.
Translate every sphere centre in +x by delta=dx/2. Also translate the entire
initial scalar profile by +delta: C=0.5+0.25 sin[2 pi(x-delta)/H]. Use the
same translated phase for reported streamwise Fourier moments. In continuum
physics these are exactly translated copies of the unshifted problem.

The grid coordinates remain (i+0.5)dx. Translation therefore changes the
discrete solid/fluid classification and staircase wall location. Record
actual voxel porosity, hydraulic resistance, maximum speed and inventory
rather than imposing the old values. Do not interpret their differences as
material improvement. OpenLB stock IndicatorSphere3D alone builds geometry;
the separate analytical distance utility only inspects exported points.

Independently converge the shifted water solve from rest using the original
stock D3Q27 TRT/Kupershtokh forcing, BounceBack and official-example lifecycle.
At each pressure-root iteration require the existing stationarity, plane
mass-flux, Mach and density checks; converge superficial flow to 1e-4 relative
using the unchanged bounded root helper, at most 8 native iterations. The
first pressure-gradient estimate comes from the accepted unshifted N32
point. Only shifted same-grid converged fields initialize later root steps.
Never accept a translated/interpolated old field as the final water solution.
Preserve the exact stock-population water balance observation at convergence.

Use the separately converged shifted field for the stock
CombinedAdvectionDiffusionRLBdynamics + AdvectionDiffusionTRTdynamics scalar
case. Keep D=1e-9 m2/s, physical dt=1.9017455796408345e-5 s and requested
end=0.2147547125969688 s, matching the parent instrumented N32 case. The
actual ceil-rounded observation time is recorded and must match the parent.
Both water and scalar applications use serial CPU_SISD builds in distinct
project folders; no stock OpenLB or pre-existing project source is modified.

Save all-step extrema and liquid/wall inventories, peak and final fields,
Fourier moments, and exact scalar population audits at steps 1/2. Once the
peak step is known, run one same-input prefix through that step with the
read-only exact population observation enabled there. Require its peak
field and common inventory prefix to reproduce the full run before using
its population audit. This prefix is instrumentation, not a new calibration.

Compare against the guarded unshifted instrumented N32 case: extrema and
their time, volume-average moments, outside-bound fractions, actual porosity
and pressure gradient. At each concentration extremum report the nearest
analytical sphere-surface gap in metres and grid cells and outward normal
relative to +x. Verify exported fluid masks against analytical sphere
classification and exact water/scalar population balances, including
full-way BounceBack population storage. These are actual population fluxes,
not a centred velocity-divergence proxy.

A strong change establishes grid-phase sensitivity, not correctness of the
better-aligned case. Passing global conservation does not establish bounds
or transport accuracy. No dispersion coefficient, adsorption service result
or DWSIM prediction is qualified. No clipping, extra diffusion, dynamics
change, fitted capacity or new experimental input is permitted.

Archive source snapshots, parent summary hashes, inputs, commands, binary
hashes, native outputs and independent reanalysis. Old numerical evidence
and the completed nine-run shear study remain unchanged.
