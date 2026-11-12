# Engineering accuracy assessment, declared before new native execution

10 November 2026. The user accepts assessing numerical error by its effect on
engineering outputs. This is a new assessment, not a retrospective pass of the
column_transport protocol. Its one-part-per-million preservation check and all
its failed runs remain unchanged. The following tolerances are chosen project
accuracy targets, not experimental uncertainties or universal standards.

## Question and fixed physical inputs

Can the stock combined regularized TRT D3Q27 transport case support a bounded,
converged comparison of usable adsorption and pumping energy at the source flow
for both particle sizes? First quantify local and average concentration errors,
then measure an inert step response. No adsorption curves calibrate transport.

Use the completed column_transport hydraulic study: synthetic periodic FCC
spheres, nominal voidage .6, diameters .80 and .45 mm, 43.8 mL/min in an 11 mm ID
column, actual DWSIM water properties at 295.15 K. Use the converged final N32
and N48 fields and the documented kinematic-viscosity scaling. Retain molecular
diffusivity 1e-9 m2/s as the declared scenario. No new physical diffusion,
capacity, velocity, inactive fraction, density or per-column fit is introduced.
The batch-qualified particle model remains frozen. These are synthetic material
transfer scenarios, not validation of the Alhamed column experiment.

## Native execution and prospective holdouts

Reuse the existing, unchanged scalar application and stock operators. Build a
separate model-1 executable outside the vendor tree. N32 base dt gives maximum
lattice speed .02. Each diameter has its own actual native field and time scale.

1. Screening: one periodic uniform-concentration and one finite step-inlet /
   zero-gradient-outlet case per diameter at N32. Run for 16 nominal pore-volume
   transits, much longer than the previous two-transit preservation check.
   For a finite case use its actual liquid control volume and Q to define time;
   exclude inlet/outlet reservoir planes. Initial step concentration is zero;
   feed is 1 kg/m3. No reactions or adsorption.
2. Refinement, if screening is a usable bounded response: repeat both tests at
   N32 half dt and N48 with the SAME physical dt as N32 half dt. Thus the grid
   comparison with N32 half dt isolates space; the N32 comparison isolates time.
   Keep elapsed physical time equal to screening, even if voxel volume changes.
3. Longer-segment holdout, if refinement supports a closure: two repeated unit
   cells, same physical flow and Dm; same 16-transit observation window, no new
   calibration. If required, extend a short/long observation window to 32 or 64
   transits in a new run. Do not append an invented/extrapolated tail.

Each stage is explicitly executed, has its own unique run directory, records
actual commands/source and input hashes, and may finish with a failed result.
A failed upstream engineering screen prevents downstream qualification. Saved
checks execute no native or DWSIM simulator. No concentration clipping,
normalization of the physical inventory or hidden mass correction is allowed.

## Accuracy screens and interpretation

The uniform diagnostic should have max local departure <=2% of the initial
concentration; maximum average physical-liquid inventory change <=.1%; final
quarter drift of min/max <=.2 percentage points. These are candidate-selection
screens, not an error bound on breakthrough. Record wall-population storage
separately from physical liquid and global kinetic mass. Finite values and global
kinetic mass closure <=1e-8 are mandatory. Water density/compressibility and
voxel conventions remain separate possible error sources.

For an initially clean step: require concentration >=-1e-8 and <=1.02, and
absolute exit plateau within 1% of the independently computed water Q times
feed concentration. Report the unnormalized outlet flux and both cumulative
boundary fluxes. Mean and variance are finite-window EXIT-FLUX diagnostics;
normalization by the observed plateau does not waive the absolute-plateau check.
The final fraction must be >=.999; the last half-window contributions to first
and second raw moments must be <=1%. These are tail screens, not a mathematical
upper bound on all unobserved future solute. No transport coefficient is accepted
from incomplete/unstable tails or negative/overshooting responses outside limits.

For a candidate effective one-dimensional closure, preserve actual Q, physical
liquid volume and length. Fit only one dispersion from the short inert response,
with the correct fixed-concentration-inlet / zero-gradient-exit reference;
predict the other moment and the longer segment without refitting. Do not call
its exit mean a reduced physical pore volume. Require mean within 5%, variance
within 10%, profile RMSE <=.03, absolute plateau within 1%, and refinement changes
<=5% in mean / <=10% in variance. A scalar D is conditional on these screens;
paper-fitted breakthrough dispersion is not an independent oracle.

Only after these requirements pass, test the unchanged Julia film-particle
column and existing actual DWSIM pump duties. Target <=5% numerical sensitivity
in 5%-feed breakthrough service time, usable uptake and energy per on-spec volume
for BOTH diameters, including axial/radial/time refinement and the native closure
alternatives. Preserve equilibrium, batch diffusivity and full inventory. Record
the service threshold and initial conditions. A tested 5% budget is an engineering
acceptance for this synthetic scenario, not experimental validation. If the
native response fails, service metrics remain unavailable; exploratory scenario
sensitivity must not be presented as a pore-derived prediction.

## Evidence and scope

See the previous [strict report](../column_transport/REPORT.md) and the
[frozen particle report](../film_particle/REPORT.md). The assessment follows the
principle that numerical accuracy concerns the intended engineering outputs
([NASA verification overview](https://www.grc.nasa.gov/WWW/wind/valid/tutorial/overview.html)).
No new laboratory data, unavailable raw datasets or specimen reconstruction is
required. The original OpenLB tree, earlier source guards, process duties,
failed experimental comparisons and interface remain untouched.
