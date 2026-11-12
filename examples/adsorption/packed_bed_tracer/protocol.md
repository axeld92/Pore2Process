# Packed-bed passive-tracer protocol, revision 3

Synthetic ordered FCC geometry, no adsorption, heat or reaction. Reuse accepted
native D3Q27/TRT flow fields at N=24/32. Generate and verify a new N=16 field
with the unchanged existing native bed application. All additions, builds and
results stay outside the vendor tree. New scalar cases use stock D3Q27<VELOCITY>,
AdvectionDiffusionTRTdynamics, sphere indicators, BounceBack and Dirichlet
boundaries. No custom collision, descriptor, geometry or boundary processor.

## Current qualification scope, recorded before revision-3 study

A finite segment cut from periodic fully developed flow has transverse
periodicity, inlet reservoir C=1 kg/m3, absorbing outlet C=0, clean initial
state. The signal is net EXIT TRACER FLUX, not concentration at the absorbing
plane or a real column RTD. Diagnose both external boundaries in uniform flow.
The kinetic control volume excludes both reservoir-node planes. Its surfaces
are the two streaming-link cuts; length=(nx-2)*dx. Reservoir-node separation
is cells*side=(nx-1)*dx and remains separately recorded. Julia uses the control-
surface length consistently. No fitted grid-dependent boundary offset.

Primary free diffusivity is a chosen 2e-6 m2/s; sensitivity is 4e-6. These are
verification scenarios, not phenol properties. Primary tau_odd=1 and magic=0.3
produce UNEQUAL tau_even=1.1; time sensitivity tau_odd=0.75 halves dt. Neither
channel is replaced by BGK. Physical D=cs2*(tau_odd-0.5)*dx^2/dt; cs2=1/3.
The positivity-favouring relaxation region is chosen after the failed trials
below, without changing the numerical acceptance thresholds or fitting data.

Run three bed resolutions N=16,24,32 with one unit-cell length; compare 24 to32.
Run N=32 at half dt and twice the bed length. Run N=24 at twice D. Uniform free-
flow references use N=16 and32, one cell, speed=0.0005 m/s. End windows are35 s
per unit-cell length. Save all cases, including any failed checks.

Measure postcollision link flux with the stock public collide/AndStream split.
Use stock lbm<D>::computeRho for KINETIC inventories, since BounceBack's
computeRho() reports prescribed wall density, not population mass. Include and
separately report temporary wall-population storage in the exact kinetic mass
balance. It is neither physical pore volume nor adsorption. No mass correction.

Screens: normalized kinetic mass residual <=1e-8; bulk concentration >=-1e-8
and <=1+1e-6; final normalized exit flux >=0.999; late-half contributions to
first and second RTD moments <1%. Uniform-flow analytical mean error <=2%,
variance error <=5%, Julia reference profile RMSE <=0.02. Final bed grid/time
changes: mean <=5%, variance <=10%. Changing bed length is a physical boundary
sensitivity, never numerical grid convergence. Export raw quantities as well
as normalized comparisons.

Independently derive the absorbing-end ADE moments and verify a conservative
Julia finite-volume reference. Derive common effective ADE coefficients from
ONLY the N=32 short-segment INERT TRACER moments. Predict the twice-long N=32
segment with those same coefficients, no refit. Reduced-model screens are
profile RMSE <=0.03, mean error <=10%, variance error <=20%. Keep closure and
native checks separate. An out-of-family moment ratio yields unavailable
closure; do not add an extra fit factor. Keep effective velocity/volume separate
from physical pore velocity/inventory. Never substitute them silently into the
batch-qualified particle model or the printed bed density.

Unique runs: results/adsorption/packed_bed_tracer/runs/<id>. Record native
requests, source/binary/input hashes, archived source bytes, fields, fluxes,
inventories, diagnostics, versions and timings. Checker validates saved native
data and replays Julia; it launches zero native simulations. No new DWSIM or
browser result is claimed. The experimental adsorption comparisons remain failed.

## Retained exploratory history

Revision1 used D3Q7, D=1e-7 and magic=3/16 at N=24/32/48 and two/four cells.
The first eight-case run is retained at20261110T125302Z-7f3c689f. Its wall
inventory diagnostic incorrectly used the prescribed wall density, and the
high-Peclet regime failed local concentration boundedness. Those results are
not an accepted native qualification despite favourable aggregate exit moments.

Separate D3Q27 trials corrected kinetic inventory and conserved mass, but
D=1e-7 still overshot the inlet concentration near curved walls. D=1e-6 removes
that upper-bound failure; the original scalar magic/rates gave tiny negative
startup excursions. Revision2 selects unequal, positivity-favouring TRT rates
and a deliberately more diffusive bounded regime. Initial halo and Dirichlet-
density initialization mistakes were corrected in project code; vendor files
were never changed. These changes do not qualify the original higher-Peclet
regime. Future extensions must revalidate it explicitly.

Revision2 run20261110T130633Z-beb0a96b passed short-segment bounds, grid/time
checks and the common-coefficient longer-response comparison, but its longer
D=1e-6 case exceeded the concentration upper bound (1.001329) and its late
second-moment contribution was1.326%. It remains unqualified. Revision3 moves
the whole bounded verification regime to D=2e-6, with D=4e-6 sensitivity; no
acceptance limit is relaxed and no long-response coefficient is fitted. This
is regime qualification, not a blind material-selection experiment.
