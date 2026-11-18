# Fixed-mass tests for the phenol breakthrough discrepancy

Use the archived N60 film, transfer rules and adsorption predictions as the
periodic reference. At 25 cm keep 17.581137888 g dry carbon, the original bed
liquid volume, 150 mg/L feed, 0.80 mm particle diameter, batch model and
additive access resistance. For other depths use the same reported bulk
density and column cross section. Compare with marker data at the original
conditions and through the same reported windows.

First test a parallel outer annulus and core at 62.5 mL/min with zero exchange.
Use an annulus one half particle diameter thick. Try wall voidages 0.47, 0.50
and 0.55 as scenarios for a higher local voidage, conserving total solid and
void volume. Solve for the common axial Ergun pressure gradient that produces
the target flow. Allocate carbon and external area by local solid volume.
Use the same N60 particle film in both zones. The two independent branches
are a zero-exchange sensitivity. The identical-voidage limit reproduces the
original prediction. The purpose is to check the direction and useful range
of the effect before conductive exchange or resolved geometry is applied.
[Radial voidage measurements](https://doi.org/10.1002/aic.690080319)
· [Ergun pressure law](https://www.industrydocuments.ucsf.edu/tid/efk76a99/pdf)

The confined-flow check uses stock OpenLB TRT and a nonadsorbing tube
boundary. Compare the radial water-flux distribution with the independent
published MRI benchmark on its own geometry. For the literature-based
11 mm tube, derive zone flow, void and carbon contact from a short synthetic
packing. Verify solvent conservation, inert tracer recovery and clean
carbon contact before using empirical interzone exchange.

Compare periodic flow with confined flow and a shared carbon history with
separate histories for the two flow pathways. Use matched transit/contact
inputs and retain both direct and slow carbon cohorts. The parcel partition
for periodic flow uses only saved tracer exposure. If late uptake still
misses, test particle-kinetic alternatives qualified independently on
withheld batch sizes. Report the density-derived dry mass separately from
the fitted equilibrium and bath-kinetics data.

Freeze every forecast before scoring. Use 10, 15 and 20 cm at 23.3 mL/min
to check how the rise widens with depth. The 10 and 15 cm series are reserved
for a shared calibration only if used. The 20 cm and both 25 cm cases
remain prediction checks. Count the common 25 cm low-flow curve once.
Figure 7 is reserved for a two-size transfer test with size-matched inputs.
Retain the existing screens of ±25% for the 5% crossing, ±20% for the midpoint,
0.10 outlet/feed RMSE and retention within the digitization reading range.
Also report area under the input-minus-output curve after the experimental
midpoint. The reading ranges characterize figure digitization.
[Alhamed Figures 4, 5 and 7](https://doi.org/10.1016/j.jhazmat.2009.05.002)
