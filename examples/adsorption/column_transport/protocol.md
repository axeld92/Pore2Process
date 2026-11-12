# Column transport and adsorption: pre-evaluation protocol

Prepared 10 November 2026, before the new native trials. Preserve the earlier
absorbing-end benchmark, failed column comparisons and batch-qualified particle.

## Question and boundaries

At fixed adsorbent inventory, flow and outlet limit, how does particle size
change service volume and pump energy? Qualify native transport before using a
coefficient for adsorption. Use stock OpenLB facilities, at least TRT, external
application/build folders and the installed examples as the API reference.

First diagnose uniform-tracer preservation in the existing FCC velocity field
with periodic ends and curved bounce-back walls. Then test the stock
zero-gradient outlet used by the upstream microMixer3d. Distinguish molecular
diffusion, effective longitudinal dispersion and particle loading diffusivity.
No clipped native concentrations, invented tails or fitted inventory/velocity.

## Numerical screens

Native kinetic mass residual <=1e-8 of transferred/initial tracer; finite values;
nonnegative concentration >=-1e-8 and upper bound <=1+1e-6. For an initially
uniform periodic tracer, max absolute departure <=1e-6. Hold physical parameters
fixed during mesh/time checks. For a qualified nonuniform breakthrough, mean
change <=5%, variance <=10%, tail contributions <=1%, plateau >=.999. Validate
uniform-flow inlet/outlet against an independent conservative Julia ADE.

Application free-liquid diffusivity starts at the previously declared phenol
scenario 1e-9 m2/s (not a new measurement). Audit the published flow/particle
range and any slower demonstration regime separately. Do not present a more
diffusive artificial tracer as qualification at the paper's operating conditions.
If a stock method fails, retain the failure and investigate another existing
OpenLB facility; do not silently relax these screens.

## Adsorption and process comparison

Use the existing FilmParticle implementation unchanged, source Langmuir
equilibrium and batch-only fitted particle diffusivity. Keep printed bulk
density, derived bed mass and scenario FCC inventory separate. No independent
fit to each column curve. The experimental comparison uses both Figure 7
diameters and the same calibration boundary. A failed comparison is a result,
not experimental validation.

A process demonstration is admissible only for a qualified native transport
regime. Preserve physical Q, pore volume and dry mass. Derive a common closure
from short native tracer data; test a longer segment without fitting it. Check
Julia mass balance, numerical refinement and limiting cases. Actual DWSIM
solve/save/reload must use native hydraulic resistance; Julia owns transient
adsorption. Report service volume, usable uptake, shaft power and energy per
treated volume. Label synthetic packing and transferred particle properties.
Do not claim an experimental process recommendation from a failed column model.

## Independent hydraulic work, before its new execution

While scalar qualification is unresolved, complete the independently testable
hydraulic branch at the paper's flow and both particle diameters. The geometry
remains synthetic FCC at nominal voidage .6; it is not the paper's specimen.
Use source dry mass 7.032455 g as a declared scenario inventory, particle-envelope
density 740/(1-.427) kg/m3, 11 mm column ID, and adjust bed length to preserve
that mass at FCC voidage .6. Do not pretend FCC voidage is the paper's .427.
Use DWSIM water properties at 295.15 K, 43.8 mL/min and chosen pump efficiency .75.
Match Reynolds number when translating native nu=1e-6 m2/s flow to that carrier.
Native pressure is adjusted to impose the same flow, not fitted to adsorption.

Native N=32/48 for both sizes; pressure at fixed flow changes <=5% on refinement;
half time step <=2%; target flow relative error <=1e-4. Retain established native
mass <=1e-6, Mach <.05, density departure <1e-4 and stationarity <=1.0001e-7.
Report operating-point resistance, not a flow-independent permeability at this
Reynolds number. Actual DWSIM pump work must match Q*deltaP/eta within 1e-8,
mass within 1e-10 relative and enthalpy balance within .001 of shaft work;
save/reload/re-solve must pass. These pumping results can be published separately
if scalar qualification fails. Service time, treated volume and cumulative
service energy then stay unavailable; no failed scalar coefficient reaches
the particle model or process calculation.
