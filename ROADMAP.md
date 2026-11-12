# Pore2Process roadmap

This roadmap separates locally demonstrated capabilities from work still to be
published or investigated. Status reflects the development workspace on
**10 November 2026**.

## 1. Establish the numerical and software reference

**Implemented and locally verified; source and checks are included.**

- Solve a synthetic capillary network in Julia from Python inputs.
- Check pressures, flows, conservation and scaling against independent
  analytical expressions.
- Return explicit SI arrays and permeability to Python; use the permeability in
  a separate Julia Darcy-column model.
- Inspect results in a local Panel interface and export completed-run snapshots.

This first reference deliberately excludes OpenLB and DWSIM. It establishes the
bridge, data contracts, numerical checks and interface conventions.

## 2. Connect pore-scale hydraulics to actual process equipment

**Implemented and locally verified; source and checks are included.**

- Run native OpenLB hydraulic cases using stock components and separate project
  folders, with at least TRT for porous-media simulations.
- Extract effective properties in Julia and calculate a process pressure load.
- Apply that load to actual DWSIM equipment, then solve, save, reload and check
  the flowsheet.
- Demonstrate that a changed pore-scale input changes the process result.

The first connected technical criterion is met locally for idealized hydraulic
cases. That establishes a working software connection; it does not validate a
real material.

## 3. Compare membrane-support designs at fixed filtration duty

**Locally demonstrated as a bounded synthetic/literature-informed study.**

The selected question is whether changing support geometry can reduce total
pumping energy while filtration flux and the effective resistance of the other
layers are held fixed.

The current work includes controlled support geometries, published-data
comparisons, explicit filtration equipment, crossflow recycle and separate feed
and circulation duties. An explicit synthetic star channel removes ambiguity
in the model's own area, perimeter and velocity definitions.

Evidence limits remain part of the result: effective non-support resistance is
not uniquely intrinsic selective-layer resistance; transferred retention is a
scenario; sustainable flux and physical fouling time at the new geometry are
unavailable.

## 4. Publish the reproducible implementation

**Initial source snapshot published.**

The repository now contains source, tests, dependency manifests, curated
evidence records and setup/launch/check commands. Full article files, local
environments, binaries and routine generated outputs remain outside version
control. Numerical, integration and browser verification are recorded separately.

See the [reproduction guide](docs/guides/reproduction.md) and
[publication record](docs/verification/publication.md). The initial source publication
preserves the original GitHub overview history; earlier development commits
remain in the local repository. Other operating systems and unattended full
native reproduction on a new machine remain to be checked.

## 5. Close the filtration hydraulic showcase

**Completed and closed for now; retained as a regression benchmark.**

The higher-Re extension checks 16 native D3Q27/TRT cases on the fully developed
laminar branch near Re=1000–2000. Spatial, temporal, domain, initial-state,
pressure-law, Mach and conservation checks precede scale translation. The
conditional operating map includes 91 actual DWSIM solve/save/reload cases.

At 12 LMH and Re=1388, widening the support from 4 to 6 µm lowers mean TMP
from 1.2211 to 1.0831 kPa, but saves 0.18166% of total shaft energy because
circulation accounts for 98.72% of the baseline duty. Membrane resistance
and process-energy dominance are therefore distinct.

Published capacity data, interpolation/monotonic policies and transfer
sensitivities remain explicit in [the final report](docs/verification/filtration/operating_map.md).
Further filtration physics is an optional future extension, not the active
next milestone.

## 6. Develop the isothermal adsorption case

**In progress; current source and comparison reports are included.**

The next research question is: at fixed adsorbent inventory and required outlet
quality, how does particle size affect usable treatment capacity and pumping
energy? This adds mass transfer and breakthrough modelling to the established
OpenLB → Python/Julia → DWSIM workflow. Source selection, model qualification
and process coupling belong to that case's own verification and publication.

Runnable cases now live under `examples/hydraulics/`, `examples/filtration/`
and `examples/adsorption/`. See the [layout guide](docs/architecture/project_layout.md) for
the updated commands and preservation of existing evidence.

The [native packed-bed tracer benchmark](examples/adsorption/packed_bed_tracer/REPORT.md)
now qualifies a bounded synthetic pore-to-macro contact handoff. One new native
flow case and eight scalar TRT cases pass; a single tracer-derived coefficient
predicts a longer segment while preserving actual Q and physical pore volume.
Higher-Peclet failures remain explicit. Next executable milestones are:

1. Qualify uniform-tracer preservation, curved-wall transport and a convective
   outlet at the selected application's Peclet range using stock OpenLB facilities.
2. Couple that qualified contact closure to the unchanged batch particle model,
   preserving source equilibrium and independently declared bed inventory.
3. Propagate changed geometry/particle inputs through service-time and resistance
   calculations into an actual DWSIM process case, with the process metric checked.

The current absorbing-end tracer result is not experimental adsorption validation
and is not directly substituted for a convective column's residence time.

The [column-transport follow-up](examples/adsorption/column_transport/REPORT.md)
has now executed 19 formal native scalar cases and three additional preservation
diagnostics. Free-flow convective-outlet controls pass. The tested porous
application configurations fail uniform preservation, including the best stock
regularized TRT candidate under smaller time steps and a finer mesh. A guarded
handoff preserves source equilibrium, batch diffusivity and declared inventory;
the failed actual native packet is rejected. The remaining adsorption/service
milestone requires resolving this transport qualification, not additional
per-curve fitted coefficients. Independent fixed-flow native hydraulics and
actual pumping duties are recorded separately in that case.

That hydraulic branch now passes 20 actual native flow/domain/initial-state
executions and five actual DWSIM case executions, including the selected
OpenLB → Julia → DWSIM route. At fixed source flow and declared dry inventory,
0.45 mm particles require 3.103 times the bed-friction pumping work of 0.80 mm
particles in the synthetic FCC scenario. This is a completed hydraulic result;
the adsorption service-time/energy tradeoff remains unavailable because native
solute transport is unqualified. Reacting and thermal extensions remain deferred.

The [engineering-accuracy assessment](examples/adsorption/transport_accuracy/REPORT.md)
has now completed 14 additional native controls. The previously reported
approximately 1% uniform-concentration drift passes the new engineering screen
over sixteen void times for both particle sizes. The obstacle is larger
nonuniform transport error: sharp/gradual fronts and a smooth periodic porous
profile violate their physical bounds. The latter removes inlet/outlet ambiguity;
paired free-flow analytical profiles pass. Next use this periodic failing case
to qualify local water/solute-flux and wall consistency before fitting a closure.
The proposed 5% service/energy budget is not yet assessed. Earlier strict
failures and the completed native/Julia/DWSIM hydraulic result remain unchanged.

A [bounded follow-up investigation](examples/adsorption/transport_accuracy/DIAGNOSIS.md)
adds four native controls. Diffusion without advection passes; half/tenth dt
and a matched-dt N48 case still fail the periodic profile bounds. Coordinates,
fluid masks, velocity conversion and TRT diffusivity assignment pass their audit.
No simple setup fix was identified. In accordance with the requested stopping
rule, **defer this high-Peclet pore-tracer branch** rather than expand the current
example into transport-method research. Failed data and the missing 5% service
accuracy qualification remain explicit.

The next practical isothermal study is a **declared column-contact sensitivity**
using the frozen batch particle component and verified pore-derived hydraulic
duties. Apply common, predeclared contact scenarios (for example 1, 4, 16 and
64 equal mixed sections) to both particle sizes, preserving equilibrium, full
inventory and the same film-correlation boundary. Do not fit each curve or
call those assumed contact patterns pore-measured dispersion. This can exercise
the process workflow and expose model dependence; it cannot supply experimental
column validation or a material-design recommendation. No new such scenario
study is claimed in the bounded diagnosis.

The user subsequently authorized a larger, finite
[transport root-cause investigation](examples/adsorption/transport_root_cause/REPORT.md).
It supersedes the earlier stopping decision for this diagnosis while preserving
that record. Matched native particle meshes and peak fields, smooth shear and
flat-wall continuum comparisons, actual population balances, grid alignment
and initial compatibility distinguish resolution from wall/coupling effects.
The stock half-way wall candidate made the porous excursions larger and was
rejected. The completed N96 particle run reduces the largest excursion by
15.6% from N64 but still exceeds the allowance 5.42-fold; final liquid mixing
also changes. This finite diagnosis qualifies no additional adsorption closure
or DWSIM prediction. The [independent follow-up](examples/adsorption/transport_verification/REPORT.md)
now reproduces both late wall-free stretching error and separate distortion
of uniform concentration during circular-wall rotation. The subsequent
[stationary wall check](examples/adsorption/stationary_curved_wall/REPORT.md)
reduces uniform drift by 86% but still fails uniform and varying tracer
accuracy. Screen any further stock setting against the independent controls
before another bed. The validated hydraulic and batch components are retained.

The [balanced OpenLB-to-Julia handoff](examples/adsorption/porous_face_handoff/REPORT.md)
preserved throughput and uniform tracer, and the
[one-cell outlet refinement](examples/adsorption/porous_outlet_tracer/N384.md)
qualified the short pulse at the archived operating point. A variance-fit
dispersion coefficient failed the [initial two-cell profile test](examples/adsorption/dispersion_holdout/REPORT.md).
The [complete two-cell holdout](examples/adsorption/residence_time_holdout/REPORT.md)
then verified the sharp rise with a frozen one-cell path model. The model
overpredicted exit-time variance by 19.78% and missed the 99% recovery
time limit. The [flow-derived exchange model](examples/adsorption/path_exchange/REPORT.md)
passed the complete two-cell comparison using the same one-cell calibration
and the saved TRT carrier flow. It reduced those errors to 2.67% and
0.61% while preserving the sharp rise. The zero-diffusion exchange control
still failed the two late metrics. Qualification is for the synthetic
two-cell inert-tracer case.

The [unfitted Figure 7 column comparison](examples/adsorption/exchange_column/REPORT.md)
now couples that exchange model to the batch-qualified particles. Both
predictions pass numerical refinement and conserve phenol. They fail the
breakthrough comparisons for 0.80 and 0.45 mm. Midpoints are 78.84 and
105.83 min versus 30.75 and 55.36 min in the digitized source. The usable
treatment capacity and energy per service cycle remain unqualified.

The [fixed-carbon attribution](examples/adsorption/exchange_attribution/REPORT.md)
found little midpoint change from Péclet matching at the original
allocation. Strong slow-path carbon enrichment and zero exchange together
reduced the midpoint errors to 7.98 and 4.16 min. The 0.80 mm midpoint
still misses its screen. The scenario identifies an interaction between
allocation and exchange, with the experimental contact pattern unresolved.

The [carbon contact study](examples/adsorption/carbon_contact/REPORT.md)
uses the saved TRT geometry to quantify external grain-wall exposure.
The slowest tenth by advection-only transit accounts for 77–88% of the
contact-weighted dry-carbon allocation. The slowest realized group
changes membership by 24–32% in one repeat. The synthetic geometry
supports strongly uneven contact with exchange. The carbon capacity
accessible to each path in the experimental bed remains unresolved.

The [unfitted contact-allocation rerun](examples/adsorption/contact_column/REPORT.md)
used the same dry mass and Péclet-matched exchange for both particle
sizes. Both show early phenol leakage within seconds. The 0.80 mm
midpoint is numerically sensitive near a plateau and remains much
earlier than Figure 7. The 0.45 mm midpoint moves later to 120.52 min.
The experimental breakthrough and adsorption service predictions
remain unqualified.

The [shared-carbon comparison](examples/adsorption/shared_contact/REPORT.md)
uses wall contact to weight film transfer area to a common particle
inventory in each axial slice. At fixed mass and Péclet-matched liquid
exchange, both diameters still show premature leakage and fail Figure 7.
The flow-weighted shared-carbon control also fails. The adsorption
service predictions remain unqualified.

The [absorbing-wall test](examples/adsorption/wall_uptake/REPORT.md)
estimated path-specific films on the saved TRT field. The wall-uptake
estimates were stable to halving the marker time step, but constant
path films with the fixed Péclet exchange overpredicted capture over
two repeats. Both unfitted Figure 7 predictions still showed premature
phenol leakage with the shared carbon inventory.

The [survivor-exposure holdout](examples/adsorption/survivor_exposure/REPORT.md)
reused the saved trajectories and held the aggregate Péclet-matched
exchange fixed. It predicts the lower second-cell uptake of first-cell
survivors for both particle sizes. The unfitted column predictions
still reach 5% of feed within five seconds, and both fail the
experimental breakthrough screens at fixed carbon mass.

The [eight-cell tracer holdout](examples/adsorption/repeated_geometry/REPORT.md)
passes for both Péclet-matched conditions. Roughly 94% and 89% of
initially fast tracer clear the aligned eight-cell field without
wall contact. Redistributing survivors within the original flow
classes slightly increases bypass at the same carrier exchange.
The effect of a physically different pore geometry remains open.

The [disordered packing comparison](examples/adsorption/disordered_packing/REPORT.md)
used a continuous eight-section TRT field at the same carbon inventory
and Péclet numbers. Initially fast wall-free passage falls by 17–19
percentage points at matched N64. The direction holds across the
tested grids and the magnitude remains grid-sensitive. The comparison
covers tracer transport and wall contact.

The [N96 disordered refinement](examples/adsorption/disordered_refinement/REPORT.md)
preserved the same packing and matched Péclet numbers. The maximum
N64-to-N96 tracer changes were 3.93 and 3.21 percentage points at
80 µs and exceeded the 2-point grid screen at both tracer time steps.
The refined disordered geometry still showed more wall contact
than the archived ordered N96 field.

## 7. Select a reacting or non-isothermal showcase

**Research direction; application and physical model not yet selected.**

Choose the material, reaction or biological system and process metric together.
Identify the quantities that must pass between scales and the published kinetic,
transport and thermal information needed to close the model. Supporting tools
such as Cantera can be evaluated when appropriate; they do not replace OpenLB
or DWSIM.

Build and validate the smallest connected case before adding parameter sweeps,
sensitivity analysis or optimization. Candidate areas include catalytic porous
media, bioprocess separations and electrochemical flow systems.

## Acceptance criterion for every end-to-end showcase

An output actually derived from OpenLB must influence an actual DWSIM process
calculation through the Python/Julia workflow. Changing a pore-scale input must
produce an explainable process-scale effect, supported by validated intermediate
results, explicit units, reproducible inputs and traceable assumptions.

Unknown experimental inputs must stay unavailable or become clearly labelled
scenario assumptions with sensitivity checks. Neither new experiments nor
reconstruction of an inaccessible specimen geometry is a gate for this project.
