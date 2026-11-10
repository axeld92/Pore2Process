# Isothermal adsorption: first connected synthetic example

First connected synthetic release verified on 10 November 2026. Numerical and
software checks pass; experimental column validation remains pending. The new
[Alhamed PDF comparison](alhamed2009_comparison.md) now supplies a separate
date-stone-carbon calibration/holdout test: evidence and numerical replay
pass, but the current reduced model fails the physical comparison.

The later [finite-bath and film/particle alternative](film_particle_adsorption.md)
now passes the withheld batch-size loading screen while still failing the
column transfer. It is a separate numerical module; the synthetic loop below
and its source snapshots remain preserved.

The hydraulic filtration example is retained as a completed regression case.
This new example adds time-dependent solute transport and adsorption, while
keeping reaction and non-isothermal additions deferred.

## Engineering question

At fixed dry adsorbent inventory, feed concentration and flow, and a declared
outlet concentration limit, how does particle size affect usable treatment
capacity and pumping energy? The intended comparison holds chemistry and bed
volume/porosity fixed. Smaller particles increase the assumed intraparticle
uptake rate and increase measured bed resistance. No continuous or commercial
optimum is inferred from three tested sizes.

## Evidence and its limits

[The ledger](../literature/adsorption_evidence.json) separates measurements,
published fitted parameters, unresolved source conflicts and chosen scenarios.
The first chemical pair is phenol/activated carbon. An author-hosted 2017
undergraduate dissertation supplies the coconut-shell-carbon starting values;
related journal article [Meneghini et al. (2021)](https://doi.org/10.1590/S1413-415220180046)
is not used as independent validation here. Original PDF SHA256, repository URL
and visually checked tables are retained. Figure 7 contains error bars; their
numerical intervals and definition have not yet been digitized/audited. The
related journal's duplicate-test/95%-confidence statement is not converted into
invented pointwise standard deviations.

The Table 2 batch fit at 25 °C gives qmax=41.69 mg/g and b=0.094 L/mg;
their SI values are 0.04169 kg/kg and 94 m³/kg solute. The Table 4 effective
diffusivity at 50 mg/L is 6.5973e-8 cm²/s = 6.5973e-12 m²/s. These are
reported fitted parameters, not independent predictions. The author source's
particle diameter is 2.3 mm. Its column uses 20 g, 10 cm bed height, 2.04 cm
diameter, reported porosity 0.44 and 20 mL/min.

Those quantities imply bulk density 611.899 kg/m³ and particle-envelope dry
density 1092.677 kg/m³. The printed density of 574.33 g/L cannot be used
unambiguously as the particle density with all the other quantities. The
Table 4 affinity unit also conflicts with Table 2, and the related journal
tables retain different fits. The workflow does not silently reconcile them.
Figure 7 requires digitization and an initial/dead-volume audit before a
clean-bed column comparison. Its radial particle model also differs from the
initial lumped LDF approximation. Experimental validation is pending.

The particle-size study therefore uses a separate synthetic FCC bed at
porosity 0.6, assumed particle-envelope dry density 1000 kg/m³, length 1 m
and area 0.01 m²: 4 kg dry adsorbent in every design. The three diameters are
0.575, 1.15 and 2.3 mm. The published equilibrium and effective diffusivity
are explicitly transferred scenario inputs. They do not identify this bed
as the source's experimental specimen.

For the later particle-size validation, the strongest candidate found is
[Alhamed (2009), “Adsorption kinetics and performance of packed bed adsorber
for phenol removal using activated carbon from dates’ stones”](https://doi.org/10.1016/j.jhazmat.2009.05.002).
Its publisher abstract describes batch and column particle-size series.
The supplied full PDF has now been inspected in the separate comparison above;
this is a separate date-stone carbon. The original release's frozen ledger
retains its earlier abstract-only access status.
Iheanacho et al. (2021), DOI 10.1016/j.sajce.2021.02.003, is an alternative
corn-cob carbon series, also not pooled into the selected model.

## Model and scale boundary

OpenLB resolves the periodic external flow around impermeable FCC spheres.
It measures directional permeability; intraparticle nanoporosity remains in
the separate adsorption law. The native free-fluid scalar-TRT benchmark is a
component verification, not a resolved packed-bed mass-transfer closure.
All new source, configuration, build products and results remain outside the
original OpenLB tree. Stock indicators, dynamics, descriptors and functors
are reused with the official example lifecycle.

Julia uses the cell-centred conservative equation

```text
eps*dC/dt + d(u*C - eps*Dax*dC/dz)/dz = -rho_bulk*dq/dt
dq/dt = k_LDF * (qeq(C) - q)
qeq(C) = qmax*b*C/(1 + b*C)
rho_bulk = (1 - eps)*rho_particle
k_LDF = 15*Deff/R_particle²
```

The intraparticle LDF model is an approximation to spherical diffusion, not
the source's complete radial diffusion/film model. The inlet specifies total
solute flux u*Cin (Danckwerts); the outlet has zero dispersive flux. An
implicit coupled concentration/loading update conserves solute. Inlet mass,
effluent mass, liquid inventory and adsorbed inventory are stored separately.
The chosen breakthrough limit is Cout/Cin=0.05; it is a scenario quality
threshold, not a legal limit. Missing crossing times remain null.

The baseline neglects axial dispersion. Factors 0.1 and 1 in
`Dax=factor*(u/eps)*particle_diameter` test a chosen positive-dispersion range.
Those factors are not fitted correlations or uncertainty bounds. Film transfer
is not separately resolved. This limits material-performance interpretation
even after all numerical and software checks pass.

Actual DWSIM water properties set the hydraulic translation. The selected
flow preserves the reference native particle Reynolds number after the carrier
viscosity/density conversion. Each actual measured native permeability supplies
the corresponding Darcy pressure load to a real stock DWSIM pump. Julia
integrates that shaft power over time to breakthrough. Phenol remains an
external dilute solute; the DWSIM fluids contain Water only. This is not a
full adsorption/regeneration process flowsheet or a reactor model.

## Verification and reproduction

The standalone Julia tests include the non-adsorbing transport limit, linear
closed-batch analytical decay, nonlinear batch equilibrium from an independent
quadratic mass balance, retarded linear-equilibrium transport, conservation,
unavailable breakthrough and invalid inputs. Spatial and temporal refinement
of the actual column case are separate checks. Python tests exercise real
Julia and actual saved simulator records. The independent checker replays
native fields, source/input contracts, column calculations, CSVs, saved DWSIM
connectivity and work/mass/enthalpy balances without launching new external
solves. The original run performs the actual solves and reloads.

The first native bed scaling failed the density-departure screen and is saved
separately. Its measured velocity can initialize the smaller-time-step runs;
every new run still evolves and passes convergence, conservation, density and
refinement screens. Scientific tolerances are not relaxed. The scalar test's
first error calculation mistakenly used an unnormalised norm; that analysis
was corrected before acceptance, with the failed log retained.

A low-head water-pump pilot exposed insufficient relative enthalpy accuracy
in the original flash route at an 8 Pa rise. The separate new adapter uses
stock forced-liquid calculation, checks freezing/saturation margins and
reloads the saved phase setting. The pilot's relative enthalpy residual
falls from 0.00743 to 4.48e-7 with the same pump work and unchanged 1e-3
acceptance limit. Calling the saturation API initially needed an explicit
stock stream-context binding; that diagnostic attempt is retained. Original
simulator sources and earlier adapters are untouched. See the official
[material-stream API](https://dwsim.org/help/api/streams/material-stream.html).

[Exact launch and test commands](../adsorption/README.md). No new dependency
installation or manifest change is needed. The existing Panel interface is
preserved; no new adsorption UI or browser pass is claimed.

Accepted results follow below.

## Accepted synthetic results

Run: `/Users/dula/Documents/ChatGPT/Pore2Process/results/adsorption/process/runs/20261110T003326Z-8ff89948/summary.json`.

**8 native FCC hydraulic cases, 7 native scalar-TRT transport cases and 4 actual DWSIM solve/save/reload cases pass.** The latter count includes three design pumps and one property anchor. The independent checker replays 11 Julia column cases and checks all four saved pump topologies and 1,178 artifacts; it launches no external solvers.

New standalone Julia tests pass **46 assertions**; the preserved Julia package passes **902 assertions**. The full Python suite passes **86 tests**. The scientific figure was visually inspected, and the actual native VTP is checked for SI fields and zero-based IDs. No new browser interface was added.

The 48/32 grid permeability change is 2.8783%. Maximum lattice Mach is 0.003824; maximum density departure is 8.328316e-05. The separate scalar benchmark has measured convergence orders 2.0022 and 2.0005.

Every design has 4 kg dry adsorbent, 50 mg/L feed, a 2.5 mg/L chosen outlet limit and 155.014850 mL/min inlet-reference flow. Solid utilisation below is relative to equilibrium loading at the feed concentration, not the infinite-concentration Langmuir capacity.

| Particle diameter (mm) | Treated volume (m³) | Feed-equilibrium capacity used | Bed pressure loss (Pa) | Actual shaft power (mW) | Bed-only shaft energy (Wh/m³) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.575 | 2.748661 | 99.726% | 127.085497 | 0.437781 | 0.04706870 |
| 1.15 | 2.738479 | 99.353% | 31.773851 | 0.109454 | 0.01176809 |
| 2.3 | 2.697734 | 97.863% | 7.945930 | 0.027372 | 0.00294294 |

Reducing diameter from 2.3 to 0.575 mm increases treated volume by 1.8878% while bed-attributable pump power rises 15.9938 times. At this long contact time, even the largest particles use most of their feed-equilibrium capacity, so faster uptake adds relatively little treatment volume.

These absolute pump powers are small because the benchmark uses low flow through a clean, open synthetic bed. They include only the bed pressure requirement. Static head, connecting pipes/valves, adsorbent manufacture/replacement, regeneration utilities and full plant cost are excluded. A large relative increase in this bed duty is not a demonstrated large change in total plant energy. No preferred industrial particle size is established.

| Particle diameter (mm) | Baseline treated volume (m³), Dax=0 | Chosen dispersion factor 0.1 | Chosen dispersion factor 1 |
| --- | ---: | ---: | ---: |
| 0.575 | 2.748661 | 2.748558 | 2.747519 |
| 1.15 | 2.738479 | 2.738157 | 2.735324 |
| 2.3 | 2.697734 | 2.696838 | 2.689586 |

These are model-form scenarios, not measured bounds or confidence intervals. The scalar OpenLB benchmark does not supply their dispersion coefficients.

Doubling the column cells changes the largest-particle breakthrough time by 0.069936%; halving the time step changes it by 0.003826%. Maximum nominal normalised solute imbalance is 5.167518e-15; maximum actual-pump enthalpy residual relative to shaft work is 4.991439e-07.

The first connected attempt exposed an overly restrictive Newton feasibility rule in the nearly zero-concentration tail for the smallest particles. A bound-constrained line search now re-evaluates the original conservative residual at each trial iterate; completed states are not clipped or mass-corrected. A dedicated stiff, fast-uptake regression passes along with the analytical tests. The failed attempt remains saved.

[Scientific figure](adsorption.png) · [Compact real-simulator reference](adsorption_reference.json) · [Commands](../adsorption/README.md)

## Next executable milestone

Obtain and audit Alhamed (2009), DOI 10.1016/j.jhazmat.2009.05.002, for a distinct particle-size calibration/holdout study. Choose its own matched equilibrium, kinetic and column conditions; do not combine its date-stone carbon with the coconut-shell fit. Retain missing conditions as unknown or explicit sensitivities. No new experiment, unavailable raw dataset or specimen reconstruction is required.

In parallel with that source work, the next numerical extension is a stock-OpenLB packed-bed tracer/transport benchmark to qualify bed-scale dispersion and film-transfer information. This first release verifies the free-fluid transport component and bed hydraulics separately; it does not claim resolved native adsorption/mass-transfer closure. A later process extension can represent bed service/switching beyond the current actual-pump load. Reaction, regeneration heat and non-isothermal coupling remain deferred.
