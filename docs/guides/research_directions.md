# Pore2Process research questions, materials and showcase options

Prepared: 8 November 2026. Status: research shortlist and discussion brief;
these applications have not been implemented or selected as the next project.

This document preserves the chemical-engineering shortlist discussed in the
project and adds bioprocess applications. The priorities now include reacting,
non-isothermal systems, a connection to the user's PhD work on electrochemical
flow cells, and relevance to the bioprocess work publicly described by DTU's
PROSYS centre.

The verified starting point is the synthetic hydraulic workflow described in
[connected_hydraulics.md](../verification/hydraulics/darcy_pump.md): OpenLB → Python/Julia →
DWSIM. It establishes the software connection and measured-permeability
translation. Reaction, biological kinetics, thermal coupling and the new
application-specific models below remain future work.

## 1. What would make a useful showcase?

A strong showcase links a controllable material feature to a measurable process
result: the same treatment, conversion or production target with less energy,
smaller equipment, greater productivity or longer operating life. A useful
question includes the constraints: product purity, conversion, enzyme inventory,
allowable temperature, water recovery or electrode loading.

European relevance is supported by priorities around water efficiency,
industrial decarbonisation and circular resources in the
[European Water Resilience Strategy](https://commission.europa.eu/topics/environment/water-resilience-strategy_en)
and [Clean Industrial Deal](https://commission.europa.eu/topics/competitiveness/clean-industrial-deal_en).
The showcase rankings here are engineering judgments about technical fit and
communication potential. Commercial demand would require customer validation.

The material must have a defined geometry, properties and reference data.
At the level of our present single-phase, isothermal flow model, changing a
material label from alumina to SiC would not change the answer at identical
geometry and boundary conditions. Material-specific effects enter through
measured structure, thermal properties, wettability, surface chemistry,
conductivity, kinetics or biological compatibility, as appropriate.

## 2. Established chemical-engineering starting points

| Application | Candidate materials | Research question | Process metric |
| --- | --- | --- | --- |
| Liquid filtration | Sintered alumina, SiC, stainless-steel filters | How do constrictions, tortuosity and graded porosity change hydraulic resistance? | Pumping energy per m³ at specified flux and separation requirements |
| Packed-bed catalysis | Ni/alumina pellets, Cu/ZnO/alumina pellets, ceramic supports | How do pellet shape, bed structure and internal transport affect conversion and pressure loss? | Reactor productivity, compressor duty and temperature limits |
| Adsorption and ion exchange | Activated carbon, zeolites, polymeric resins | How do particle size and pore accessibility affect breakthrough and regeneration? | Treated volume per cycle, purity and regeneration consumption |
| Porous heat exchangers and thermal storage | Alumina particles, ceramic honeycombs, SiC or metallic foams | Does improved heat transfer justify additional flow resistance? | Useful recovered heat, charging/discharging time and blower power |

Filtration is closest to the validated hydraulic example. Adsorption adds
transport and transient behaviour; catalysis adds reaction and heat effects.
Thermal storage offers an energy-focused extension without immediately adding
complex chemistry. Packed-bed flow should first be compared with appropriate
pressure-drop references; reactive pellet models additionally need independent
kinetic and diffusion data.

## 3. European showcase opportunities from the original shortlist

### 3.1 Structured CO₂ adsorbents

**Materials:** amine-functionalised silica, zeolite-based monoliths, or another
well-characterised sorbent formed into honeycombs, lattices or foams.

**Question:** Can contactor geometry reduce pressure loss while maintaining CO₂
capture productivity?

**Metrics:** electricity and regeneration heat per tonne captured, CO₂
recovery/purity and contactor productivity. Compare structures at explicit
constraints on sorbent chemistry, inventory, bed volume and operating cycle.

This gives an attractive comparison between packed particles and engineered
structures. Relevant European research includes
[SINTEF's 3D-CAPS work](https://www.sintef.no/en/publications/publication/0198cc7ed0b9-54ddef79-009d-456b-b4b9-d9d297f7908d/)
and [TU Delft's monolith research with Skytree](https://www.tudelft.nl/me/over/afdelingen/process-energy/research/research-sections/engineering-thermodynamics/projects/experimental-study-and-modelling-of-monolith-structures-for-direct-air-capture).

An industrial CO₂-containing gas is a practical initial target before a full
DAC cycle. OpenLB-derived flow and transport information must be combined with
adsorption, heat and humidity data. Julia could represent the transient
adsorption cycle, with DWSIM calculating connected gas handling, heat duties and
product conditioning through a verified interface. A cyclic adsorber is not
already implemented in our workflow.

### 3.2 Water reuse and micropollutant treatment

**Materials:** ceramic filtration supports, granular activated carbon or
ion-exchange resin beads.

Two distinct questions are possible:

- Can a graded ceramic support reduce filtration energy with the same selective
  layer and filtration flux?
- Can an adsorbent architecture extend operation before breakthrough without
  excessive pumping, media replacement or regeneration consumption?

**Metrics:** kWh/m³, treated volume before breakthrough, media consumption,
water recovery and outlet contaminant concentration.

The revised EU wastewater framework strengthens the relevance of micropollutant
removal; see the [Council's explanation](https://www.consilium.europa.eu/en/press/press-releases/2024/11/05/urban-wastewater-council-adopts-new-rules-for-more-efficient-treatment/).
Concrete research examples include
[Wageningen's activated-carbon/PFAS work](https://www.wur.nl/en/activity/carbon-against-pfas-adsorption-low-concentration)
and [DECHEMA's PFAStrea project](https://dechema-dfi.de/kwi/en/Divisions/Chemical%2BTechnology/Current%2BProjects%2BChemical%2BTechnology/2027%2BPFAStrea%2BF962.html).

For a PFAS showcase, adsorption capacity, competing substances and regeneration
are central. Capture and destruction are separate objectives. A hydraulic
calculation alone cannot establish dissolved-contaminant removal, and an
ordinary microfiltration support is not a molecular PFAS separation model.

### 3.3 Flow-battery electrodes

**Materials:** carbon felt, graphite felt, carbon paper or engineered carbon-fibre
structures.

**Question:** What compression or fibre arrangement gives the best electrical
performance after pumping losses are included?

**Metrics:** net power, round-trip energy efficiency, limiting current and
pumping consumption at defined operating conditions.

Compression and structure affect electrolyte transport, electrical contact and
cell performance. Relevant European work includes
[KIT's study of compression and transport properties](https://publikationen.bibliothek.kit.edu/1000099205/45123221)
and [carbon-electrode research in BatCAT](https://www.nmbu.no/en/faculties/faculty-science-and-technology/optimizing-redox-flow-batteries-through-role-carbon).

OpenLB would supply pore-scale transport information; Julia would represent a
validated reduced electrochemical model; DWSIM would calculate circulation,
pumping and thermal management. Electrolyte properties must reflect composition,
temperature and state of charge. This is a distinctive showcase and a natural
connection to electrochemical flow-cell research.

### 3.4 PEM-electrolyser porous transport layers

**Materials:** porous titanium fibres, sintered titanium powder or structured
Ti foils.

**Question:** How should pore structure and wettability balance water supply,
oxygen removal and electrical contact?

**Metrics:** cell voltage at specified current density, hydrogen-specific
electricity consumption and operating limits.

Examples of European component research include
[Fraunhofer's controlled-porosity titanium work](https://www.iwm.fraunhofer.de/en/research-projects/CapS-PTL.html)
and [ZBT's structured-foil PTLs](https://zbt.de/en/projects/metalfoil/).

A credible efficiency result requires gas–liquid behaviour and an
electrochemical/contact-resistance model. Our present permeability calculation
would establish only the initial transport stage. Pumping energy alone would
not establish a meaningful improvement in electrolyser electrical efficiency.

### 3.5 Electrified structured catalytic reactors

**Materials:** electrically heated SiC structures carrying a catalyst, for
example a supported nickel reforming catalyst.

**Question:** Can structure improve heat delivery and reaction utilisation while
keeping pressure loss acceptable?

**Metrics:** electrical energy per unit product, conversion, reactor productivity
and maximum temperature gradients.

There is experimental precedent for
[electrically heated SiC structured reforming catalysts](https://www.sciencedirect.com/science/article/pii/S0920586120307860).
Industrial heat electrification also has a current European policy context,
including the Commission's
[industrial heat auction](https://climate.ec.europa.eu/news-other-reads/news/commission-publishes-terms-and-conditions-eur1-billion-industrial-heat-decarbonisation-auction-2026-09-24_en?prefLang=pl).

This is a strong route to non-isothermal reacting systems. Cantera could supply
validated reaction/transport models where suitable mechanisms exist, Julia
could translate local behaviour into reactor performance, and DWSIM could
calculate the surrounding plant. Electrified fossil-feed reforming still needs
a complete carbon balance; electrification by itself is not a zero-carbon claim.

### 3.6 Selective lithium recovery

**Materials:** lithium-selective manganese-oxide or titanium-based sorbents,
composite beads or porous electrochemical electrodes.

**Question:** How does transport architecture influence lithium recovery and
selectivity in a competing-ion brine?

**Metrics:** lithium recovery and purity, water/reagent consumption, energy per
kilogram recovered and sorbent-cycle productivity.

Fraunhofer's [M-BAT project](https://www.ise.fraunhofer.de/en/research-projects/m-bat.html),
starting in 2026, provides a European example linking geothermal-brine recovery
to downstream product production.

Competitive uptake, regeneration and brine thermodynamics would need substantial
experimental input. Geometry alone cannot establish selectivity. This is an
ambitious later showcase with a strong resource-security story.

### 3.7 Membrane-assisted solvent recovery

**Materials:** surface-modified ceramic nanofiltration layers and their porous
supports.

**Question:** When does membrane preconcentration reduce the distillation duty
needed to recover solvent at the required purity?

**Metrics:** heat/electricity per kilogram of recovered solvent, recovery,
purity and membrane area.

The Dutch [ReSolve project](https://ispt.eu/projects/resolve-re-use-of-industrial-solvent-process-streams/)
illustrates the direction. DWSIM would have a substantial role in the
distillation/recycle flowsheet. OpenLB could resolve suitable support-scale
transport, while selective-layer behaviour requires an appropriate model and
measured data. Nanometre-scale selectivity cannot simply be inferred from the
present continuum capillary-flow model.

## 4. Bioprocess connections to DTU PROSYS

PROSYS publicly describes research spanning biocatalysis, fermentation,
resource recovery and mathematical modelling across scales. Its biocatalysis
work includes enzyme kinetics, thermodynamics and stability. These are strong
connections for a reacting Pore2Process application.
[Source: PROSYS bio-based transformations](https://www.kt.dtu.dk/research/prosys/research-areas/bio-based-transformations).

Its research overview also identifies oxygen supply for bio-oxidations,
multistep biocatalysis and in-situ product removal.
[Source: PROSYS research areas](https://www.kt.dtu.dk/research/prosys/research-areas).

Two useful examples from the public record are:

- **FLUIZYME:** the PROSYS projects page describes research on enzyme stability
  under dynamic mixing and fluid-interface exposure, connecting small-scale
  experiments to larger-scale prediction. This motivates questions about
  activity loss under realistic flow histories.
  [Source: PROSYS projects](https://www.kt.dtu.dk/research/prosys/projects).
- **MORE-GAS:** DTU's June 2025 announcement describes a project beginning in
  January 2026 on model-based assessment of heterogeneity in gas fermentation,
  involving PROSYS/Biosustain researchers and pilot-scale data collection.
  [Source: DTU announcement](https://www.dtu.dk/english/news/all-news/dtu-project-receives-dkk-40-million-grant-from-the-novo-nordisk-foundation?id=af885a71-bf96-41ad-9263-61dc18234fd4).

The candidates below are proposed Pore2Process applications inspired by those
areas. A theme match does not mean a specific colleague is currently pursuing
our proposed case or that its data are available. Historical DTU studies are
identified as precedents. No contacts have been approached.

## 5. Bioprocess showcase candidates

### B1. Continuous immobilised-lipase esterification

**Materials/system:** lipase on a macroporous polymer or silica carrier.
A candidate chemistry is a flavour/fragrance esterification, with the actual
substrate pair selected for accessible kinetic, equilibrium and activity data.

**Question:** How do carrier pore accessibility and reactor temperature affect
ester productivity per unit enzyme, pressure loss and catalyst lifetime?

OpenLB could resolve flow and transport around or within a suitable resolved
carrier geometry. Julia would represent reaction, diffusion/effectiveness,
water-related inhibition or equilibrium effects where supported, and enzyme
deactivation. DWSIM could calculate feed conditioning, separation and recycle.

**Thermal role:** begin with imposed wall/inlet temperatures and a coupled
energy balance. Temperature-dependent kinetics and stability are meaningful
even if self-heating is small. Do not manufacture a hotspot by assigning an
unsupported heat source.

**Metrics:** conversion/selectivity at fixed throughput, product per mass of
enzyme over its operating life, separation duty and pumping energy.

**Fit:** a manageable first reacting biocatalysis case. A published
[continuous packed-bed immobilised-lipase study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7751923/)
provides a reactor-modelling precedent, although its specific reaction and
isothermal experiment do not validate our proposed esterification or heat model.

### B2. Porous enzyme reactor for glucose oxidation

**Materials/system:** glucose oxidase with catalase immobilised on a porous
silica/polymer carrier or an appropriate membrane support. Buffer conditions
and enzyme loading must be explicit.

**Question:** Can carrier geometry improve oxygen access and glucose conversion
while limiting pressure loss, peroxide exposure and loss of enzyme activity?

The chemistry offers coupled substrate/oxygen transport and a reaction cascade.
Julia would model the measured enzyme kinetics, catalase coupling, activity
loss and energy balance. OpenLB would supply validated local transport and,
in a later stage, thermal information. DWSIM would calculate the surrounding
liquid/gas handling and heat duties.

**Thermal role:** evaluate reaction heat and its temperature feedback alongside
oxygen limitation. Calculate the possible adiabatic temperature rise before
claiming strong thermal gradients. Include buffer/speciation contributions in
measured or calculated heat where relevant.

**Metrics:** space-time yield, enzyme utilisation, oxygen consumption/supply
energy, maximum temperature and cooling duty at a specified conversion.

**Fit:** particularly good for PROSYS's bio-oxidation and enzyme-stability themes.
DTU researchers have experimentally studied glucose oxidation in reactor
cascades, including the role of catalase and oxygen availability.
[Source: Lindeque and Woodley study](https://backend.orbit.dtu.dk/ws/portalfiles/portal/275410319/Lindeque_and_Woodley_2022_.pdf).
Our immobilised, non-isothermal configuration would be a new model requiring
its own validation.

### B3. Biofilm trickle-bed biomethanation

**Materials/system:** structured polymeric or ceramic packing supporting a
methanogenic biofilm. A thermal comparison could consider a more conductive
carrier, subject to actual wetting and biological-compatibility data.

**Question:** How does packing geometry balance gas–liquid transfer, biofilm
transport, pressure loss and heat removal during H₂/CO₂ conversion to methane?

The idealised net conversion is `CO₂ + 4 H₂ → CH₄ + 2 H₂O`. A biological model
also needs biomass/yield terms, uptake rates, temperature dependence and
appropriate pH/inhibition assumptions. Equilibrium conversion is not a
substitute for those kinetics.

**Thermal role:** this is a strong candidate for reaction-heat and cooling
feedback. Water phase, evaporation/condensation and enthalpy reference states
must be consistent. Start with prescribed wetting/biofilm thickness; add
multiphase dynamics or growth after separate benchmarks.

**Metrics:** methane productivity and quality, H₂ conversion, recycle/compression
energy, cooling/recoverable heat and the allowable operating-temperature range.

Julia would translate local mass/heat transfer and biological kinetics into a
reactor model. DWSIM could solve gas recycle, compression, water separation and
heat recovery. An appropriate reactor interface must be verified.

**Fit:** strongest overall candidate for a genuinely non-isothermal biological
flagship. A [2019 PROSYS PhD thesis](https://orbit.dtu.dk/en/publications/biomethanation-of-synthesis-gas-in-trickle-bed-reactors/)
provides direct precedent on trickle-bed design, wetting, temperature and scale-up.
Its syngas results are a research reference, not automatic validation for a
simplified H₂/CO₂ feed or a different microbial community.

### B4. Porous gas distributors for aerobic gas fermentation

**Materials/system:** sintered stainless-steel or ceramic gas distributors in
a fermentation producing biomass/protein or another specified product. Select
one organism, gas substrate and kinetic dataset before implementation.

**Question:** Does distributor pore architecture improve gas utilisation and
productivity enough to compensate for compression, circulation and cooling?

This places a real porous component inside a recognisable bioprocess. OpenLB
could study local bubble formation or transfer near a prescribed interface.
Julia would need a validated bridge from local information to bubble populations,
gas transfer and biological uptake at reactor scale. Distributor permeability
alone does not determine whole-reactor `k_L a`.

**Thermal role:** include metabolic heat, oxygen limitation and cooling demand
where supported by the selected culture's data. Temperature and dissolved-gas
heterogeneity can affect uptake and productivity.

**Metrics:** product per reactor volume/time, gas utilisation, electricity per
kilogram product and cooling duty at specified product quality.

**Fit:** close to the scale-heterogeneity questions in
[MORE-GAS](https://www.dtu.dk/english/news/all-news/dtu-project-receives-dkk-40-million-grant-from-the-novo-nordisk-foundation?id=af885a71-bf96-41ad-9263-61dc18234fd4).
The porous-distributor showcase is our proposal, not a claim about that project's
chosen hardware. The multiphase-to-reactor scale transition makes this a
substantial extension.

### B5. Microbial electrosynthesis in a porous flow-through cathode

**Materials/system:** carbon felt, carbon cloth or a structured carbon cathode
with an acetogenic biofilm; CO₂-to-acetate is a candidate reaction system.

**Question:** How do electrode architecture and biofilm thickness affect product
formation and net electricity consumption once pumping and product recovery
are included?

This is the clearest bridge between electrochemical flow-cell research and
bioprocess engineering. OpenLB could resolve pore flow and species transport.
Julia would couple electrochemical response, biofilm uptake, diffusion and
appropriate acid–base behaviour. The model must distinguish direct electron
uptake from H₂-mediated conversion according to the selected experimental
system. DWSIM could calculate circulation, gas handling and product recovery.

**Thermal role:** quantify electrical, reaction and biological heat consistently.
Avoid counting all electrical input as heat and also adding the same energy
again as a reaction contribution. Thermal gradients may be secondary at the
chosen scale; they must be demonstrated rather than assumed.

**Metrics:** acetate productivity/titre, coulombic efficiency, total electricity
per kilogram product, pumping losses and recovery duty.

**Fit:** high personal/research fit, with more model coupling than an enzyme bed.
DTU research on modified carbon-felt biocathodes provides a precedent for
linking electrode structure with microbial conversion.
[Source: original DTU-hosted study](https://backend.orbit.dtu.dk/ws/files/132273509/c7se00073a.pdf).
This is a wider DTU precedent, not a claim that it is a current PROSYS project.

### B6. Packed-bed solid-state fermentation

**Materials/system:** a porous bed of agricultural residues, such as bran and
cellulosic material, supporting a selected industrial fungus for enzyme or
other bioproduct production. Initial geometry may be synthetic but must be
labelled as such.

**Question:** How do bed architecture and aeration remove metabolic heat while
maintaining oxygen access and avoiding excessive drying?

OpenLB could provide local airflow, heat-transfer and transport information for
representative structures. Julia would model growth, substrate consumption,
metabolic heat and moisture at bed scale. DWSIM could represent humid-air
conditioning, blowers, heating/cooling and exhaust/recycle handling.

**Thermal role:** central. Solids/gas temperature differences, evaporation and
moisture retention can matter, so a single equilibrium temperature or a dry
fixed-geometry model must not be assumed sufficient.

**Metrics:** temperature uniformity, oxygen-limited fraction, product yield,
water loss, batch time and energy per kilogram product.

**Fit:** excellent for demonstrating biological heat generation in a porous
medium. It is a general bioprocess option; no specific current PROSYS project
is asserted here. A [pilot-scale experimental transfer study](https://www.sciencedirect.com/science/article/pii/S1385894720333714)
and [Wageningen's packed-bed modelling thesis](https://research.wur.nl/en/publications/macroscopic-modelling-of-solid-state-fermentation/)
provide starting references. Changing geometry, growth and moisture create
substantial validation work.

### B7. Enzyme-retaining membrane reactor with product removal

**Materials/system:** ceramic or functionalised polymer membranes retaining or
immobilising enzymes, coupled to a selected enzyme reaction and separation.

**Question:** Can controlled pore architecture and product removal increase
useful catalyst productivity without excessive membrane resistance or activity
loss?

**Metrics:** product per enzyme inventory/time, enzyme retention/lifetime,
product purity, membrane area and energy per kilogram recovered product.

OpenLB could describe suitable support-scale flow/transport; Julia would couple
kinetics, retention and product inhibition where supported; DWSIM would solve
the connected separation/recycle process. Thermal coupling would address
activity, stability and external heat duties, with self-heating assessed from
data. Molecular retention or membrane selectivity requires separate closure.

**Fit:** direct thematic relevance to PROSYS's
[separation technology research](https://www.kt.dtu.dk/research/prosys/research-areas/separation-technology).
The historical DTU [CerEnzyme project](https://ssc-energy.dtu.dk/projects/cerenzyme)
also describes inorganic supports designed around enzyme immobilisation.
These establish topic relevance, not availability of a ready-made thermal model.

Protein chromatography on porous resin beads is another useful downstream
showcase: bead transport/binding → breakthrough → purity/recovery and cycle
productivity. It would normally be a weaker choice for demonstrating substantial
reaction heat, so it is not a leading non-isothermal candidate.

## 6. Ranking with the new interests included

These are qualitative planning judgments, not measured performance rankings.

| Goal | Leading option | Reason |
| --- | --- | --- |
| Fastest next hydraulic material case | Graded alumina/SiC filtration support | Reuses the existing measured-resistance workflow |
| Manageable first bioreaction | Immobilised-lipase bed | Adds kinetics, temperature dependence and catalyst lifetime with a bounded initial scope |
| Strong PROSYS biocatalysis connection | Immobilised oxidase/catalase reactor | Couples pore transport, oxygen supply and enzyme stability |
| Strongest reacting/non-isothermal bio flagship | Biofilm trickle-bed biomethanation | A porous carrier, conversion, gas transfer and a meaningful heat-removal problem |
| Most direct bridge to the PhD | Microbial electrosynthesis in a porous cathode | Connects electrochemical transport to biology and product recovery |
| Strong biological heat/moisture demonstration | Solid-state fermentation | Material structure affects airflow, drying, oxygen access and metabolic-heat removal |
| Close gas-fermentation scale-up theme | Porous distributor plus gas-fermentation model | Connects component design to reactor heterogeneity and productivity |
| Broad separation/energy flagship | Structured CO₂ adsorption contactor | Clear material-to-cycle-to-energy connection |
| Most natural substantial Cantera role | Electrified structured catalytic reactor | Validated chemical kinetics and thermal transport can be central |

The earlier recommendation for a quick hydraulic continuation remains:

**Selected on 8 November 2026:** the filtration question below is now the
active research direction. See [the bounded first study](../verification/filtration/layered_support.md)
and [the article/data workspace](../../literature/README.md). Biological and
bioelectrochemical options remain future candidates.

> At fixed filtration flux and unchanged selective layer, how much can support
> geometry reduce total pumping energy, and when does the selective layer
> become the dominant resistance?

With the additional interests, three biological questions deserve discussion:

1. **Biocatalysis:** At fixed enzyme inventory and conversion, can carrier
   architecture reduce oxygen-supply and pumping demand while preserving
   activity under the measured temperature history?
2. **Thermal biofilm reactor:** At fixed methane productivity/quality, which
   packing structure limits temperature excursions with the lowest combined
   recycle, compression and cooling demand?
3. **Bioelectrochemistry:** At fixed acetate production and recovery targets,
   which porous cathode/biofilm structure minimises total electricity use?

A sensible portfolio would keep one tractable enzyme case as the reacting
learning example and one longer-term flagship: biomethanation for thermal
bioprocess engineering or microbial electrosynthesis for the electrochemical
connection. A final selection should follow the availability of kinetics,
thermal measurements and a suitable experimental reference.

## 7. Route to reacting, non-isothermal systems

The thermal model should answer a physical question, not simply colour a
geometry with temperature. A staged implementation would be:

1. **Nonreacting transport:** benchmark flow, passive-species transport and
   heat transfer separately, using the appropriate existing OpenLB examples.
2. **Known heat source:** validate the energy balance and boundary conditions
   with a prescribed, independently checkable source before introducing biology.
3. **Isothermal reaction:** validate a well-mixed or one-dimensional kinetic
   reference, then test resolved transport against its limiting cases.
4. **Temperature-dependent reaction:** couple measured kinetics and heat of
   reaction; check isothermal and adiabatic limits and overall conservation.
5. **Biological complexity:** add enzyme deactivation, biofilm diffusion,
   growth, inhibition, pH, moisture or phase changes one at a time as required.
6. **Process effect:** pass the validated response/model into the same traceable
   Python/Julia/DWSIM boundary and compare actual process outputs.

Before expecting significant self-heating, estimate the adiabatic temperature
rise using the selected stoichiometry, conversion, enthalpy and heat capacity.
Dilute enzyme systems may show small temperature rises even with appreciable
reaction enthalpy. Temperature-dependent rates and deactivation can still make
thermal operation important. Stronger thermal demonstrations may be found in
biomethanation, aerobic metabolism or solid-state fermentation, depending on
actual throughput and heat removal.

For biofilm/gas systems, begin with a defined geometry and prescribed interfaces
where physically justified. A prescribed film can validate local transport; it
does not validate film formation, wetting dynamics or whole-reactor transfer.
For gels, membranes and microbial structures below the resolved scale, use
independently justified effective transport rather than claiming direct resolution.

## 8. Roles, contracts and validation

The required architecture is unchanged:

- **OpenLB:** the detailed pore/component-scale flow, species and thermal stage,
  using existing suitable components and official example structure. The
  [OpenLB code policy](../architecture/openlb_policy.md) applies. Composition of existing
  facilities for a new application still needs a fresh benchmark.
- **Julia:** custom reduced kinetics/biology, scale translation, analysis and
  reactor-level models, independently usable without the UI.
- **Python:** coordination, interfaces, data, visualisation and reproducible
  input/output snapshots.
- **DWSIM:** actual process-stream, equipment, recycle and utility calculations.
  Reuse a suitable supported unit/interface where available; investigate the
  target-platform coupling before choosing a custom implementation.
- **Cantera, optional:** validated chemical kinetics, thermochemistry and
  transport where suitable mechanisms/data exist. Enzyme regulation, microbial
  growth, biofilm behaviour and deactivation still require validated biological
  constitutive models; adding Cantera does not supply them automatically.

An explicit pressure-loss relation is useful for beds, carriers and membrane
supports. Reacting showcases additionally require a validated reaction/thermal
response. DWSIM must receive outputs that actually change its material or
energy calculation; a biological reactor picture or an equilibrium-only
substitution would not establish the intended coupling.

Record, as applicable:

| Quantity | Contract issue |
| --- | --- |
| Permeability and pressure loss | m² and Pa; direction, representative area, validity range, Darcy/Forchheimer basis |
| Species transport | m²/s diffusivity or explicitly defined transfer coefficient; species order and concentration basis |
| Reaction rates | mol/(m³ s), mol/(m² s) or per enzyme/biomass mass, with the relevant volume/area/mass stated |
| Enzyme/biomass inventory | Dry/wet basis, loading, activity convention and changing inventory |
| Heat transport | W/(m K) conductivity, W/(m² K) interfacial coefficient, actual exchange area and phase temperatures |
| Heat generation | W/m³ on a declared volume basis; reaction sign, enthalpy reference and water phase |
| Electrochemical response | Potential/reference electrode, current sign, area basis, electron balance and heat accounting |
| Biological chemistry | pH/speciation, buffer assumptions, maintenance, growth/yield, inhibition and deactivation |
| Process outputs | Product flow/purity, conversion, utilities, recovery, recycle and equipment conditions |

Validation should include mass/element/energy conservation, and charge/electron
balance where relevant; spatial/time-step convergence; independent kinetic and
thermal limiting cases; and comparison with an experiment at defined conditions.
Thermal, species and reaction coupling have not yet been verified in this project.

## 9. Information to seek before choosing a bio case

A short discussion with relevant colleagues could establish:

- An enzyme/organism and reaction for which measured kinetics and temperature
  dependence can be used and shared.
- A defined carrier, electrode, membrane, packing or distributor geometry, with
  porosity, dimensions and material properties.
- Measurements of heat release, cooling demand or temperature profiles; whether
  the experiment was thermostated and what uncertainty was reported.
- A reference for diffusion, uptake, deactivation, pH and gas transfer as needed.
- A feasible process target: conversion/productivity/purity, oxygen or H₂
  utilisation, enzyme lifetime, energy consumption or recovery duty.
- Data availability, confidentiality and whether a synthetic demonstrator can
  be checked against a publishable reference.

The immediate executable planning step is to choose one enzyme/organism and
reference dataset, estimate the relevant transport and heat scales, and write
one bounded hypothesis. This document expands the shortlist; it does not
initiate a new simulation or change the validated hydraulic implementation.
