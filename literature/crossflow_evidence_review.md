# Crossflow performance and circulation energy: evidence search

**Receipt update, 9 November 2026:** the four requested articles have now been
received as full Markdown, metadata and figures. See [the received-paper
review](crossflow_papers_review.md) and [numerical extraction](crossflow_extraction.json).
The acquisition requests and access descriptions below document the initial
search and are superseded for García 2008, Chiu 2005, Chiu 2006 and Frederic
2018. No further article is required before the bounded reference benchmark.

Screened 9 November 2026. This is a literature selection and model plan;
no new crossflow correlation has been fitted and no simulation has changed.
The completed [hydraulic recycle case](../docs/crossflow_recycle.md) remains
the reference. Search provenance, access levels, source locations and explicit
unknowns are recorded in [the screening ledger](crossflow_screening.json).

## Decision

Use a published **critical-flux versus crossflow-velocity experiment** as the
next bounded benchmark. At a prescribed filtration flux, it can identify
circulation rates that meet the paper's fouling criterion. Combine that
constraint with the module's hydraulic cost before attempting an optimum.
The present clean-water model assigns circulation a cost but no filtration
benefit, so it cannot answer this question by itself.

The strongest direct lead is **García and Chiu (2008)**: ceramic channels,
critical flux, axial pressure loss and specific energy in one study. Obtain
its full article first, with **Chiu and James (2005)** as the companion
experimental-method lead. **Frederic et al. (2018)** is the most relevant
separate model lead for support/selective-layer contrast and spatial fouling.
None has yet supplied accepted numerical calibration data for this extension.

This extends the application from clean water to a declared suspension feed.
It does not turn the existing Zhu-calibrated support scenario into a validated
prediction for the new suspension or membrane.

## What the supplied papers can supply

| Supplied source | Relevant evidence found | Role and limit for the next step |
| --- | --- | --- |
| Mo et al. (2022), main article and Word supplement | Fixed-layer structure/resistance model. The introduction and bibliography point to Frederic et al. (2018). | Retain the layer framework and follow that reference. No suitable velocity-performance series identified in the supplied material. |
| Biesheuvel and Verweij (1999) | Support permeability, processing, porosity and mechanical strength. | Retain the existing support comparison and its failed reduced-model screen. No crossflow-performance calibration identified. |
| Zhu et al. (2021), §2.3 and Fig. 6 | Crossflow test mode; DI-water permeance at 1 bar; matched M1 support/whole-membrane series. | Keep the accepted bounded impedance extraction. Testing in crossflow does not establish a measured velocity-benefit curve. |
| Jan et al. (2025), §§2.1, 2.3 and references [28], [30] | Single-channel dimensions, constant-flux crossflow setup, coating/permeance/rejection comparison. References lead to Chen (2022) and Jan (2024). | Useful module/setup and ceramic-material leads. No suitable velocity sweep identified. Changing the SiC coating changes the selective layer; this is not our fixed-layer design comparison. |

Jan's reported tube is 7 mm inside diameter, 10 mm outside diameter and
100 mm long. Those dimensions belong to Jan's specimen, not to Zhu's M1
membrane or a future star-channel module. They must not be combined into a
supposed experimentally matched membrane. The source hashes and exact local
locations are in the ledger; prior extractions remain unchanged.

## New sources and their distinct jobs

### 1. Main circulation-benefit benchmark: García and Chiu (2008)

F. García and T. Chiu, **“Economic aspects of critical flux operability in star
shaped microfiltration membranes: Influence of some operating conditions,”**
*Journal of Membrane Science* 325(2), 641–646.
[DOI 10.1016/j.memsci.2008.08.033](https://doi.org/10.1016/j.memsci.2008.08.033).

The publisher abstract links increasing critical flux to increasing axial
pressure loss; its reported energy minimum is at 0.9 m/s for that study.
Indexed publisher sections identify TiO₂ suspensions and **Table 1** containing
Reynolds number, pressure loss and specific energy. This directly addresses
the missing performance/circulation tradeoff.

**Request:** full PDF and Markdown, original figures, Table 1, equations and
any available supplement. Extract the velocity–critical-flux series at fixed
feed conditions, channel cross-section and wetted perimeter, active length,
membrane area, pore size/material, temperature/viscosity if reported, flux-step
criterion and duration, and the energy equation/system boundary. The figure
number for the critical-flux series has not been verified.

**Limit:** publisher abstract and indexed sections screened; full tables,
figures and equations not inspected. The reported minimum is a source-specific
lead, not a recommended operating speed or an optimum for Pore2Process.
Do not copy its energy values into our shaft-energy comparison until the
normalization, included duties and pump efficiency convention are checked.

### 2. Companion experimental protocol: Chiu and James (2005)

T. Y. Chiu and A. E. James, **“Critical flux determination of non-circular
multi-channel ceramic membranes using TiO₂ suspensions,”**
*Journal of Membrane Science* 254, 295–301.
[DOI 10.1016/j.memsci.2005.01.037](https://doi.org/10.1016/j.memsci.2005.01.037).

This is a cited companion lead for the same ceramic-channel/particle family.
Its DOI metadata and bibliographic trail were checked; the full experimental
article was not inspected. Retrieve it with the 2008 article to check the
critical-flux protocol, channel definition and velocity/chemistry conditions.
Do not assume the two papers use identical modules, specimens or conditions,
or count overlapping measurements as independent validation.

### 3. Geometry-to-fouling model: Frederic et al. (2018)

E. Frederic, C. Guigui, M. Jacob, C. Machinal, A. Krifi, A. Line and P. Schmitz,
**“Modelling of fluid flow distribution in multichannel ceramic membrane:
Application to the filtration of produced water,”**
*Journal of Membrane Science* 567, 290–302.
[DOI 10.1016/j.memsci.2018.09.021](https://doi.org/10.1016/j.memsci.2018.09.021).

Found through Mo's bibliography. The publisher describes a theoretical
Darcy/resistance model with **2D membrane cross-sections**, support/selective
permeability contrast, fouling and circular/square channels. It could connect
our original support question to spatial flux redistribution.

**Request:** full article, equations, geometry/boundary-condition figures and
available supplements. Check where the fouling closure comes from and what
experimental or numerical comparisons support it. This is a model reference,
not independent experimental validation, and its cross-sectional formulation
must not be presented as an already verified axial crossflow-pressure model.

### 4. Independent alternative: Tanudjaja et al. (2017)

H. J. Tanudjaja, V. V. Tarabara, A. G. Fane and J. W. Chew,
**“Effect of cross-flow velocity, oil concentration and salinity on the critical
flux of an oil-in-water emulsion in microfiltration,”**
*Journal of Membrane Science* 530, 11–19.
[DOI 10.1016/j.memsci.2017.02.011](https://doi.org/10.1016/j.memsci.2017.02.011).

Publisher sections identify **Fig. 4: critical flux versus crossflow velocity**.
The study uses direct observation through the membrane and varies oil loading,
velocity and salt. Its reported domain includes 0.1–0.4 m/s and 50–500 ppm oil.
Request the full article if the TiO₂ benchmark lacks a usable series or if we
select oily-water treatment. The membrane, droplet physics and critical-flux
criterion require their own benchmark; its coefficients cannot be transferred
to a ceramic/TiO₂ case as validation.

### 5. Time-dependent ceramic fouling: Chen et al. (2022)

M. Chen, S. G. J. Heijman, M. W. J. Luiten-Olieman and L. C. Rietveld,
**“Oil-in-water emulsion separation: Fouling of alumina membranes with and
without a silicon carbide deposition in constant flux filtration mode,”**
*Water Research* 216, 118267.
[DOI 10.1016/j.watres.2022.118267](https://doi.org/10.1016/j.watres.2022.118267),
[public university PDF](https://ris.utwente.nl/ws/portalfiles/portal/297029948/1_s2.0_S0043135422002305_main.pdf).

This is Jan's reference [28]. Public PDF text was inspected. **§2.5.1 fixes
crossflow at 0.44 m/s**, so it cannot identify velocity dependence. **Fig. 4**
compares normalized permeance over filtration cycles for D0/D90 at
90/100/110 LMH. The protocol includes backwashing and flushing, with room
temperature reported as **22 ± 3 °C**; that is a reported operating-temperature
range, not an invented flux error bar.

Useful later for cycle/fouling/cleaning duty or for a separate ceramic oily-water
benchmark. Retain the full cycle protocol when reading the curves. PDF text
and captions were screened, but curve pixels were not visually verified or
digitized; direct local download returned HTTP 403. No additional numeric
fouling dataset has been accepted. The abstract's low-fouling observations
must not be converted into a universal critical-flux threshold.

### 6. Biological-feed alternatives

**T. Y. Chiu, M. V. Lara Dominguez and A. E. James (2006), “Critical flux and
rejection behaviour of non-circular-channelled membranes: Influence of some
operating conditions,”** *Separation and Purification Technology* 50(2), 212–219.
[DOI 10.1016/j.seppur.2005.11.028](https://doi.org/10.1016/j.seppur.2005.11.028).
Publisher sections describe synthetic activated-sludge-process effluent,
star-channel ceramics, hydrodynamics, critical flux and solids retention.
This is a concrete biological-water-treatment option. Full text is needed
before selecting a curve or carrying its retention law into a process case.

**R. W. Field, D. Wu, J. A. Howell and B. B. Gupta (1995), “Critical flux concept
for microfiltration fouling,”** *Journal of Membrane Science* 100(3), 259–272.
[DOI 10.1016/0376-7388(94)00265-Z](https://doi.org/10.1016/0376-7388(94)00265-Z).
The publisher abstract reports constant-flux experiments with yeast, yeast
debris and dodecane emulsions in flat-sheet/tubular systems. It supplies both
the foundational method and a potential yeast-processing benchmark. Read
the full paper before treating any cited yeast data as newly measured there.

These are optional separate feed cases, not extra prerequisites for the
current TiO₂/ceramic lead. No bioelectrochemical application is implied.

### 7. Definitions and limitations

P. Bacchin, P. Aimar and R. W. Field (2006), **“Critical and sustainable fluxes:
Theory, experiments and applications,”** *Journal of Membrane Science* 281, 42–69.
[DOI 10.1016/j.memsci.2006.04.014](https://doi.org/10.1016/j.memsci.2006.04.014),
[author-hosted review](https://www.patricebacchin.fr/images/doc/paper/bacchin_jms_2006.pdf).
The PDF was downloaded to the ignored literature inbox and its first page
visually checked; the review has not been read in full.
Use this review for terminology and protocol limitations, not as an independent
experimental calibration set. A short-step critical flux is not proof of
indefinite fouling-free operation; sustainable flux also depends on feed,
time horizon and the tolerated cost of fouling/cleaning.

## Next executable model and acceptance plan

1. Read the main/companion articles and select **one matched feed, membrane and
   operating-condition series**. Transcribe tables first; digitize plotted
   points from original PDF figures when needed. Preserve paper/specimen,
   figure/table, units, conditions, method and reported uncertainty. Unknown
   conditions remain unknown; reading sensitivity is not experimental error.
2. Define the physical module. Inlet superficial channel velocity is
   `u_in = Q_module,in / A_open`; hydraulic diameter is `4*A_open/P_wetted`.
   Membrane filtration area is a different quantity. Star-channel inner/outer
   dimensions are not a circular bore diameter. State inlet versus mean
   velocity and account for permeate removal where significant.
3. Implement a bounded Julia interpolation or justified low-parameter relation
   `J_crit(u, feed conditions)` within the measured domain. The initial
   constraint is `J_target <= J_crit` under the source protocol. Any operating
   margin is an explicit scenario with sensitivity, not an experimental datum.
   Use a held-out subset if the number of independent points permits it;
   otherwise report a reproduction check without calling it validation.
4. Bound `delta_p_channel(u)` using that module's reported measurements or
   a geometry-appropriate correlation. Check Reynolds number and entrance/
   turbulent effects. Do not transplant the present assumed linear channel
   resistance over a turbulent range without evidence. If OpenLB is extended,
   use stock facilities in separate project folders and recheck convergence.
5. Pass the admissible operating points and hydraulic load through the saved
   Python/Julia-to-DWSIM boundary. Retain fixed target flux and unchanged
   selective layer within each support-design comparison. The process result
   is shaft energy per volume of permeate **subject to the declared filtration
   criterion**. Reproduce the paper's energy boundary separately from our
   two-pump boundary; distinguish shaft, electrical and cooling duties.
6. Introduce feed concentration and rejection balances before claiming a
   suspension-recycle prediction: concentration generally changes with recovery.
   Water-only DWSIM properties are, at most, a declared dilute-carrier
   approximation with sensitivity. Do not silently label the old pure-water
   recycle case a validated slurry or emulsion process.

In this first extension, support geometry may alter clean hydraulic resistance
while the same empirical fouling criterion is held fixed as a scenario. A claim
that support geometry itself changes fouling requires additional justified
local-flux modelling or evidence. A time-dependent fouling law, cleaning
optimization, reaction and non-isothermal feed physics remain later steps.

The immediate acquisition request is the **2008 main article and 2005
companion**, followed by **Frederic 2018** for the geometry link. Add the 2006
biological-feed paper if exploring that alternative in parallel. Supplements
are useful when available; their absence, unavailable raw data and new lab
measurements are not acceptance gates.
