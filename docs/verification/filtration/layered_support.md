# Fixed-flux filtration: support geometry and the energy limit

> **Dated verification record.** This retains the execution/audit history.
> Use the [maintained scientific report](../../../examples/filtration/layered_support/REPORT.md) for the current explanation
> and the example README for current commands.

Selected on 8 November 2026. Research question:

> At a fixed filtration flux and unchanged selective layer, how much can support
> geometry reduce total pumping energy, and when does the selective layer become
> the dominant resistance?

The first material family to investigate is asymmetric alumina/SiC ceramic
filtration membranes. Alumina is the initial literature candidate. Zhu's M1
series now supplies a bounded effective-impedance reference; the next support
geometry is synthetic and literature-informed, rather than an exact specimen
reconstruction. No new experiments or inaccessible raw data are required;
see [verification scope](../../architecture/verification_scope.md). Bioelectrochemistry remains an option in
[the broader portfolio](../../guides/research_directions.md).

The [controlled fixed-opening comparison](graded_support.md) is now
complete. The 4 → 6 µm interior widening reduces raw support resistance by
35.88% and actual DWSIM transmembrane shaft energy by 11.30% in the declared
literature-anchored scenario. Both terminal necks remain fixed in this family;
it does not cover every possible asymmetric graded support.

## Model and the answer it can support

Begin with fully wetted, clean, Newtonian water at constant temperature and
normal flux J through layers of the same membrane area. The selective layer is
a fixed measured or explicitly assumed hydraulic resistance, not another copy
of the OpenLB support geometry. Define

```text
Rs = Ls / ks                         [m^-1]
Rtotal = Rs + Rselective + Rfixed     [m^-1]
delta_p = mu * J * Rtotal            [Pa]
Q = J * Amembrane                    [m^3/s]
Pshaft = Q * delta_p / eta           [W]
eshaft = delta_p / (eta * 3.6e6)      [kWh/m^3]
```

Rfixed can represent explicitly supplied intermediate/interface resistance.
In the current synthetic sweep it is zero. It must not hide unknown measured
contributions when fitting an actual membrane. There is no osmotic term in
this pure-water model.

Let f=ks,new/ks,baseline at unchanged support thickness, and let
alpha=Rselective/Rs,baseline. With Rfixed=0, independent algebra gives

```text
fractional shaft-energy saving = (1 - 1/f) / (1 + alpha)
maximum saving if support resistance vanished = 1 / (1 + alpha)
selective resistance fraction = Rselective / Rtotal
```

This is the controlling limit: improving the support can only remove its own
share of the pressure load. A tenfold improvement in support permeability does
not imply a tenfold reduction in total pumping duty. At equal efficiency and
fixed flux/area, fractional savings in pressure, shaft power and shaft energy
per filtrate volume are identical.

We define **dominant** as greater than half of the total resistance:
Rselective > Rs+Rfixed. A **strongly dominant** selective layer has at least 90%:
Rselective >= 9(Rs+Rfixed). Equality at 50% is a crossover, not dominance.
These are declared reporting thresholds, not a material-specific law.

## Executable first study

The first application calculation extends the completed hydraulic benchmark
with **fresh native OpenLB simulations**, reusing its existing case header and
stock numerical facilities. No original OpenLB files are edited.

Two synthetic straight-pore support cells have radii 4 and 5 µm, identical
20 µm square representative area and 40 µm axial simulation length. They are
not reconstructed alumina/SiC samples. Enlarging the pore also changes porosity
and the support's surface opening. Thus this is not an equal-porosity or
fixed-interface grading experiment. Holding selective resistance constant is
an explicit model constraint; unchanged selectivity has not been demonstrated.

OpenLB solves both at N24 and N40, with 10 Pa imposed pressure drop. It measures
flow and interior tap pressure drop. Julia infers axial permeability using the
representative cell area. Lower-pressure (1 Pa), doubled-pressure (20 Pa) and
doubled-domain cases at N24 check flow scaling and the stability of that
property. The lower and upper pressure gradients bracket the support's
operating gradient in the upscaled example. Numerical acceptance
requires analytical permeability/profile errors <1%, mass-flow spread <0.1%,
N24→N40 k change <0.5%, and pressure/domain k changes <0.1%. Pressure taps
exclude the ends: the inferred bulk k does not include reservoir/entrance losses.

The scale translation treats the cells as aligned, continuing pores through
the thickness; it makes no isotropic-permeability or real-network-connectivity
claim. It uses a synthetic 1 mm support, 1 m² membrane area,
100 L m^-2 h^-1 flux, water at 298.15 K and 100000 Pa inlet absolute pressure,
and 75% shaft-to-hydraulic pump efficiency. These are illustrative choices,
not extracted experimental conditions. The same water density/viscosity come
from actual DWSIM IAPWS-IF97 properties and configure OpenLB and Julia.

Four separately defined selective layers span alpha=0.25, 1, 10, 100. Within
each pair the **same absolute selective resistance** is retained when support
k changes. Alpha always uses the baseline 4 µm support; it is never recomputed
against the improved support. These are sensitivity scenarios, not four
measurements. Units and snapshots are saved in JSON and CSV.

Julia calculates each pressure load. Python supplies that load to an actual
DWSIM pump as its discharge requirement and sets mass flow from rho*J*A.
Every pump is solved, saved, reloaded and solved again. There is no membrane
separation unit or retentate/permeate split in this flowsheet. A separate
pressure-loss unit is unnecessary for the present externally specified load;
it becomes useful when integrating a feed/retentate/permeate system or
crossflow circulation. That future interface must be investigated separately.

**Energy scope:** total here means all modelled transmembrane resistance
contributions. It is pump **shaft** energy, not electricity at the meter.
Motor losses, crossflow circulation, backwash, fouling, concentration
polarisation, solute rejection, cleaning, gravity and auxiliaries are absent.
No plant-wide energy or lifetime/material optimum is claimed.

## Verified numerical results on this Mac

The study passed on 8 November 2026. The final evidence snapshot is
`results/filtration/runs/20261108T192257Z-2d24de62/summary.json`.
It reuses verified native measurements from the first complete run
`20261108T191117Z-142f753e`, adds the 1 Pa operating-range check and performs
fresh DWSIM solves. The default launch command performs all native solves
afresh. A compact tracked record is [filtration_reference.json](../../../examples/filtration/layered_support/reference/filtration_reference.json).

Fine-grid measured permeability is **2.511215119e-13 m²** for 4 µm pores and
**6.130870215e-13 m²** for 5 µm pores. The support resistance decreases by
**59.0398%** at unchanged thickness. This improvement produces the following
actual DWSIM shaft-energy changes at fixed flux:

| Selective resistance / baseline support resistance, alpha | Selective share before → after | Total modelled shaft-energy saving |
| --- | --- | --- |
| 0.25 | 20.00% → 37.90% | 47.2319% |
| 1 | 50.00% → 70.94% | 29.5199% |
| 10 | 90.91% → 96.07% | 5.3673% |
| 100 | 99.01% → 99.59% | 0.5846% |

For alpha=1, specific shaft energy decreases from **7.293065001e-5** to
**5.140158564e-5 kWh/m³**, while required pressure decreases from
**196.9128** to **138.7843 Pa**. These small absolute values belong to this
synthetic resistance choice. They are not a claim about commercial membrane
energy consumption. Other selective layers change the absolute values and
the achievable percentage benefit.

The selective layer becomes dominant when alpha>1 for the baseline and
alpha>0.409602 for the improved support. It reaches 90% at alpha>=9 and
alpha>=3.686416, respectively. Support improvement therefore moves the system
toward selective-layer limitation; this is expected, not a simulation failure.

Seven native cases cover the two supports/two grids and three pressure/domain
checks. The N24→N40 permeability changes are below 0.5%; both fine-grid
permeability and velocity-profile errors are below 0.1%. Across 1–20 Pa,
inferred k changes by less than 0.003%; doubling the domain changes it by
less than 1e-6%. The doubled-domain flow ratio is 0.5000833, within the stated
0.1% criterion. The upscaled estimated maximum pore Re is below 0.004.

All eight DWSIM configurations pass mass/energy, pressure, shaft/energy-stream
power and save/reload/solve checks. Shaft power agrees with Q*delta_p/eta
within 1e-9 relative; the largest actual flash energy residual is **0.0916%**
of shaft power, below the declared 0.1% criterion. Outlet density and viscosity
changes remain below 0.01%, supporting the constant-property approximation.

Direct Julia tests pass **280 assertions**, including **65** new layer-model
assertions. Python integration/UI regression tests pass **13/13**, and seven
native invalid/unconverged-input cases fail cleanly. Stored-result verification
checks **71 artifacts** and recomputes eight Julia layer loads; an altered CSV
copy is rejected. That replay check does not rerun DWSIM. All 2,070 checked
original OpenLB files still match the pinned archive.

The comparison figure was visually inspected and regenerated with points/lines
on logarithmic axes; its renderer/input/image hashes are in the run's
`figure_provenance.json`. The CLI has no new browser interface. Existing Panel
tests passed, but the previously documented manual rotation-gesture check
remains separate from this study. No installation or coupling blocker remains.
The synthetic graded-geometry comparison is the next implementation milestone;
experimental reconstruction of a particular material is outside current scope.

## Experimental evidence and next geometry

The current analytical and numerical checks validate the model implementation.
They do not validate its prediction for a manufactured membrane. The
[article workspace](../../../literature/README.md) records the supplied papers and
the extraction contract. A separate, limited support closure has
now been calibrated to Biesheuvel measurements; its primary held-out screen
fails. See [the experimental report](calibration.md).
Whole-membrane permeance cannot be assigned to the selective layer
unless the other layer contributions are independently identified.

The three requested main articles and available figures have now been reviewed.
[The combined evidence review](../../../literature/reviews/filtration/evidence_review.md)
retains nine Biesheuvel experimental support-property cases and the Jan
whole-membrane comparison, separately from Mo's model results. They describe
different specimens and must not be subtracted to manufacture a selective-layer
resistance. The Biesheuvel table/formula/powder sizes are now audited against
the same study in the author thesis. Mo's Word supplement is now fully
reviewed; its reference Zhu et al. (2021), DOI 10.1016/j.ceramint.2021.01.088,
has now supplied [a bounded matched M1 series check refined from the PDF](zhu_pdf.md).
Both the nominal and reading-sensitivity screens pass with the clearer reading;
a low-flux model discrepancy remains. Unreported conditions and uncertainty
limit the strength of literature comparisons, but do not block the computational
geometry study. A permeability match alone does not identify a unique structure.

Mo, Li and Yang (2022) directly study matching support structure to a fixed top
layer. The user-supplied main article and all seven figures have now been
reviewed; see [the evidence review](../../../literature/sources/mo2022/review.md). It supports
the hydraulic plateau and highlights the intermediate-layer contribution.
Its reported layer resistances are model quantities, not accepted experimental
inputs. Supporting Information Section S1 and Tables S2/S3 have now been read
from the supplied Word file. They resolve the references but do not provide
all numerical reproduction settings or a matched experimental layer series.
[Article record](https://pubmed.ncbi.nlm.nih.gov/35640511/).

The public Jan et al. (2025) paper provides whole-membrane context, but changes
the selective coating. Its uncoated membrane already contains a selective
layer. Those values cannot isolate a support-only benefit.
[Full paper](https://ris.utwente.nl/ws/portalfiles/portal/501138486/1-s2.0-S0011916425005612-main.pdf).

That geometry milestone is now implemented using Zhu's traceable effective
non-support impedance for scale translation, retaining the same-study holdout
and its limitations. Its OpenLB geometry keeps the selective-side mouth, interface region,
support thickness and representative area fixed, while varying only the
interior widening/constriction profile. The case composes existing OpenLB cylinder/cone
indicators, preserving the official example lifecycle and external code layout.
It checks geometry connectivity, mesh resolution at the narrowest feature,
domain/tap placement, conservation, pressure range and linearity. Mechanical
strength and coating feasibility remain literature-based scenario constraints,
without claiming a manufactured or certified design.

The result quantifies the conditional support benefit. The subsequent
[explicit DWSIM hydraulic flowsheet](once_through.md) now passes 18 cases
with recovery and feed/permeate/retentate balances. It accounts for fresh-feed
pressurization, prescribed upstream loss and retentate letdown. The
[crossflow recycle extension](crossflow_recycle.md) also passes: separate
pressurization/circulation pumps, a reference channel-resistance law, purge,
temperature control and 25 actual DWSIM cases. At the nominal loop point,
circulation dominates and widest-support shaft-energy saving is 1.91%.
Next, select published module/channel geometry and a transport or fouling law
that gives circulation a performance role. Compare support benefit with uncertainty in the fixed
layer/interface, operating variability and mechanical constraints. Neither
the original transmembrane-only metric nor the once-through/recycle scenarios
are total plant electricity demand. Reacting, non-isothermal
bioprocess showcases remain later milestones rather than additions to this
bounded hydraulic study.
