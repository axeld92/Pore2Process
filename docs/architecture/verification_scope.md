# Verification scope for the computational tool

Clarified by the user on 8 November 2026. Pore2Process is a computational
workflow, using published articles and their available supporting information
as its experimental evidence. New experiments, raw-data acquisition, author
contact, private datasets or external repositories are not project prerequisites.
The available evidence is determined paper by paper, without assumptions based
on the authors' country or institution.

## What we will establish

1. **Numerical verification:** analytical capillary/Poiseuille and Darcy
   references, conservation, grid/domain/pressure-range checks and meaningful
   software tests. These assess whether the implemented equations are solved
   correctly and the numerical result is sufficiently resolved.
2. **Literature-based model evaluation:** compare with the published tables
   and figures we actually have; record digitization, specimen scope and
   unreported conditions. Reserve points for prediction checks where useful,
   without calling a same-study split a new independent experiment.
3. **End-to-end workflow verification:** a change in a synthetic or
   literature-informed OpenLB geometry changes a derived property, Julia's
   translated load and an actual DWSIM process result. Preserve units, input
   snapshots, assumptions and saved outputs.

These are sufficient to build and assess a credible computational showcase.
OpenLB remains the required detailed pore solver; Julia, Python and DWSIM retain
their defined roles. All OpenLB additions stay outside its original codebase,
using its existing components and example lifecycle.

## Meaning of pore geometry validation

That phrase previously referred to checking whether a simulated connected 3D
pore structure represents a particular experimental specimen. Matching its
permeability alone cannot establish that: multiple geometries can have the same
effective resistance. That stronger specimen-specific claim is **outside the
current acceptance scope**, rather than unfinished work that blocks progress.

Our geometries may be synthetic design hypotheses informed by published pore
sizes, porosity or processing trends. We can study a controlled geometric
change and test its numerical behaviour without claiming to reconstruct the
actual ceramic specimen. Tuning a model parameter to match a published
permeance is calibration of effective behaviour, not proof of a unique structure.

## Handling information the papers do not report

- Retain unknown temperature, viscosity, replicate count and measurement
  uncertainty as unreported; do not repeatedly request inaccessible data.
- Prefer viscosity-independent measured permeance/impedance comparisons when
  appropriate, as in the Zhu M1 series.
- When a simulator needs missing physical inputs, choose and document scenario
  conditions, use consistent properties across the connected calculation and
  examine sensitivity where the assumption matters. Those conditions are not
  attributed to the experiment.
- Treat digitization bounds as reading sensitivity, not experimental error
  bars or a confidence interval. Missing experimental uncertainty limits the
  strength of a claim; it does not prevent software implementation.
- Preserve genuine model discrepancies and previously declared screens. A
  failed comparison is a result to understand, not a reason to invent data,
  relax a threshold after fitting, or abandon a verified software connection.

Mechanical-strength observations from the literature may constrain a scenario,
but device certification and experimental validation of unchanged selectivity
are outside this software demonstration. Filtration pressure-load energy must
also remain distinct from total plant or crossflow circulation energy.

## Next executable milestone

The first synthetic fixed-opening OpenLB comparison is now
[implemented and verified through Julia and DWSIM](../verification/filtration/graded_support.md),
using Zhu M1 impedance as a bounded scale-translation reference. Its raw
geometry result and common baseline calibration are reported separately.

The [explicit hydraulic filtration stage](../verification/filtration/once_through.md) is now verified
for 18 actual DWSIM once-through cases, with Julia's closure, declared recovery,
feed/permeate/retentate streams, all unit/overall balances and saved/reloaded
flowsheets. Original OpenLB data remain intact and are replayed before reuse.
Its prescribed upstream-loss sensitivity does not predict an axial crossflow
channel pressure profile or solute rejection.

The [crossflow recycle extension](../verification/filtration/crossflow_recycle.md) now passes 25 actual
DWSIM cases, with declared circulation ratios, a separate channel-resistance
law, purge and temperature control. Fresh-feed and circulation work, membrane
loss, pressure letdown and removed heat remain distinct. Independent stream
balances, convergence/seed/tolerance checks, saved connectors and reload solves
are verified. This is a lumped hydraulic scenario, not resolved axial crossflow
or evidence of improved separation/fouling performance.

Next, use published channel/module geometry and an application-specific
transport/fouling relation to evaluate a circulation-cost/performance tradeoff.
Continue to distinguish analytical verification, literature support, calibration
and scenario assumptions. No inaccessible raw data or new experiments are prerequisites.
