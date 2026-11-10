# Isothermal adsorption: first synthetic connected example

[Scientific report: purpose, theory, results and conclusions](REPORT.md).

Question: at fixed dry adsorbent inventory, feed flow and outlet-quality limit,
how does particle size change treated volume before breakthrough and pumping
energy? Phenol on activated carbon is the starting literature family. The
initial bed is a synthetic periodic FCC sphere arrangement, not a reconstructed
experimental specimen. Reaction, heat coupling and regeneration are deferred.

First release verified on 10 November 2026: 8 native bed cases, 7 native
scalar-TRT cases, 4 actual DWSIM cases, 46 new Julia assertions and 86 Python
tests. The existing Julia package still passes 902 assertions. Independent
replay checks 1,178 artifacts and 11 column calculations. This is a synthetic
prototype; experimental column validation and native bed transport closures
remain pending.

A separate [finite-bath and film/particle module](film_particle/README.md)
now tests Alhamed's date-stone carbon. It passes the declared withheld batch
loading screen but fails column transfer. This does not replace the accepted
synthetic loop or promote it to experimental validation.

The following [contact and inventory audit](column_contact/README.md) keeps
those particle parameters frozen and tests a common mixed-section column.
Its analytical/conservation checks pass; both diameter comparisons still fail.
The [report](../../docs/verification/adsorption/column_contact.md) distinguishes density-derived bed
mass, unavailable run masses and published fitted breakthrough descriptors.

From the project root on the existing macOS ARM64 setup:

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/preflight.py --output results/adsorption/preflight.json
uv run python -m examples.adsorption.evidence
uv run python -m examples.adsorption.transport
uv run python -m examples.adsorption.bed
julia --project=julia julia/adsorption/runtests.jl
uv run python -m examples.adsorption.run
uv run python -m examples.adsorption.check
uv run python -m examples.adsorption.reference
uv run pytest -q
```

No new Python/Julia dependencies are required. The selected Julia executable,
shared project and JuliaCall configuration remain in the existing bridge.
The standalone numerical module is `julia/adsorption/AdsorptionColumn.jl`, loaded
once per process; the frozen filtration package files are unchanged. Calls
stay on the main Python thread. Native applications use stock OpenLB 1.9r0
facilities and the existing external build helper/GCC 15.2 configuration.

The two native applications retain the helper's executable basename
`capillary3d`, in separate `.local/adsorption/transport/` and
`.local/adsorption/bed/` directories. The first is a prescribed-velocity D2Q5
scalar-TRT analytical test. The second is a D3Q27/TRT sphere-bed flow solve.
There are no new descriptors, collision models, templated geometry classes or
OpenLB data structures. Source, builds and output stay outside original OpenLB.

The scalar diffusion rate belongs to the odd TRT channel; the even OMEGA is
derived with magic=3/16. The benchmark checks a periodic Fourier mode against
its independently calculated decay and advection. It is not a measured bed
dispersion coefficient. The bed uses stock sphere indicators, halfway
bounce-back, exact-difference forcing, ValueTracer and velocity/density
functors. It checks grid refinement, time step, pressure scaling, particle
scaling, conservation, Mach and density variation. The initial acoustic scaling
was rejected. Smaller saved time steps retain the same 1e-4 density screen;
optional native velocity estimates still undergo new convergence checks.
The revised suite uses even tau=0.5+0.0009*N; its separate temporal check uses
tau=0.525. All use magic=3/16. Initial fields are stock equilibrium estimates,
not exact population checkpoints. A separate cold run checks initial-state
invariance.

Bed outputs can be reused only when their complete command, native initial-field
hash and binary match. Changed sources/requests require a fresh output folder;
archive `results/adsorption/bed/` before a changed study. Accepted process runs
have unique IDs and save complete native evidence, inputs, implementation and
dependency-lock snapshots. They do not depend on a later mutable native run.

Python owns coordination, exports and plotting. Julia owns a conservative
cell-centred, implicit column balance with Langmuir equilibrium and a spherical
intraparticle linear-driving-force approximation. Actual OpenLB permeability
sets Julia's pressure requirement, which drives a real DWSIM water-carrier pump.
Every pump is solved, saved, reloaded and solved. Phenol uptake and transient
inventory are external Julia calculations; no native DWSIM adsorption unit is
claimed. This initial process case has no bed switching, regeneration,
separation train or full economic objective.

Results include SI JSON, comparison and sensitivity CSVs, concentration/loading
profiles, full breakthrough/mass-inventory histories, real `.dwxmz` pump cases,
native field CSVs, an actual native-field VTP, hashes and a scientific figure. Missing breakthrough is JSON
null. CSV node indices are zero-based; x/y/z and axial cell centres are metres.
Positive flow is +x in OpenLB and inlet-to-outlet in the column. Concentration
is kg solute/m³ liquid; loading is kg solute/kg dry solid. Particle-envelope
dry density differs from bed bulk density, which is `(1-eps)*rho_particle`.

Baseline axial dispersion is deliberately zero. Chosen positive-dispersion
sensitivities use `Dax = factor * interstitial_velocity * particle_diameter`
with factors 0.1 and 1; these are scenarios, not measured correlations or
experimental confidence intervals. Film resistance and actual bed-scale
transport closures remain future work. Published fitted intraparticle
diffusivity is transferred unchanged as a scenario, not re-fitted to the column.

Read [the source audit and scientific status](../../docs/verification/adsorption/synthetic_bed.md)
before interpreting this synthetic loop as a material prediction. The supplied
Alhamed (2009) PDF has now been reviewed in a [separate date-stone-carbon
study](../../docs/verification/adsorption/alhamed.md). Its evidence and numerical checks
pass, while the reduced model fails calibration and particle-size transfer.
The original synthetic example remains unchanged. Reproduce the new comparison
with `uv run python -m examples.adsorption.alhamed.source`,
`uv run python -m examples.adsorption.alhamed.study` and
`uv run python -m examples.adsorption.alhamed.check`.
