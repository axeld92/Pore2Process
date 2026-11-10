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

See the [reproduction guide](docs/reproduction.md) and
[publication record](docs/code_publication.md). The initial source publication
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
sensitivities remain explicit in [the final report](docs/operating_window.md).
Further filtration physics is an optional future extension, not the active
next milestone.

## 6. Develop the isothermal adsorption case

**In progress separately; implementation publication pending.**

The next research question is: at fixed adsorbent inventory and required outlet
quality, how does particle size affect usable treatment capacity and pumping
energy? This adds mass transfer and breakthrough modelling to the established
OpenLB → Python/Julia → DWSIM workflow. Source selection, model qualification
and process coupling belong to that case's own verification and publication.

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
