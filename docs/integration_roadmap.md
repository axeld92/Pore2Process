# Toward the research workflow

The required architecture remains OpenLB for detailed pore-scale physics,
Julia for custom models/analysis/scale translation, Python for coordination,
data/visualization/interface, and DWSIM for process simulation. Cantera may
support future chemistry; it does not replace either mandatory simulator.

| Milestone | Verified state on 10 November 2026 |
| --- | --- |
| Python–Julia capillary-network reference | Complete numerics, bridge, exports and Panel app; canvas rotation also verified in the browser on 9 November. |
| Native tool compatibility | Complete: stock OpenLB channel and actual DWSIM pump cases, saved/reloaded and solved. Builds now remain outside OpenLB. |
| Detailed hydraulic reference | Complete: stock-component 3D OpenLB case, four-grid study, pressure/time-step/stop checks, conservation and SI exports. |
| Julia scale translation | Complete: measured permeability, Darcy flow and inverse pressure-load functions, usable directly from Julia. |
| Property-driven DWSIM case | Complete: actual DWSIM water properties configure the pore and translation stages; measured k sets a real pump discharge requirement. |
| End-to-end design comparison | Complete for idealized 50/40 µm straight pores: actual OpenLB → Julia → DWSIM results, intermediate validation, source/artifact hashes and saved flowsheets. |
| First research material/application | Selected: fixed-flux clean-water filtration, unchanged selective layer; alumina/SiC ceramic supports are the initial literature family. Synthetic layer study passes with seven native cases and eight actual DWSIM configurations. |
| Experimental support reference | Implemented: audited Biesheuvel Table 2, Julia reduced-model fit and held-out check. Error improves but the primary 25% screen fails; discrepancy is retained. |
| Matched support/whole-membrane reference | Implemented: Zhu M1 PDF reading, 11-point calibration/10-point holdout in Julia; nominal max error 4.56%, conservative reading-corner max 7.88% against the unchanged 10% screen. Bounded effective impedance; earlier JPEG result preserved. |
| Controlled fixed-opening geometry comparison | Complete: 4/5/6-µm interiors, 24 verified native cases, fixed literature-informed non-support impedance and nine actual DWSIM pump cases. Widest design reduces native body resistance 35.88% and conditional transmembrane shaft energy 11.30%. |
| Explicit hydraulic filtration equipment | Complete: 18 actual DWSIM once-through cases with feed/permeate/retentate, declared recovery, stock pump/splitter/valves, per-unit/overall balances, actual saved connector checks and reload solves. A separate assumed 20 kPa upstream loss reduces widest-design savings to 3.81%. |
| Crossflow recycle and circulation duty | Complete: 25 actual stock-DWSIM recycle cases, separate pumps, purge, temperature control, channel-resistance law, balances, seed/tolerance checks and saved/reloaded solves. Circulation is 83.09% of nominal shaft duty; widest-support saving is 1.91%, falling to 0.25% at threefold module flow. |
| Literature-constrained crossflow benchmark | Complete as a conditional scenario: real Julia source constraints and passive-solute balances, 32 solved/saved/reloaded DWSIM cases, two flow conventions and two policies. Selected tested speeds are 0.9/1.4 m/s; widest-support savings span 7.11–10.99%. Absolute source flow convention and physical fouling time remain unvalidated. |
| Explicit star-channel hydraulic definition | Complete: declared six-point polygon, unambiguous area/perimeter/velocity convention, 17 actual stock-TRT OpenLB cases, independent Julia duct reference and checked Re envelope 4.94–34.59. |
| TRT support update and biological-feed sensitivity | Complete: 24 revalidated support cases, 28 actual DWSIM recycle cases including seed validation; distributed Julia pressure/flux and conservative COD balance. At 12 LMH and Re=20 the widest support saves 11.20% of total pump shaft energy. Retention is a transfer scenario; sustainable flux remains unavailable. |
| Local connected interface | Complete: actual Julia/DWSIM runs, support controls, native 3D inspection, labelled profiles, per-session snapshots/downloads, invalid recovery and preserved original capillary tab. Browser and downloaded-file checks pass. |
| Higher-Re hydraulics, source audit and operating map | Complete: 16 stock D3Q27/TRT cases, 91 actual DWSIM solve/save/reload cases and a bounded critical-flux policy comparison. Native laminar-branch envelope is approximately Re=1000–2000; physical transition and sustainable operation are not established. |

The first connected technical acceptance criterion is met for the synthetic
hydraulic demonstration. See [the results](connected_hydraulics.md) and
[reproduction commands](../coupling/README.md). The separately calculated
porous load drives a DWSIM pump; there is no porous unit operation in that
flowsheet. The original reference remains available in its own app tab; the
new filtration tab runs the verified Julia/DWSIM workflow with per-session
inputs and downloads.

The selected research step is [the filtration energy question](filtration_research_case.md).
Axial support permeability supplies a layers-in-series pressure requirement at
fixed flux. The process metric is pump shaft energy per m³ of filtrate, with
the selective layer unchanged within each design pair. A synthetic straight-pore
screen establishes the resistance partition before introducing a graded body
with the selective-side opening held fixed. The [literature workspace](../literature/README.md)
records the published evidence and its limits. The [computational verification scope](verification_scope.md)
uses the papers and available supplements; unavailable raw data, new experiments
and exact specimen reconstruction are not acceptance gates.
All three requested main articles have now been reviewed. The
[evidence review](../literature/filtration_evidence_review.md) records nine
Biesheuvel experimental support-property cases, Jan's whole-membrane
comparison and Mo's fixed-top-layer model. The Biesheuvel table is now visually
corroborated using the author thesis; [the calibration result](experimental_calibration.md)
records the failed primary screen. Mo's Word supplement has now been reviewed;
[its evidence report](../literature/mo2022_supplement_review.md) identifies
Zhu et al. (2021), DOI 10.1016/j.ceramint.2021.01.088, as the matched-layer
source. Its supplied Markdown/figures and original PDF support [the refined M1 study](zhu_pdf_refinement.md).
Unreported hydraulic conditions and uncertainty remain explicit limits, rather
than requests for inaccessible data. Intermediate-layer
resistance must be identified separately from the selective layer in a
three-layer membrane; these papers cannot be combined as a matched specimen.

The **synthetic fixed-opening comparison is now complete**; see
[the numerical/process result](graded_support_study.md). Its common baseline
calibration exposes unresolved structural scaling rather than claiming an
absolute specimen reconstruction. All original OpenLB files are preserved.

The [explicit hydraulic filtration milestone](process_filtration.md) now passes.
The [installed-interface investigation](dwsim_filtration_interfaces.md) selects
stock splitter/pressure-loss operations while preserving Julia's hydraulic law.
Recoveries 25/50/75% and separate prescribed upstream losses 0/20 kPa are tested;
no chemical rejection or total-plant energy prediction is claimed.

The [crossflow recycle milestone](crossflow_recycle.md) now passes with a stock
mixer, cooler, generic recycle block and two separate pumps. The independent
channel resistance changes pressure loss with module flow. All 25 cases pass
mass/enthalpy, pressure, property, convergence and saved-archive checks, and an
insufficient-iteration case is rejected. No unavailable dataset was required.

The [crossflow literature search](../literature/crossflow_evidence_review.md)
has now led to receipt and review of four further articles. The
[received-paper review](../literature/crossflow_papers_review.md) retains six
García (2008) table rows and separate Chiu (2005/2006) figure readings. Source
energy arithmetic agrees within 0.912%, but the absolute bundle-flow convention
and some borrowed table/figure values remain unresolved. Frederic's 2D model
has no calibrated physical fouling time. No new simulation model is fitted.

The [bounded benchmark](crossflow_benchmark.md) now passes 846 Julia assertions,
50 Python tests and 32 actual DWSIM cases, with a separately rejected
nonconvergent solve. At 100 LMH, the table policy selects 0.9 m/s and the
cross-source monotonic policy selects 1.4 m/s among the tested points.
Two flow conventions and 0.8/1.0 capacity fractions remain explicit scenarios.
Finite-pass passive-solute balances hold module-inlet concentration fixed by
adjusting fresh concentration. DWSIM retains Water properties; this is not a
validated suspension-property or geometry-dependent fouling model.

The [explicit geometry and lower-flux biological-feed case](star_channel_showcase.md)
now complete that software extension. Native source, builds and results remain
outside OpenLB, and porous-media runs require at least TRT. The declared
geometry removes ambiguity in our own flow convention without identifying the
paper's unreported convention. Its low-Re hydraulic case does not inherit the
paper's Re=1388 critical-flux limit.

The [higher-Re continuation](operating_window.md) now completes that numerical
step, the supplied-paper audit and a material–operation map. Spatial, time-step,
periodic-domain, initial-state and pressure-law checks pass before the measured
coefficient enters Julia and actual DWSIM cases. The capacity policies retain
equal-Re/equal-velocity transfer, source reading bounds and shape-free monotonic
alternatives. These support conditional discrete-grid choices; they do not
validate physical transition, indefinitely sustainable flux or a fouling-time law.

The current hydraulic milestone is complete and is retained as a regression
example. A concrete optional hydraulic
verification, if selected, is a permeating-wall channel benchmark against the
existing distributed Julia model, using stock OpenLB boundary facilities in a
new project folder. It would test the present slowly varying approximation;
it would not resolve the separate experimental capacity-transfer question.
Reacting and non-isothermal additions are deferred at the user's request and
will require a selected application plus published kinetic/thermal inputs.
No inaccessible raw dataset, specimen reconstruction or new laboratory
measurement is a gate for these computational milestones.

The filtration hydraulic example is **closed for now** and remains a regression
benchmark. The active next example is **isothermal adsorption**, being developed
separately. Its research question is how particle size changes usable treatment
capacity and pumping energy at fixed adsorbent inventory and required outlet
quality. Its implementation and validation will be published in a separate
update after that work is ready.

Implementation choices for that milestone:

- Synthetic design parameters and declared clean-water scenario conditions;
  literature comparisons do not imply an exact manufactured specimen.
- Representative volume and boundary conditions; directional permeability or
  another property; analytical/convergence checks and applicable published comparisons.
- The validity range of a Darcy relation, including inertia, anisotropy,
  heterogeneous geometry, temperature/composition changes or reaction.
- Whether DWSIM should receive a pressure-load requirement, a resistance
  correlation or a dedicated pressure-loss unit operation for that application.
- Whether the next workflow needs a connected interface, parameter sweeps or
  sensitivity analysis. The validated CLI precedes optimization and additional
  process complexity.

The OpenLB code policy in [openlb_code_policy.md](openlb_code_policy.md) applies
to every extension: use existing components, follow official example structure,
and keep every project addition outside the original source tree.
