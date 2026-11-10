# Implementation and acceptance — 7 November 2026

Latest model/evidence step, 10 November: the [common column contact and
inventory audit](column_contact_audit.md) preserves the batch-qualified
particle component unchanged. The printed bulk density implies a derived
7.032455 g bed mass; no fraction-specific weighed mass is reported. The
observed breakthrough window does not establish saturation or prove that
density wrong. Table 5 DL/Pe remain fitted breakthrough-shape descriptors,
not independently measured fluid dispersion. A common equal-section model
is selected on 0.8 mm only, then transferred to 0.45 mm with fixed particles,
film correlation, isotherm and full active inventory. **Both comparisons fail**:
t50 is 87.49/118.32 min versus 30.75/55.36 min in Figure 7. Conservative
implementation and selected time/radial refinements pass. No simulator or
dependency changes were made; no new OpenLB/DWSIM solve or browser check is
claimed. Next is a synthetic stock-TRT packed-bed inert-tracer RTD benchmark,
without claiming specimen reconstruction or experimental adsorption validation.
Fifty new Julia assertions, 51 preserved particle assertions and 128 Python
tests pass. Independent replay checks 24 cases, 66 calibration evaluations,
83 artifacts and 30 source files. The existing PyVista UI regression required
host graphics access after an AppKit abort in the restricted runner; no code
change or new adsorption browser acceptance is implied.

Previous model/evidence step, 10 November: the [separate finite-bath and
film/particle model](film_particle_adsorption.md) is implemented. Source
equilibrium matches all seven published batch endpoints within 2.03%; shared
transport fitted only to the 0.8/1.5 mm curves predicts the withheld 0.45/0.225
mm loading curves with RMSE 1.33/0.42 mg/g. **Both the fully transferred and
one-anchor column routes still fail Figure 7.** No equilibrium or physical
inventory refit is hidden in that comparison.

The new Julia suite passes 51 assertions and the full Python suite passes
112 tests. Independent replay verifies 57 saved Julia cases, seven analytical
endpoints, 141 artifacts and 20 source files. Conservation and selected
radial/axial/time refinements pass. The original synthetic example and the
earlier failed LDF study both reverify, including original OpenLB source and
saved DWSIM evidence. No dependency, lockfile or simulator-source change was
made. Next is a common column contact/flow-distribution audit; reactive and
thermal additions remain deferred.

Previous evidence step, 10 November: the [supplied Alhamed PDF
comparison](alhamed2009_comparison.md) retains 48 digitized Figure 7 markers,
seven published batch fits and all 18 column-fit rows for separate date-stone
carbon. Ten actual Julia column calculations pass source/artifact replay,
conservation and the declared refinement checks. **The reduced model fails
calibration and held-out particle-size prediction.** This is retained as a
failed physical comparison; it does not invalidate the completed synthetic
example or establish experimental hydraulic validation. No new OpenLB or
DWSIM solves, simulator modifications, dependency changes or reactive/thermal
additions were made. Next is a finite-bath equilibrium/uptake audit followed
by a tested film-plus-intraparticle alternative using the same paper. The
full Python suite passes 94 tests; the preserved synthetic-example checker
still passes all source, artifact and saved-flow checks.

Previous completed milestone, 10 November: the [first synthetic isothermal adsorption
release](adsorption_example.md) passes **8 native FCC bed cases, 7 native
scalar-TRT reference cases, 46 new Julia assertions and 86 full-suite Python
tests**. The existing Julia package also passes all 902 assertions. Three
particle-size cases and a property anchor complete actual DWSIM solve/save/
reload checks. Independent replay covers 1,178 artifacts, 11 Julia column
calculations and all four saved pump topologies. All 2,070 checked original
OpenLB files remain unchanged; Python/Julia locks and the completed examples
are preserved.

At fixed 4 kg adsorbent inventory and 155.015 mL/min flow, reducing particle
diameter from 2.3 to 0.575 mm increases conditional treatment volume before
the chosen 5% breakthrough limit by 1.8878%, while bed-attributable pump power
rises 15.9938 times. This is a low-flow, clean synthetic bed with long contact
time; the absolute bed duty is small and excludes the rest of the plant.
The chemical inputs are literature-informed transfer scenarios. Experimental
column validation, resolved native bed mass-transfer closures and a full
adsorption-service/regeneration flowsheet are not completed or claimed.

The evidence audit retains an ambiguous printed density, conflicting affinity
units/fits and undigitized error bars. No unsupported density is silently
substituted. A rejected native acoustic scaling, a low-head flash diagnostic
and an initial stiff-uptake solver failure are retained. Revised scaling,
stock forced-liquid DWSIM and a conservative residual line search pass the
unchanged acceptance screens. The figure was visually inspected; no new
adsorption Panel UI/browser verification is claimed. See
[commands](../adsorption/README.md) and [the compact reference](adsorption_reference.json).

Earlier completed milestone, 9–10 November: the [higher-Re operating map](operating_window.md)
completes the three requested hydraulic/evidence/design-map steps. **16 native
D3Q27/TRT cases pass 144 numerical screens; 91 actual DWSIM cases solve, save
and reload; 902 Julia assertions and 74 Python tests pass.** A deliberately
one-iteration recycle is rejected. Independent replay checks 682 local
artifacts, native/support ancestry, 1,620 capacity comparisons and 1,458
discrete selections. All 2,070 checked original OpenLB files remain unchanged.

The measured fully developed laminar-branch envelope is approximately
Re=1000–2000. Spatial, temporal, domain, initial-state and pressure-law checks
pass before the native coefficient is used by Julia. This does not validate
physical transition, entrance loss or permeating-wall CFD. The longer-domain
case required a retained native-field restart and larger wall-time budget;
scientific tolerances were unchanged. The DWSIM adapter input guard now allows
20,000 iterations; the independent perturbed-seed case takes 8,796 iterations
and reproduces the nominal pump/cooling duties within 1e-8 relative.

At 12 LMH and Re=1388, widening the support interior from 4 to 6 µm reduces
mean TMP from 1.2211 to 1.0831 kPa, but total pump shaft energy falls only
**0.18166%**, from 0.056703822 to 0.056600813 kWh/m³. Circulation supplies
98.7234% of baseline shaft duty. Fixed effective non-support resistance rises
from 68.5026% to 77.2330% of membrane resistance. This separates membrane
resistance dominance from total process-energy dominance.

The supplied-paper audit preserves Chiu's biological-feed critical-flux
protocol, TiO2/shared-data caveats, reading sensitivities and unknown physical
fouling time. Equal-Re/equal-velocity transfer and linear/monotonic capacity
policies can select different grid points; they remain conditional scenarios.
The figure was visually inspected. This report/CSV extension adds no new Panel
controls or browser-pass claim. Existing examples, the live interface, locks
and unrelated files are preserved. Reaction and non-isothermal work remain
deferred. [Commands](../operating_map/README.md) and
[the numerical reference](operating_map_reference.json) support reproduction.

Earlier completed milestone, 9 November: the [explicit-channel and biological-feed
showcase](star_channel_showcase.md) passes **17 native channel cases, 24 freshly
verified TRT support cases, 28 actual DWSIM crossflow solve/save/reload cases,
878 Julia assertions and 61 Python tests**. A separate eight-case capillary
TRT regression, three geometry audits and 25 native invalid/unconverged checks
also pass. All 2,070 checked original OpenLB files remain unchanged.

At 12 LMH and inlet Re=20, widening the support interior from 4 to 6 µm with
fixed 4 µm openings reduces total pump shaft energy by **11.20%**; fixed
effective non-support resistance rises from 68.5% to 77.2% of membrane load.
The Chiu biological-effluent case adds a conservative COD balance and an
explicit transferred retention sensitivity. Sustainable flux and fouling time
at the new geometry/speed remain unavailable. DWSIM uses a pure-water carrier
with its stock forced-liquid flash route and an independently checked phase
assumption. The new recycle mass tolerance is tighter, not relaxed.

The Panel app now has capillary and filtration tabs. Scientific/interface
verification is distinguished in [the report](star_channel_showcase.md) and
[browser record](browser_checks.md). Earlier runs and their source assumptions
remain preserved; earlier BGK results are not relabelled as TRT.
Actual browser runs, geometry changes, native 3D rotation, invalid recovery,
session isolation, repeated tab use and a complete Firefox-downloaded ZIP all
pass. The four requested software extensions are complete. The next research
step and retained physical unknowns are explicit in the report and roadmap.

Earlier milestone, 9 November: [crossflow recycle and circulation energy](crossflow_recycle.md)
passes 25 actual DWSIM recycle cases and a deliberate nonconvergence rejection.
Separate feed/circulation pumps, purge, channel resistance and temperature
control close all unit/overall mass and energy balances. Circulation is 83.09%
of nominal shaft duty; widest-support savings are 1.91%, falling to 0.25% at
threefold module flow. Current regressions pass **518 Julia assertions and
43 Python tests**. Published channel geometry and an application-specific
transport/fouling law are the next research extension.

The preceding milestone, 9 November: [explicit DWSIM hydraulic filtration](process_filtration.md)
passes 18 actual once-through geometry/recovery/channel-loss cases, including
all unit/overall balances, actual saved connector ports and save/reload solves.
The earlier OpenLB results remain unchanged and are reverified. At 50% recovery,
the widest support saves 11.30% with zero separate upstream loss, or 3.81% with
the assumed 20 kPa upstream loss. Its regression counts were **422 Julia
assertions and 33 Python tests**. The crossflow extension is completed above
and recorded in detail below.

The preceding milestone, 8 November: [the fixed-opening interior-widening
study](graded_support_study.md) passes all 24 native cases and nine actual
DWSIM design/scenario cases. Native support resistance falls by 35.88% and
conditional transmembrane shaft energy by 11.30% for the widest interior.
Its regression counts at completion were **372 Julia assertions and 25 Python tests**.
The sections below preserve the earlier milestones.

The preceding literature extension, 8 November: [the Zhu PDF refinement](zhu_pdf_refinement.md)
passes the unchanged 10% nominal and reading-corner screens (4.56% and 7.88%
maximum error). The earlier JPEG reading and failed Biesheuvel 25% screen are
preserved below. These are bounded literature comparisons, not exact specimen
reconstructions or experimental confidence intervals.

[Scope clarified by the user](verification_scope.md): this is a computational
tool using papers and available supplements. New experiments, inaccessible raw
data and actual-specimen 3D geometry validation are not project prerequisites.
Numerical verification, literature-based evaluation, explicit assumptions and
traceable end-to-end coupling remain required.

The Julia solver, Python bridge, CLI, exports, and Panel application are
implemented locally. Numerical acceptance passed on this Mac. Browser checks
are recorded separately below; HTTP/import checks are not visual evidence.

## Reproducible environments

The standard-library preflight ran before package installation, then ran again
in `.venv` with actual import probes. Python 3.12.9 and Julia 1.11.6 were reused.
uv 0.13.7 was installed in `.tools/`; all direct Python dependencies and pytest
were resolved in `uv.lock`. `uv sync --locked --offline` succeeded after setup.
The `julia/` project contains the package, exact PythonCall 0.9.37 compatibility,
and generated manifest. JuliaCall is also 0.9.37. Repeating bootstrap verified
the project without changing the lock/manifest hashes. The normal user Julia
depot is used consistently. No runtime installation occurs in the app.

## Numerical and integration acceptance

Direct Julia tests passed **191/191 assertions** (161 scientific/geometry/scaling
and 30 zero-drop/validation checks), including a final run with PythonCall in
the shared manifest. Tolerances are those in the supplied handoff: 1e-9 relative
for nonzero reference flow/permeability, 1e-9 times the pressure scale for
pressure, and 1e-9 times reference total flow for absolute conservation errors.
The tests cover pressure offsets/reversal, radius/viscosity/length scaling,
single-capillary and asymmetric networks, deterministic geometry/orientation,
closed transverse boundaries, zero drop, and invalid/non-finite inputs.

The final Python suite passed **11/11 tests** with macOS graphics access:
direct Julia/bridge equality; array shapes/dtypes/index conversion; analytical
values and SI/display conversions; changed pressure/radius/viscosity; no shared
mutable arrays; JSON/CSV/VTP/PNG/ZIP round trips; zero-drop nulls and validation;
rejection of Julia calls from worker threads; and UI state/snapshot behavior.
The UI test attaches the view to a Bokeh Document and checks that the same
mounted VTK model survives repeated runs. The full suite cannot be run inside
the restricted terminal sandbox because PyVista's Cocoa initialization aborts
there; it passes with normal host graphics access. This restriction is separate
from the numerical solver.

The generated CLI result in `results/default/` contains all five required
files. JSON uses standard nulls for unavailable values, CSV endpoints are
zero-based and signed from-to flows are labelled, and VTP data round-trips
correctly. The documented CLI was rerun after implementation stabilized.

| Quantity | Observed Python → Julia default result |
| --- | --- |
| Nodes / throats | 96 / 224 |
| Length / representative area | 0.005 m / 1.6e-5 m² |
| Inlet flow | 7.853981633974476e-10 m³/s |
| Permeability | 2.454369260617024e-12 m² |
| Default column flow | 1.227184630308512e-8 m³/s |
| Relative flow error | 1.22e-15 |
| Maximum pressure error | 6.39e-14 Pa |
| Maximum interior imbalance | 3.67e-25 m³/s |
| Inlet/outlet mismatch | 1.55e-24 m³/s |
| Normalized conservation error | 1.97e-15 |

These are observed values; analytical reference values and portable thresholds
are documented independently in the README. Small roundoff differences across
processes are expected. `docs/verified_reference.json` preserves the observation
and versions. One cold initialization took 2.25 s, first bridge call including
compilation 1.75 s, warmed default core 0.128 ms, and the largest UI network
(1280 nodes / 3456 throats) core 3.52 ms. These are individual local timings
after package precompilation, exclude plotting/serialization/export, and do not
establish a performance benchmark or speedup.

## Browser acceptance

The initial in-app-browser check displayed the colored 3D network, numeric
default results, disabled Run state, and pressure plot. Reload/control calls
then timed out, so subsequent checks used a separate Firefox tab. Firefox
displayed the model explanation and controls, labelled pressure profile,
pressure colorbar, physical-scale tubes, and axis ticks in mm.

A live second-run check caught Bokeh 3.10 moving synchronous callbacks onto
workers and a resulting Julia GC/GIL deadlock. Async event-loop wrappers and
a main-thread bridge guard fixed it; repeated Run then completed. The browser
showed 47.1239 µL/min at 100 Pa and 94.2478 µL/min at 200 Pa with unchanged
permeability. Changing column length from .02 to .04 m changed column flow from
0.736311 to 0.368155 mL/min without changing network run ID `cd92d4ad`.
The browser download was saved as `results/pore2process-cd92d4ad.zip` and opened
independently: all five files were present, inlet pressure was 200 Pa, column
length was .04 m, and exported flow matched the changed simulation.

Further final browser checks and any remaining limitations are recorded in
`docs/browser_checks.md`.

The final build also passed live invalid-diameter validation with a concise
message, recovered for another successful default run, and displayed a zero-drop
run with unavailable permeability/column flow. The 3D pane remains mounted across
updates and browser wheel zoom was visually verified. The original unverified
drag-to-rotate gesture was closed by an actual canvas drag on 9 November;
`browser_checks.md` records the screenshot. There are no remaining installation, numerical, bridge, or export
blockers on this machine.

## Scope and next step

This is an idealized model/software-connection validation, without empirical
validation for a real porous material. The app deliberately uses bounded
synchronous solves and one persistent Julia runtime, with per-session inputs
and completed snapshots. No background worker or cancellation is implemented.
The VTK client needs browser graphics support and Panel's vtk.js assets on first
load; offline browser asset availability has not been established.

The project-wide architecture now explicitly requires OpenLB for the detailed
pore stage and DWSIM for the process stage, with Julia scale translation and
Python coordination. Native compatibility checks completed on 8 November 2026;
see [compatibility_check.md](compatibility_check.md). The next executable
integration step and unresolved application choices are recorded in
[integration_roadmap.md](integration_roadmap.md). Heterogeneous network radii and
imported geometry remain possible later research extensions. No remote
repository or deployment was required.

The 8 November stock-component hydraulic extension and connected workflow
now pass. OpenLB's existing components replace the previous custom cylinder
helper, all project additions/build outputs stay outside its original tree,
and 2,070 original files match the archive. The four-grid hydraulic study and
its sensitivity checks pass; direct Julia now has **215 passing assertions**,
the original Python suite has **11/11**, and **25** native invalid/unconverged
cases fail cleanly. See [the current benchmark](hydraulic_benchmark.md).

The first complete **OpenLB → Python/Julia → DWSIM** comparison also passes:
actual 50/40 µm pore simulations supply measured permeability, Julia derives
a fixed-flow pressure requirement, and actual saved/reloaded DWSIM pump cases
show the expected power change. All stages share DWSIM-derived water properties.
See [verified connected results and limitations](connected_hydraulics.md) and
[launch commands](../coupling/README.md). The original app, bridge and dependency
locks/manifests remain unchanged. This adds a CLI demonstration and scientific
figures, without a new connected UI or a real-material validation claim.

## Selected filtration application — 8 November 2026

The [fixed-flux, fixed-selective-layer study](filtration_research_case.md) now
passes. Seven fresh native support/validation cases, including 1–20 Pa and a
doubled domain, supply measured permeabilities to a separate Julia layer model.
Eight actual DWSIM pump configurations are saved, reloaded and solved. The
final snapshot reuses the verified original support exports while adding the
lower-pressure check and rerunning DWSIM; source run IDs are retained.

Changing synthetic support radius from 4 to 5 µm reduces support resistance by
59.04%. Actual modelled total shaft-energy savings range from 47.23% at a 20%
baseline selective share to 0.585% at a 99.01% share. No experimental calibration,
fixed-interface grading, mechanical-strength validation or selectivity claim
is implied. The absolute energy is conditional on the supplied layer resistance.

Direct Julia now passes 280 assertions (65 new), Python has 13 passing tests,
seven native invalid/unconverged cases fail cleanly, and replay verifies 71
artifacts plus rejection of an altered copy. The original 33-artifact connected
demonstration still verifies, the original native build still succeeds and
the OpenLB archive comparison still checks 2,070 files. Dependency locks and
the bridge runtime configuration remain unchanged. The scientific figure was
visually reviewed; no new browser interface was built or claimed.

Experimental inputs will come from articles. The first full-text/supplement
request and extraction contract are in [literature/README.md](../literature/README.md);
local article Markdown/PDFs can go in ignored `literature/inbox/`. At this
milestone no experimental calibration dataset had been accepted. The next step was to
review Mo et al. (2022) and its experimental validation sources, select a matched
ceramic-membrane dataset, then introduce a stock-indicator graded support with
the selective-side opening held fixed. Bioelectrochemistry and the reacting,
non-isothermal bioprocess ideas remain later options.

## Main-article evidence reviewed — 8 November 2026

Mo (2022), Biesheuvel and Verweij (1999), and Jan et al. (2025) have now been
read from the supplied Markdown with metadata and available figures. The
[combined review](../literature/filtration_evidence_review.md) and structured
extraction records preserve locations, units, source hashes and eligibility.
Nine experimental support-property rows from Biesheuvel and two published
whole-membrane permeances from Jan were transcribed as candidates. No production
simulation parameters were changed or experimental calibration claimed.

At this review the outstanding checks were Mo's Supporting Information
(Section S1 and Tables S2/S3), the Biesheuvel Table 2 audit, and Jan's supplement
for interface and test-condition evidence. The Biesheuvel audit is now complete
as recorded below. Intermediate resistance
must remain distinct from selective resistance; values from these different
specimens cannot be subtracted to infer a measured selective-layer property.

## Experimental support closure — 8 November 2026

The public author thesis reproduces the Biesheuvel study. Its Table 2, CK
formula and nominal powder means were rendered and visually inspected;
[the audit](../literature/biesheuvel_source_audit.json) records all locations
and hashes. All nine extracted cases and the permeability multiplier agree.
This resolves the table audit without creating another independent dataset.

The new [fixed protocol](../calibration/protocol.json) fits one Julia closure
multiplier to AKP-30 temperature-series rows 2–4, then evaluates held-out
hold-time rows 5–7. Fitted b is 2.98438. Primary held-out mean absolute relative
error improves from 66.66% to 21.11%, but row 5 has +30.89% error and fails the
25% project screen. Training, transfer and out-of-domain errors are all retained.
No threshold or split was tuned after fitting; no confidence intervals were
invented. This is a retrospective same-study check of a reduced closure, not
experimental validation of an OpenLB geometry or the complete membrane.

The existing Julia energy model separately propagates measured versus fitted
support k with identical assumed fixed-layer resistance. It shows material
model error materially changes the conditional savings. Native OpenLB and
DWSIM cases were preserved; no new native run or browser check is claimed.

Final direct Julia tests pass **315 assertions**, including **35** calibration
checks. The complete Python suite passes **16/16**, including three new tests
of real-Julia fit/oracle agreement, held-out isolation, source/split guards,
conditional energy and replay/tamper rejection. The final saved nine-case run
replays successfully. Enabling `--require-screen-pass` correctly rejects its
scientific outcome while the normal integrity checker passes. Plots were
visually reviewed. Locks, Julia environment/depot and bridge configuration
remain unchanged. [Commands](../calibration/README.md), [result report](experimental_calibration.md)
and [compact reference](experimental_calibration_reference.json) preserve the outcome.

At this stage the next evidence request was Mo's S1/S2 and original PDF,
followed by its matched experimental sources. The Word supplement has since
been received and reviewed as recorded below. The subsequent scope clarification
removes raw-data acquisition and specimen reconstruction as project gates;
unknown conditions and mechanical feasibility limit claims instead.

## Mo Word supplement reviewed — 8 November 2026

The supplied supplement was read directly from its Word XML and rendered with
the bundled document tool. All eight pages, three tables and the flowchart
were visually reviewed. [The extraction](../literature/mo2022_supplement_extraction.json)
preserves source/render hashes, model settings, ten S2 rows, nine S3 rows,
19 references and explicit unavailable values. All structured table entries
were checked against the original XML, and all 15 supplied permeance values
were checked against the existing actual Julia SI-conversion function.

S2 contains six permeance rows and four pore-size rows, with no row containing
both. Model temperature/pressure are not measurement conditions for all cited
experiments. The printed density 3.7 kg/m³ is retained and flagged rather than
silently corrected. Primary-source checks also find a tubular-thickness
convention ambiguity and distinguish bare-support from coated-membrane values.
No new calibration or simulation acceptance is inferred from this compilation.

The priority experimental paper at that stage was **Zhu et al. (2021), DOI
10.1016/j.ceramint.2021.01.088**. Its publisher abstract describes similar
transfer-coated layers over supports with different fluxes. Its main article
and original PDF have since been supplied and audited below; unreported data
are not further acquisition requirements. [The review](../literature/mo2022_supplement_review.md)
preserves the earlier request and extraction protocol. Mo's Word supplement
does not need resubmission. Simulation code, locks, original OpenLB files and
the previously failed Biesheuvel prediction screen remain unchanged.

## Matched Zhu M1 hydraulic series — 8 November 2026

The supplied full Markdown/metadata and all ten scientific figures were read.
No supplement was identified; the main article suffices for this bounded study.
The direct publisher check returned 403, so absence of a supplement is not proved.
All 21 Figure 6 M1 black diamonds were digitized with raw pixels/axis anchors
retained and a numbered source overlay inspected. Four touching pairs were
deblended explicitly. Source hashes and deterministic reproduction are saved.

A frozen retrospective protocol fits a positive effective non-support
impedance in Julia using 11 odd-rank points; ten even-rank points remain held
out. The support coefficient is fixed at one. Nominal held-out mean/max flux
errors are **2.44%/4.53%**, passing the 10% engineering screen. Conservative
marker/axis reading corners reach **10.78%**, failing that screen; the lowest
training point has an **11.04%** residual. These limits are retained rather
than changing the split/model. The fitted impedance is **2.51170e8 Pa s/m**,
corresponding to **1433.29 LMH** limiting flux at 1 bar. No unreported viscosity
or temperature is assigned and no intrinsic support k is inferred.

The conditional fixed-100-LMH energy calculation places 90% non-support
dominance at a support flux of **12.90 m/h** at 1 bar. It estimates **30.47%**
endpoint energy savings, compared with **42.86%** translated from the paper's
reported M1 flux gain. This disagreement and mechanical constraints prevent
an experimentally validated energy-performance claim. No new native OpenLB,
DWSIM case or UI is claimed. The earlier failed Biesheuvel result is preserved.

Final direct Julia tests pass **347 assertions** and the complete Python suite
passes **19 tests**, including fit/oracle agreement, held-out isolation,
digitization reproduction, impedance-versus-μR equivalence, reading sensitivity,
energy/threshold checks and export replay/tamper rejection. [The report](zhu_matched_validation.md)
and [commands](../calibration/README.md) explain reproduction and limits.
Dependency locks/runtime/bridge configuration and original OpenLB files remain
unchanged; the Julia package adds only the independent measured-impedance API.

The later original-PDF refinement is recorded separately. Next, preserve the
effective impedance range as a scale-translation reference and specify a
controlled synthetic graded-support stock-component OpenLB case. Use the papers
as available, declare scenario assumptions, and verify grid/domain/pressure
behaviour. No raw-data acquisition or actual-specimen reconstruction is required.

## Original Zhu PDF and clarified computational scope — 8 November 2026

The original nine-page PDF has now been inspected at Figure 6 and its full
relevant page. Its native 890 × 906 image is larger than the prior 623 × 624
conversion. A separate, reproducible 21-point reading retains the same model,
11/10 split, pixel bounds and 10% screen. Nominal held-out mean/max errors are
**2.66%/4.56%**, while conservative reading corners reach **7.88%**: both screens
pass. Fitted effective impedance is **2.50947e8 Pa s/m**. The lowest-flux
training residual is now **13.56%**, and all nominal held-out predictions remain
low. These limitations and the earlier JPEG/Biesheuvel failures are preserved.
[The PDF report](zhu_pdf_refinement.md) and [reference](zhu_pdf_reference.json)
record source hashes, saved outputs and exact reproduction commands.

The complete Python regression suite passes **20/20 tests** after adding a
PDF-reading Julia-oracle/source/replay check. Both PDF-v2 and JPEG-v1 saved runs
pass integrity and numerical replay, and both digitizations regenerate exactly.
The Julia numerical implementation is unchanged from the preceding 347-assertion
direct test run; no new full direct Julia suite is claimed. No dependencies,
locks, runtime/bridge setup, original OpenLB files, DWSIM cases or UI changed.
No new native simulation or browser check was needed for this reading update.

The user clarified that this is a computational-tool project based on papers
and available supporting information. [The accepted scope](verification_scope.md)
requires numerical checks and bounded literature comparisons, without new
experiments, inaccessible raw-data acquisition or exact-specimen 3D pore
reconstruction. Unreported conditions/uncertainty stay unknown or become
explicit scenario assumptions where a simulator needs inputs. A synthetic
graded-support comparison can now proceed through the existing OpenLB → Julia
→ saved-data → DWSIM path without another experimental dataset.

## Fixed-opening interior-widening study — 8 November 2026

Three synthetic pores retain both 4-µm mouths/terminal necks, a 40-µm body
and a 20-µm representative cell. Only the interior radius changes to 4, 5 or
6 µm. All 24 actual native measurements pass analytical baseline, grid,
conservation, pressure/reversal, domain and time-step checks. Fine-grid
resistance changes are at most 0.026%; baseline analytical permeability and
profile errors are 0.0820% and 0.0780%. Three finest-grid material exports
match the declared profile on every axial plane and each form one connected
fluid component. All 2,070 checked original OpenLB files match the archive.

Julia carries measured body-resistance ratios through a single common
baseline calibration to Zhu's 3.12 m/h support response, retaining fixed
effective non-support impedance. Actual DWSIM energy savings are **7.96%**
and **11.30%**, versus native resistance reductions **25.29%** and **35.88%**.
Reading-bound scenarios remain separate from experimental uncertainty. The
non-support share rises from 68.50% to 77.23%; 90% is not reached in this family.
No unique experimental structure, intrinsic selective resistance or total-plant
energy performance is inferred.

The final run `20261108T222740Z-26e9f37b` retains 24 hash-verified native
measurements and nine freshly solved/saved/reloaded DWSIM design cases, plus a
property anchor. Replay verifies 270 artifacts and recomputes numerical analysis
and Julia loads without new external solves. An altered full-run CSV is rejected.
The independent Julia suite passes **372 assertions**, the complete Python
suite **25 tests**, and **18** native input/nonconvergence guards pass. Geometry
and comparison figures were visually inspected; SI VTP files reload with correct
bounds, shapes, float64 fields and finite pressure. No new browser interface
was built or checked for this CLI study.

Locks, dependencies, bridge and DWSIM adapter are preserved. The existing
native launch helper gains only an optional timeout; earlier defaults remain
600 seconds and this study requests 1,800. The 647-second widest solve justified
that execution allowance without changing numerical thresholds. Interrupted
pilots/partial runs and earlier failed literature screens remain preserved.
[Results and limits](graded_support_study.md), [compact reference](graded_support_reference.json)
and [commands](../graded_support/README.md) describe exact reproduction.

The subsequent [explicit hydraulic filtration milestone](process_filtration.md)
is completed below; circulation/recycle duty remains the following extension.

## Explicit hydraulic filtration — 9 November 2026

The [target-version interface investigation](dwsim_filtration_interfaces.md)
selected a stock pump, splitter, upstream loss valve, permeate pressure-loss
valve and retentate letdown valve. Julia supplies the same calibrated native
resistance-derived TMP. DWSIM solves water properties, flashed stream states
and pump energy. The available CrossflowUF model prescribes flux and does not
calculate our pressure–flux closure; no new custom membrane unit was created.

Final run `20261108T231207Z-316d61cb` stores 18 actual solved/saved/reloaded
flowsheets across three designs, recoveries 25/50/75% and prescribed upstream
losses 0/20 kPa. Replay checks 99 local artifacts plus the unchanged upstream
270-artifact study, reanalyzes 24 native cases and audits three geometries.
All seven material streams, twelve directed connectors and every unit/overall
balance pass. Max relative mass residual is zero at exported precision and
max enthalpy residual 2.64e-5 of shaft power versus a 1e-3 limit. Viscosity and
density departures remain below the declared 0.1% screens.

At 50% recovery the full fresh-feed rate is twice the permeate rate; its shaft
energy doubles relative to the earlier permeate-only pump. The assumed 20 kPa
upstream loss reduces widest-design savings from 11.30% to 3.81%. Actual DWSIM
renders and scientific energy figures were visually inspected. Stock rounded
energy labels are supplemented by exact W annotations and JSON values.

The direct Julia suite passes **422 assertions**, the full Python suite
**33 tests**. New tests call actual Julia and replay the actual saved DWSIM
record; changed flows/enthalpies/properties, artifact bytes and splitter ports
are rejected. These saved-record tests are distinguished from the fresh
external solves performed by the study command. Numerical/load replay does
not rerun OpenLB or DWSIM. Cold Julia initialization and DWSIM setup/solve
timings are recorded separately, without speedup claims.

One pilot failed to compile an optional annotation because of an incorrect
overload; the installed signature was used and a fresh complete run passed.
Intermediate completed runs are retained; the final run also measures cold
Julia initialization before upstream replay and checks derived-output overflow.
No numerical thresholds were
relaxed. Original simulations, installed distributions, prior failures and
unrelated files are preserved.

[Commands](../process_filtration/README.md), [results/model limits](process_filtration.md)
and [the compact reference](process_filtration_reference.json) give reproduction.
That subsequent pressure-consistent mixer/recycle/circulation-pump case is now
complete, as recorded below. No new experimental dataset was required.

## Crossflow recycle and circulation energy — 9 November 2026

Final run `20261108T233356Z-b884df5e` retains 25 actual solved/saved/reloaded
DWSIM recycle cases, plus a deliberately rejected one-iteration solve.
The 21 design cases vary three pore interiors, module/permeate flow ratio,
overall water recovery and independent channel resistance. Four variants test
initial guesses and tenfold tighter tolerances. Original OpenLB measurements,
the common literature scaling and effective non-support impedance are unchanged.

The installed generic recycle interface and official example were inspected
before selecting stock mixer, cooler, recycle, two pumps, two splitters and
three valves. Every case checks 13 material streams, three energy streams,
all ten unit balances plus overall mass/enthalpy, and 26 saved connections.
Convergence takes 16–185 direct-substitution iterations; the nominal case takes
48. Different seeds/tolerances reproduce work and heat duties within 1e-8
relative. Maximum relative mass and energy residuals are 2.12e-11 and 2.45e-5.

Early pilots exposed an approximately 1.9e-8 K pure-water tear re-flash residual.
The documented stock TP override on the temperature-controlled tear streams
resolves it with the original 1e-8 K tolerance retained. The diagnostic pilot
logs and failed states remain under `.local/crossflow_recycle/`; stock code and
convergence algorithms are unchanged. A compile-time namespace ambiguity was
resolved with an explicit FluentAPI alias. Two intended overflow tests were
corrected to use inputs that actually overflow; no scientific screen changed.

At 50% recovery, module/permeate ratio 5 and nominal channel resistance, actual
feed/circulation powers are 0.753768/3.703687 W, with 2.454001 W removed as heat.
Circulation contributes 83.09% of shaft work. The widest support saves 1.91%
of total shaft energy, falling to 0.25% at ratio 15. Cooling heat remains
separate from electricity. This is a lumped clean-water hydraulic cost study;
it does not infer optimal circulation, resolved axial pressure/shear, solute
rejection, fouling or complete plant utility demand.

The checker verifies 138 local artifacts, replays real Julia loads and can
reverify the earlier 18 filtration cases and original 24 native cases. It
performs no new external solves. Direct Julia tests pass **518 assertions**,
including 96 recycle checks; the full Python suite passes **43 tests**.
Altered mass, density, cooling, convergence flags, connector ports and artifact
bytes are rejected. Actual DWSIM renders and the scientific energy figure were
visually inspected. All 2,070 checked original OpenLB files still match the
pinned archive, with the pre-existing config exception recorded separately.
Implementation/runtime hashes match the accepted run. No dependencies, locks,
original simulator sources or unrelated `assets/` were changed.

[Results and assumptions](crossflow_recycle.md), [compact reference](crossflow_recycle_reference.json)
and [commands](../crossflow_recycle/README.md) document reproduction. Next is an
application-specific published channel geometry and transport/fouling relation
that gives circulation a performance benefit to compare against its cost.

## Literature-constrained crossflow benchmark — 9 November 2026

Run `20261109T005324Z-230bef0c` completes the bounded source benchmark and
conditional process comparison: **32 actual solved/saved/reloaded DWSIM cases**,
**846 passing Julia assertions** and **50 passing Python tests**, including
the existing VTK tests. Ancestor evidence is reverified and the deliberately
one-iteration case is rejected. The checker covers 174 artifacts and all
32 saved archive topologies without launching external solvers.

At fixed 100 LMH, the table-based policy selects 0.9 m/s among tested points;
the separate monotonic figure-bound policy selects 1.4 m/s. The unresolved
source flow convention is represented by literal and factor-seven scenarios.
The widest support saves **7.11–10.99%** of total shaft energy at those selected
points. These are conditional water-carrier scenarios, not absolute module
calibration or a continuous optimum. The original non-support impedance and
OpenLB support measurements remain fixed.

Julia also balances a passive retained solute through the finite-pass loop.
Source module-inlet concentration stays fixed by adjusting fresh composition;
DWSIM still contains Water only. High circulation requires 88–2511 stock
direct-substitution iterations. Our adapter's input ceiling was raised to
5000, with its old default and stock algorithm preserved. Seed/tolerance
checks pass; no numerical acceptance screen was relaxed.

The host preflight and locked environment synchronization pass with the same
Python/Julia/.NET/DWSIM installations and unchanged dependency manifests.
The initial sandbox-only Julia import failure came from compilation-cache
permissions; actual host tests and imports pass. The scientific figure was
visually inspected. No new OpenLB solve or original OpenLB source change was
required. Exact commands, results, source limits and the next explicit
star-channel hydraulic reference are documented in
[the report](crossflow_benchmark.md) and
[launch instructions](../crossflow_benchmark/README.md).
