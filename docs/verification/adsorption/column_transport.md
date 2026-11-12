# Column transport and process qualification, 10 November 2026

The maintained [report](../../../examples/adsorption/column_transport/REPORT.md)
and [commands/contracts](../../../examples/adsorption/column_transport/README.md)
separate accepted numerical components from failed application qualification.

## Native scalar evidence

- Formal 19-case run: `20261110T151006Z-1aaf231f`.
- Additional three-case preservation diagnosis: `20261110T153541Z-923381eb`.
- Stock free-flow zero-gradient and Robin exit controls pass; no original
  OpenLB source was modified. The legacy convection setter was not patched.
- Application uniform tracer fails all four tested stock collision choices.
  Combined regularized TRT is best; smaller dt improves its error but the
  declared preservation/bounds screens still fail. The N48 follow-up also fails.
- The saved-native checker verifies data and recomputes the Julia reference;
  this does not promote failed scientific screens. It runs no native simulation.

## Numerical components and preservation

Direct Julia boundary-reference tests cover independent inlet/outlet moments
and a conservative uniform profile. The unchanged particle suite is rerun;
its batch qualification and previous failed column comparisons remain intact.
New Python tests call real Julia in nonadsorbing limits, exercise a hypothetical
qualified packet for the frozen-parameter contract, and reject the actual failed
native packet. Hypothetical contract fixtures are not native qualification.

`results/adsorption/column_transport/preservation.json` verifies the unchanged
source guards of the synthetic process (26), batch particle (20), earlier column
contact (30) and absorbing tracer (38) studies. Vendor checks compare 2,070
original distribution files to the pinned archive; the documented pre-existing
config exception remains unchanged. Existing dependency locks and unrelated
branding files are preserved.

## Native hydraulic execution and resumption

The fixed-flow branch uses source flow/diameters and a separately declared FCC
packing/inventory scenario. Its native numerical kernel is copied unchanged
from the established application. N32/N48 and half-dt cases independently solve
the pressure needed for the common flow. Actual DWSIM pumping is separate from
the unqualified adsorption/service branch.

Serial coordinator `20261110T150020Z-bac62eb7` preserves its completed native
cases and carrier-property anchor. Resume attempt `20261110T152946Z-e510f586`
had an argv-path bookkeeping error and produced no new native case. The next
resume `20261110T153300Z-d6fd3c16` preserved additional native cases; an import
helper name was corrected. Completed source/argv/output bytes were retained.
The active corrected study is `20261110T154704Z-d14d13ed`; copied original
coordinators and input hashes remain archived. Resumption changes scheduling,
not dynamics, geometry, discretization, physical conditions or acceptance limits.
Stopped coordinators were removed after their native children finished.

The final Julia-to-DWSIM handoff is governed by the separate
[process protocol](../../../examples/adsorption/column_transport/process_protocol.md):
Julia must calculate the load actually supplied to DWSIM, with a matched-Reynolds
guard. It does not count merely agreeing with the earlier Python arithmetic as
that connection. It leaves service time/volume unavailable while native solute
transport is unqualified. The completed process run is
`20261110T163606Z-3e5ef256`: both actual Julia-supplied DWSIM loads pass
solve/save/reload and work/mass/enthalpy checks. Native pressure drops are
308.809/958.188 Pa for .80/.45 mm; actual shaft duties are .300574/.932636 mW,
or .114374/.354884 Wh per m3 pumped in the friction-only scenario. The ratio is
3.10285. No usable-capacity tradeoff is inferred.

The final hydraulic study contains 18 actual native executions (12 reused
completed cases plus six new in its final coordinator). Two additional native
executions qualify cold initialization and a perturbed two-cell periodic domain.
The diagnostic wrapper was subsequently made portable: identical saved native
inputs, binary and data were hash-checked/reanalyzed, with zero new native runs.
The current diagnostic reference points to `20261110T163145Z-691455f4`, retaining
its actual native source `20261110T161859Z-d67ce326`. In total this milestone
executes 20 hydraulic and 22 formal/diagnostic scalar cases. Earlier exploratory
scalar probes are retained separately as unqualified diagnostics.

Three DWSIM cases belong to hydraulic verification (including the carrier anchor);
two additional actual cases demonstrate Julia supplying the pump load. These
are five executions, not five independent designs. The saved checkers execute
zero new native/DWSIM simulations. They verify 1,790 scalar artifacts/29 sources,
432 sensitivity artifacts/24 sources, 2,483 hydraulic artifacts/28 sources,
362 cold/domain artifacts/23 sources and 47 final process artifacts/36 sources.

Final host-graphics Python regression: **174 passed** in 63.77 s. Direct Julia:
26 boundary-moment + five profile + ten operating-point assertions pass; the
unchanged particle component passes 51 assertions. Repository checks pass for
18 examples and 13 paper folders. Figure QA is separate from browser acceptance.

## Environment and scope

Actual machine preflight is saved at
`results/adsorption/column_transport/preflight.json`: macOS 26.6.2 ARM64,
Python 3.12.9, Julia 1.11.6, JuliaCall/PythonCall 0.9.37, uv 0.13.7, NumPy 2.5.4,
Panel 1.10.0, PyVista 0.49.2 and Matplotlib 3.11.2. Native GCC is 15.2.0;
the existing private .NET SDK is 10.0.403. Actual DWSIM assembly version is
10.2.15.0. No new runtime installation occurred.

No new laboratory evidence, reconstructed specimen, experimentally measured
fluid dispersion, interactive browser pass, reactive model or non-isothermal
adsorption model is claimed. The current blocker is numerical transport
qualification, not access to more experimental data.
