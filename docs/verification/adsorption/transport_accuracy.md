# Engineering transport accuracy — 10 November 2026

The maintained [scientific report](../../../examples/adsorption/transport_accuracy/REPORT.md)
and [run guide](../../../examples/adsorption/transport_accuracy/README.md) distinguish
the accepted uniform drift from the larger nonuniform transport failures.

Actual native executions:

| Study | Native run ID | Cases |
| --- | --- | ---: |
| Sixteen-transit uniform and sharp front | `20261110T180254Z-e4d0abff` | 4 |
| Two known smooth-input durations | `20261110T180944Z-096133a1` | 4 |
| Short stock regularized-inlet diagnosis | `20261110T181253Z-eecee644` | 2 |
| Periodic smooth-profile isolation | `20261110T182020Z-41c6cae7` | 4 |

These are 14 actual OpenLB executions, not independent experimental material
designs. Native fields and the existing DWSIM carrier properties are reused;
there are no new hydraulic, DWSIM, adsorption-fit or browser runs.

Both uniform cases pass the new engineering screen. Sharp fronts produce
−.3565..1.3370 and −.3402..1.3280 kg/m3 for .80/.45 mm particles with 1 kg/m3
feed. Gradual inputs and a stock regularized inlet do not resolve those bounds.
Known inlet mean/variance are independently audited by Julia. In periodic
controls, the free-flow profile errors are .2215/.2217% of initial amplitude;
porous profiles depart substantially from the initial .25–.75 bounds despite
having no inlet/outlet. No qualified dispersion or service result is produced.

Saved checks bind native/analysis source snapshots, actual requests, build and
artifact hashes. The three front/uniform assessments are replayed separately
from the periodic-profile checker. Their checks run zero new native or DWSIM
simulations. The old column_transport checker was also replayed successfully,
preserving its strict failures, passed hydraulics and source guards. OpenLB's
archive verifier still checks 2,070 original files. The pre-existing vendor
config exclusion is unchanged; no original OpenLB file was edited.

Environment: macOS 26.6.2 ARM64, Python 3.12.9, Julia 1.11.6,
JuliaCall/PythonCall .9.37, GCC 15.2, pinned OpenLB 1.9r0. No installation or
dependency lock changed. Exact inventory, simulator and test logs are retained
under `results/adsorption/transport_accuracy/`.

The full host-graphics Python suite passes 180 tests in 56.34 s. Repository
checks cover 19 examples, 13 paper folders and 161 Markdown files with no
errors. All three saved-data figures were visually inspected. Their exact
PNG/source hashes are in `figure_sources.json`. No new browser was exercised.

The report records a timing-convention deviation: finite sharp fronts run for
15.899 actual void times because the nominal runner uses full-cell porosity
when excluding reservoir planes. Periodic uniform windows are 16.000. Raw
durations and inventories remain unchanged; the deviation does not affect the
observed large concentration-bound failures.

## Bounded follow-up and stopping decision

Native run `20261110T183934Z-2a431270` adds four actual executions, reusing the
unchanged profile binary and original OpenLB sources. A no-advection control
passes; N32 half/tenth dt and N48 at the same dt as N32 half-step fail the
concentration-range screen. The [diagnosis](../../../examples/adsorption/transport_accuracy/DIAGNOSIS.md)
records results, parameter/coordinate/velocity audits, source references and the
decision to defer the high-Peclet branch. No unique root cause or correction is
claimed. No new Julia adsorption or DWSIM case was executed.

The saved checker verifies 27 source files and 573 artifacts and independently
replays the field audit, with zero new native/DWSIM runs. Ten targeted tests
pass in 3.09 s; the prior full-suite 180-test run is separate evidence and was
not repeated for this isolated diagnosis. The new saved-data figure was visually
inspected. The original vendor archive check still covers 2,070 files.
