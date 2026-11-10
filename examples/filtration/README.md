# Filtration hydraulics

This completed hydraulic study is closed for now. Each stage has a scientific
report covering its purpose, theory/equations, results, discussion and conclusions,
and a separate run guide. The
[final operating-map report](operating_map/REPORT.md) answers the support-design
and total-pumping-energy question at the closing process conditions.

| Stage | Scientific report | Run guide |
| --- | --- | --- |
| Layers-in-series demonstration | [Support benefit and dominance](layered_support/REPORT.md) | [Commands](layered_support/README.md) |
| Published support/membrane calibration | [Calibration and held-out results](calibration/REPORT.md) | [Commands](calibration/README.md) |
| Fixed-opening support geometry | [Native resistance and scale translation](graded_support/REPORT.md) | [Commands](graded_support/README.md) |
| Once-through circuit | [Recovery and upstream losses](once_through/REPORT.md) | [Commands](once_through/README.md) |
| Crossflow recycle | [Feed versus circulation duties](crossflow_recycle/REPORT.md) | [Commands](crossflow_recycle/README.md) |
| Literature capacity policies | [Capacity and source-flow ambiguity](literature_benchmark/REPORT.md) | [Commands](literature_benchmark/README.md) |
| Explicit low-Re channel | [Distributed hydraulics and COD scenario](star_channel/REPORT.md) | [Commands](star_channel/README.md) |
| Final higher-Re map | [Process comparison and conclusions](operating_map/REPORT.md) | [Commands](operating_map/README.md) |

Follow [ordered reproduction](../../docs/reproduction.md) through the ancestor
studies before launching a later stage. Existing immutable runs remain under
their original `results/` directories. Historical BGK evidence and subsequent
TRT verification are labelled separately in the reports. Source measurements,
reading sensitivities and model assumptions are not conflated.

Detailed execution/audit records and compact numerical references remain in
`docs/`; the scientific reports link to them without requiring them for the
basic explanation of each example.
