# Compact evidence

`transport_reference.json` records the executed 19-case native study and its
failed application qualification, separate from artifact/source verification.
The full arrays, native logs and source snapshots stay under ignored results.
Unavailable dispersion is JSON null and the failed packet is rejected by the
particle-model bridge.

`hydraulic_reference.json` records 18 actual fixed-flow native executions and
three actual DWSIM cases. `hydraulic_diagnostic_reference.json` records two
additional actual cold/domain runs; portable-wrapper replay is counted separately.
`sensitivity_reference.json` retains three additional failed scalar diagnoses.
`process_reference.json` records the selected actual OpenLB → Julia → DWSIM
handoff and two additional DWSIM executions, retaining unavailable service
quantities. Read [the report](../REPORT.md).
