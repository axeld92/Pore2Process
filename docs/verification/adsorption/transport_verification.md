# Independent transport checks, 11 November 2026

The [report](../../../examples/adsorption/transport_verification/REPORT.md)
compares full three-dimensional stretching with a qualified continuum
reference, then tests one circular wall with an analytical reference.
The [run guide](../../../examples/adsorption/transport_verification/README.md)
records commands, environment, native binary data and interpretation.

- Stretching: `20261111T004319Z-61b04382`, five native executions. N96 passes
  midpoint RMS (0.474%) and fails final RMS (7.46%) and peak excursion
  (0.023213 kg/m³). Both controls pass their applicable screens.
- Curved wall: `20261111T011153Z-1a03a0b9`, eight native executions. Fine
  no-flow diffusion passes; fine rotation fails with 47.1% RMS in the
  outer shell. Under matched rotation, uniform initial concentration of 0.5
  gives a whole-window range
  of 0.103112–0.979615 kg/m³. Four/full-depth disk fields agree exactly.
- Four independent reference tests pass. Read-only replay verifies 26 source
  files for each stage, 122 stretching artifacts and 188 curved artifacts.
  Both use the same archived stock TRT executable and pass the 2,070-file
  vendor integrity check. No full-bed, adsorption or DWSIM run was launched.
