# Native tracer/contact verification - 10 November 2026

Maintained explanation: [scientific report](../../../examples/adsorption/packed_bed_tracer/REPORT.md).
Commands: [run guide](../../../examples/adsorption/packed_bed_tracer/README.md).

- Native run: `20261110T131217Z-7d3c3c65`; one new steady FCC flow and eight
  passive-scalar cases. All active collisions use stock TRT, with stock D3Q27,
  geometry, bounce-back and scalar Dirichlet boundaries. All 2,070 checked vendor
  files remain unchanged, excluding the documented pre-existing config exception.
- N=16/24/32, half dt, doubled length and doubled diffusivity are checked. The
  selected regime uses Dm=2e-6 m2/s, scalar magic=0.3 and unequal TRT times.
  Uniform-flow analytical, concentration bounds, kinetic mass and tail checks
  pass. Final grid/time moment changes pass the predeclared screens.
- Native saved-data checker: 8 scalar cases, 1,099 artifacts, 38 source files,
  reconstructed inventories/link transfers and two Julia diagnostic predictions.
  It launches no new native simulation.
- Selected translation: `20261110T132745Z-1528e655`. One effective coefficient
  comes from the short native tracer mean, with measured Q, physical pore volume
  and therefore v=QL/V fixed. Longer-segment mean/variance/absolute-flux errors
  are 0.0668%/0.979%/0.512%; normalized-profile RMSE is 0.000339.
- Conservative checker independently rederives D, replays six Julia calculations,
  checks 11 artifacts and 5 source files, and verifies physical flow, volume and
  mass identities. No new native/process solve is run by this stage.
- Direct Julia: 36 tracer/reference assertions and 51 unchanged particle
  assertions pass. Full host-graphics Python suite: **150 passed**. Figures were
  visually inspected. Exact versions and source/input hashes are saved.
- The earlier failed native regimes remain retained, including upper-bound and
  tail failures. A historical diagnostic transcription is corrected explicitly
  in the report; historical snapshots are unchanged. No threshold was relaxed.
- No new DWSIM solve or browser interface is claimed. The earlier adsorption
  comparisons and completed hydraulic/process examples are preserved. This is
  synthetic absorbing-end contact qualification, not an experimental column or
  realistic-Peclet adsorption prediction.
