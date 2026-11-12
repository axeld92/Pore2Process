# Actual Julia-to-DWSIM operating-point handoff

Declared before this additional process execution, 10 November 2026. Consume
only the completed native fixed-flow study whose grid/time/Mach/density/mass
and stationarity screens passed. Julia validates Reynolds matching (relative
flow mismatch <=1e-4), translates the measured operating-point resistance into
the declared sample load and returns SI Q, deltaP and reference work. It reuses
the existing independent Julia pressure-load function. The apparent K is not
used to extrapolate to another Reynolds number. Actual DWSIM pump inputs must
come from those Julia outputs. Solve, save, reload and re-solve both sizes.

This is a bed-friction circuit with zero elevation and zero ancillary pressure
load. It is not reconstruction of the source downflow apparatus. Use source
43.8 mL/min and diameters; synthetic FCC .6 porosity; full source-density-derived
dry mass; actual water carrier properties at 295.15 K; chosen pump efficiency
.75. Retain all native and Julia input snapshots and outputs. DWSIM work must
agree with Julia within 1e-8 relative, mass within 1e-10 and enthalpy within .001
relative to shaft work. Keep source, scope and simulator execution counts.

The higher-Re operating point must also pass the separately declared native
cold-initial-state and perturbed two-cell-domain checks before this handoff.

Optional algebraic sensitivities use efficiencies .6/.75/.85 and common extra
pressure duties 0/1000/10000 Pa. These are explicitly chosen process scenarios,
not uncertainty bounds and not additional DWSIM executions. They use the pump
work law verified by the actual base cases. Current failed scalar evidence must
be rejected; no new native-derived adsorption/service prediction is produced.
