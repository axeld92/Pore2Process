# One-cell axial dispersion fit and two-cell holdout

The [report](../../../examples/adsorption/dispersion_holdout/REPORT.md)
records cohort `20261112T033234Z-32e27f5e`. A coefficient of
3.24979e-6 m²/s fitted to the N384 one-cell exit variance was frozen
before the two-cell pore calculation. The holdout reused the saved
TRT flow with N288 tracer resolution per repeat and an 80 µs step.

The frozen model predicted 5% arrival 45.42 ms too early. The run
stopped at 205.76 ms after the cumulative-profile RMS error lower
bound exceeded 0.030. The outlet profile failed, so full-window
pore mean, variance and tail metrics are unavailable. The pore mass
residual stayed below 2.40e-15 of the injection. Read-only replay
verified 25 sources and 41 artifacts. No full-bed, OpenLB hydraulic
or DWSIM run was launched.
