# Flow-derived path exchange and the complete two-cell tracer curve

Cohort `20261112T073700Z-e27daca0` used the same saved N96 TRT flow
to estimate path exchange from Julia marker pairs, preserving the
previous one-cell marginal calibration. The predictions were frozen
before scoring the existing complete two-cell curve. The primary model
passes all predeclared screens, with 0.61% error at 99% recovery and
2.67% in variance. The rise-duration error is 1.08%. The matched
advection-only pairing still fails 99% recovery and variance. Half-step,
independent-seed and raw one-cell checks pass. Read-only replay verified
26 sources and 45 artifacts. Adsorption remains on hold.
[Report and data contract](../../../examples/adsorption/path_exchange/REPORT.md)
