# Stationary-wall time-step refinement, 11 November 2026

The [completed comparison](../../../examples/adsorption/stationary_wall_time_step/REPORT.md)
halved the time step of the stock regularized TRT half-way setup on the
same 128-cell disk and preserved the physical flow, diffusivity and duration.
Maximum uniform departure fell from 0.06714 to 0.04936 kg/m³. Varying-tracer
RMS across the liquid increased from 4.83% to 6.74% at the midpoint and
from 11.82% to 12.80% at the final time against the same Julia reference.

Both tracer tests failed their accuracy limits. The exact-midpoint repeat
reproduced the full-duration varying case's history through that step. The
cohort `20261111T144732Z-c53cd814` contains three disk runs and zero bed
or process runs. Read-only replay checked 20 source hashes and 88 artifacts.
[Run guide](../../../examples/adsorption/stationary_wall_time_step/README.md)
