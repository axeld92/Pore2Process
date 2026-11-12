# Fixed-flow tracer refinement, 11 November 2026

The [report](../../../examples/adsorption/porous_tracer_refinement/REPORT.md)
records three Julia runs using the archived N96 TRT flow on one FCC cell.
Halving the N96 step changed the final tracer profile by 1.079% of its
initial amplitude, above the 0.1% criterion. Subdividing each voxel into
eight cells preserved the original throughput, and the matched-step
N192 run reduced the estimated numerical share of final variance loss
from 30.27% to 5.63%. Final variance and wave phase still missed the
spatial limits.

Cohort `20261111T160638Z-156b2680` passed read-only replay of 25 sources
and 50 artifacts. It produced zero new OpenLB, full-bed or process runs.
[Run guide](../../../examples/adsorption/porous_tracer_refinement/README.md)
