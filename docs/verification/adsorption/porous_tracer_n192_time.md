# N192 time-step refinement, 11 November 2026

The [report](../../../examples/adsorption/porous_tracer_refinement/N192_TIME.md)
compares the exact archived N192 tracer mesh and face flows at half
the original time step. Final concentration variance increased 0.1255%
and fitted wave phase shifted by +4.15e-6 rad. Final profile RMS change
was 0.217% of initial amplitude, above the 0.1% time-step criterion.
The estimated numerical share of final variance loss fell to 5.42%.
The study generated one Julia tracer run and zero new mesh, OpenLB,
full-bed or process runs. Read-only replay checked 26 sources and
33 artifacts.
