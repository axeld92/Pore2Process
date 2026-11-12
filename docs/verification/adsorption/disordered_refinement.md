# N96 disordered-packing tracer refinement

Cohort `20261112T171431Z-193be93e` used the archived sphere
centres and inlet cut, stock OpenLB 1.9r0 D3Q27 TRT and the same
carbon inventory and Péclet numbers. The N96 native and
conservative flow checks passed. All four tracer cases recovered
the full inert dose. The primary maximum N64-to-N96 changes
were 3.93 points for 0.80 mm and 3.21 for 0.45 mm. Both
also exceed the 2-point screen at 40 µs.
The time-step screens passed after expanding the 0.80 mm
ensemble. Read-only replay verified 40 sources and 277 artifacts.
[Report and reproduction](../../../examples/adsorption/disordered_refinement/REPORT.md)
