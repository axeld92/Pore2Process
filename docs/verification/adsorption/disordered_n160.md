# N160 disordered-packing tracer refinement

Cohort `20261113T130058Z-5ba9e799` reused the N128 packing,
inlet cut and Péclet targets. The stock OpenLB 1.9r0 TRT
flow and conservative handoff passed. All four runs recovered
the full inert dose and passed the time-step checks. Maximum
N128-to-N160 differences were 0.98 and 0.79 points at 80 µs
and 0.76 and 0.71 points at 40 µs for 0.80 and 0.45 mm,
respectively. All four comparisons meet the 2-point grid limit.
[Report and reproduction](../../../examples/adsorption/disordered_n160/REPORT.md)
