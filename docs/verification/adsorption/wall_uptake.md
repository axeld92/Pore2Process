# Path-specific film transfer from wall uptake

Cohort `20261112T143948Z-126af769` reused the saved TRT carrier flow
and the Péclet-matched diffusivities and exchange matrices. It preserved
7.032455155 g of shared dry carbon in both column cases. The 80 and
40 µs perfect-sink marker tests passed their time-step checks. The
constant-film reduction overpredicted two-period uptake by 1.09 and
1.42 percentage points at 80 µs for the 0.80 and 0.45 mm particles.
Path-level second-period errors remained above their limit at both
time steps.

The frozen column predictions passed numerical refinement and failed
both Figure 7 comparisons. The 5% crossing occurred within five seconds
for both sizes. The 50% crossings were 5.84 s and 122.61 min against
30.75 min and 55.36 min, respectively. Read-only replay verified 36
sources and 59 artifacts. [Report and reproduction](../../../examples/adsorption/wall_uptake/REPORT.md)
