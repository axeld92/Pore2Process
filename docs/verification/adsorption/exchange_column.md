# Unfitted exchange-column comparison with Alhamed Figure 7

Cohort `20261112T080827Z-1b52e000` coupled the 40-class synthetic tracer
exchange model to the batch-qualified spherical particle model. Both
diameters used the source geometry and reported voidage, the dry mass
derived from printed bulk density, and the same unfitted film correlation.
Six Julia forward runs were frozen before reading the Figure 7 markers
for scoring. Both diameters passed numerical refinement and conservation.

The exchange model predicted midpoints of 78.84 and 105.83 min for
0.80 and 0.45 mm, versus digitized values of 30.75 and 55.36 min.
Both failed the early, midpoint and curve-error checks. Matched
controls with zero exchange also failed. Read-only replay verified
32 sources and 41 artifacts. This comparison qualifies no adsorption
service time or energy prediction.
[Report and reproduction](../../../examples/adsorption/exchange_column/REPORT.md)
