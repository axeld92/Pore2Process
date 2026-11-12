# Stationary curved wall, 11 November 2026

The [completed report](../../../examples/adsorption/stationary_curved_wall/REPORT.md)
compares stock OpenLB TRT with a refined Julia radial continuum reference.
All seven native cases in `20261111T122659Z-56628d90` completed. The
first reference-only attempt `20261111T121850Z-40561f07` missed its spatial
screen; it launched no native cases. One additional radial doubling passed
without changing tolerance or physics. Eight Julia reference assertions passed.

At N128 the zero-wall-speed case reduces the largest uniform drift from
0.479615 to 0.065761 kg/m³ at the same peak velocity as archived rotation.
The varying case still gives 12.29% final RMS error globally and 13.10%
near the wall, with a 0.086006 kg/m³ transient excursion. Closed inventory,
depth equivalence and vendor integrity pass. Read-only replay checked 30
sources and 187 artifacts. No new full bed or process prediction was run.
