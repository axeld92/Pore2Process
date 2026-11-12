# N128 wall tracer in the disordered packing

The [N96 study](../disordered_refinement/REPORT.md) found that the
N64-to-N96 profile differences exceeded two percentage points for
both particle-matched conditions. This refinement repeats the
eight-section test on the same disordered geometry at N128.

The 32 sphere centres, 0.8 mm sphere diameter, axial origin,
11.113 mm length and transverse dimensions come from the archived
N96 input. The continuous sphere volume is 8.57864234 × 10⁻⁹ m³,
with 7.032455155 g as the column-scale dry carbon basis. The
0.80 and 0.45 mm conditions match Péclet numbers 24,985.45 and
14,054.32 using the accepted mean interstitial velocity.
[Protocol](protocol.md)

The external OpenLB application reuses the stock 1.9r0 D3Q27
TRT collision, Kupershtokh forcing, BounceBack and sphere
indicators. A 1024 × 128 × 128 lattice resolves the same full
disordered domain. The interpolated N96 velocity and density
initialize the N128 solve. The relaxation offset preserves the
N96 ratio of lattice time step to voxel width, and the starting
driving gradient extrapolates the two accepted grid values
linearly in inverse resolution. The new lattice evolves to the
same native acceptance criteria.
[OpenLB 1.9 release](https://www.openlb.net/download/)

The Julia face projection preserves the native directional
throughputs and balances each pore voxel. The eight-section
perfect-sink wall tracer uses 8,192 flux-weighted markers at
each of 80 and 40 µs, with the same inlet seeds as N96. The
initially fast cohort is the 90% with the shortest first-section
inert transit times. The grid test uses the largest absolute
difference over sections 3–8 in cumulative wall-free survival,
conditional capture, and initially fast wall-free survival.
