# N160 wall tracer in the disordered packing

The [N128 study](../disordered_n128/REPORT.md) exceeded the 2-percentage-point
grid limit in three of four matched comparisons. This study repeats the
eight-section test with the same disordered packing at N160.

The sphere-centre file, 32 spheres, axial inlet cut and dry carbon basis
of 7.032455155 g come from N128. The two particle-matched cases use the
same target Péclet numbers, 24,985.45 and 14,054.32, using the accepted
carrier speed for the diffusivity match.
[Protocol](protocol.md)

The external stock OpenLB 1.9r0 TRT application uses D3Q27,
Kupershtokh forcing, BounceBack walls and sphere indicators.
The 1280 × 160 × 160 grid is initialized from the archived N128 field,
with the same lattice time step to voxel-width ratio. The initial pressure
gradient extrapolates the accepted N96 and N128 gradients in inverse
resolution. The native acceptance criteria are carried over from N128.

The conservative face projection supplies Julia with locally balanced
rates and preserves the native directional throughputs. The eight-section
tracker uses 8,192 flux-weighted markers, the same inlet seeds as N128,
and 80 and 40 µs time steps. A solid-wall encounter counts as capture.
The grid test takes the largest absolute change over sections 3–8 in
total wall-free survival, conditional capture among survivors and
initially fast wall-free survival.
