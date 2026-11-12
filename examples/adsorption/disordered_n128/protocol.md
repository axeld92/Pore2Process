# N128 grid test of the disordered eight-section packing

Reuse the archived N96 sphere centres, carbon mass, inlet cut,
and both particle-matched Péclet numbers. Solve a 1024 × 128 × 128
periodic lattice with the same stock OpenLB 1.9r0 D3Q27 TRT,
Kupershtokh forcing, BounceBack walls and sphere indicators.
The N96 velocity and density initialize the run by periodic
cell-centred interpolation. Scale the relaxation offset so the
initial lattice velocity conversion matches N96. Start from
the driving-gradient trend of the accepted N64 and N96 fields.
Use the unchanged native steady-flow, Mach, density and
velocity-match criteria before accepting the N128 field.

Use the conservative Julia face projection qualified against
the N64 archive and used at N96. Keep the same perfect-sink
eight-section tracker, 8,192 inlet markers and seed, 12 s
observation horizon, and matched 80 and 40 µs tracer time
steps. Require at least 99.9% inert recovery and a voxel
solid-volume difference of at most 0.5% from N96. Match Pe
with the actual accepted N128 mean interstitial velocity.

Compare N96 and N128 at each matched time step. For each
particle-matched condition, require at most 2 percentage
points of maximum absolute difference over sections 3–8 in
total wall-free survival, capture among entering survivors,
and initially fast wall-free survival. Use 80 µs as primary.
Check the 1-point time-step screen separately. A supplemental
marker ensemble, if needed for a borderline result, is kept
separate from the original grid comparison.

The comparison covers tracer transport and external wall
contact. Adsorption and Figure 7 are outside this test.
