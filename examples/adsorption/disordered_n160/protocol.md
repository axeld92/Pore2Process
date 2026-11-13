# N160 grid test of the eight-section disordered packing

Reuse the N128 sphere centres, carrier conditions, tracer inlet and
7.032455155 g column-scale dry carbon basis. Solve a 1280 × 160 × 160
periodic lattice with the stock OpenLB 1.9r0 D3Q27 TRT application,
Kupershtokh forcing, BounceBack walls and sphere indicators.
Use periodic cell-centred interpolation of the N128 moments as the
initial condition. Maintain the previous ratio of lattice time step
to voxel width and extrapolate the N96–N128 pressure gradient in
inverse resolution. Apply the original native acceptance criteria.

Use the conservative face projection qualified against the N64
reference and used at N96 and N128. Match both archived Péclet
numbers using the accepted N160 axial interstitial speed. Use the
same eight-section perfect-sink wall tracker with 8,192 markers,
the same inlet seeds, 12 s horizon and both 80 and 40 µs steps.
Require at least 99.9% inert recovery and a voxel solid-volume
difference of at most 0.5% from N128.

Check the maximum absolute N128-to-N160 change over sections 3–8
for total wall-free survival, conditional capture among incoming
survivors and initially fast wall-free survival. Use the 2-point
screen at both matched time steps, with 80 µs as the primary
comparison. Check the 1-point time-step screen separately and
retain any supplemental ensemble as a separate diagnostic.

The milestone covers tracer transport and external wall contact.
Adsorption and Figure 7 are outside the test.
