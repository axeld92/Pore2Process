# Transfer the N160 memory model to the fixed-mass phenol column

Use the Alhamed Figure 7 feed, flow, bed dimensions and porosity. Preserve
7.032455155 g of dry carbon, the printed-density inventory basis, the
batch-qualified spherical loading diffusivity and the source Langmuir
isotherm. Use the same superficial-velocity Wakao film scenario with
multiplier one. Figure 7 observations enter only the post-freeze score.

Reconstruct the inlet through section 2 from the tracer training markers.
For sections 3–8, reuse the saved 40-state transition and delay moments
and verify them against independent probes. Give each allowed transition
a two-point, positive transit law matching its first and second delay
moments. With mean m and variance v, the lower delay is
max(m/2, m − sqrt(v)), the upper delay is m + v/(m − lower), and weights
give mean m. A zero-variance branch has one delay.

Repeat the eight-section pattern to cover the actual bed. Choose the
nearest number of whole patterns using the archived N160 width for
0.80 mm and scaling it by particle diameter for 0.45 mm. Scale all
delays by one common factor per bed so their flow-weighted total
hold-up equals the reported liquid volume. At joins between patterns,
mix the outgoing streams into one feed to the next pattern. This is
a closure beyond the eight verified sections. Record the actual
section length and Péclet number for both curves.

Each section has one shared radial loading profile with equal dry
mass per section. Delays hold the advecting liquid, and the outlet of
each section transfers phenol through separate films to its shared
particles. Assign the section's geometric external particle area to
its receiving states by the share of independent probes with a solid
wall hit in that section. This is a contact proxy. Use no inferred
adsorption probability or fitted dispersion. Run an area-by-flow
control with the same transport and a one-age-group control with the
same path-level contact shares. Use 80 µs probe inputs for the primary
curves and check the 40 µs inputs at the selected column resolution.

Verify inert mean and variance first on sections 3–8 at the original
probe scale, within the existing 5% and 10% limits. For the column,
check physical bounds, dry/void inventory, film balance and overall
phenol balance. Refine the loading-shell mesh and time steps through
18000 s. Use the established 2% crossing and 0.005 curve-difference
limits. Freeze all curves before scoring the open Figure 7 markers
with the established RMSE, 5% crossing and midpoint screens.
