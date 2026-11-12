# Absorbing-wall film transfer and shared-carbon comparison

Keep the saved, balanced N96 TRT velocity, solid voxel geometry,
two-period inlet and outlet treatment, flux-weighted injections
and advection-only inlet twins. Reuse the two Péclet-matched diffusivities
and the 80 µs and 40 µs marker settings. Treat a sampled first hit on a
carbon face as irreversible uptake in a perfect-sink tracer test.
Continue the same marker without uptake to record its counterfactual
transit and first contact in the second period. Use the first hit
only for the absorbing mass balance. The second-period local test
resets marker weight at the first seam. Two-period survival uses
the original dose without resetting it. The outer axial reflections
are excluded from carbon uptake.

Rank paths using the same diffusion-free inlet twins and 40 flow
quantile classes as the preceding shared-carbon test. For each
class in the first period, measure the captured fraction p and
the ratio a of surviving occupation time to unreacted occupation
time. Here surviving time for a passage ends at first hit or exit.
The perfect-sink balance gives NTU = p/a and per-period conductance
K_i = Q_i*NTU in the target column. This ratio remains finite when
every sampled marker in a class is captured. A class with no
sampled hits has zero estimated conductance. Retain the 80 µs
conductances and check the half-step sensitivities. Use the second
period and full two-period survival as diagnostic comparisons for
constant per-period conductance, with the fixed Péclet exchange
matrix. Report any failure of this transfer separately from the
Figure 7 results.

Check that the added observer preserves each archived transit
time and the hit/no-hit flags. Test the absorbing occupation
calculation against a constant-hazard oracle. For time-step
sensitivity require a change in first-period total capture at
most 0.005 and in total inferred conductance at most 15% of the
finer estimate. For transfer to the second period require a
flow-weighted absolute class capture error at most 0.01. Require
the two-period survival error at most 0.01. These thresholds apply
to fractions of the inlet tracer dose. Continue the requested
column comparison even if these diagnostics reveal limits of the
constant-conductance approximation.

Use one shared particle radial profile per slice, with the same
film area shares as the preceding contact-based shared-carbon
case. For each path divide its inferred conductance by the
corresponding area per period to give its effective film
coefficient. An area with zero conductance has zero radial
boundary flux. Leave the external area, total dry carbon mass
of 7.032455155 g, particle diffusivity, isotherm, carrier volume
and Péclet-matched exchange fixed. Retain the previously declared
correlation-based shared-carbon cases as controls. The perfect-sink
coefficients describe transfer on the synthetic carrier field.

Complete both adsorption runs through 18,000 s. Use the existing
conservation, bounds, refinement and Figure 7 screening criteria.
Use independent finer startup or midpoint windows if necessary,
and flag crossings that remain sensitive. Freeze the physical
inputs and resulting curves before scoring. No Figure 7 data fit
any model parameter. No new OpenLB hydraulic or DWSIM run is part
of this comparison.
