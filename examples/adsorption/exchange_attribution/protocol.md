# Fixed-inventory attribution of the Figure 7 midpoint discrepancy

Cross three dry-carbon allocations with three path-exchange laws for both
particle sizes. Hold the total dry mass, flow, liquid volume, path flow
shares, path delay distribution, source equilibrium, film correlation and
batch-qualified particle diffusivity fixed at their earlier values. Use the
same inputs for every factor level and give neither Figure 7 curve a fit.

For class flow share (w_i) and relative liquid delay (s_i), use carbon
fractions proportional to (w_i s_i) (the original equal carbon per void
volume), (w_i) (equal carbon per unit throughput), or (w_i s_i^2)
(a fixed sensitivity that concentrates more carbon in the slow paths).
Normalize each to one. The last two allocations are scenarios. The
quadratic rule sets a deliberately stronger dependence on contact path,
with no parameter estimated from a column. The liquid allocations always
remain (w_i s_i). Report the allocation to the fastest and slowest
tenth by flow and the implied local voidage range.

For each allocation use no exchange, the original synthetic marker law,
and a Péclet-matched marker law. For the latter keep the original
dimensionless carrier velocity and solid mask, and set the marker
diffusivity in the reference geometry to (D_{ref}Pe_{ref}/Pe_{column}).
Use the two source column Péclet numbers based on liquid velocity and
scaled period length. Only the ranking transitions change. The frozen
one-cell delay marginals, overall residence time and normalized flow
shares are shared. Reuse the same flux-weighted marker launch, seed,
65,536 markers, and 80 µs and 40 µs steps. Require at least 99.9% pairs,
compare their binned two-period cumulative responses and 99% time and
variance using the previous 0.005/5%/5% step-sensitivity criteria. Use the
80 µs law. These controls isolate the modeled effect of diffusion on
exchange within the same synthetic field. The experimental packing,
flow pattern and local solid allocation remain unmeasured.

Reuse the four frozen baseline runs where the allocation and exchange
match exactly. For new cases check the predicted midpoint at two
discretizations. Start with 1/2 spatial subdivisions, 24/48 radial
shells, 0.08/0.04 s startup steps for 24 s, and 8/4 s subsequent
steps. Require a midpoint change at most 0.5 min and a maximum
normalized curve change at most 0.005 after startup. Refine to four
subdivisions, 96 shells, 0.02 s startup steps and 2 s main steps if
needed. Use the finest result. Keep all phenol balances and the prior
physical bounds. The midpoint is the primary attribution metric. Do
not interpret differences smaller than the resolution of this check.

Freeze the allocation and exchange laws, all column inputs and the
completed predictions before scoring either source curve. Report the
conditional changes in midpoint when changing one factor at fixed
other factor, interactions when both change, and the best remaining
midpoint gap across the declared scenarios. Use the same source RMSE
and midpoint screens for context. An apparent improvement is a
sensitivity, with no inference of the actual carbon or exchange law
from these data. No new hydraulic, OpenLB, or DWSIM simulation is
part of the comparison.
