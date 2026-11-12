# Unfitted contact-allocation rerun of both Figure 7 columns

Use the qualified carbon-wall study and the Péclet-matched laws from
the fixed-carbon attribution study. Keep the same 40 transit classes,
one-period delay distribution, column void volume, source operating
conditions and frozen batch particle model. Keep 7.032455155 g dry
carbon in every case, with no inactive fraction and no change to
the film correlation.

The primary allocation is the previously reported advection-ranked
carbon proxy. Recompute it over the existing 40 probability bins using
the same inlet-paired diffusion-free trajectories. For each solid
octant, split its mass among the ranks in proportion to the sum of
the observed reflection lengths on its surface. Use the 80 µs case
at each diameter's Péclet match and record the corresponding 40 µs
allocation sensitivity. The rank label of each adsorption class is
assigned by its flow quantile. This maps the advection-ranked mass
shares onto the equal-quantile diffusive tracer classes of the
Péclet-matched exchange law. The two rank labels can disagree for an
individual marker.

Run a rank-label sensitivity using the measured contact shares
conditioned directly on realized diffusive rank. Use the same 80 µs
events, total mass and exchange matrix. This variant uses the exact
class labels from which the exchange matrix was estimated. Freeze
both choices ahead of scoring. They are scenario closures of
grain-envelope contact and have no independent calibration against
experimental carbon capacity.

Use the unchanged conservative adsorption kernel and start with
2/4 spatial subdivisions per path period, 48/96 radial shells,
0.04/0.02 s startup steps for 24 s, and 4/2 s main steps through
18,000 s. Require physical bounds and normalized total mass and
film-flux residuals below 1e-8. Require a change in the predicted 5%
and 50% crossings of at most 2% and a maximum normalized curve change
after startup of at most 0.005. Resolve crossings in the first 12 s
with a separate short-run refinement when needed. Retain both complete
and short simulations. Use the finest available short curve up to
12 s for early scoring and the full curve thereafter. Record the
concentration difference at the join and qualify it to 0.005 of the
feed. Refine the full duration if needed for later crossings or the
late curve. No calibration may use the observations or comparison
outcomes.

Use the established Figure 7 marker comparison, with fractional RMSE
at most 0.10, 50% crossing error at most 20%, and resolved 5% crossing
error at most 25%. Retain the archived void-proportional,
Péclet-matched cases as the control for the carbon-allocation change.
This is a forward comparison using synthetic packing. No new OpenLB
hydraulics or DWSIM runs are part of this test.
