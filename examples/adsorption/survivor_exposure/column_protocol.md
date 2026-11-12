# Shared-carbon column transfer after the two-cell holdout

Run both diameters only after the frozen quartile model passes the
holdout at both time steps. Use the 80 µs quartile model for each
primary curve. Retain the 40-class residence-rank model as a control
and transfer the 40 µs quartile model at the selected column
resolution to check sensitivity to the wall-tracking time step.

Use the same Péclet-matched flux between parent classes. Within each
parent class use its frozen column transit delay for all subgroups.
This keeps the carrier volume, flow-weighted delays and the joint
inert transit distribution through two repeats. The refined flux
propagates the differing subgroup concentrations at the repeat
boundaries. Divide the archived fluxes by the receiving marginal to
recover the receiving matrix and verify recovery of the original
aggregate Péclet matrix within numerical balancing tolerance.

Map the first-cell empirical survival to a plug-flow transfer number
as `NTU = -log(survival)`, so that a perfectly absorbing surface
recovers that group's measured outlet survival in the plug limit.
For an all-hit group, use half a survivor divided by its marker count
to obtain a finite coefficient. Record its flow share and the effect
of this endpoint regularization on the two-cell prediction. A group
without hits has zero film conductance. The per-period conductance is
the group's liquid flow times NTU.

Use solid octants as in the previous contact study. Distribute each
octant's share of the fixed total geometric particle area among the
groups using their first-period reflection-distance shares. This is
an area-allocation proxy. Its aggregate equals the original external
area, and there is still only one particle radial loading profile per
axial slice. Keep the dry carbon mass, particle size and density,
batch diffusivity, isotherm, liquid volume and operating conditions
from the shared-carbon study. Neither holdout-cell wall events nor
Figure 7 influence the film, area, exchange or loading inputs.

Use the existing shared-carbon variable-film Julia solver and the
previous refinement, startup and conservation checks. Complete
18,000 s for both diameters. Save the selected curves and checks
before scoring against Alhamed's Figure 7. The experiment tests
transfer of the synthetic closure to both column cases.
