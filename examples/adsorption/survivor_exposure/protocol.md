# Wall-exposure memory in the two-cell holdout

Reuse the four saved perfect-sink marker ensembles on the balanced
N96 TRT carrier, with both particle diameters and the 80 and 40 µs
time steps. Retain the common dry carbon inventory, the batch particle
model and the Péclet-matched exchange at all 40 original path classes.
This holdout uses only first-hit events at solid voxel faces.

## Candidate fixed before scoring

Use the transit ranks of the inert diffusing markers for the 40 parent
classes, consistent with the archived exchange matrices. Within each
parent class, sort by the advection-only twin transit time and split
into four equal-flow groups. Assign quartiles separately in each period
using only inert times. The twin rank is a proxy for the advective
route's persistent exposure to the wall. A baseline uses just the 40
parent classes.

Use the empirical first-period survival in each group as the local
perfect-sink transfer rule. Estimate the flow between subgroups from
the paired inert transit times. Add one total pseudocount split
uniformly inside each allowed coarse pair, and scale rows, columns and
parent-pair blocks iteratively to the specified fluxes. Zero parent
pairs stay zero. Subgroup inlet and outlet flow must match, and summed
flows between parent pairs must reproduce the frozen Péclet-matched
exchange. Apply first-cell survival before transporting surviving
phenol to the second cell. Apply the same subgroup survival rule
there. This retains a different subgroup composition among survivors.

Freeze the subgroup rules, inert exchanges and predictions before
reading second-period hit values for scoring. The first-hit time in
the second period is a response only. For diagnosis, also apply the
frozen second-cell subgroup rule to the observed survivor population
at the seam, using first-cell hits and second-cell inert times. This
isolates local uptake error from errors in the estimated survivor
composition. This conditional calculation is a diagnostic.

## Gate

Retain the original 1-percentage-point tolerances on absolute
two-period survival error and on flow-weighted absolute error in the
fresh, unit-weight second-period capture across the original
advection-only classes. Also require second-period capture error at
most 1 percentage point among first-cell survivors, normalized by
surviving dose. Report fresh capture and class errors using the
diffusive parent classes as additional diagnostics. All three gates
must pass for each diameter at both archived time steps. Report the
baseline alongside the candidate. Evaluate Figure 7 only if the
candidate passes. A column extension would use the same shared carbon
inventory, batch particle inputs and aggregate Péclet exchange.

The test probes the ability of this coarse state to transfer sampled
wall uptake across two repeats of a synthetic pore geometry. The
perfect-sink first hit provides no new measured carbon uptake data.
[Encounter-based diffusion and escape](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.107.044105)
