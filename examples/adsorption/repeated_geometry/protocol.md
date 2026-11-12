# Eight-cell tracer holdout at fixed carrier exchange

Tile the archived balanced N96 OpenLB TRT field eight times. Join
the same axial faces as in the qualified two-period geometry. Use
the same Péclet-matched diffusivities at 80 and 40 µs for both
particle sizes, 16,384 flux-weighted markers per case, and a new
inlet seed shared across the four cases. Track each marker's
first passage through every seam. A sampled solid reflection
counts as perfect-sink wall contact. Later reflections and
passage times describe the counterfactual unabsorbed trajectory.
Exclude the outer reflecting inlet and outlet from wall contact.

Verify the longer tracker first against a deterministic straight
chain and a same-seed subset of the archived two-cell outputs.
The holdout has a 12 s synthetic observation horizon. Require
at least 99.9% inert outlet recovery and retain explicit bounds
for trajectories censored by the observation horizon.

Freeze the prediction before running the new tracer cases.
Propagate the same first-cell survival probabilities and two-cell
quartile exchange from the previous study, without recalibration.
At each period apply the absorbing survival rule, then transport
the remaining dose through the frozen exchange. Report cumulative
survival and the newly captured share of survivors entering each
cell. Report the same values for markers initially in the fastest
90% of inert transit times. That rank boundary already exists in
the 40-class closure.

The matched comparison resets the subgroup composition within
each receiving parent class after cells 2, 4 and 6, distributing
survivors in proportion to the subgroup's frozen flow share.
This changes only the allocation of surviving tracer within each
class. It keeps the same liquid flux between original classes,
the same perfect-sink wall probabilities and the same parent
class tracer mass at each reset. It represents loss of fine path
memory. The two-cell data alone cannot establish a decorrelation
length for an actual disordered bed.

For each particle size and time step, require a maximum error
of one percentage point over cells 3 through 8 for total
cumulative survival, conditional capture among all incoming
survivors, and cumulative survival in the initially fast cohort.
Compare the 80 and 40 µs results at the same seed and require
the same limits for time-step sensitivity. Report the model
contrast with reset as a sensitivity of the frozen closure.

The tracer geometry uses the same perfect-sink surface in both
predictions. The pinned 7.032455155 g dry carbon inventory
and Péclet-matched carrier exchange are carried forward from the
column case. The tracer calculation uses a perfect sink and
does not consume carbon capacity. This milestone ends with the
tracer comparison. No adsorption column or Figure 7 fitting is
part of this study.
[Encounter-based diffusion and escape](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.107.044105)
