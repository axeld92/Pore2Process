# Travel-time memory in the section-specific N160 model

Reuse the already sampled N160 stock TRT trajectories. The archived
reference ensemble and the independent operator-probe ensemble retain
their separate seeds. Keep the packing, carbon mass, Péclet numbers,
two tracer time steps and ten local transit-time path classes.

Initialize tracer mass and the first two elapsed-time moments at section
2 using only original sections 1–2. Within each path class, split
markers by elapsed-time rank, from youngest to oldest. Use four groups
per path for the primary model and two and eight groups for numerical
sensitivity. Group sizes differ by at most one training marker within
each path. Each probe subset uses proportional group sizes.

For each downstream seam pair, use elapsed times from the independent
probes to **rank** markers within the donor and receiver path classes.
Tally the transition probabilities and the first two weighted moments of
the receiving section's local delay for each pair of refined states.
No cumulative probe arrival moment is used as a prediction or delay.
Collapsing the refined states must recover the original section-specific
path exchange and local delay moments. One age group per path must
reproduce the existing model.

Julia advances mass and both elapsed-time moments through the refined
states. Within a donor group, the residual variation of elapsed time is
assumed independent of the transition and new delay. The two independent
probe halves give a sampling sensitivity check with four groups per path.
The resulting rules and predictions are frozen before the original
section 3–8 arrivals are accessed by the scorer.

The acceptance test is a maximum absolute variance error of 10% at
every section 3–8 for each of the four particle size and time-step
conditions, using the four-group model with the pooled probes. The
5% mean-arrival limit is retained. This checks the reduced description
on the same fully probed packing against independent trajectories.
Adsorption and Figure 7 are outside the scope.

[Joint transport of tracer age and velocity](https://doi.org/10.1002/2017WR020578)
