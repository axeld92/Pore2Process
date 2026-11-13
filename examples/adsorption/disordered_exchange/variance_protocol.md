# Attribute the held-out travel-time variance error

Keep the frozen two-section Julia prediction and its ten transit-time
classes. Use all four saved N160 cases. Sections 3–8 supply diagnostic
observations only. No coefficient, class count, or acceptance threshold
is updated from those sections.

For a marker arriving at a seam, let A be its elapsed time and D its
time through the next section. The variance gained across that section
is Var(D) + 2 Cov(A,D). Use observed donor-class elapsed-time means to
ask what this gain would be under four combinations of transition
frequencies and conditional delay moments. Each input can come from
the frozen training seam or the observed seam. Within each combination
the elapsed time depends only on the donor class. A zero-count pair
uses the receiver class's pooled delay moments when an alternate
transition table gives that pair positive weight. Report the affected
probability.
Repeat the audit using the nearest sampled donor rank in that receiver
class for missing pairs to check the sensitivity to this choice.

Attribute the two input changes by averaging each change with the other
input frozen and with it observed. This shares their interaction equally.
Measure missing travel-time memory as the remaining gap between the
fully observed, memoryless gain and the measured gain. Split this into
sorting of elapsed time by the receiving class within a donor class,
and elapsed-time/delay covariance within a donor–receiver pair. Include
the difference between the frozen increment and its counterpart
initialized with measured donor ages as a separate propagation term.

Check that the four contributions sum to the observed minus frozen
variance increment at every seam, and that their cumulative sum gives
the variance error at every held-out section. Report signed cumulative
terms at section 7, the largest-error section identified by the prior
frozen holdout, and at the outlet. Also report the maximum across seams
and show both tracer time steps. The observed-input combinations are
diagnostic substitutions. The original prediction and its 10% variance
limit remain the qualification result.

This attribution concerns inert tracer travel times in one synthetic
packing. It does not use wall-hit outcomes, adsorption or Figure 7.
[Spatial Markov transport in heterogeneous porous media](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.78.026308)
