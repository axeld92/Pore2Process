# Two-section calibration and disordered tracer holdout

Use the four archived N160 marker ensembles for the 0.80 and 0.45 mm
Péclet-matched conditions at 80 and 40 µs. Use only the first two
first-passage times and first two wall-hit clocks to estimate the model.
Freeze the coefficients and eight-section predictions before scoring the
remaining six sections. The ensembles share the original geometry,
7.032455155 g column-scale carbon basis and inlet sampling. Each case
has its own measured exchange and contact counts.

Define ten local transit-time classes using equal-flow deciles, with
stable ranking and the archived 90% rank boundary for initially fast
tracer. Count joint inert transitions from section 1 to 2. Normalize
each donor row to obtain a transition probability and verify that it
preserves the carrier class fractions. For each donor–receiver pair,
count the markers avoiding wall contact in section 2. Their proportion
within that pair is the wall-free probability conditional on the
transition. The inert exchange multiplies this probability to obtain
a substochastic transition for the wall-free dose. The missing mass
is perfect-sink wall contact. An empty pair has zero exchange, zero
wall-free transfer, and no inferred conditional wall probability.

Initialize at the end of section 2 using the measured wall-free
inventory by receiving class for both the full and initially fast
doses. Repeat the same inert and wall-free transitions for sections
3–8. This tests whether one exchange and wall-contact rule transfers
along this disordered packing. A control uses the same inert exchange
and section-2 initial state with wall survival depending only on the
receiving class. Neither model uses further geometry or outcomes.

For inert timing, start with the empirical first and second raw moments
of cumulative arrival time by class at section 2. Add the measured
section-2 transit-time moments conditional on each transition at every
subsequent step. Given the donor and receiver classes, treat the new
transit time as independent of earlier times. Compare the predicted
cumulative mean and variance against sections 3–8 using the earlier
5% and 10% relative limits. Also report the change in the observed
joint transition distribution from seam to seam, in total variation,
and the predicted fraction of initially fast tracer that remains fast
at every section. Use the first two observed classes to initialize
that persistence diagnostic.

Score maximum absolute errors over sections 3–8 for total wall-free
survival, conditional capture among incoming survivors, and cumulative
wall-free survival of the initially fast cohort. Use the existing
2-percentage-point engineering screen for these holdout observables.
Report each metric separately for all four cases, with 80 µs primary.
Verify inert and wall-free dose balances, probability bounds and
monotonic capture. This is a spatial holdout on one synthetic packing
using the same marker identities. It evaluates tracer transport and
perfect-sink wall contact.
[Wall encounters and first passage](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.107.044105)
