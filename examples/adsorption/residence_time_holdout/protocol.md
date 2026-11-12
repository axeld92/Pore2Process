# One-cell path residence times and a complete two-cell tracer holdout

Freeze one parallel-path model using only the qualified N384 one-cell
outlet pulse. Represent the impulse residence times by nonnegative
weights on the 80 µs grid. Remove the effect of the known 4 ms
rectangular injection using a nonnegative fit to the outlet mass in
each interval. Normalize the observed one-cell recovery to unity
for this fit and report the omitted fraction. Correct the mean to the
independently known liquid volume divided by throughput by adjusting
delays beyond one nominal transit. The adjustment rises linearly
to a constant amplitude over one transit and tapers to zero over
the last transit of the measured support. For a path of time T
across one cell, predict 2T across two repeats. This tests complete
persistence of fast and slow paths, a hypothesis for repeated
geometry. It includes no refit, path exchange, adsorption, or
unmeasured tail extrapolation. Apply the same physical 4 ms input
pulse to both predictions.

The one-cell reconstruction must conserve probability and match the
archived cumulative response within 0.005 absolute. Lock the fitted
path distribution, predicted curves, and their hashes before running
or using the two-cell comparison. A prefix of the two-cell response
was already published for the previous dispersion-only test. It is
excluded from calibration, and the late part is a new holdout.

Use the same archived N96 stock TRT carrier flow, balanced shared
flows, two-cell geometry, N288 tracer mesh per repeat, 80 µs step,
molecular diffusivity of 1e-9 m²/s and open boundaries as the prior
test. Run all 21,710 steps through 1.7368 s. The established early
rejection rule is disabled so the full reference curve is available.
Check the reproduced early prefix and retain the final concentration.

Require the complete pore curve to conserve mass to 1e-8 of the
injected dose, stay within 1e-8 of the 0–1 kg/m³ bounds, recover
at least 99.9% and add no more than 1% to either raw time moment
over the last nominal transit. At the matched times use fractions
of the injected dose. Require cumulative-profile RMSE at most
0.030, mean exit time error at most 5% and exit-time variance
error at most 10%. Separately require at most 5% relative error in
5% arrival, 10% relative error in the rise duration from 10% to
75% recovered, 10% relative error in time to 99% recovery, and
25% relative error in the small fraction still in the bed at
0.60 s. That fixed tail time is twice the 0.30 s diagnostic in
the archived one-cell test. All criteria must pass. The known
end of the old early prefix does not set any threshold.
