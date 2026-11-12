# N192 time-step check on the saved TRT flow

The [protocol](n192_time_protocol.md) uses the exact N192 subcell faces
and initial concentration field from the [matched-mesh study](REPORT.md).
The bounded Julia transport solver, molecular diffusivity and both
physical observation times remain fixed.

| Interval endpoint | Original step | Half step | Steps in the interval |
| --- | ---: | ---: | ---: |
| 0.1073725554 s | 163.4285 µs | 81.7143 µs | 657 → 1,314 |
| 0.2147641283 s | 163.4575 µs | 81.7287 µs | 657 → 1,314 |

The final comparison uses the liquid concentration variance and phase of
the fitted longitudinal fundamental at the same physical time. Relative
variance changes use the original N192 result as the denominator. Phase
changes use the shortest signed angular difference, and profile RMS is
normalized by the initial 0.25 kg/m³ amplitude.

## Results

| Observation | Relative variance change | Phase change (rad) | Profile RMS change |
| --- | ---: | ---: | ---: |
| 0.1073725554 s | +0.0417% | −4.2804e-4 | 0.178% |
| 0.2147641283 s | +0.1255% | +4.1543e-6 | 0.217% |

At the final time, variance increased from 0.0129861738 to
0.0130024674 kg²/m⁶, a difference of 1.6294e-5 kg²/m⁶. The phase
increased from 2.2464397941 to 2.2464439484 rad. The fitted amplitude
increased by 0.0259%. These time-step effects on final variance and
phase were small relative to the preceding N96–N192 mesh differences
of 8.94% and 0.0976 rad. The complete profile still missed the 0.1%
time-step RMS criterion at both observations.

The numerical share of final half-variance loss decreased from 5.63%
to 5.42%. The smaller-step solution stayed in 0.2500335–0.7499665 kg/m³,
and its maximum relative liquid inventory drift was 2.23e-16.
Every second new timestamp reproduced the archived time exactly.

## Verification record

One Julia tracer execution completed in cohort
`20261111T162208Z-42f74625`. Read-only replay checked 26 source hashes and
33 artifacts. The five existing bounded-solver analytical checks passed.
The calculation used the archived N192 faces and initial field without
rebuilding the mesh or launching an OpenLB, full-bed or DWSIM run.
[Run and replay commands](README.md#n192-time-step-check)
