# Effect of the measured water-tracer spread on conversion

**Bottom line.** The published water-tracer spread reduced predicted
conversion by at most 2.24 percentage points under a nominal rectangular
inlet pulse. The model still overpredicted all 20 reactor measurements.
Even the broadest moment match allowed by a 7 min inlet window reduced
conversion by at most 4.08 points and overpredicted all 20 runs.
[Comparison](reference/bed_comparison.csv)

## Deriving the mixing input

Figure A1a gives the blue-dextran outlet concentrations at a water flow of
1.2 mL/min. Dividing the 41-point digitized curve by its trapezoidal integral
gives a mean exit time of 19.0906 min and a variance of 27.7799 min². Table A1
gives 19.07 min and 27.83 min². Normalizing the curve uses its timing and
shape without requiring the absolute dye dose and recovery to agree.
[Figure A1 extraction](../../../../literature/sources/lorenzo2025/README.md#catalyst-and-water-tracer-data)

Assume the 8.4 mL dose enters at the stated water flow as an independent 7
min rectangular pulse starting at zero on the printed time axis. Actual
inlet timing and concentration are unavailable. Under independent
convolution the means and variances add. Subtracting the pulse mean of
3.5 min and variance of 4.0833 min² leaves 15.5906 min and 23.6966 min²
for the bed. The pulse-free relative variance is 0.09749. The Julia model
matches it to an axial Péclet number using the Danckwerts inlet and
zero-gradient outlet expression

```text
variance / mean² = 2/Pe - 2(1 - exp(-Pe))/Pe²
```

This gives Pe = 19.46 from the digitized curve and 19.36 from the printed
Table A1 moments after the same pulse subtraction. The paper tabulates a
separate correlation-derived Pe for the tracer experiment.
[Water-tracer source record](../../../../literature/sources/lorenzo2025/extraction.json)

Two pulse sensitivities use an impulse at time zero and an impulse at
7 min. If the inlet is supported on that window, the latest impulse gives
the largest possible pulse-free relative variance under independent
convolution. Its inlet variance is zero and its inlet mean of 7 min gives
the smallest residual bed mean. It is a conditional bound on the moment
match, requiring all dye to enter at the final time.

## Conversion check

The calculations use the original 2:1 batch-derived rate coefficients,
3.3 mol/L acid feed, 6 g catalyst and 20 combinations of temperature and
flow. Only the variance-based Pe enters the dispersion scenarios.
Transferring the same dimensionless spread from water to the reaction liquid
at every reacting flow is an assumption. The runner writes the predictions
before reading any of the column conversion measurements.
[Frozen benchmark](../../pentyl_propanoate/REPORT.md)

| Scenario | Danckwerts Pe | Largest conversion decrease vs. plug flow (pp) | Held-out RMSE (pp) | Mean overprediction (pp) |
| --- | ---: | ---: | ---: | ---: |
| Plug flow | ∞ | 0.00 | 13.35 | 11.86 |
| Nominal rectangular pulse | 19.46 | 2.24 | 11.84 | 10.32 |
| Impulse at the start | 25.20 | 1.78 | 12.15 | 10.64 |
| Maximum spread on the nominal inlet window | 9.41 | 4.08 | 10.60 | 8.99 |

At 60 °C and 3 mL/min, the nominal pulse case predicts 69.1% versus
49.7% measured. The maximum spread case predicts 67.3%. Doubling the
reacting mesh from 2,048 to 4,096 cells changed the calculated conversion
by less than 2.0 × 10⁻⁷ as a fraction, and the largest normalized acid
balance residual was below 6.7 × 10⁻¹².
[Run record](reference/summary.json)

The conversion gap persists with substantially more spreading than the
synthetic OpenLB water bed gave. This test supports inspecting the
transfer of the batch-derived catalyst rate before further refining that
packing. It leaves dispersion in the reacting liquid and the cause of the
water-tracer discrepancy unresolved.

## Reproduce

Run from the project root with the established JuliaCall environment.
Outputs include the frozen predictions, their held-out comparison and a
summary with source hashes.

```sh
source scripts/env.sh
uv run --locked python -m examples.biocatalysis.pentyl_propanoate_transfer.mixing.run \
  --output results/biocatalysis/pentyl_propanoate_transfer/mixing/reproduction
```
