# Why does the batch model overpredict packed-bed esterification?

**Bottom line.** The measured water-tracer spread reduced predicted conversion
by at most 2.24 percentage points. The catalyst concentration basis changes
the 2025 model's result by up to 13.61 points and explains much of the
difference between the original prediction and the published curves.
Direct transfer from measured batch progress still overpredicts conversion.
The cause of the lower apparent activity in the reacting bed remains open.

## Measured mixing has a small effect in the reaction model

The normalized Figure A1a blue-dextran outlet curve gives a mean of
19.0906 min and variance of 27.7799 min². Subtracting the moments of an
assumed 7 min rectangular inlet pulse gives a pulse-free mean of
15.5906 min and variance of 23.6966 min². These moments give a
Danckwerts axial Péclet number of 19.46. The original synthetic OpenLB
water result gave Pe = 243.37 in the same boundary convention.

Using the measured water spread as a reacting-liquid sensitivity reduces
RMSE from 13.35 to **11.84 percentage points**. All 20 predictions are
still high. Even the largest relative spread compatible with an
independent inlet pulse in the assumed 7 min window gives **10.60
points RMSE** and all 20 still high. These results apply to this
axial-dispersion model and the assumed transfer of dimensionless
spreading from water to the reaction liquid.
[Measured-tracer sensitivity](mixing/README.md)

## The published curves are sensitive to the catalyst concentration basis

The 2025 article prints a positive acid-disappearance rate with a
positive source sign in its acid balance. The reconstruction uses
acid consumption, its printed fixed inlet concentration and zero
outlet gradient. It uses the Bi–Bi Ping–Pong coefficients in SI Table
S3 and internal effectiveness factors in Table B1. The apparent
minute and g/L rate basis reproduces 68 digitized 2025 batch markers
with 4.39 points RMSE. Interpreting the coefficients on an hourly
basis gives 62.43 points RMSE. The authors leave the rate units and
the mapping of the weighed catalyst to the bed rate coefficient
unspecified.
[Equation and source audit](../../../literature/sources/lorenzo2025/author_model.md)

| Catalyst convention in the liquid balance | Catalyst concentration (g/L) | RMSE to authors' Figure 6c lines (pp) | RMSE to measurements (pp) |
| --- | ---: | ---: | ---: |
| 6 g divided by total bed volume | 213.15 | 1.05 | 3.81 |
| 6 g divided by interbead liquid volume | 333.05 | 11.77 | 10.46 |

Both use the printed internal effectiveness factors and the interstitial
residence time of 18.02 min at 1 mL/min. That time is based on 18.02 mL
of liquid in a 28.15 mL bed. Using the bulk concentration in the
liquid-volume balance gives **64% of the full mass exposure**, equivalent
to 3.84 g on the liquid basis before the effectiveness factor. This
reconstruction comes close to the plotted curves. Its largest deviation
of 1.66 points exceeds the plot-reading half-range of 0.5 points.
The reported catalyst loading remains 6 g. The authors' exact
implementation and the fraction of catalyst active in the bed remain
open questions.
[Reconstruction and comparisons](published/README.md)

## Direct transfer tests the batch-to-bed assumption

We then evaluated the measured batch curves at the equivalent catalyst
exposure without imposing a fitted rate law. In plug flow, divide 6 g
by the inlet acid molar flow. In batch, multiply elapsed time by
the initial catalyst-to-acid ratio. At 60 °C and 3 mL/min the
bed exposure is 0.6061 kg min/mol acid, equivalent to 79.05 min
at the 2025 printed batch loading or 81.81 min at the 2026
nominal 10% loading. Linear interpolation of measured points
gives 71.2% and 72.8% conversion respectively, compared to
49.7% in the packed bed.
[2025 batch marker extraction](../../../literature/sources/lorenzo2025/batch_2to1_markers.md)

| Batch transfer on the common 15 conditions | RMSE (pp) | Mean overprediction (pp) |
| --- | ---: | ---: |
| Original fitted plug-flow model | 14.93 | 13.64 |
| 2025 batch measurements | 15.59 | 14.97 |
| 2026 printed batch conversion | 16.92 | 15.36 |
| 2026 conversion recalculated from acid mole fractions | 18.40 | 16.97 |

All four routes overpredict all 15 conditions. The 2025 data also
cover three cases at 0.7 mL/min, and all 18 covered predictions are
high. The remaining low-flow cases exceed the available batch times.
These tests retain the assumption of proportional catalyst exposure
and the plug-flow comparison. Their result makes the choice of
fitted rate equation an insufficient explanation of the gap.
[Direct comparison](reference/bed_comparison.csv)

## Next steps

First check the conversion basis for the reacting liquid. The paper's
outlet measurement is a concentration, so the reproduction uses a
constant liquid-volume assumption to compare it with molar conversion.
Available phase-equilibrium and composition data could bound the effect
of reaction-generated water and liquid-volume change under the reported
conditions. Then test independently supported catalyst contact or
internal transport assumptions at the full 6 g inventory. The
20 conversion measurements can continue to serve as the prediction check.
