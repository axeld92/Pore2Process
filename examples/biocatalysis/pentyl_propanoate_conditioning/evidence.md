# Independent catalyst conditioning evidence

The [numeric source record](inputs.json) retains source conditions,
unreported quantities, and PDF hashes for the newly extracted materials.
The prior [bead evidence](../pentyl_propanoate_transport/particle/evidence.json)
and [water-activity evidence](../pentyl_propanoate_water_balance/evidence.md)
supply the geometry, moisture, diffusivity and liquid-property inputs.
The independent studies use the same commercial immobilised lipase or
its reported support; none measures the actual initial state of the
reactive bed in this example.

## First liquid contact

The 2026 stirred-batch procedure heats the immobilised lipase together
with 29.7 mL pentan-1-ol, then introduces 10.3 mL propanoic acid at the
reaction temperature. The catalyst has had alcohol-only contact, of
unreported duration. The 2025 column procedure starts pumping the
premixed acid and pentanol after heating the jacketed column to the
setpoint. The paper does not specify an alcohol-only preconditioning
step or report the initial pore-liquid composition and catalyst
history for each run. Thus neat-acid exposure in the bed must not be
assumed. [Batch procedure](https://doi.org/10.1021/acs.iecr.6c02375) ·
[Bed procedure](https://doi.org/10.1021/acs.iecr.5c02780)

## Swelling and liquid exchange

Chaussard et al. measured Novozym 435 bead-size distributions before
and after one hour of stirring in equimolar ethyl butyrate and
racemic 2-pentanol. Wet-to-dry mean diameter ratios for the three
sieved fractions are 1.208, 1.169 and 1.091. The dry and wet
distributions are separate samples, with dispersions retained in
the [previous particle evidence](../pentyl_propanoate_transport/particle/evidence.json).
The observations establish solvent-dependent swelling; they do
not give rates of swelling in either target feed, nor a measured
change in the active fraction. The original reactor paper gives
593 µm mean diameter and two separate porosity estimates (0.40
from aqueous tracers and 0.22 from dry nitrogen sorption). Its
tortuosity of 5.3 was adopted from another reaction.
[Chaussard et al. (2024), SI Section 1, Tables S1–S2](https://doi.org/10.1021/acsengineeringau.4c00030)

The 2025 paper's water tracer experiment at 1.2 mL/min gives mean
exit times and variances of 19.07 min and 27.83 min² for the larger
blue dextran, and 24.27 min and 56.84 min² for sodium azide. Use the
separately [digitised curves](../../../literature/sources/lorenzo2025/water_tracer.csv)
and [source record](../../../literature/sources/lorenzo2025/extraction.json)
for independent checks of these moments. The nominal pulse volume
is 8.4 mL. No measured injection shape or water temperature is
available. The integrated tracer doses and printed doses do not
agree; only the normalized shapes enter the conditional exchange
calculation. [Lorenzo et al. (2025), Appendix A Figure A1, Table A1](https://doi.org/10.1021/acs.iecr.5c02780)

## Water uptake

Kuroiwa et al. equilibrated water with Novozym 435 in separate
organic solvents at 40 °C, with the liquid sampled after 100 h
in closed vials and analyzed by Karl Fischer titration. Water
in the enzyme preparation is inferred from the water balance.
Figure 6 and Table 1 give the enzyme's fitted Langmuir parameters:

| Liquid | Capacity (g water/g catalyst) | Half-capacity concentration (g/L) |
| --- | ---: | ---: |
| 2-methyl-2-butanol | 0.025 | 0.42 |
| Acetone | 0.045 | 0.42 |

The acetone observations originate in the authors' previous study
and are shown with a fitted curve. The 100 h incubation is a
measurement endpoint. No time-dependent uptake law or fitted
parameter uncertainty is reported. A common range up to 1 g/L
is used for equilibrium interpolation and the bead relaxation
screen. The capacity alone supplies a deliberately generous
retention scenario for more concentrated liquids. These
solvents and their isotherms do not specify the actual capacity
or water activity in the reacting acid/pentanol/ester mixture.
[Kuroiwa et al. (2015), Section 2.5, Eq. 4, Figure 6 and Table 1](https://doi.org/10.5539/jfr.v4n5p133)

The existing fresh-catalyst moisture source reports 0.7 ± 0.2%
water for a separate commercial lot. Starting completely dry
maximizes potential additional uptake in the batch capacity
scenario. Subtract the high moisture inventory in a separate
sensitivity. [Melgosa et al. (2015), Section 2.1](https://doi.org/10.1016/j.supflu.2014.11.003)

## Repeated reaction and acid contact

Lorenzo et al. studied reuse with the same alcohol and preparation,
replacing propanoic acid with acetic acid. They used a 40 mL
batch at 40 °C, 2:1 pentanol:acid, 3.6 mol/L initial acid, 10%
catalyst relative to acid and 100 rpm. Their Table 1 gives
84–87% final conversion and 7.8 h to the stated final-conversion
criterion in every one of ten cycles. The lower panel of Figure 2
shows cycles 1, 5 and 10. It has sufficient resolution for
interpolating the time to a common 70% conversion (about 4.54,
4.91 and 4.96 h). The coordinates, axis calibration and reading
range are in [inputs.json](inputs.json); the modest apparent
rate ratio is inferred from equal progress and is independent
of our packed-bed data. The 2023 experiment does not measure
intrinsic catalyst activity or use propanoic acid.
[Lorenzo et al. (2023), Section 2.3, Figure 2, Table 1](https://doi.org/10.3390/pr11061640)

Meneses et al. tested Novozym 435 with propanoic acid and benzyl
alcohol. Table 1 of the journal article (the equivalent table in
the author's thesis is Table 7, Chapter 3) reports residual
esterification activities on a separate lauric acid/n-propanol
assay (1:1, 5% catalyst, 60 °C, 40 min, triplicate):

| Exposure before assay | Activity (U/g) |
| --- | ---: |
| Fresh | 54.6 ± 0.7 |
| Benzyl alcohol contact | 46.2 ± 0.4 |
| Propanoic acid contact | 1.7 ± 0.7 |
| Second fed-batch cycle, acid:alcohol 1:1 | 1.0 ± 0.2 |
| Second fed-batch cycle, acid:alcohol 1:3 | 1.1 ± 0.4 |
| Second fed-batch cycle, acid:alcohol 1:5 | 18.5 ± 1.2 |

The source does not identify the statistic represented by the
plus/minus values. The two fed-batch cycles lasted 24 h each
at 50 °C and 150 rpm, with acid fed in the first 5 h, molecular
sieves, and acetone washing and drying between cycles. The
separate-solvent exposure duration and temperature are unavailable
in the cited propanoic-acid study. The endpoints demonstrate
susceptibility to propanoic acid and solvent history. They
cannot supply the rate of damage during 2:1 mixed-feed
esterification in pentanol.
[Meneses et al. (2019)](https://doi.org/10.1007/s00449-019-02159-w) ·
[Author thesis, Ch. 3, printed p. 53](https://repositorio.ufsc.br/xmlui/bitstream/handle/123456789/219184/PENQ0862-T.pdf?isAllowed=y&sequence=-1)

## Separate check on the duration of conditioning

The filled-circle series in the 2025 Supporting Figure S2 shows
pentyl propanoate at alcohol:acid 2, with a sample at 10 min
and subsequent measurements every 30 min. Flow and temperature
are not identified for this series. Digitised conversions at
30–180 min vary over about 1.22 percentage points. This is a
comparison to published stability, not an activity measurement,
and the precision of marker reading is not experimental
uncertainty. In particular, stable conversion cannot establish
the initial activity of the beads or rule out loss before the
first steady sample.
[Lorenzo et al. (2025), SI Figure S2](https://doi.org/10.1021/acs.iecr.5c02780.s001)
