# Does the way the lipase first contacts the liquid explain lower bed conversion?

**Bottom line.** The manner of first contact is a credible difference.
Published studies show swelling of the immobilised beads, reversible
water uptake and a strong loss of residual activity following propanoic
acid exposure. Ordinary diffusion into liquid-filled pores is fast.
The water-uptake scenarios give appreciable startup moisture memory,
especially near a dry inlet, but much less water retention relative to
reaction-generated water at substantial conversion. A simple large
gradual loss (or recovery) during the three-hour run is inconsistent
with the scale of the plotted conversion stability. An early persistent
loss, a distinct intrabead liquid environment, or a different
catalyst-mass dependence remains possible. No independently measured
conditioning factor for the actual reacting bed was available.

## The procedures really differ

In the 2026 stirred-batch procedure, the catalyst is heated with
pentanol before the acid is added. The 2025 packed-bed procedure
starts pumping the premixed substrates once the cartridge reaches
its operating temperature. The details of the prefill and bead
conditioning history in that procedure are not reported.
The latter exposes the catalyst to a mixture, not necessarily to
neat acid. [Batch study](https://doi.org/10.1021/acs.iecr.6c02375)
· [Bed study](https://doi.org/10.1021/acs.iecr.5c02780)

## Physical entry into the pores is fast in the water experiment

Chaussard et al. observed 9–21% increases in the mean diameters of
three sieve fractions of Novozym 435 after an hour in a reacting
organic-liquid analogue. That confirms swelling; it does not
establish the swelling time in either target feed. In the earlier
[transport audit](../pentyl_propanoate_transport/REPORT.md), measured
porosity and size contrasts reduced predicted conversion, without
closing the gap. [Wetting studies](evidence.md#swelling-and-liquid-exchange)

In the actual reactor's water test, sodium azide has access to
additional storage within the beads whereas blue dextran mostly
follows the interbead pathways. Under a shared aqueous pulse and
simple mobile-to-immobile exchange model, their **normalized**
tracer time distributions imply a mean individual return time
of 0.92 min (1.01 min using the printed moments). Allowing any
shared inlet pulse constrained to the nominal seven-minute
window gives a conditional upper envelope of 1.77–1.86 min.
The additional mean storage corresponds to about 6.24 mL.
Recovered tracer masses (6.72 and 2.48 mg) disagree with the
printed injected masses (8.4 and 0.84 mg); these estimates use the
time distributions, not absolute tracer mass. The inferred return
time concerns exchange in a water-wetted bed, not wetting or
sorption in the esterification liquid.

## Water uptake can preserve startup differences

At 40 °C the separate Kuroiwa et al. study gives fitted
equilibrium capacities of **25 and 45 mg of water per gram**
of immobilised catalyst in two analogue solvents. If applied to
the nominal 1.018 g batch dose, these would retain at most
25.4–45.8 mg. At 50% acid conversion the batch has produced
1,237.5 mg of water. These full-capacity scenarios with initially
dry beads would therefore remove 2.1–3.7% of the generated water
from the bulk liquid. At 2% conversion, in contrast, the
within-range equilibrium scenarios retain 16.8–26.5 mg out of
49.5 mg produced. Fresh catalyst moisture would reduce the net
capacity. The actual reacting-liquid isotherm remains unavailable.
[Kuroiwa et al. (2015), Figure 6 and Table 1](https://doi.org/10.5539/jfr.v4n5p133)

At equal conversion, with the **maximum** analogue uptake in
the finite batch and no ongoing net retention in the steady bed,
the earlier liquid-activity projections give water-activity rises
of up to 0.040 at 20% conversion and 0.027 at 50%.
Evaluating the separate controlled-water-activity initial-rate
data only within their measured activity range corresponds to
relative bed/batch rates as low as 0.856 and 0.958, respectively.
The large early difference in water inventory is interesting,
but the independent activity-rate data do not cover the driest
part of the reaction. These solvent-specific sensitivities are
not a target-liquid activity model or a transferred rate correction.

Julia then solved reversible water storage and radial diffusion
within a 593 µm bead at 40 °C with the separately sourced
diffusivity and a liquid film (see [protocol](protocol.md)).
Imagine replacing the surrounding liquid by another liquid with
a prescribed water concentration. How much of the initial-to-final
water-inventory difference still remains during equilibration?
For the slowest bed flow, including an illustrative displacement
age along the bed, the following ranges span the two isotherms
and the 0.42 and 1 g/L concentration steps:

| Particle porosity scenario | Moisture change | Whole bead after 30 min | Whole bead after 3 h | Outer 80 µm after 3 h |
| --- | --- | ---: | ---: | ---: |
| 0.40, aqueous basis | Wetting | 28–57% | 0–2.3% | 0–1.1% |
| 0.40, aqueous basis | Drying | 49–65% | 4.2–14.2% | 2.2–7.9% |
| 0.22, nitrogen basis | Wetting | 41–66% | 0–11.8% | 0–4.7% |
| 0.22, nitrogen basis | Drying | 58–72% | 11–25.6% | 6.0–14.8% |

Here 0% means equilibrated water inventory; 100% means no change
from the initial state. The numbers describe water, not active
enzyme fraction. Without sorption, the corresponding bed-average
water difference is already only 0.03% after 30 min. Adsorption
can therefore be important for startup moisture history even when
liquid diffusion is fast. Drying toward a nearly water-free inlet
retains the largest memory in these cases. The actual inlet water
content, transient sorption, reactive pore-liquid diffusivity and
enzyme activity at these moisture states were not measured.
The largest refinement difference on this calculated inventory
was 0.253 percentage points, with normalized water conservation
error below 4.0e-8.

## Published acid contact effects are large enough in principle

In the same research group's pentanol/acetic-acid study, ten reuse
cycles retained 84–87% final conversion. The plotted times to
70% conversion for cycles 1, 5 and 10 were about 4.54, 4.91
and 4.96 h. Their ratio gives a mild apparent rate sensitivity of
0.914 (marker-reading range 0.878–0.952). This is evidence from
the same alcohol and preparation, but a different acid.
[Lorenzo et al. (2023), Figure 2 and Table 1](https://doi.org/10.3390/pr11061640)

Meneses et al. report residual enzyme activity of **54.6 ± 0.7 U/g**
fresh, **46.2 ± 0.4 U/g** after benzyl-alcohol contact and only
**1.7 ± 0.7 U/g** following separate propanoic-acid contact.
After two fed-batch cycles even the alcohol-rich 1:5 case had
18.5 ± 1.2 U/g. These independent observations support the
possibility of substantial chemical/conditioning damage.
They use a different alcohol, exposure history and activity
assay. The separate neat-acid exposure duration was not reported
and its effect cannot be turned into a time law for the mixed feed.
[Meneses et al. (2019), residual activity table](https://doi.org/10.1007/s00449-019-02159-w)

To check magnitude, apply each independent apparent activity
ratio separately to the prior central-film bed calculation.
The weighed inventory remains **6.0 g** in every run; its activity
per gram, not its mass, is changed. No new kinetics or factor is
fitted to the 20 reactor measurements.

| Relative-rate scenario | Relative rate | Conversion at 60 °C, 3 mL/min | RMSE across 20 runs (pp) |
| --- | ---: | ---: | ---: |
| Prior central film | 1.000 | 59.00% | 9.58 |
| Acetate reuse time ratio | 0.914 | 56.33% | 7.83 |
| Benzyl alcohol contact / fresh | 0.846 | 54.03% | 6.36 |
| Propanoic acid, two 1:5 reuse cycles / fresh | 0.339 | 29.61% | 18.04 |
| Separate neat acid / benzyl alcohol contact | 0.037 | 4.15% | 58.20 |
| Measured | | 49.70% | |

The mild changes move the prediction in the correct direction;
the large residual-activity losses make it much too low. Neither
experiment gives a transferable specific rate factor for this bed.
As a magnitude illustration, the previous concentration-basis
comparison corresponds to 64% of full rate exposure. If a fraction
of the catalyst suffered the separate neat-acid loss, the same
aggregate activity would correspond to about 37% of the mass
being affected. This is an algebraic example using an unrelated
severe endpoint, not an estimate of damaged catalyst.

## Can conditioning keep evolving over each three-hour run?

The target paper's Supporting Figure S2 propionate time trace
(2:1 feed, flow and temperature unspecified) is approximately
74% conversion after 10 min and 77–78% from 30 to 180 min.
The span over the latter period is about 1.22 percentage points.
A simple hypothetical first-order activity decay reaching 0.64
after 30 min would continue to 0.069 after 180 min; the alternative
initial wetting deficit reaching 0.64 activity after 30 min would
recover to 0.998 after 180 min. These two illustrations use the
same 0.64 rate-basis magnitude, without fitting Figure S2.
Steady-state snapshots across the operating grid then predict
conversion falls of 31–59 points or rises of 2.9–13.2 points,
respectively. This is a scale check, not a transient simulation
of the unidentified time-trace conditions.

The observed stability weighs against simple large, ongoing
changes of activity over that interval. It does not rule out
rapid conditioning or irreversible loss before 30 min, slow
water changes with little effect on rate, or a distinct steady
intrabead environment. The original catalyst sequence is a
plausible reason for such differences, but the experiment does
not identify which reason applies.

## Next executable step

Recheck proportionality between the apparent batch rate and the
weighed immobilised catalyst mass using the independent 2026
batch curves at four different loadings at 40 °C. The current
scaling takes the low-loading batch rate to the very high local
solids concentration of the packed bed. Compare mass-proportional
and concentration-dependent rate models on held-out batch loadings
first; then transfer a separately supported model to the fixed
6 g bed. In parallel, seek published time-resolved CALB residual
activity measurements with propanoic acid in primary alcohols,
including alcohol-only preconditioning followed by acid contact
if available. Only an independently supported exposure/rate rule
should be used to predict a persistent loss. Keep the twenty
reactor conversions exclusively for scoring.
