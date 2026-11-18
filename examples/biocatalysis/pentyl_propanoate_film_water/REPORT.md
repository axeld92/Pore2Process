# Can the liquid film or water activity explain the lower bed conversion?

**Bottom line.** The independently sourced packed-bed film estimate
lowers predicted conversion by at most **1.40 percentage points**.
Even the slower diffusivity scenario lowers it by at most **1.65
points**. All 20 held-out conversions stay overpredicted at the full
6.0 g loading. Independent batch experiments show that water activity
can change the lipase rate substantially. The size of the effect
depends on the liquid and preparation, and the water activity the
enzyme experiences in this particular bed remains unknown.

## How much does the film slow the reaction?

The film is the liquid next to the bead surface. Acid has to cross
it before it reaches the immobilised enzyme. With the central inputs,
the corrected Schmidt number ranges from about **702 to 4,723**.
Measured pentan-1-ol viscosity and a propanoic-acid diffusion estimate
give a film coefficient of **6.80–22.0 µm/s** over the 20 conditions.
At the inlet, the film lowers the local calculated rate to **91–95%**
of its film-free value. Pure pentanol data and analogues supply the
transport inputs, with the actual reacting mixture's properties
still unresolved. [Property sources and correlation](film/evidence.md)

| Case | Conversion at 60 °C and 3 mL/min | Error across 20 runs (RMSE, pp) |
| --- | ---: | ---: |
| Full 6 g, film-free reference | 60.23% | 10.46 |
| Wilson–Geankoplis, central diffusivity | 59.00% | 9.58 |
| Wilson–Geankoplis, slower diffusivity | 58.80% | 9.44 |
| Wilson–Geankoplis, faster diffusivity | 59.15% | 9.69 |
| Reactor Equation B.2, corrected Schmidt number and central diffusivity | 56.95% | 7.70 |
| Reactor Equation B.2, corrected Schmidt number and slower diffusivity | 56.45% | 7.27 |
| Isolated stationary sphere, central diffusivity | 55.15% | 7.34 |
| Measured | 49.70% | |

The two low-transfer alternatives lower predictions by up to 5.44
points for Equation B.2 and 6.07 points for the isolated sphere.
Equation B.2 gives lower Sherwood numbers in this low-Reynolds,
high-voidage bed than the independent packed-bed correlation. The
isolated-sphere scenario represents stagnant liquid. Even these
cases overpredict every observed conversion. The comparison keeps
the printed internal effectiveness and the earlier rate law fixed.
[Predictions and full comparison](film/reference/bed_comparison.csv)

Refining from 512 to 1,024 axial cells changed conversion by
at most **0.00117 percentage points**. The largest relative axial
acid balance error was **$3.71\times10^{-10}$**, and the largest
inlet film balance residual was **$4.27\times10^{-13}$ mol/(L min)**.
[Calculation record](film/reference/summary.json)

## What the batch results say about water

Increasing initial water activity from about 0.12 to 0.97 lowered
the initial rate of the same commercial lipase by **58–83%** in
MTBE in two ethyl butyrate studies. In the 2023 study, the decrease
was **62%** in hexane and **50%** in acetonitrile. These are changes in
initial rate. The studies measured ester formation early in the run.
[Tjørnelund et al.](https://doi.org/10.1002/cctc.202201207) ·
[Dutta Banik et al.](https://doi.org/10.3390/catal7080227)

An endpoint-based log-linear response gave maximum held-out
intermediate errors of **8–12% of the driest initial rate** within the
separate series. Transferring the 2023 MTBE slope to the 2017 MTBE
batch gave **21%** of the anchor rate at water activity 0.97, where the
other batch measured **49%**. Both proportions use that batch's own
rate at water activity 0.23 as the anchor. This spread gives the
scale of the error from transferring a water-activity response even
between nominally similar experiments.
[Batch comparison](water_activity/reference/comparison.csv)

A closer solvent-free analogue with the same acid and isobutanol
slowed to **19%** of its no-added-water rate at the largest water
dose. Water activity was unreported in that experiment. A separate
same-acid, 1-octanol study found a rate maximum at water activity
**0.45**, so added hydration can also help on the dry side of an
optimum. [Kuperkar et al.](https://doi.org/10.1016/j.molcatb.2013.10.024) ·
[Nordblad and Adlercreutz](https://doi.org/10.1002/bit.21706)

## Next step

Constrain whether the packed-bed beads see a different water
activity or hydration history from the original pentanol batch
at the same conversion. Both reactions generate the same amount
of water per mole of ester. The batch-fitted rate can already
reflect water generated in that batch. Published activity or
phase data for the reacting mixture and controlled-water batches
of the same chemistry would support a bed-specific test. With
analogue data alone, an explicit sensitivity bracket can show how
large a bed-to-batch difference would have to be. The 20 reactor
observations remain reserved for checking the resulting predictions.
[Evidence and limits](water_activity/evidence.md)
