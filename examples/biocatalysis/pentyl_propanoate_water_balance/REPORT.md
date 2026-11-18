# Do packed-bed beads see more water than stirred-batch beads?

**Bottom line.** At the same conversion and temperature, both reactors
have the same **bulk water activity** if they start with the same feed
and retain the same products. The nominal batch headspace is too small
to account for an appreciable difference in that water inventory.
The film around a bead can raise its local water activity, especially
early in the reaction. With the ordinary packed-bed transfer estimate,
the increase is small at 50% conversion. The actual activity inside
the beads remains unresolved.

## What the water and phase balances establish

At 50% acid conversion, both reactors have made half a mole of water
per mole of acid fed. On the dry 2:1 feed basis, that is **0.167 mole
fraction water** in the total liquid. This composition gives one
equilibrium water activity at fixed temperature and pressure. If the
liquid separates, both phases have the same water activity at
equilibrium. The published ternary tie lines illustrate the difference
between concentration and activity. Organic phases contain about
0.176–0.232 mole fraction water and share their activity with the paired
water-rich phases at 0.953–0.990 mole fraction water.
[Amamehi et al., Tables 1–2](https://doi.org/10.12962/j2964710X.v5i1.19035)

The published stirred reactor had 40 mL of charge in a 50 mL vessel and
a condenser. The nominal **10 mL** above the liquid can hold **0.51–2.91
mg** of saturated water vapor from 40 to 80 °C. Even the 80 °C value is
under **0.6% of the water produced by 20% conversion**. In the activity
scenarios this headspace inventory changes bulk activity by at most
**0.0018**. The estimate covers the nominal vapor volume. Condenser
hold-up and any sustained water loss remain unquantified. The 75
selected published batch records give identical water and ester mole
fractions at the printed precision. The authors obtained the
batch water values by mass balance. [Lorenzo and Ortega, Figure 1 and Table 1](https://doi.org/10.1021/acs.iecr.6c02375) ·
[IAPWS vapor pressure](https://www.iapws.org/relguide/Supp-sat.html) ·
[Inventory calculations](reference/headspace.csv)

A separate Lipozyme 435 lot contained 0.7 ± 0.2% water by mass. Using
0.9% and complete release gives about **9 mg** in the 10% stirred-batch
charge and **54 mg** from the packed bed's 6 g. The latter would be
equivalent to the water made at about **5% conversion** in one bed
void volume if released all at once. This is a possible startup
inventory. The actual starting moisture, release time and ongoing
water exchange with the solid are unresolved.
[Melgosa et al., Section 2.1](https://doi.org/10.1016/j.supflu.2014.11.003)

## How much wetter can the bead surface get?

Water made by the reaction has to leave the bead. The comparison
uses the previously checked packed-bed film correlation, a measured
water diffusivity in pentanol scaled with temperature, and the prior
film-free reaction rate as the source for a screening calculation.
The reference is a stirred-batch bead whose surface sees the bulk
liquid. If the stirred-batch film has the same resistance as the
bed film and the bulk liquids match, the two surfaces also match.
[Transport evidence](evidence.md#water-transport-and-rate-sensitivity)

| Water-transfer scenario | Largest surface activity rise over available points at 20–75% conversion |
| --- | ---: |
| Packed-bed correlation and central water diffusivity | 0.069 |
| Packed-bed correlation and slower water diffusivity | 0.077 |
| Isolated bead in stagnant liquid | 0.189 |

The two packed-bed peaks occur at the low-flow, 40 °C end near 20%
conversion. Including the approach from zero conversion raises the
packed-bed maxima to 0.099–0.110. These are conditional maxima on
the tested activity curves and use an ideal stirred reference.
The stagnant case shows how much stronger water buildup could be
with much weaker transfer. At **60 °C, 3 mL/min and 50% conversion**,
the central film gives **0.034 mol/L** of local water enrichment.
The two neat-pentanol activity curves translate this into **0.0045–0.0071**
additional activity. At 50 °C and the same flow and conversion,
the projections that also use measured ternary phase data give
roughly **0.006–0.015**. Across the liquid projections at 50 °C,
the absolute bulk estimates at 50% conversion span about 0.58–0.91.
The reactor-temperature four-component liquid activity surface is
unavailable. [Representative comparisons](reference/representative.csv) ·
[Measured binary activities](https://core.ac.uk/download/pdf/188223613.pdf)

Using the independent MTBE batch response to illustrate scale,
the central film and saturated headspace together could lower the
**local rate** by up to **23%** at supported points beyond 20%
conversion. The lower-diffusivity case raises this to **25%**.
At the matched 60 °C, 50% conversion example the analogue gives
about **1%**. These are rate sensitivities evaluated at the two
activities and the kinetic transfer to pentyl propanoate is
unverified. The batch law can already reflect the effect of the
water generated along its own reaction path. The 20 packed-bed
conversions stayed held out.
[Independent batch rates](https://doi.org/10.1002/cctc.202201207) ·
[Retained numerical results](reference/summary.json)

## Next step

Audit catalyst conditioning as a possible source of a persistent
bead difference. The stirred-batch experiment heated the beads in
pentanol before adding acid. The column hydration and priming history
are less well constrained. Published bead water uptake or
conditioning and reuse tests could support a targeted comparison
at the same bulk composition. The present results provide a scale
for how much activity the ordinary external film could contribute.
