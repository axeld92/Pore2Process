# Batch evidence for the water-activity response

Two published batch studies independently set the initial water activity of
commercial immobilised CALB and measured initial ester formation. In both
studies butyric acid reacted with ethanol to make ethyl butyrate in methyl
*tert*-butyl ether (MTBE). From the lowest to the highest water activity, the
rate fell by 83% in one study and 58% in the other. The solvent, alcohol, acid,
and concentration differ from the solvent-free pentyl propanoate reaction.
[Tjørnelund et al. (2023)](https://doi.org/10.1002/cctc.202201207) · [Dutta Banik et al. (2017)](https://doi.org/10.3390/catal7080227)

## Controlled batch rates

[Batch rate table](batch_initial_rates.csv) preserves the published
activity numbers, units and uncertainty where tabulated. The water activity
column gives nominal salt setpoints on a dimensionless scale. Relative activity
is the mean specific activity divided by the mean at the lowest setpoint in
the same study and solvent. The ratios avoid transferring absolute rates across
catalyst normalizations and substrates. Blank entries indicate unavailable
tabulated data or no separate reading from a plot.

### Tjørnelund et al. (2023)

- Table 1 and Figure 4a give duplicate-run specific activities for Novozym 435.
  The CSV includes all three solvents for this immobilised preparation. Keep
  the solvent series separate. The response at $a_w=0.97$ relative to 0.12 is
  0.167 in MTBE, 0.385 in hexane and 0.500 in acetonitrile.
- Initially 50 mM butyric acid and 50 mM ethanol in each solvent. Table 2
  gives 3 mL and 1 mg enzyme protein in either MTBE or acetonitrile, and
  5 mL and 0.5 mg protein in hexane. The authors assumed 10% protein by
  mass for the immobilised formulation. Reaction stirring was 800 rpm at
  room temperature. An exact reaction temperature is unavailable.
- Pre-equilibration was at room temperature for 24 h with stirring at
  100 rpm and gas-phase humidity tracked with a hygrometer. The nominal
  salt setpoints were 0.12 (LiCl), 0.33 (MgCl₂), 0.53 (Mg(NO₃)₂·6H₂O),
  and 0.97 (K₂SO₄). Ethanol was pre-equilibrated separately. The authors
  replaced evaporated solvent with pre-equilibrated solvent before adding
  ethanol. Karl Fischer data demonstrate that changing the solvent changes
  how much water is present at a given activity.
- Specific activities come from a linear regression of ester formation
  at up to 10% conversion, scaled by volume and divided by enzyme protein
  mass. The reported ± values are preserved. The main text gives duplicate
  experiments without specifying how the plotted error bars were calculated.
  Use the rates for within-series ratios and retain the absolute units as
  printed. For example, a literal 26 mmol/min/mg with 0.5 mg protein would
  consume 10% of the 0.25 mmol starting acid within 0.12 s, which makes the
  printed absolute scale difficult to reconcile with the sampling method.
  [Article, Tables 1–2 and Experimental Section](https://doi.org/10.1002/cctc.202201207) · [DTU copy of the published PDF](https://backend.orbit.dtu.dk/ws/portalfiles/portal/304333369/ChemCatChem_2022_Tj_rnelund_Optimised_Protocol_for_Drying_Aqueous_Enzyme_Solutions_in_Organic_Solvents_Comparison_1_.pdf)

### Dutta Banik et al. (2017)

- Initially 25 mM each of butyric acid and ethanol in 3 mL MTBE, approximately
  10 mg Novozym 435, 30 °C, 800 rpm and two independent sets at each condition.
  Medium and immobilised catalyst were separately pre-equilibrated overnight
  above saturated salt solutions. The listed salts were LiCl, potassium
  acetate, MgCl₂, Mg(NO₃)₂·6H₂O and K₂SO₄. Initial slopes were measured by
  GC-FID at below 10% conversion with samples at one-minute intervals.
- Table 1 lists the five experimental rates alongside a separate series of
  molecular simulation settings. Rate values in the CSV are from that table.
  Their ordering and correspondence to experimental setpoints were checked
  against Figure 1 and Section 3.1.2. The table's $a_w$ column belongs to
  the simulation conditions in Section 3.2.
- We rendered the published Figure 1 from the DTU PDF and read marker
  centres along the horizontal axis, using the ticks from 0 to 1, to
  approximately 0.01. The plotted first and third activities are 0.11 and
  0.38, compared with Methods setpoints of 0.12 and 0.33. Preserve both
  coordinates for sensitivity to this unresolved source discrepancy. Figure 1
  shows error bars. Their approximate vertical half-spans in rate units at the
  displayed resolution are 0.010, unresolved, 0.009, 0.019 and 0.010 in order.
  These are readings from the graphic and the paper does not specify the
  spread statistic. Table 1 and Figure 1 report `mM/min/mg`, while the Methods
  paragraph gives `mmol/min/g`. The within-series ratios avoid the discrepancy
  in the printed dimensions. [Article, Figure 1, Table 1 and Methods](https://doi.org/10.3390/catal7080227) · [DTU copy of the published PDF](https://backend.orbit.dtu.dk/ws/portalfiles/portal/137397759/catalysts_07_00227.pdf)

## Water addition in a solvent-free propanoate analogue

Kuperkar et al. used Novozym 435 to make isobutyl propanoate from propanoic
acid and isobutanol. Section 2.2 gives a 50 mL baffled stirred vessel, an
equimolar base charge of 0.1 mol each reactant (7.5 mL acid and 9.3 mL
alcohol), 200 rpm and 5% w/w catalyst. Section 3.1.1 gives 40 °C as the
base test temperature. Figure 4 and Section 3.1.4 do not restate the fixed
conditions for the water series. The experiment used water doses, with no
inert organic solvent. Water activity and residual moisture at zero added
water are unreported. [Kuperkar et al. (2014), Sections 2 and 3.1.4](https://doi.org/10.1016/j.molcatb.2013.10.024)

[Water-addition data](batch_water_addition.csv) record the Figure 4B initial
rate response. The rates at zero and the largest dose have relative values
of 1 and 0.19. The 0.2 entry is approximately 0.027 mol/min/g as read from
the 0 to 0.07 y axis on the rendered article page 147. The other three rates
are printed in the text and agree with the plot. Figure 4 supplies no
rate error bars. Section 2.3 describes low-conversion slopes and about ±3%
agreement between its two analytical methods. That comparison provides no
replicate uncertainty for the rates. [Kuperkar et al. (2014), Figure 4 and Section 2.3](https://doi.org/10.1016/j.molcatb.2013.10.024)

The authors label the water axis in M, but also equate 0.1 and 0.3 on that
axis to added volumes of 1.8 and 5.4 mL. Those doses correspond to approximately
0.1 and 0.3 mol of water. That amount would require a reference volume of 1 L
to match the listed molarities, while the reactor capacity is 50 mL. Keep the
water-axis labels and specified dose volumes as reported. Added water can
change reactant concentrations and phase behavior along with enzyme
hydration. The series supports a qualitative rate sensitivity to a wetter
initial mixture. It cannot supply a concentration-normalized inhibition
constant or a water activity curve for the pentanol bed.
[Kuperkar et al. (2014), Sections 2.2 and 3.1.4](https://doi.org/10.1016/j.molcatb.2013.10.024)

## Implications for pentyl propanoate

The initial-rate design measures an apparent kinetic response to water activity
before appreciable ester accumulates. It leaves the relative contributions of
active-site hydration, support wetting and local substrate availability
unresolved. A reduction in eventual equilibrium conversion is a separate
effect. The between-solvent differences show that even one enzyme formulation
and one reaction can have several rate responses to the same change in water
activity. [Tjørnelund et al. (2023), Figure 4](https://doi.org/10.1002/cctc.202201207)

The same commercial lipase gave a rate maximum at $a_w=0.45$ for propanoic acid
and 1-octanol. This closer acid analogue supports a possible rate increase with
hydration on the low-activity side. The accessible abstract reports the
optimum without a numerical rate series. It supplies a constraint on assuming
monotone inhibition and no coefficient for the pentanol reaction.
[Nordblad and Adlercreutz (2008), pp. 1518–1524](https://doi.org/10.1002/bit.21706)

The pentyl propanoate batch study varied loading, temperature and feed ratio.
Water in those runs was estimated by mass balance and changes with conversion.
The packed-bed study dried the feed but did not publish its residual water or
bead water activity. The studies therefore leave the relationship between
local water mole fraction and activity, $a_w=\gamma_w x_w$, unmeasured for the
reacting liquid and immobilised enzyme. The independent batch ratios can be
used for a sensitivity study over their measured initial activity ranges.
At a given conversion, a batch and an ideal plug-flow bed starting from the
same feed share the stoichiometrically generated water per mole fed. A rate
law already calibrated on the batch can absorb its own water response. A
bed-specific rate correction needs a difference in the local water activity
or catalyst hydration histories at the same conversion.
A quantitative kinetic transfer would need published matched pentanol and
propanoic-acid rate measurements at controlled initial activity and a water
activity model for the reacting mixture at bed conditions.
[Lorenzo and Ortega (2026), Table 1 and SI Table S2](https://doi.org/10.1021/acs.iecr.6c02375) · [Lorenzo et al. (2025), Sections 2.1–2.3](https://doi.org/10.1021/acs.iecr.5c02780)
