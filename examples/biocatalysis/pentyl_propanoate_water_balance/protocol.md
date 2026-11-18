# Matched-conversion water protocol

Compare at equal temperature, pressure and molar acid conversion
$X$ from one mole of acid and two moles of pentanol. For dry feed
and retained products, the moles of acid, alcohol, water and ester
are $(1-X,2-X,X,X)$. Both bulk liquids then have water mole fraction
$X/3$ and the same equilibrium water activity, independent of the
activity-coefficient model. In two liquid phases at equilibrium,
the water chemical potential is common across the phases.

Use $X=0$–0.75 on a 0.005 grid across the frozen five temperatures,
four flows and 6.0 g bed. Compare film fluxes at the same bulk
conversion. Use the reported conditions and batch rate law only.
The 20 measured bed conversions stay held out.

## Bulk water inventory

The nominal finite batch headspace can hold at most
$n_{w,g}^{\rm sat}=V_g p_w^*(T)/(RT)$ of water. This uses a closed
vapor space at saturation with all volume at the liquid temperature.
Subtract at most the smaller of this inventory and the water present
to screen a batch bulk activity depression. Water retained in the
external condenser and sustained loss through a vent belong to
separate balances. Add the independent fresh-catalyst moisture
inventory as a separate complete-release sensitivity. For the
packed bed, compare fresh-bead moisture with one fluid-filled
void volume to show the scale of a possible startup pulse.
A steady liquid inventory balance excludes ongoing accumulation.

## Activity models

Use the hot binary VLE as a composition surrogate for a liquid
with all organics represented by pentanol. Evaluate the full and
lower-repeat-spread curves only inside available mole fractions.
For 40 and 50 °C also compare projections of measured ternary
tie lines at 30 and 50 °C, interpolating activity at 40 °C.

For each ternary liquid below the organic binodal, intersect the
fixed dry acid-to-alcohol ratio with the organic endpoints.
Interpolate between zero water activity and saturation using both
a linear curve and the normalized shape of the hot binary curve.
Use the water mole fraction of the paired aqueous phase as the
approximate saturation activity. In the measured two-phase region,
reuse the prior tie-line and lever-rule interpolation and use the
paired water-rich composition as the common activity estimate.
Run each method with the ester excluded and with the ester
counted as pentanol. These are sensitivity calculations based on
measured lower-order systems. A measured four-component activity
surface is unavailable.

## Water produced at the bead surface

At each matched $X$, use the frozen film-free Julia rate with the
reported internal effectiveness to screen the water production
rate. Verify that this rate decreases over the matched progress
grid. Treat the full nominal bead exterior as wetted. Let $a_l$
be this area per liquid volume and $k_w$ the water film coefficient.
On the mol/L/min rate basis,
the steady water concentration rise in mol/L is

$$
C_{w,s}-C_{w,b}=\frac{r_{\rm water}}{60 k_w a_l}.
$$

Take the stoichiometric bulk composition and add this increment
to water, holding the three organic amounts fixed. This isolates
the water buildup signal. The comparison uses the ideal stirred
reference with its bead surface at the bulk liquid activity.
Equal films with equal bulk inventories give zero difference.
The water-specific calculation leaves relative organic gradients,
batch film strength, and intrabead hydration history unresolved.

For an order-of-magnitude kinetic comparison, take the ratio of
independent batch initial rates interpolated in log rate at the
two computed activities. Apply the interpolation only where
both values fall inside that batch data's water-activity range.
Use it as a differential sensitivity of the apparent rate.
Check the headspace calculation against a pure-water reference,
and test a same-tie-line total composition and the zero-offset
rate ratio.
