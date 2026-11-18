# Film and water-activity protocol

## External liquid film

Use the 2025 full-inventory reconstruction of the 2:1 alcohol-to-acid
feed. The five temperatures span 40–80 °C. The flow rates are 0.7, 1, 2
and 3 mL/min. Retain 6.0 g of immobilised lipase, the 593 µm spherical
bead approximation, void fraction 0.64, the Table B1 internal effectiveness,
and the prior axial dispersion and boundary conditions.
[Frozen source conditions and kinetics](../../../literature/sources/lorenzo2025/author_model.json)

Use superficial liquid velocity $u_s=Q/A$, bead diameter $d_p$ and
liquid density $\rho$, dynamic viscosity $\mu$ and molecular diffusivity
$D$ in SI units. The ratio of momentum to mass diffusivity is
$Sc=\mu/(\rho D)$. The particle Reynolds number is
$Re=\rho u_s d_p/\mu$, so $Re\,Sc=u_s d_p/D$. Get the film
coefficient $k_f$ in m/s from $Sh\,D/d_p$.

Use Wilson–Geankoplis as the primary low-Reynolds packed-bed correlation.
Recompute the film using the published reactor correlation with the
corrected Schmidt number and use the analytic isolated-sphere
$Sh=2$ value as another low-transfer scenario. The Wilson–Geankoplis
property cases use the central, slower and faster diffusivity proxies
from the [source record](film/properties.json).

The nominal surface per liquid volume is $a_l=6(1-\varepsilon)/(d_p
\varepsilon)$. The film uses a common effective mobility for the
stoichiometric change in all four concentrations. At each axial
location solve

$$
60\, k_f a_l C_{A,0}(X_s-X_b)=\eta_I r(X_s).
$$

$X_b$ and $X_s$ are reaction progress in the bulk and at the bead
surface, respectively. Let $C_{A,0}$ and the rates use mol/L and
mol/(L min). The factor 60 puts film flux and reaction on the same
time basis. The internal factor $\eta_I$ acts once on the rate, and
the reversible reaction keeps its bulk equilibrium. This is a
steady, spherical-bead film estimate with fixed average properties
at each temperature and the whole nominal exterior wetted. Mixture
diffusivities, separate species mobilities and wetting remain uncertain.

Check a linear-reaction film solution, the high-transfer limit, local
and integrated acid balances, and the 512-to-1,024-cell change.
Require an absolute conversion change at most $10^{-4}$, relative
axial acid balance error at most $10^{-9}$, and the maximum inlet
film balance residual at most $10^{-10}$ mol/(L min). Write the forward
results before opening the measured conversions for scoring. Compare the
film penalty with the size and direction of the original conversion gap.

## Independent batch water activity

Use initial ester-formation rates from controlled water-activity batches.
Keep each study and solvent separate. Normalize by a measured rate
in the same series to preserve the published relative response. Use
endpoints for the simple interpolation
$\ln(r/r_{\rm low})=b(a_w-a_{w,\rm low})$ and score all intermediate
points. Treat this fit as descriptive over the measured interval.

Use the Tjørnelund MTBE endpoint slope to predict the changes in the
separate Dutta Banik MTBE series. Normalize the receiving series at its
own measured rate at $a_w=0.23$ and score the three higher activities
inside the donor range. Evaluate both the nominal and approximate
figure coordinates from the receiving paper. Report residuals as
fractions of the anchor rate. Alongside the controlled-activity
results, retain the Kuperkar solvent-free water-dose series as evidence
about added water with its original dose labels.
[Source conditions, locations and uncertainty](water_activity/evidence.md)

The target batch and bed water activities are unavailable. The 20
reactor conversions provide the film comparison only. The independent
batch analysis estimates the relative initial-rate response in each
published solution.
