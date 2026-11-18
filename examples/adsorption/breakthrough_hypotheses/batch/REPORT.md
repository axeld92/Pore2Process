# Particle uptake and equilibrium in the phenol breakthrough test

Varying the particle diffusivity within the tested batch-qualified range
did not fix the two 25 cm breakthrough curves. Slower particle uptake moved
the first leakage in the right direction at high flow, while reducing the
amount of phenol captured. Faster particle uptake captured more phenol but
delayed first leakage further. Both alternatives still captured too little
phenol after the published curve crossed 50% outlet/feed.

## Published equilibrium and carbon inventory

Alhamed's Table 1 reports a packed bulk density of 740 kg/m³, a voidage
of 0.427, and two equilibrium fits. The nominal tube has an 11 mm inside
diameter in Section 2.3 and Table 2. The corresponding dry masses are
7.032, 10.549, 14.065 and 17.581 g for the 10, 15, 20 and 25 cm beds.
These are density times bed volume. The paper does not list separately
weighed masses or voidages for the plotted runs. The derived dry density
of the particle envelope, 1291.45 kg/m³, uses the same reported voidage
and is kept fixed in the bath and column tests. [Alhamed, Tables 1–2 and §2.3](https://doi.org/10.1016/j.jhazmat.2009.05.002)

The equilibrium laws use dissolved concentration `C` in mg/L and
equilibrium loading `q` in mg phenol/g dry carbon.

```text
Langmuir    q(C) = 111.5 × (0.139 C) / (1 + 0.139 C)
Freundlich  q(C) = 24.6 × C^(1/2.74)
```

Each law was tested with the same finite-bath balance
`V C_end + m q(C_end) = V C_initial`. All seven Table 3 fitted endpoint
loadings are within 2.03% with Langmuir and 1.12% with Freundlich.
The calculated residual liquid concentrations span 5.56–33.28 mg/L
with Langmuir and 6.39–36.02 mg/L with Freundlich. At 150 mg/L column
feed the laws give 106.40 and 153.15 mg/g, respectively. Table 3 therefore
provides an equilibrium check at low residual liquid concentration. Its
200 and 400 mg/L values are starting concentrations and cannot distinguish
the equilibrium loadings at 150 mg/L. A column marker at the end of a
finite observation window does not establish a saturation capacity.
For the two 1 g, 250 mL bath concentrations, the printed Freundlich
fit and finite-bath balance give 48.40 and 90.99 mg/g. Section 3.1.1
quotes 48.7 and 97.0 mg/g from its Freundlich calculation. The audit
uses the printed coefficients and checks the independent mass balance
against the Table 3 fitted endpoints.
[Alhamed, Tables 1 and 3 and §3.1.1](https://doi.org/10.1016/j.jhazmat.2009.05.002)

Table 1 cites Y. A. Alhamed, *Phenol removal using granular activated
carbon from dates’ stones*, Bulgarian Chemical Communications, pp. 26–35
(2008). Its original concentration range and isotherm points were
unavailable in the sources inspected. The 2009 reference prints volume
41 and the journal indexes the title under volume 40. Neither equation
was refitted for the column. [Journal author index](https://bcc.bas.bg/BCC_Volumes/Volume_40_Number_4_2008/Volume_40_Number_4_2008_PDF/Index-vol40.pdf)

### Capacity sensitivity using the reported GAC4 plateau

Pasupulety et al. (2019, DOI 10.34049/bcc.51.3.5070) report sample GAC4
with a ZnCl₂-to-dry-stone mass ratio of 0.5 and carbonization at 700 °C
for 3 h. Its BET area of 951 m²/g, total pore volume of 0.456 cm³/g and
micropore volume of 0.355 cm³/g match Alhamed (2009) Table 1. Both papers
cite Alhamed (2008) for the preparation or characterization. The
independence of the carbon batches is unconfirmed.
[Pasupulety et al., experimental section and Table 1](https://bcc.bas.bg/BCC_Volumes/Volume_51_Number_3_2019/BCC-51-3-2019.pdf)

Table 2 gives a fitted Langmuir phenol capacity of **103 mg/g** for GAC4.
The 1.50 mm phenol baths used 0.1 g of carbon in 25 mL with initial
concentrations of 10–150 mg/L. This supplies at most 37.5 mg phenol per g
carbon. The 103 mg/g plateau is therefore an extrapolation of the fit.
The paper reports 25 °C and 24 and 48 h of shaking for the control.
The published tables supply the capacity. The paired equilibrium
concentrations and loadings, phenol affinity coefficient, and fit
uncertainty for GAC4 remain unavailable. The individual article and
full-issue indexed texts agree on the capacity.
[Pasupulety et al., adsorption measurements and Table 2](https://bcc.bas.bg/BCC_Volumes/Volume_51_Number_3_2019/BCC-51-3-2019.pdf)

Fixing the plateau at 103 mg/g, the original Langmuir maximum relative
Table 3 endpoint difference of 2.032% gives a common affinity interval
of 0.13409–0.21347 L/mg. For each endpoint, the two loading limits give
residual concentrations by the finite-bath mass balance, then affinity
limits through `b = q/[C_end(q_max - q)]`. This interval is a
consistency check against the 2009 fitted endpoint values. Its midpoint,
0.173782 L/mg, gives a maximum endpoint difference of 0.492% and a loading
of 99.19 mg/g at 150 mg/L. With the original diffusivity and stirred film
coefficient, Figure 2 loading RMSEs are 2.34, 1.56, 1.40 and 0.51 mg/g
from largest to smallest particle size. The largest change with shell
or time refinement is 0.031 mg/g.

The capacity sensitivity used the saved fine confined flows, separate
wall/core carbon histories, and original carbon inventories and column
film coefficients. All three forecasts passed the numerical refinement
checks at 48 shells and a 4 s main step. They were sealed before scoring
Figures 4–5. The two low-flow cases passed timing and shape and retained
too little phenol. The high-flow 5% arrival remained close to the paper
and its midpoint arrived too early. All three final retention checks
missed. The capacity sensitivity yields another batch-consistent
equilibrium curve with a different column prediction.

| Held-out bed and flow | 5% crossing (min) | 50% crossing (min) | Calculated retention (mg/g) | Published retention (mg/g) |
| --- | ---: | ---: | ---: | ---: |
| 20 cm, 23.3 mL/min | 58.00 | 169.19 | 51.03 | 56.37 |
| 25 cm, 23.3 mL/min | 102.06 | 223.21 | 50.56 | 57.13 |
| 25 cm, 62.5 mL/min | 16.00 | 54.26 | 47.34 | 55.67 |

Replay the saved sensitivity from the checkout root. Its summary SHA-256
is `ee6b40fc9c2217892f9b88589bcbc887476c36b5bb7b120b82c313b14db68756`.

```sh
source scripts/env.sh
uv run python -m examples.adsorption.breakthrough_hypotheses.run --check \
  results/adsorption/breakthrough_hypotheses/gac4_sensitivity/20261118T010430Z-60cb2086
```

## Langmuir batch check

The Figure 2 bath started at 200 mg/L, with 1 g carbon in 250 mL liquid,
and was shaken at 200 rpm at 22 ± 1 °C. Only 1.50 and 0.80 mm markers informed
the kinetic searches. The 0.45 and 0.225 mm curves remained held out.
All cases use the same inverse-diameter rule for the stirred film.
Figure 2 starts at 30 mg/g, so neither small-particle curve resolves early
uptake well. The first digitized points for 0.45 and 0.225 mm are
45.87 mg/g at 71.31 min and 47.98 mg/g at 76.40 min, about 94% and 99%
of the independently computed Langmuir bath loading. The complete
1 mL sampling schedule and pointwise measurement uncertainty
are unavailable. Table 3's `h` values are fitted extrapolations.
[Alhamed, §2.2, Table 3 and Figure 2](https://doi.org/10.1016/j.jhazmat.2009.05.002)

| Profile | Loading diffusivity (10⁻¹² m²/s) | Stirred film at 0.80 mm (µm/s) | 1.50 mm RMSE (mg/g) | 0.80 mm RMSE (mg/g) | 0.45 mm RMSE (mg/g) | 0.225 mm RMSE (mg/g) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Archived fit | 1.8460 | 77.32 | 2.100 | 1.621 | 1.327 | 0.420 |
| Lower diffusivity | 1.5824 | 1000 | 2.188 | 1.693 | 1.313 | 0.420 |
| Higher diffusivity | 2.4577 | 45 | 2.700 | 1.572 | 1.347 | 0.420 |
| Lower-film check | 2.9193 | 40 | 3.083 | 1.563 | 1.355 | 0.420 |

The per-size center-marker screen is 3 mg/g. The 40 µm/s profile
failed on 1.50 mm and was excluded. The higher-diffusivity profile
passes at the nominal marker coordinates. Moving the marked positions
by the declared two-pixel sensitivity produces a 1.50 mm RMSE range
of 1.96–3.87 mg/g. The two extrema are deterministic sensitivities under
the chosen model and threshold. The reading range describes image
digitization. Doubling radial shells or halving the step changed the
scored loading by at most 0.39% of reference equilibrium loading over
all four profiles. [Batch protocol](protocol.md)

## Transfer to the two 25 cm beds

Each forecast retains 17.581 g of carbon, the original liquid volume,
the saved N60 OpenLB TRT flow and contact rules, and its separately
estimated liquid film and slow-cohort access resistance. Changing
the loading diffusivity leaves the packed-column film fixed.
The stirred-bath film coefficients do not set the column film.
The curves were saved and hashed before the Figure 4–5 comparison.
[N60 column reference](../../porosity_matched_column/REPORT.md)

The crossing times are when the outlet reaches 5% or 50% of the feed.
Last-marker retention includes liquid and carbon inventory. The last
column integrates additional phenol retained from the **published**
50% crossing to that curve's last marker.

| Flow (mL/min) | Curve | 5% time (min) | 50% time (min) | Outlet/feed RMSE | Retained (mg/g) | Retained after published midpoint (mg/g) |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 23.3 | Figure 4 | 125.74 | 274.40 | — | 57.13 | 10.71 |
| 23.3 | Archived N60 | 143.50 | 236.51 | 0.060 | 53.88 | 8.43 |
| 23.3 | Lower diffusivity | 132.83 | 228.23 | 0.071 | 52.93 | 8.71 |
| 23.3 | Higher diffusivity | 162.00 | 249.08 | 0.066 | 55.20 | 7.77 |
| 62.5 | Figure 4 | 16.16 | 71.63 | — | 55.67 | 26.08 |
| 62.5 | Archived N60 | 23.16 | 59.23 | 0.057 | 50.51 | 21.18 |
| 62.5 | Lower diffusivity | 20.52 | 54.35 | 0.073 | 49.02 | 21.00 |
| 62.5 | Higher diffusivity | 28.62 | 68.19 | 0.067 | 52.79 | 21.12 |

The lower diffusivity still gives a 27.0% late 5% crossing at high
flow against the 25% screen. Its midpoint is 24.1% early against the
20% screen. The higher diffusivity improves the midpoint and total
retention at high flow, but delays the 5% crossing by 77.1% and captures
no additional phenol after the published midpoint. The Figure 5
duplicate gives the same low-flow conclusions. Every transfer passed the
24-to-48-shell and 8-to-4-second numerical comparison, with all nominal
crossings moving by less than 0.74%. [Alhamed, Figures 4–5](https://doi.org/10.1016/j.jhazmat.2009.05.002)

## Published Freundlich sensitivity

With the printed Freundlich law, requalifying the same film and spherical
particle model gives a loading diffusivity of 1.23760 × 10⁻¹² m²/s.
The stirred-film estimate is 1000 µm/s at the upper search bound, so
the batch test does not establish a unique film coefficient.
The four batch errors are 0.765, 0.773, 1.036 and 0.241 mg/g in order
of decreasing particle size. Both withheld size curves pass, and
the largest numerical refinement change is 0.152% of the reference
finite-bath equilibrium loading. No equilibrium coefficients were
fitted. [Freundlich protocol](freundlich_protocol.md)

| Flow (mL/min) | Curve | 5% time (min) | 50% time (min) | Outlet/feed RMSE | Retained (mg/g) | Retained after published midpoint (mg/g) |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 23.3 | Figure 4 | 125.74 | 274.40 | — | 57.13 | 10.71 |
| 23.3 | Freundlich prediction | 151.87 | 323.33 | 0.091 | 64.40 | 15.48 |
| 62.5 | Figure 4 | 16.16 | 71.63 | — | 55.67 | 26.08 |
| 62.5 | Freundlich prediction | 21.64 | 76.06 | 0.037 | 61.59 | 31.08 |

At high flow this change improves the midpoint and the overall curve
error, yet the 5% crossing is still 34% late. It overshoots total
retention by 7.27 and 5.93 mg/g, and uptake after the published
midpoint by 4.78 and 5.00 mg/g at low and high flow. The two published
isotherms with independently qualified kinetics bracket the
observed end-window and late uptake. The equilibrium choice
therefore matters to the discrepancy. The tested batch conditions
provide weak discrimination at the higher concentrations in the
bed, and the two breakthrough predictions each miss part of the
measured behavior. The Freundlich runs pass the same 48-shell,
4-second refinement, with dosage and packed-bed inputs unchanged.
[Alhamed, Figures 4–5](https://doi.org/10.1016/j.jhazmat.2009.05.002)

[Commands](README.md) · [Compact verification record](reference_summary.json) · [Conditional shared depth transfer](conditional/README.md)
