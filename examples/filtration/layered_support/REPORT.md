# Layers in series: the limit of improving a support

## Research question

At fixed filtration flux and an unchanged selective layer, how much can a
support-permeability improvement reduce pumping energy? When does the fixed
layer dominate the membrane's hydraulic resistance? This first application
isolates that mechanism with synthetic straight pores and assumed layer ratios.

## Theory and assumptions

Clean Newtonian water crosses fully wetted layers of common membrane area.
Flow is steady, isothermal and linear; osmosis, fouling and gravity are absent.
Support thickness $L_s$ and axial permeability $k_s$ define an area-specific
resistance $R_s=L_s/k_s$ in m⁻¹. The model is

$$
R_{tot}=R_s+R_{sel}+R_{other},\qquad
\Delta p=\mu J R_{tot},\qquad Q=JA_m.
$$

$J$ is superficial membrane flux in m/s, $A_m$ is membrane area in m² and
$\mu$ is viscosity in Pa s. $R_{sel}$ denotes the stipulated selective-layer
resistance; $R_{other}=0$ in this initial sweep. An actual pump is assigned
the Julia pressure load and water mass flow in DWSIM. Its independent work
and specific-energy references are

$$
W=\frac{Q\Delta p}{\eta},\qquad
e=\frac{\Delta p}{\eta\,3.6\times10^6}
\quad[\mathrm{kWh/m^3\ filtrate}].
$$

Let $f=k_{s,new}/k_{s,0}$ at unchanged thickness and
$\alpha=R_{sel}/R_{s,0}$. The same absolute selective resistance is retained
for both supports; $\alpha$ is always defined against the baseline. Then

$$
1-\frac{W_{new}}{W_0}=\frac{1-1/f}{1+\alpha},\qquad
\lim_{f\to\infty}\left(1-\frac{W_{new}}{W_0}\right)=\frac{1}{1+\alpha}.
$$

The selective share is $\alpha/(\alpha+1/f)$. It exceeds 50% when
$\alpha>1/f$ and reaches 90% when $\alpha\ge9/f$. These are reporting
definitions of dominance, not material-specific empirical laws.

## Geometry and operating conditions

The two OpenLB support cells have radii 4/5 µm, a 20 µm square sample area
and 40 µm axial length. Julia infers permeability from native measured flow
and internal pressure taps. The analytical comparison is
$k=\pi r^4/(8A_{cell})$; it is not the simulation's production answer.

Upscaling assumes a 1 mm aligned-pore support at 100 LMH
($100\times10^{-3}/3600$ m/s), 1 m² area, 298.15 K and 75% pump efficiency.
DWSIM provides common water properties. Four assumed selective resistances
give $\alpha=0.25,1,10,100$. Neither the geometry nor these ratios is an
extracted experimental specimen. Enlarging the entire pore also enlarges its
surface opening and changes porosity; this is not yet a fixed-interface study.

## Results

The retained original study is run `20261108T192257Z-2d24de62`, using the
historical BGK-era native evidence. Its fine-grid permeabilities are
$2.511215119\times10^{-13}$ and $6.130870215\times10^{-13}$ m². At unchanged
thickness, support resistance falls by 59.0398%.

| $\alpha$ | Selective share, before → after | Actual DWSIM shaft-energy saving |
| --- | --- | ---: |
| 0.25 | 20.00% → 37.90% | 47.2319% |
| 1 | 50.00% → 70.94% | 29.5199% |
| 10 | 90.91% → 96.07% | 5.3673% |
| 100 | 99.01% → 99.59% | 0.5846% |

At $\alpha=1$, required pressure falls from 196.9128 to 138.7843 Pa and
specific shaft energy from $7.293065001\times10^{-5}$ to
$5.140158564\times10^{-5}$ kWh/m³. These small absolute values arise from
the selected synthetic resistance, not commercial membrane measurements.

Seven native cases check two grids and pressure/domain changes; eight actual
DWSIM configurations check work, balances and save/reload. Fine-grid
permeability/profile errors are below 0.1%, and inferred permeability changes
by less than 0.003% across the tested 1–20 Pa range. The largest DWSIM
enthalpy residual is 0.0916% of shaft power, below its 0.1% screen.

## Discussion

Improving support permeability can remove only the support's contribution.
When the fixed layer already supplies almost all resistance, a large pore-level
change yields a small energy saving. The original 59% resistance reduction
is therefore compatible with savings ranging from 47% to below 1%.

Here “total” covers the modelled transmembrane load. DWSIM contains a pump and
streams, without retentate, recycle or circulation. Shaft work excludes motor
losses and other plant utilities. Holding selective resistance fixed is a
model constraint, not evidence of unchanged selectivity after a larger opening.

The table preserves the original run; later porous calculations require TRT.
The subsequent [fixed-opening comparison](../graded_support/REPORT.md) controls
the mouths and was later revalidated with TRT for the explicit-channel study.
That later evidence is identified separately rather than relabelling this run.

## Conclusions

Series resistance explains both the diminishing return and the dominance
threshold. A defensible material-design comparison next needs fixed interface
geometry and an explicitly supported effective fixed resistance.

## Reproduction and evidence

[Commands](README.md) · [Original scientific record](../../../docs/filtration_research_case.md)
· [Compact reference](../../../docs/filtration_reference.json).
This synthetic sweep uses analytical verification, not experimental calibration;
the [literature-calibration example](../calibration/REPORT.md) treats that separately.
