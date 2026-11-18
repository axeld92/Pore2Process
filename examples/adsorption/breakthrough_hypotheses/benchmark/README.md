# Radial water flow in a packed tube

Water flowed preferentially through the outer part of Eghbalmanesh et al.'s
packed bed. At a particle Reynolds number of 10, the outer half particle
diameter carried 37.04% of the flow through 26.53% of the cross-sectional
area. Two independent Ergun zones predicted 55.23% through this region when
given its measured voidage. Resolving the radial voidage and including
no-slip wall friction gave 35.26%. The benchmark tests the hydraulic
assumption behind H1. Adsorption and interchange of liquid between the
annulus and core require separate evidence. [Eghbalmanesh et al. (2024), Fig. 2](https://doi.org/10.1002/aic.18322)

## Reference and calculation

The reported tube had a 21 mm inside diameter and a length of 61 mm. It held
consolidated polypropylene spheres of diameter 3.00 ± 0.05 mm. Water
at 20 °C flowed at 69.3, 207.8 and 346.4 mL/min, corresponding to particle
Reynolds numbers 10, 30 and 50 based on the empty-tube velocity. The authors
used 1000 kg/m³ and 1.00 mPa s in their
particle-resolved calculation and averaged in the packed interior over
15.875–40.875 mm. The experimental reference is the boundary-postprocessed
MRI curve (BPP-MRI). Section 5.1.3 quantifies uncertainty for raw velocity.
An uncertainty for the processed cumulative curve is unavailable.
[Article, §§3.1–3.3 and 5.1–5.2](https://doi.org/10.1002/aic.18322)

The archive contains `Figure_2_a.csv`, `c` and `e` with radial superficial
velocities and voidage, and `b`, `d` and `f` with normalized cumulative flow.
The coordinate $S=(R-r)/d_p$ increases from wall to centre. The archived
cumulative curve is the flow inside the radius $R-Sd_p$. The complementary
fraction is the flow within the wall annulus of thickness $Sd_p$. The
annulus occupies cross-sectional area fraction $1-(1-Sd_p/R)^2$.
[Original data, Fig. 2 CSVs, CC BY-NC 4.0](https://doi.org/10.4121/4b8d0613-5ee7-4971-b695-835f287336e8)

The calculation independently integrates each published MRI superficial
velocity sample over its cylindrical annulus to check the flow convention
and normalization. The two-zone Ergun sensitivity averages the published
MRI voidage separately inside and outside the outer half particle diameter.
The extended Brinkman equation uses the same Ergun drag, symmetry on the
axis and no slip at the wall. The effective momentum viscosity is
$\mu_\mathrm{eff}/\mu=2\exp(0.0035\,\mathrm{Re}_p)$.
The two-zone and distributed-voidage calculations assess wall friction and
radial structure in separate comparisons. All solve for a common axial
pressure gradient at the reported total flow using measured voidage and
published closure coefficients.
[Giese et al. (1998), liquid-flow experiments](https://doi.org/10.1002/aic.690440225)
· [Kwapinski et al. (2004), Eqs. 2–6](https://doi.org/10.1002/ceat.200407001)

## Verified results

Area and volume-flow shares refer to this specimen's annulus.
MRI and PR-CFD values are linearly interpolated from the archived normalized
curves. The PR-CFD curve is renormalized by its value at the wall, about
0.9996 of nominal flow.

| Particle Reynolds number | $S$ | Area share (%) | MRI (%) | PR-CFD (%) | Two-zone Ergun (%) | Brinkman–Ergun (%) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | 0.25 | 13.78 | 28.62 | 30.69 | 28.68 | 29.06 |
| 10 | 0.50 | 26.53 | 37.04 | 37.33 | 55.23 | 35.26 |
| 10 | 1.00 | 48.98 | 53.77 | 55.08 | 68.91 | 55.62 |
| 30 | 0.50 | 26.53 | 35.99 | 37.32 | 43.43 | 31.34 |
| 50 | 0.50 | 26.53 | 35.27 | 37.37 | 46.91 | 36.55 |

| Root-mean-square error against processed MRI, percentage points | Re 10 | Re 30 | Re 50 |
| --- | ---: | ---: | ---: |
| Independently integrated MRI velocity | 0.12 | 0.13 | 0.13 |
| Authors' particle-resolved simulation | 1.12 | 0.92 | 1.64 |
| Two-zone Ergun | 9.08 | 4.96 | 6.75 |
| Two-zone Brinkman–Ergun | 6.77 | 4.39 | 5.63 |
| Brinkman–Ergun | 1.10 | 2.45 | 1.44 |

The independently integrated mean velocity was 1.0013–1.0015 times the
nominal inlet velocity before normalization. The integration is a check
on two representations of the same measurements. The differences between
PR-CFD and MRI compare models with data and carry no experimental confidence
interval. The benchmark computes its root-mean-square errors at the MRI
sample locations within the actual radius. The authors' Table 3 reports
1.13, 0.94 and 1.66 percentage points for the same cumulative comparison.

The numerical checks reproduce the analytical cumulative flow in a
clear no-slip circular pipe within $8\times10^{-7}$ in flow fraction
and the uniform-voidage Ergun solution within $10^{-10}$ in normalized
velocity. Doubling the radial mesh from 160 to 320 rings per particle
diameter changes any reported Brinkman flow fraction by at most 0.011
percentage points. The force residual is below $4\times10^{-11}$ of the
pressure gradient at the finer resolution. Adding wall friction at the
two-zone voidages reduces the error at all three Reynolds numbers. The
largest improvement comes from also resolving the radial voidage.

## Geometry sensitivity at Alhamed's conditions

The 11 mm tube and 0.80 mm particle size give $D/d_p=13.75$. The outer
half diameter is 0.40 mm thick and occupies 14.02% of the cross section.
At this diameter ratio, independent spherical-bed studies establish the
presence of radial voidage variation. The two shapes below are scenario
inputs for Alhamed's carbon bed. The smooth profile follows Kwapinski
et al.'s expression for nearly spherical packings with some variation in
size and shape. The oscillatory profile follows de Klerk's correlation for
equal spheres as reproduced in Das et al.'s Eq. 15. Their Fig. 6 includes
packings at $D/d_p=12$ and 16. The sphere data show high voidage at the
wall, a trough about half a particle diameter in, and a second peak near one
particle diameter. The Alhamed carbon shape, size distribution and local
voidage are unavailable. [Benenati & Brosilow (1962)](https://doi.org/10.1002/aic.690080319)
· [Sederman et al. (2001), MRI of glass-sphere packings](https://doi.org/10.1016/S0032-5910(00)00374-0)
· [de Klerk (2003)](https://doi.org/10.1002/aic.690490812)
· [Das et al. (2017), Eq. 15 and Fig. 6](https://pure.tue.nl/ws/files/56182313/dasDNS2017.pdf)

The scenario sets the core voidage by preserving the original global void
fraction of 0.427 and solid inventory. With constant solid density, wall
solid volume fraction sets its carbon mass share. The assumed water viscosity
of 0.833 mPa s gives particle Reynolds numbers of 3.92 and 10.52 at the
two rates. Scenario inputs come from geometry, assumed fluid properties
and hydraulic closures.

| Assumed porosity shape | Mean voidage in outer 0.40 mm (%) | Carbon in outer 0.40 mm (%) | Brinkman flow share at 23.3 mL/min (%) | At 62.5 mL/min (%) |
| --- | ---: | ---: | ---: | ---: |
| Oscillatory | 55.04 | 11.00 | 20.14 | 20.26 |
| Smooth | 59.74 | 9.85 | 31.74 | 29.78 |

These scenarios give the wall annulus a larger share of water flow than of
carbon. Two-zone Ergun on the same mean voidages gives 33.14–47.83% wall
flow across the two rates and shapes. This helps interpret the high-flow
sensitivity as omission of wall friction and radial structure. The 20–32%
Brinkman range is conditional on the two adopted shapes and the shared wall
closure. Persistent bypass depends on exchange between flow paths, which a
radial average does not determine. The confined stock-TRT and tracer
experiments address that question separately.

## Reproduce

Run from the project root. The benchmark was verified with Python 3.12.9
and Julia 1.11.6.

```sh
source scripts/env.sh
uv run --locked python examples/adsorption/breakthrough_hypotheses/benchmark/fetch.py
julia --startup-file=no examples/adsorption/breakthrough_hypotheses/benchmark/radial.jl
```

The Python script downloads three bounded HTTP ranges totalling 194,439
bytes from the published archive and checks the CRC-32 of each extracted
record. Use `--archive /path/to/SupplementaryMaterial.zip` with a local copy
of the full archive. Both scripts keep inputs and generated CSVs under
`results/adsorption/breakthrough_hypotheses/benchmark/` by default. To use
another directory, pass `--output DIRECTORY` for Python and `DIRECTORY` as
the Julia script's argument.

`radial_benchmark.csv` records the cumulative flow comparisons and
`benchmark_checks.csv` records normalization, residuals and grid changes.
`near_wall_scenarios.csv` records the assumed carbon-bed profiles' effects.
The archive also includes the reconstructed particle centres and radii.
`sphere_coordinates.csv` converts them to metres using the reported 25 mm
field of view over 128 pixels. Its axial `x` uses the archived slice number.
Transverse `y` and `z` use the row and column coordinates relative to the
reported tube centre. The sheet supplies 758 particle rows (the article's
Table 4 lists 800), with a 20.9766 mm tube diameter and sphere extents from
6.525 to 58.325 mm along the archived axis. `geometry.json` records these
values and 76 small reconstructed wall intersections of at most 6.9 μm.
The comparison above uses the published radius and profiles. A separate
particle-resolved OpenLB comparison would require explicit treatment of
those geometric details and a grid check. [Archive, coordinate workbook](https://doi.org/10.4121/4b8d0613-5ee7-4971-b695-835f287336e8)
