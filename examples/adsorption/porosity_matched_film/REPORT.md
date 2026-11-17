# Clean-wall phenol transfer at 42.7% void fraction

At 62.5 mL/min, the finer run projects **0.110%** escape from the clean
25 cm bed when scaled by carbon area and **0.109%** when scaled by
liquid volume. The high-flow conclusion stays below 5% on both grids
and both finer-grid time steps. The film conductance remains
grid-sensitive.

The [60% void fraction study](../disordered_film_conductance/REPORT.md)
gave different clean-film estimates depending on how it scaled the
result to the column. The present packing matches the column's nominal
42.7% void fraction and keeps the same 32 spheres of 0.80 mm diameter.
It reduces the transverse area and preserves an approximately 11.1 mm
axial sample. The 25 cm model keeps its fixed 17.581 g inventory and
0.10210 m² of external carbon area.

Separate native runs at each interstitial velocity use stock OpenLB
1.9r0 D3Q27 TRT, Kupershtokh forcing and BounceBack. The OpenLB geometry
uses periodic orthogonal boundaries and the stock sphere indicator. The
project app derives its rectangular dimensions and liquid viscosity in
a separate source and build folder. The vendor archive stays intact.
[OpenLB geometry documentation](https://openlb.gitlab.io/release/content/03_geometry.html)

Both original N45 pilot fields reached their flow targets. They exceeded
the density limit, and the high-flow pilot also exceeded the Mach limit.
They supplied only the starting velocities and scaled density deviations
for the lower-time-step runs. Every field used for tracer calculations
passed the original screens.

| Transverse grid | Full grid | Spacing (µm) | Voxel void fraction | Maximum Mach over both flows | Largest relative density departure |
| --- | --- | ---: | ---: | ---: | ---: |
| N45 | 432 × 45 × 45 | 25.77 | 0.427046 | 0.03171 | 0.001617 |
| N60 | 576 × 60 × 60 | 19.33 | 0.426955 | 0.03222 | 0.001622 |

The target mean interstitial velocities were 9.570 and 25.670 mm/s,
respectively. The largest velocity-vector mismatch in any accepted field
was below 0.101%. Each tracer run used its own measured axial throughput and
liquid volume to match section Péclet numbers of 13,291.35 and 35,652.76.

The Julia transfer solves for locally balanced face flows while
preserving cross-sectional throughput. The maximum axial-plane or
transverse-plane throughput error was $2.41\times10^{-14}$ of the
axial flow, and the RMS face corrections were 1.61–2.90%. Flux-weighted markers provide the
fraction reaching a clean wall and the time spent in the liquid before
capture or exit. Their ratio, multiplied by the liquid volume, is the
conductance of the synthetic section. Its film coefficient divides that
conductance by the area of the modeled spheres. The projected clean-bed
result retains the earlier access resistance and recomputes the exposure
of each flow path from the denser trajectories. The 25 cm projection
uses fully mixed joins between the eight-section patterns.
[Film estimate at 60% porosity](../disordered_film_conductance/REPORT.md)

## Clean-bed estimates at both flow rates

All rows use the 25 cm column with 17.581 g of carbon and the previous
slow-access resistance. Film conductance is in $10^{-6}$ m³/s.
Both escape columns include the wall-contact distribution estimated
from the new trajectories.

| Grid and tracer step | Flow (mL/min) | Film coefficient (µm/s) | Conductance scaled by carbon area | Escape scaled by carbon area (%) | Escape scaled by liquid volume (%) |
| --- | ---: | ---: | ---: | ---: | ---: |
| N45, 20 µs | 23.3 | 127.3 | 6.520 | 0.00041 | 0.00040 |
| N60, 40 µs | 23.3 | 106.6 | 5.466 | 0.00305 | 0.00300 |
| N60, 20 µs | 23.3 | 107.8 | 5.526 | 0.00293 | 0.00288 |
| N45, 20 µs | 62.5 | 211.7 | 10.830 | 0.03155 | 0.03119 |
| N60, 40 µs | 62.5 | 176.5 | 9.032 | 0.12663 | 0.12528 |
| N60, 20 µs | 62.5 | 181.1 | 9.269 | 0.11011 | 0.10892 |

Every tracer ensemble recovered the full inert dose and its mean
residence time agreed with liquid volume divided by flow within 0.59%.
On N60 at 20 µs, 64.31% and 48.94% of the markers hit a wall in the
short packing at low and high flow, respectively. The previous 60%-void
case gave 29.10% and 16.91%. The comparison includes the different
packing arrangement and velocity fields. Matching void fraction
brings the carbon-area and liquid-volume scale-ups within 0.21% in
conductance.

## High-flow cutoff and numerical sensitivity

At 62.5 mL/min, the N60, 20 µs estimate gives an effective film
conductance of $9.269\times10^{-6}$ m³/s, including the fixed access
resistance for the slow fraction. The previous ideal-distribution
threshold was $3.147\times10^{-6}$ m³/s. The new estimate uses the
wall-contact distribution from the porosity-matched flow. Its 95%
tracer-sampling interval for escape is **0.094–0.132%**, conditional
on that distribution. The largest upper bound across the checked
grid, time-step and scaling combinations is **0.148%**.

Halving the tracer step from 40 to 20 µs on N60 changed high-flow
conductance by 2.62% and projected escape by 0.017 percentage
points. Refining from N45 to N60 at 20 µs lowered conductance by
14.41% and changed escape by 0.079 percentage points. The grid
change exceeds the 3% conductance limit. The escape change is small
relative to the 5% cutoff, so the result resolves the cutoff for
this packing without qualifying the coefficient as grid-converged.
At 23.3 mL/min, the time-step and grid conductance changes were
1.09% and 15.25%, respectively.

This comparison addresses the clean bed before carbon loading, with
one synthetic packing and the same assumed mixing between repeats.
Figure 7 was left unfitted.

Read-only replay verified 34 source files and 1,137 local artifacts.
The OpenLB vendor integrity check covered all 2,070 original files.

[Scientific protocol](protocol.md) · [Commands](README.md) · [Compact record](reference_summary.json)
