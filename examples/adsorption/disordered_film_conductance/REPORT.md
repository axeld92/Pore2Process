# Phenol transfer through the liquid film in the saved N160 flow

At 62.5 mL/min, the tracer-based coefficient supports a **62.3% increase**
in effective film conductance when applied to the column's fixed carbon
surface area. The earlier ideal-allocation test required 9.94%. Transferring
by liquid volume gives an 18.9% decrease. The two results distinguish the
liquid-film estimate from the assumption used to transfer it to the column.
[Earlier clean-film bound](../disordered_early_bypass/REPORT.md)

## How the coefficient was estimated

The study used the balanced, saved OpenLB 1.9r0 D3Q27 TRT flow through
32 synthetic 0.80 mm spheres. The packed section was 11.11 mm long, with
periodic transverse boundaries and a liquid void fraction of 0.60. The
same velocity shape was reused at liquid Péclet numbers of 13,291 and
35,653. Those numbers match the 23.3 and 62.5 mL/min conditions of the
25 cm column. The conversion back to physical units uses the column
interstitial velocity and phenol diffusivity of $10^{-9}$ m²/s.
[N160 flow and tracer verification](../disordered_n160/REPORT.md)

An inlet pulse represents a fixed dose, with each tracer weighted by
inlet flow at injection. Each tracer contributes to wall uptake at its
first carbon contact and to the liquid inventory until contact or exit.
For a continuous inlet at concentration $C_{in}$, the corresponding
uptake rate is $QC_{in}(1-S)$ and the average concentration in the liquid
is $QC_{in}\bar T_a/V_l$. Here $Q$ is the flow rate, $S$ the fraction
that exits before any wall contact, $\bar T_a$ the mean time before
capture or exit and $V_l$ the voxelized liquid volume. Dividing uptake by
average concentration gives the film conductance. Dividing by the
analytical external sphere area gives the film coefficient. This treats
the carbon surface concentration as zero and isolates transfer across
the liquid. The existing column correlation is the comparison.
[First-contact capture](https://doi.org/10.1103/PhysRevE.107.044105)
· [Packed-bed film correlation](https://doi.org/10.1016/0009-2509(78)85120-3)

The archived 80 µs runs supplied 24,576 trajectories at each rate. New
40 and 20 µs runs used a matched subset of 8,192 inlet seeds. These
time steps are in the archived flow's time units. Halving
the step checks a known sensitivity of first-contact estimates near
absorbing boundaries. Whole-trajectory bootstrap intervals quantify
tracer sampling variation. The earlier N128–N160 comparison checked
wall-contact probabilities. The film coefficient here comes from one
packing and one spatial grid. [Absorbing-boundary time discretization](https://doi.org/10.1051/ps:2001112)

| Flow (mL/min) | Film coefficient (µm/s) | 95% tracer sampling interval (µm/s) | Previous correlation (µm/s) | First and second halves (µm/s) |
| ---: | ---: | ---: | ---: | ---: |
| 23.3 | 63.74 | 61.12–66.16 | 31.88 | 61.60, 66.01 |
| 62.5 | 90.54 | 85.63–95.22 | 55.61 | 86.79, 94.51 |

At 20 µs, 29.10% and 16.91% of the inlet pulse hit a wall before exit.
The unreacted mean residence times agreed with liquid volume divided by
flow within 0.22% across all ensembles. The largest change when going from
80 to 40 µs on the matched subset was 3.06% in estimated conductance, so
the protocol called for the 20 µs check. The change from 40 to 20 µs was
0.05% and 1.96% at the low and high flows, respectively, within the 3%
limit. Both estimates remained on their respective sides of the high-flow
conductance threshold.

## How the coefficient was transferred to the column

Both comparisons retain 17.581 g of carbon, its external area of
0.10210 m² in the spherical particle idealization, and the previous
access model with half on the direct path and the other half behind
the same extra resistance. The reduced coefficient for the slow
half is recalculated with the new liquid-film coefficient and that
fixed resistance in series. The model then combines the two halves
in parallel. [Previous access assumptions](../carbon_access_transfer/REPORT.md)

The surface-area transfer applies the measured coefficient per square
metre to the column's external carbon surface area. The liquid-volume
transfer instead preserves the synthetic packing's conductance per cubic
metre of liquid before applying the same access resistance. These yield
different projections because the nominal column void fraction is 0.427
versus 0.60 in the synthetic packing. The nominal column area per cubic
metre of liquid is approximately twice the synthetic value. Agreement
between time steps tests the pore-scale calculation. The scale-up relies
on which feature of the synthetic packing transfers to the denser bed.

## 25 cm clean-film check

All conductances in the table include the previous slow-access resistance.
The minimum column conductance for 5% escape assumes the most favorable
distribution of conductance across the flow paths.

| Flow (mL/min) | Previous | Minimum for 5% | Carbon-area transfer | Liquid-volume transfer |
| ---: | ---: | ---: | ---: | ---: |
| 23.3 | 1.651 | 1.173 | 3.277 | 1.640 |
| 62.5 | 2.862 | 3.147 | 4.645 | 2.320 |

Conductances are in $10^{-6}$ m³/s. The high-flow sampling intervals
were 4.395–4.885 by carbon area and 2.195–2.439 by liquid volume, in the
same units.

Keeping the earlier wall-contact allocation gives 0.38% and 4.76%
clean-bed escape by the carbon-area transfer at 23.3 and 62.5 mL/min.
At high flow the 95% sampling interval is 4.27–5.35%, which straddles
the cutoff. With the liquid-volume transfer, the predictions are
3.61% and 16.64%. The high-flow carbon-area point estimate therefore
supports the increase, with the 5% cutoff still unresolved for the
current allocation. The synthetic-to-column area transfer is the
larger uncertainty.

The analysis reused the existing flow field and kept the batch particle
model outside the clean-film diagnostic. No adsorbing-bed rerun or Figure 7
fit was performed. A clean-wall comparison at a porosity closer to 0.427
would test the area transfer directly in a denser synthetic packing.

[Protocol](protocol.md) · [Commands](README.md) · [Compact record](reference_summary.json)
