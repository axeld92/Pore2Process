# Measured pore permeability changes actual pump duty

## Purpose

Can a quantity obtained from OpenLB actually change a DWSIM process result
through Python and Julia? The controlled input is pore radius. The process
metric is pump shaft power at fixed throughput through an idealized porous slab.

## Theory and scale translation

Two straight pores, with radii 50 and 40 µm, share a 400 µm length and a
200 µm square representative cell. Each native flow calculation measures
discharge $Q_p$ and the pressure difference between internal taps. Julia infers

$$
k=\frac{\mu L_{taps}Q_p}{A_{cell}\Delta p_{taps}}.
$$

The upscaled slab has $L_b=1$ m, $A_b=1$ m² and prescribed total water flow
$Q_b=10^{-5}$ m³/s. Repetition of the aligned-pore medium and Darcy linearity
are assumptions. The required pressure rise and independent pump-work oracle are

$$
\Delta p_b=\frac{\mu L_b Q_b}{kA_b},\qquad
W_{shaft}=\frac{Q_b\Delta p_b}{\eta}
=\frac{\dot m\Delta p_b}{\rho\eta},\qquad \eta=0.75.
$$

The actual DWSIM Steam Tables calculation at 298.15 K and 100000 Pa supplies
$\mu=0.0008900825175722$ Pa s and $\rho=997.0474354081$ kg/m³. These
properties configure OpenLB, Julia and $\dot m=\rho Q_b$. Gauge pressure load
becomes an absolute pump discharge pressure $100000+\Delta p_b$ Pa.
DWSIM solves a stock pump with connected material/energy streams, saves its
flowsheet, loads it and solves it again. The porous slab is represented by
the externally calculated load; there is no custom porous-medium unit here.

An independent cylindrical argument gives $k\propto R^4$ at fixed cell area,
so the 40/50 µm permeability ratio should be $0.8^4=0.4096$ and the fixed-flow
shaft-power ratio its inverse, 2.44140625. Measured permeability, rather than
this expression, drives the process calculation.

## Results

The completed connected evidence is run `20261108T181450Z-e1b190f9`.

| Radius | Measured permeability (m²) | Slab pressure rise (Pa) | Actual DWSIM shaft power (mW) |
| --- | ---: | ---: | ---: |
| 50 µm | $6.130870215\times10^{-11}$ | 145.18045 | 1.9357394 |
| 40 µm | $2.511215119\times10^{-11}$ | 354.44296 | 4.7259061 |

The measured permeability ratio is 0.409601742; the actual power ratio is
2.44139587. Both geometries were checked at N24/N40. Fine-grid analytical
errors are below 1%; conservation, stopping, units and save/reload checks pass.
Pump work agrees with the independent measured-property reference within
$10^{-9}$ relative. Enthalpy residuals are 0.0409% and 0.0168% of shaft work,
below the declared 0.1% screen. Estimated process pore Re is below 0.015.

## Discussion

This establishes a causal software connection: changing geometry changes a
native measured property, a Julia equipment load and actual DWSIM work. The
very small powers belong to the stated open-pore, low-throughput scenario;
they do not estimate a commercial filtration plant's utility demand.

The table is preserved evidence from the original BGK-era connected study.
The [separate capillary benchmark](../capillary_benchmark/REPORT.md) was later
revalidated with TRT. That later component verification is not presented as a
fresh execution of this entire two-design process study. The partial attempt
during folder organization did not complete and supplies none of the table's
accepted results.

Entrance losses, random connectivity, other equipment head, motor efficiency
and real-material calibration are absent. DWSIM's pump check verifies the
coupling and work accounting, not the independent accuracy of its water
property package.

## Conclusions

The end-to-end property-to-process path is demonstrated for an idealized
hydraulic load. A meaningful application must additionally identify which
resistance can change and which process duties remain fixed; that is the
purpose of the [layered filtration example](../../filtration/layered_support/REPORT.md).

## Reproduction and evidence

[Commands and contracts](README.md) · [Connected-study record](../../../docs/connected_hydraulics.md)
· [Compact measured/simulator reference](../../../docs/coupling_reference.json).
The numerical reference is analytical; no experimental material dataset is used.
