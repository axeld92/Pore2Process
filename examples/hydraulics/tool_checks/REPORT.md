# Native-engine checks: channel flow and pump automation

## Purpose

Before coupling models, can OpenLB and DWSIM execute independently on the
target machine, and can Python inspect their numerical outputs? These are
installation and component checks, rather than a porous-material design study.

## Theory and test definitions

The unmodified official OpenLB plane-channel example solves steady laminar
flow between no-slip walls separated by $H$, with length $L$ and maximum
velocity $U_{max}$. The independent reference is

```math
u_x(y)=4U_{max}\frac{y}{H}\left(1-\frac{y}{H}\right),\qquad
\Delta p=\frac{8\rho\nu U_{max}L}{H^2}.
```

The stock probe uses D2Q9/BGK, Re=10, relaxation time 0.8 and a 2 m × 1 m
domain. It is a compatibility test of an unchanged official channel example,
not a new porous-media BGK calculation. Velocity/pressure error norms,
boundary-node counts and reported stopping are inspected. A permeability
export or separate integrated conservation test is not part of this probe.

The DWSIM probe solves pure water through a stock pump, with an energy stream,
using IAPWS-IF97 properties. Inputs are 298.15 K, 100000 Pa inlet, 10 kg/s and
75% efficiency. Outlet pressures 1000000 and 1900000 Pa give pressure rises
of 900000 and 1800000 Pa. The independent work relation is

```math
W=\frac{\dot m\Delta p}{\rho_{in}\eta}.
```

Mass and enthalpy-flow balances are checked separately. Efficiency is a
fraction, pressure is absolute Pa, and DWSIM kW/kJ kg⁻¹ outputs are converted
to W/J kg⁻¹. Both flowsheets are saved, reloaded and solved again.

## Results

| OpenLB channel diagnostic | N32 | N64 |
| --- | ---: | ---: |
| Inlet / outlet boundary nodes | 31 / 31 | 63 / 63 |
| Relative velocity L2 error | 0.00279288 | 0.000699683 |
| Relative pressure L2 error | 0.00737092 | 0.00184006 |

Doubling resolution reduces each error by about four. A separate N40 attempt
has zero outlet nodes, density drift and unacceptable pressure error; it is
rejected even though its velocity stopping criterion was met.

| DWSIM pump diagnostic | Baseline | Double pressure rise |
| --- | ---: | ---: |
| Shaft work (W) | 12035.5356965 | 24071.0713931 |
| Mass imbalance (kg/s) | 0 | 0 |
| Enthalpy-balance residual (W) | 0.000783894 | 0.000753268 |

The native ARM64 OpenLB 1.9r0 and DWSIM 10.2.15 probes completed on the
development Mac. Automation uses an external C#/.NET adapter launched by
Python. No DWSIM GUI verification is inferred from the automation result.

## Discussion

The rejected mesh illustrates why a small velocity error or executable exit
code is insufficient: the intended boundary conditions must exist. Two
successful grids support this stock channel check, not every future geometry.
The pump's exact factor-two scaling establishes unit and equipment behavior
using the property-package density; it is not independent water-property
validation. The probes do not pass a pore-derived property between engines.

## Conclusions

Both native engines are executable and inspectable on the tested platform.
The next numerical stage is the [3D capillary](../capillary_benchmark/REPORT.md),
followed by the [actual property-to-pump connection](../darcy_pump/REPORT.md).

## Reproduction and evidence

[Setup/launch commands](README.md) · [Compatibility record](../../../docs/verification/hydraulics/tool_checks.md)
· [Compact reference](reference/compatibility_reference.json).
The analytical channel and pump expressions above are the test oracles;
no experimental dataset is used.
