# Three-dimensional capillary: verifying measured permeability

## Purpose

Does a native OpenLB velocity/pressure calculation recover the hydraulic
resistance of a straight cylindrical pore, and is its inferred permeability
stable under numerical and operating changes? Unlike the network reference,
this example resolves flow inside the capillary before Julia extracts a
property from measured fields.

## Theory and measurement definition

A cylinder of radius $R=50$ µm and length 400 µm lies along $+x$. The
representative sample area is $A_{cell}=4\times10^{-8}$ m², including solid
area; it is not the open pore area $\pi R^2$. Water-like reference properties
are $\mu=10^{-3}$ Pa s and $\rho=1000$ kg/m³. Gauge pressure decreases from
0.1 to 0 Pa. The flow is steady, isothermal, Newtonian and laminar, with
no-slip walls and no gravity.

The continuum reference obeys incompressible mass and momentum balance,

```math
\nabla\cdot\mathbf u=0,\qquad
\rho(\mathbf u\cdot\nabla)\mathbf u=-\nabla p+\mu\nabla^2\mathbf u.
```

For fully developed axial flow, the convective term vanishes and the radial
viscous balance produces the profile below.

For an axial pressure gradient, fully developed momentum balance yields

```math
u(s)=\frac{\Delta p}{4\mu L}(R^2-s^2),\qquad
Q=\int_0^R2\pi s u(s)\,ds=\frac{\pi R^4\Delta p}{8\mu L}.
```

Here $s$ is radial distance. Comparing with Darcy flow gives the independent
sample-area permeability:

```math
k_{exact}=\frac{\pi R^4}{8A_{cell}}
=6.13592315154\times10^{-11}\ \mathrm{m^2}.
```

OpenLB evolves its lattice populations with stock collision, streaming, wall
and pressure-boundary facilities. The accepted current porous-flow benchmark
uses D3Q19/TRT with magic parameter $3/16$. It starts from a numerical initial
state; the analytical profile is a verification oracle. Julia instead uses
the measured plane discharge and interior pressure taps:

```math
k_{measured}=\frac{\mu L_{taps}Q}{A_{cell}\Delta p_{taps}}.
```

Interior taps isolate the bulk pressure gradient from boundary effects. This
measurement is not a reservoir-to-reservoir entrance-loss coefficient. Profile
error, plane-integrated mass spread and inlet/outlet mismatch are independent
diagnostics. Reynolds number uses pore diameter and mean pore velocity;
lattice Mach and density departure screen the weak-compressibility approximation.

### What TRT relaxes

In the unforced bulk, opposite discrete velocities $i,\bar i$ define even/odd
population parts $f_i^{\pm}=(f_i\pm f_{\bar i})/2$. With the corresponding
second-order equilibrium parts, the stock TRT collision and streaming are

```math
f_i^{\star}=f_i-\omega_+(f_i^+-f_i^{eq,+})
-\omega_-(f_i^--f_i^{eq,-}),\qquad
f_i(\mathbf x+\mathbf c_i\Delta t,t+\Delta t)=f_i^{\star}(\mathbf x,t).
```

For this D3Q19 model, lattice viscosity is
$\nu_{lb}=(1/3)(1/\omega_+-1/2)$ and
$\nu=\nu_{lb}\Delta x^2/\Delta t$. The odd rate follows the chosen magic
parameter,

```math
\left(\frac{1}{\omega_+}-\frac12\right)
\left(\frac{1}{\omega_-}-\frac12\right)=\frac{3}{16}.
```

OpenLB supplies these operations, its D3Q19 velocities, equilibrium and
Bouzidi wall treatment; project code does not recreate them. This bulk
description is distinct from the stock boundary-node updates. The installed
1.9r0 `collision::TRT` and `TRTdynamics` definitions were checked for this
report; changing relaxation or discretization requires its own verification.

## Results

Current TRT evidence is run `20261109T145749Z-e4c8d041` under
`results/hydraulic_trt/`. Eight solver cases and three geometry audits pass.

| Finest-grid TRT quantity | Recorded value |
| --- | ---: |
| Cells across diameter | 40 |
| Measured permeability | $6.13304992846\times10^{-11}$ m² |
| Permeability error against analytical value | 0.04683% |
| Relative velocity-profile L2 error | 0.04990% |
| All-plane mass spread | 0.03361% |
| N32 → N40 permeability change | 0.02052% |

Pressure reversal/doubling, smaller time step and tighter steady stopping also
pass the unchanged acceptance screens. Export replay checks SI coordinates,
zero-based IDs, artifact hashes and real Julia reanalysis. The retained
[benchmark record](../../../docs/hydraulic_benchmark.md) distinguishes these
TRT results from the older BGK tables and their separate run identities.

## Discussion

Agreement in permeability alone could conceal compensating geometry/profile
errors. The additional profile, conservation, connectivity and operating
checks reduce that risk. Curved-wall voxelization can give non-monotonic errors;
the grid-pair change is evidence of stability, not a claimed formal convergence
order or error bar for arbitrary pore shapes.

The result characterizes one aligned pore and one sample-area convention. It
does not establish isotropic permeability, a connected manufactured pore
network or predictive accuracy for a ceramic specimen. Historical BGK results
are retained as historical evidence; they do not substitute for the TRT
revision. The extra partial capillary rerun made during folder organization
is not part of this accepted benchmark.

## Conclusions

Measured native flow can supply a verified axial permeability through the
Python–Julia boundary. Downstream use should consume that measured property,
with its settings and provenance, as illustrated by the
[Darcy/pump coupling](../darcy_pump/REPORT.md).

## Reproduction and evidence

[Commands and export contracts](README.md) · [Current TRT and historical results](../../../docs/hydraulic_benchmark.md)
· [Earlier compact BGK reference](../../../docs/hydraulic_reference.json).
No experimental article is needed for the cylindrical analytical solution.
