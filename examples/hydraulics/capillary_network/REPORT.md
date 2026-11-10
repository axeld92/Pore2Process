# Capillary network: from pore flow to a Darcy column

## Purpose

Can a Python interface call an independently usable Julia numerical model,
recover a checkable effective property, and use it at a larger scale? This
example establishes that reference before introducing detailed OpenLB flow or
DWSIM equipment. The material is a synthetic lattice of identical capillaries.

## Geometry, assumptions and equations

The lattice has $n_x,n_y,n_z$ pressure junctions, spacing $a$ in metres,
capillary radius $r$ in metres and Newtonian viscosity $\mu$ in Pa s. Junctions
add no storage or pressure loss. Flow is steady, incompressible, fully developed
and laminar; gravity, entrance effects and transverse periodicity are absent.

For zero-based logical coordinates $(i,j,k)$, the node ID is
$i+n_x(j+n_yk)$, with $x$ varying fastest. Its position is
$(ia,(j+1/2)a,(k+1/2)a)$. Positive-$x$, positive-$y$ and positive-$z$ neighbours
are connected once. Positive throat flow runs from its first endpoint to its
second. The inlet/outlet faces have prescribed pressures $p_{in},p_{out}$.

Hagen–Poiseuille flow gives each throat's conductance and signed discharge:

```math
g_{ij}=\frac{\pi r_{ij}^{4}}{8\mu\ell_{ij}},\qquad
q_{ij}=g_{ij}(p_i-p_j).
```

Here $g$ is in m³/(Pa s) and $q$ in m³/s. Julia assembles the sparse
conductance Laplacian $K$. Conservation at interior junctions gives

```math
K_{II}p_I=-K_{IB}p_B.
```

There is no interior solve when $n_x=2$. Summing signed nodal net outflows on
the inlet, and their negatives on the outlet, gives $Q_{in}$ and $Q_{out}$
without double-counting boundary-face edges.

The sample length is $L=(n_x-1)a$ and representative area is
$A=n_yn_za^2$. Each transverse junction represents one cell of area $a^2$;
the declared transverse domain is $[0,n_ya]\times[0,n_za]$. The area between
outer node centres is therefore not the sample area used for permeability.

```math
k=\frac{\mu LQ_{in}}{A(p_{in}-p_{out})},\qquad
Q_{column}=\frac{kA_{column}\Delta p_{column}}{\mu L_{column}}.
```

Permeability $k$ is in m². At zero pressure drop the exact constant-pressure,
zero-flow state is returned; inferred $k$ and column flow are unavailable.
The column assumes repetition of the same idealized medium at a larger scale.

An independent analytical reference follows from $n_yn_z$ parallel paths, each
containing $n_x-1$ capillaries in series. Transverse pressure differences vanish:

```math
p(x)=p_{in}-\Delta p\,x/L,\qquad
Q_{exact}=\frac{n_yn_z\pi r^4\Delta p}{8\mu L},\qquad
k_{exact}=\frac{\pi r^4}{8a^2}.
```

These expressions check the matrix solution; they do not generate its answer.

## Results and verification

Defaults are $(6,4,4)$ nodes, $a=1$ mm, $r=50$ µm, $\mu=1$ mPa s and
$p_{in}/p_{out}=100/0$ Pa. The column has length 0.02 m, area $10^{-4}$ m²
and pressure drop 1000 Pa.

| Quantity | Recorded result |
| --- | ---: |
| Nodes / throats | 96 / 224 |
| Length / representative area | 0.005 m / $1.6\times10^{-5}$ m² |
| Pressure at successive x planes | 100, 80, 60, 40, 20, 0 Pa |
| Inlet flow | $7.853981633974476\times10^{-10}$ m³/s |
| Flow in display units | 47.1238898 µL/min |
| Inferred permeability | $2.454369260617024\times10^{-12}$ m² |
| Column flow | $1.227184630308512\times10^{-8}$ m³/s = 0.7363108 mL/min |
| Maximum pressure error | $6.39\times10^{-14}$ Pa |
| Normalized conservation error | $1.97\times10^{-15}$ |

The [numerical reference](reference/verified_reference.json) preserves the
observations. Tests cover pressure/radius/viscosity scaling, pressure offsets
and reversal, single-chain and asymmetric lattices, zero flow and invalid
inputs. Nonzero flow and permeability must agree within $10^{-9}$ relative;
pressure and conservation have separately declared tolerances. Python tests
call Julia and check array copies, indexing, units and export consistency.

## Discussion

The small errors establish implementation agreement with a problem that has an
exact solution. They do not measure uncertainty in a real material. Doubling
radius multiplies flow and permeability by 16, whereas doubling viscosity
halves flow and leaves permeability unchanged. These trends distinguish
geometry from fluid properties at the software boundary.

The model contains no resolved three-dimensional velocity field or realistic
network disorder. Python exports coordinates/pressure as float64 and endpoint
IDs as zero-based int64; Julia's one-based IDs are converted once in the bridge.
All exported geometry remains in SI units despite mm/µm display controls.

## Conclusions

The Python–Julia connection, conservation accounting and permeability-to-column
translation are verified for this idealized problem. The next level of spatial
physics is the [native capillary benchmark](../capillary_benchmark/REPORT.md).

## Reproduction and evidence

[Setup and commands](README.md) · [Model/data-contract details](../../../docs/guides/setup.md)
· [Acceptance record](../../../docs/verification/first_simulation.md).
The governing analytical reference is derived above; no experimental dataset
is used. This report summarizes retained evidence and introduces no new solve.
