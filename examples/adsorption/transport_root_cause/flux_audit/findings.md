# Local flux and wall-location findings

The exact population diagnostics distinguish a conservation error from an
incorrect spatial distribution. The native N32 sine-profile peak has local
collision, streaming and signed-link balance residuals below 1e-15 in native
concentration units. Fluid-wall link exchanges cancel when full-way BounceBack
storage is included. **No unshifted post-collision population is negative even
at the concentration peak**: the minimum is positive, about 5.83e-5. Thus
concentration boundedness and individual population positivity are separate
properties in this case. Numerical mass conservation does not guarantee a
physically admissible spatial concentration distribution.

The actual two-step, initially constant N32 control confirms a distinct coupling
imperfection. Its first-step concentration pattern matches the prediction from
stock `equilibrium<D>::firstOrder`, sampled at the actual upstream velocities,
to 5.55e-16 kg/m3. The maximum relative first-step change is 1.1912e-4. This is
the scalar equilibrium streaming stencil's failure to preserve a constant
exactly; it is neither an approximate centred velocity divergence nor evidence
that OpenLB's actual water populations create mass. The two-step control alone
does not establish the long-time size of this effect or explain the much larger
sine-profile excursions.

## Where the worst concentration peaks occur

The peak snapshots place all 16 severe excursions above 0.1 kg/m3 at
wall-adjacent nodes near the **downstream poles of the specified spheres**, for
both N32 and N48. The nearest smooth-surface outward normal has x component
0.99724 at N32 and 0.99881 at N48; positive x is the mean flow direction.
Nearest-surface gaps are 0.312 and 0.696 grid spacings. These are low-speed
regions: the maximum-concentration node speeds are about 0.000755 and
0.000499 m/s, respectively.

The smaller errors extend further into the liquid. Of 672 N32 fluid nodes
outside the engineering allowance, 608 have no solid in their immediate
26-neighbour geometric neighbourhood; for N48 this is 1040 of 1136 nodes.
Peak location is an association, not proof of where the errors first arise.
The declared smooth solid surfaces are convex spheres; voxel corners and
constricted void passages should not be described as physical concave solid
surfaces. The independent no-solid shear control is needed to separate bulk
high-Peclet transport errors from wall amplification.

## What stock OpenLB facilities offer

The installed scalar TRT alias uses first-order equilibrium and explicitly
notes a small diffusion-term error. `CombinedAdvectionDiffusionRLBdynamics`
reconstructs first-moment nonequilibrium before TRT: the even nonequilibrium
modes are removed. Consequently, retuning the TRT magic parameter at fixed
diffusivity provides essentially no bulk freedom in this regularized
configuration. No parity or diffusion-conversion mistake was found.

A separately archived candidate composes stock SecondOrder equilibrium, TRT,
external velocity and CombinedRLB. Its initial populations use the stock
second-order equilibrium explicitly, because CombinedRLB's ordinary equilibrium
initialization is hardcoded to first order. Five actual native runs reject this
candidate as a correction: porous concentration extrema remain about
.1180–.8827 kg/m3, while the smooth-shear profile error remains about 38.59% of
amplitude. Uniform-flow analytical and constant-shear controls pass. This
isolates the omitted second-order equilibrium contribution as an insufficient
explanation of the dominant failures. See the [candidate report](candidate/README.md).
The existing source and gates are preserved.

`RegularizedHeatFlux3D` explicitly rejects curved walls. The installed
`BouzidiAdeDirichletPostProcessor` imposes fixed concentration and is unsuitable
for an impermeable particle. Generic reflected Bouzidi populations exist, but
using their curved interpolation for a scalar would require independently
qualifying no-flux behavior, global mass balance and represented wall location.
The file named BoolakeeNeumann concerns solid-mechanics boundary traction and
is not a scalar no-flux implementation.

OpenLB also contains an explicit finite-difference advection/diffusion module
with upwind and central schemes plus zero-Neumann boundary facilities. Its
upwind source implements pointwise `u*grad(C)`, not conservative face mass
fluxes. Conservation in spatially varying voxel velocities would need separate
verification. First-order upwind numerical diffusion is of order `u*dx/2`,
roughly 2.8e-7 m2/s at representative N32 velocity: hundreds of times the stated
molecular diffusivity. Removing oscillations that way would not, by itself,
validate tracer spreading or a column dispersion coefficient. No alternative
finite-difference production solver was implemented.

The [machine-readable findings](findings.json) preserve exact values and hashes.
Recompute the local findings with:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.flux_audit.findings
uv run pytest tests/adsorption/test_transport_flux.py -q
```

All conclusions here are limited to the actual archived native cases and the
read-only installed-source review. No adsorption calibration, DWSIM calculation
or experimental claim follows from these diagnostics.
