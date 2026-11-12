# Prospective transport verification, 11 November 2026

This is a separate continuation of the [packed-bed investigation](../transport_root_cause/REPORT.md).
Run the wall-free case first, then the isolated curved wall. No full-bed,
adsorption or DWSIM run is part of this study. These are prescribed velocity
fields, not new hydraulic solutions or measurements.

Both tests use a cell-centred mesh of width 1.389172554 mm, D=1e-9 m²/s,
dt=1.9017455796408345e-5 s and 11,293 steps (0.2147641283 s). The velocity
scale is 0.01293729518 m/s. Initial concentration is 0.5 kg/m³ plus a mode
with 0.25 kg/m³ amplitude. Use stock OpenLB 1.9r0 D3Q27, the scalar-specific
regularization with stock advection–diffusion TRT, first-order equilibrium and
magic 0.3, as in the bed. OpenMP changes execution only. The curved case uses
stock cylinder geometry and full-way BounceBack.

## 1. Three-dimensional stretching

Let k=2π/L. In the periodic cube prescribe

`u=(U[1+0.5 cos(ky)cos(kz)], 0.35 U sin(kz), 0.35 U sin(ky))`.

Each component is independent of its own coordinate, so the continuum
velocity is exactly divergence free. The transverse motions change y and z
and their two shears do not commute. Initial concentration is
`0.5 + 0.25 sin(kx)`. Use full 32³, 64³ and 96³ cubes. Include a uniform-flow
analytical control and a constant-concentration control at N32. No thin-depth
assumption is permitted for the stretching runs.

For `C=0.5+0.25 Im[exp(ikx) f(y,z,t)]`, the exact PDE reduces to
`f_t = D(f_yy+f_zz-k²f) - uy f_y - uz f_z - ik ux f`, initially f=1.
The independent reference splits exact Fourier diffusion, phase rotation and
the two exactly translated transverse subflows symmetrically. Qualify at the
actual middle and final time: spatial and temporal refinement must each change
the envelope by at most 1e-4 of initial amplitude on the comparison grid.
Compare 64→128, 128→256 or 256→512 transverse points and 256→512 or 512→1024
steps per half-window; retain actual comparisons. Verify signs independently
against backward characteristics at zero diffusion and analytical translation
at uniform velocity. A failure in the reference prevents a claimed comparison.

## 2. One curved wall

Fluid fills a disk of radius 0.4 mm, extruded periodically along z. The exterior
is solid, including the square's periodic seam; there are no inlet, outlet
or end walls. Prescribe rigid rotation `u=(-Ω(y-yc), Ω(x-xc), 0)` with
Ω=U/R. This idealizes liquid rotating with the wall; there is no continuum
normal velocity or solute flux. The stationary voxel representation can make
that same tangential velocity cross individual voxel faces. A zero-flow pair
separates the diffusive boundary check from this advective effect.

Use the first Neumann mode of angular order 2:
`C=0.5+0.25 J₂(αr/R)/J₂(α) cos[2(θ-Ωt)] exp[-D(α/R)²t]`, where
α=3.0542369282271403 solves J₂′(α)=0. The analytical derivative gives zero wall
flux. Extend the mode at t=0 to initialize solid populations. Compute the
oracle independently with a power series and check the derivative root using
[NIST's derivative identities](https://dlmf.nist.gov/10.6) and
[Bessel zero convention](https://dlmf.nist.gov/10.21). The disk diffusion-only
cases set Ω=0 with the same mode. The constant case sets amplitude zero.

Use N32/64/96/128 across the surrounding square, four periodic z layers, and
a separate N32 full-depth match. The full-depth comparison must reproduce
all saved fluid fields and step extrema, with inventory proportional to depth.
Use diffusion-only N32/N128 and constant rotation at N128. Report actual
voxel volume and compare the continuum at actual fluid centres. Count only
fluid and active reflected wall populations in the closed numerical inventory.

## Screens and finite scope

Require correct native coordinates, material, SI velocity, initial condition,
diffusion and finite data. Numerical inventory drift must stay below 1e-8
relative. For nonuniform fields compare RMS and maximum error over initial
amplitude at both middle and final times, and retain the largest transient
range violation. Declare an engineering screen of 0.02 kg/m³ excursion and
1% RMS; for the disk also apply 1% RMS in the fixed outer 50 µm liquid shell.
Report the strict physical range separately. Constant tolerance is 1e-9 kg/m³.
Keep scientific failures and run the independent disk stage after the completed
stretching stage even if the stretching error exceeds the screen. Do not tune
D, clip values, change collision or launch a full bed to rescue either result.
No general grid/time convergence, column transport or experimental prediction
is implied by a bounded test at the chosen OpenLB time step.
