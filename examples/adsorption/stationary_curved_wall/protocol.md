# Stationary circular-wall protocol, 11 November 2026

Follow the [rotating-wall benchmark](../transport_verification/REPORT.md)
with the same 0.4 mm radius, 1.389172554 mm square, Bessel initial mode,
D=1e-9 m²/s, dt=1.9017455796408345e-5 s and 11,293 steps. The disk has periodic
depth and no ends. Prescribe steady incompressible angular shear
`u=(-Ω(r)(y-yc), Ω(r)(x-xc), 0)`, with `Ω(r)=Ω0[1-(r/R)²]` in the disk and
`Ω0=(3√3/2)U/R`. The speed is exactly zero at the stationary smooth wall,
and its peak is U=0.01293729518 m/s at r=R/√3, matching the old maximum speed.
This is a prescribed scalar-transport control, not a new hydraulic solve.
Differential rotation deforms the wave; the rigid-rotation result is a
matched comparison for the *constant* tracer, not the wave's reference.

Use stock OpenLB 1.9r0 cylinder indicator, D3Q27 with the scalar-specific
regularization and stock advection–diffusion TRT, magic 0.3, and full-way
BounceBack. No original library files, dynamics or boundary operators are
changed. Run N32/64/96/128 with the original m=2 nonuniform mode, constant
controls at N64/N128, and one N32 full-depth equivalent. Thin cases have
four axial layers; the N32 comparison must reproduce all saved fields and
every-step extrema, with the correct depth scaling for both inventories.
The archived no-flow and rigid-rotation cases remain separate. No clipping,
added diffusion, fitting, full bed, adsorption or DWSIM runs.

## Independent continuum reference

Set s=r/R and write `C=0.5+A Re(exp(2iθ) f(s,t))`, with A=0.25 kg/m³.
Then `f_t=(D/R²)(f_ss+f_s/s-4f/s²)-2iΩ(s)f`, with f regular at the origin
and f_s=0 at the wall. Initial `f=J₂(αs)/J₂(α)`, α=3.0542369282271403.
This uses [Bessel's equation](https://dlmf.nist.gov/10.2) and its
[derivative identity](https://dlmf.nist.gov/10.6), independent of OpenLB.
The substitution f=s²g gives weighted diffusion `s^-5 (s^5 g_s)_s` with
zero origin flux and wall condition `g_s(1)=-2g(1)`. Julia uses radial
finite volumes and Crank–Nicolson with the standard tridiagonal solver.
Compare spatial and temporal refinements separately at native midpoint and
final times. Require maximum radial mode change ≤1e-4 of amplitude on 4,097
radii, including the wall. Finite budget: 256→512→1024→2048 radial cells,
2048→4096→8192 time subdivisions per observation interval. Check pure
diffusion, rigid rotation and diffusion-free shear against exact modes.

## Interpretation and screens

Verify initialized fields, wall geometry, SI velocities, time and TRT
diffusivity. Require finite fields and ≤1e-8 relative closed numerical
inventory drift, counting fluid plus only active reflected wall populations.
Report physical liquid inventory separately. Check every-step bounds and
the exact saved peak field. Use the existing 0.020 kg/m³ engineering excursion
allowance, 1% RMS of initial amplitude globally and in the fixed outer 50 µm
shell at midpoint and final, and a strict 1e-9 kg/m³ constant tolerance.
Also report strict physical bounds and immediate wall-link error. Track
velocities at wall-adjacent cell centres; these need not be zero even though
the prescribed continuum wall velocity is zero. Retain failed results.
The matched physical timestep permits spatial comparison; this is not a
native temporal-convergence study or validation of a fixed particle bed.

## Reference extension recorded before any native run

The initial reference-only attempt (`20261111T121850Z-40561f07`) missed the
radial screen at midpoint: 1.181e-4 for 1024→2048 cells. Each spatial doubling
reduced the change by about four. Allow one further radial doubling to 4096,
keeping the original tolerance, physics, time budget and all native cases.
The stopped attempt and its original source snapshot remain in results.
