# Balanced face flows on the stationary disk

Keep the archived smooth disk, 0.4 mm radius, periodic axial depth equal to
four OpenLB cells at N=128, steady quadratic angular shear, peak speed
0.012937295182166382 m/s, and molecular diffusivity 1e-9 m²/s. Solve the
axially invariant transport equation over the full disk in Julia, using
annular sectors that end at the smooth wall. Initialize the same constant
0.5 kg/m³ or Bessel mode of amplitude 0.25 kg/m³ using cell averages.

Prescribe volume flow by integrating the azimuthal speed across each radial
face of an annular sector. Both azimuthal faces of a cell have the same
volume flow, with opposite signs in its balance. Radial and axial carrier
flows are zero. Use one flux per internal face for tracer transport and
zero tracer flux at the origin and wall. Use centered face concentration
for advection and two-point diffusion with polar face metrics. Advance the
full real concentration field with a symmetric radial–azimuthal–radial
split and Crank–Nicolson within each part. Use no clipping, fitted diffusion
or mode-dependent solver.

Reuse the archived independently converged Julia radial mode at both original
observation times. Compare volume-weighted errors over the annular sectors
using averaged reference values and compare interpolated solutions at the
original OpenLB liquid cell centres. Keep the 1% RMS limit (relative to the
initial 0.25 kg/m³ amplitude) globally and in the outer 50 µm at both times,
the 1e-9 kg/m³ uniform maximum departure limit, the 0.02 kg/m³ engineering
bound, and the relative liquid inventory limit of 1e-8. Report the 1e-8 kg/m³
physical-bound check separately. Check volume and local carrier-flow closure.

Run the varying case at 40×256 and 80×512 radial by azimuthal cells with a
2e-4 s time step cap. At the finer mesh halve that cap, keeping the exact
reference times, and also run the uniform tracer at that refined setting.
Require spatial changes below 0.5% and time-step changes below 0.1% of the
initial amplitude as RMS on the fixed comparison points in both regions
and at both times. If the varying solution or spatial check misses its
screen, permit one further mesh doubling with the two time caps and repeat
the uniform check there. Retain every executed case.

This test verifies tracer transport with known balanced face flows on a
smooth disk. OpenLB flow-field transfer and full-bed qualification are
separate stages.
