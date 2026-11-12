# Balanced face flow handoff on one periodic FCC unit cell

Reuse the converged stock OpenLB 1.9r0 D3Q27 TRT fields at N=32, 64 and
96 for the 0.8 mm spheres and 1.389172554 mm periodic cell. Keep the archived
water-to-carrier velocity conversion and molecular diffusivity 1e-9 m²/s.
The comparison uses one pore geometry cell per resolution, with no inlet,
outlet, adsorption, or new hydraulic, full-bed or process execution.

Reconstruct one signed volume flow for each Cartesian fluid-to-fluid face
from the mean of the two adjacent normal velocities times face area. Solid
faces carry zero flow. Use a minimum Euclidean correction of the shared
face flows, subject to zero net flow in every fluid cell and preservation
of the native mean flow in each coordinate direction. Audit imported
coordinates, mask, viscosity conversion and native throughput, then
report unbalanced and balanced plane flows and cell residuals. The axial
plane flow must agree with the archived volume-throughput value to 1e-9
relative, and the local residual must be within 1e-10 of maximum face flow.
Report RMS and maximum corrections relative to face scales and flag RMS
corrections above 5% of initial face-flow RMS as a fidelity failure.

Integrate conservative centered advection and two-point molecular
diffusion on those faces in Julia at the stationary-disk observation
times, 0.10737255542652152 and 0.21476412830883945 s. Use classical RK4.
Check constant 0.5 kg/m³ throughout against the 1e-9 kg/m³ limit. Track
the varying profile initialized by 0.5+0.25 sin(2πx/L) kg/m³. Require
1e-8 relative liquid inventory drift and the prior 0.02 kg/m³ excursion
allowance beyond 0.25–0.75 kg/m³. Report strict physical bounds separately.
Verify the discrete scalar variance budget using the integrated physical
diffusion dissipation and a normalized residual of at most 1e-5 of initial
variance. Use no clipping, fitted dispersion or added artificial diffusion.

Check a periodic straight-flow sine against its analytical translation
and diffusion. At N64 compare the FCC varying field to a separately
evaluated matrix exponential action on the same discrete face operator.
Require <=0.1% RMS of initial amplitude at both times. Halve the integration
step at N96, with the same RMS tolerance. For the different saved hydraulic
meshes, use global liquid variance and the amplitude of the fitted
longitudinal fundamental as spreading measures. Compare N64 to N96 at
each time. Require both measures to change by at most 5% relative to N96
and the fundamental phase to change by at most 0.05 rad. This tests
spatial sensitivity of this periodic spreading control. It supplies
neither a continuum field oracle in FCC nor a column dispersion coefficient.
Keep negative results and distinguish passed checks from an unqualified
spreading calculation.
