# Prospective second-order-equilibrium diagnostic

10 November 2026. Test whether the first-order scalar equilibrium contributes
materially to the demonstrated high-Peclet excursions. This is a diagnostic,
not a change to the accepted hydraulic, particle or process models.

Compose installed stock `equilibria::SecondOrder`, `collision::TRT`,
`AdvectionDiffusionExternalVelocityCollision` and
`CombinedAdvectionDiffusionRLBdynamics`. No dynamics implementation or descriptor
is recreated. Use stock `equilibrium<D>::secondOrder` for initial populations,
because CombinedRLB's inherited initialization explicitly uses first order.
The periodic cases have no inlet/outlet boundary whose equilibrium is affected.
No clipping, diffusivity increase, concentration source, fitted factor or
velocity modification is permitted. Keep physical D=1e-9 m2/s, geometry and
velocity, mesh N32, common physical dt and observation window unchanged.

Run five actual native cases: porous sine and uniform concentration; smooth
periodic shear sine and uniform concentration; and free uniform-velocity sine
with an analytical solution. Use the existing parent controls for comparison.
Preserve all actual command lines, initialization, build/source/artifact hashes.

Report the whole-history concentration range, conservation, final shear-profile
error against the independent converged Fourier reference, and uniform-flow
analytical error. Improvement in just one case does not qualify a correction.
The engineering bound allowance remains 0.02 kg/m3; existing 1% smooth-profile
and constant-preservation screens remain unchanged. A failed candidate remains
archived as failed; no adsorption or DWSIM prediction follows from it.
