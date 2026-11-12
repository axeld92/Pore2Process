# Prospective stock regularized-inlet diagnostic

10 November 2026. The sharp inlet and two ramp durations with the original
Dirichlet inlet already show substantial transient negative concentrations.
The uniform .80 mm case passes the new engineering screen. Changing the ramp
alone is therefore not demonstrated to resolve the front problem.

Test ONE additional stock boundary choice: replace the scalar inlet's stock
AdvectionDiffusionDirichlet with stock RegularizedTemperature, explicitly mixed
with AdvectionDiffusionTRTdynamics. The installed 1.9r0 regularizedTemperature3D.h
uses CombinedAdvectionDiffusionRLBdynamics and regularized boundary momenta on
flat normals; corners follow its stock fallback, and unsupported normal types
retain the installed behavior. Do not write a boundary algorithm. This API is
used by the installed porousPlate3d example and other examples.

Keep the ramp-1-transit input, stock combined regularized TRT bulk, bounce-back
walls, zero-gradient outlet, native flow, physical Dm and geometry unchanged.
First run each particle size at N32 for two control-volume transits AFTER the
ramp, as an early boundedness diagnosis. These short trials cannot qualify a
tail, dispersion or process prediction. Record concentration-extrema locations
using existing geometry node indices; that diagnostic does not alter updates.
If the front is bounded, execute a new sixteen-transit trial and the full
prospective refinement/holdout protocol. If it remains unphysical, retain the
failed candidate and keep adsorption service quantities unavailable.

This tests an inlet contribution, not a claim that the root cause is known.
Existing OpenLB sources and all earlier cases remain unchanged.
