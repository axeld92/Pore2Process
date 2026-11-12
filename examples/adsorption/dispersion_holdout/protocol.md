# One-cell dispersion fit and two-cell inert-tracer holdout

Use the qualified archived N384 one-cell pulse as the only calibration
response. Fix liquid flow, volume, and physical length from the original
saved TRT reconstruction. Use a one-dimensional constant-coefficient
advection–dispersion equation with prescribed total inlet solute flux
and zero diffusive outlet flux. Fit one effective axial dispersion
coefficient to the measured one-cell exit-time variance. The rectangular
4 ms inlet pulse contributes its known variance, pulse width squared
divided by 12, to the modeled exit-time variance. Keep the modeled mean
fixed by liquid volume divided by flow, plus the pulse mean. Check the
finite-volume model against the existing analytical moments and refine its
mesh and time discretization.

Before opening any two-cell tracer results, save the fitted coefficient,
one-cell comparison and complete two-cell prediction, with hashes.
Repeat the N96 periodic geometry twice axially and restore the original
paired seam flows between the repeats. Keep the open boundaries only at
the two outer faces, the same TRT carrier flow, N288 tracer cells per
period, the same molecular diffusivity, inlet pulse and 80 µs time step.
Use the verified Julia brick kernel and a planned comparison window of
1.7368 s, twice the one-cell window. Keep hydraulic inputs and adsorption
parameters fixed. Run no full bed.

Score the held-out cumulative outlet mass fraction normalized by the
known injected mass on the fixed time grid. Use the existing 0.03 profile
RMSE, 5% mean and 10% variance limits. A complete passing response must
also conserve mass within 1e-8 of injection, obey the physical bounds,
recover at least 99.9% and contribute no more than 1% of either raw
moment in the last nominal transit. A profile failure is conclusive
once the square root of accumulated squared error divided by the full
planned sample count exceeds 0.03. All remaining squared errors are
nonnegative. After at least 10% of the tracer has exited, the solver may
stop at that point and report the verified lower bound, early arrival,
and mass balance. In that case full-window pore moments and tail
metrics are unavailable. No extrapolation, refit, adsorption prediction,
or DWSIM execution is included.
