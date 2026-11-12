# Prospective smooth-input diagnosis

10 November 2026, after observing the in-progress sharp-front d800 case but
before building/running this candidate. The sharp-front baseline is retained:
its recorded local extrema already include negative and >1.3-feed values.
This addendum does not change its tolerances or promote that response.

Test whether an under-resolved imposed front is an important contributor. Keep
the same stock bulk/wall/inlet/outlet operators, physical geometry, flow and
Dm=1e-9 m2/s. Replace only the inlet time history with the installed OpenLB
PolynomialStartScale, clamped at its endpoint by the application in the same
way as official example start-ups. Record the actual imposed concentration at
each interval midpoint. Use ramp durations .5 and 1 physical unit-cell transit
for BOTH particle sizes, independently of adsorption. This is a controlled inert
tracer input, not a smoothed adsorption prediction or an increased diffusivity.

Run to 16 control-volume transits AFTER the ramp. Compare two measured responses
after removing the KNOWN inlet mean/variance: for a linear time-invariant system,
mean_response=mean_transport+mean_input and variance_response=variance_transport
+variance_input. Evaluate the actual sampled input with independent Julia
moment quadrature. The stock quintic ramp has derivative Beta(3,3), hence input
mean=T/2 and variance=T^2/28; check sampled values against these expressions.
Do not divide Fourier signals or differentiate noisy data, clip responses,
invent a tail or fit any adsorption curve. Outlet absolute flux and its tail
remain independently screened. Baseline sharp-front behavior is not erased.

Keep the main protocol's 2% concentration, positivity, mass and plateau screens.
Adjusted transport moments from the two ramp durations must agree within 5%
(mean) and 10% (variance); otherwise the ramp has not produced a robust tracer
measurement. A small tail contribution alone does not bound an unobserved tail.
Before any closure/process acceptance, retain spatial/time refinement and an
unfitted longer-segment prediction. Report a failed candidate as failed.

Source reference: installed OpenLB 1.9r0 analyticalF.hh PolynomialStartScale and
the ordinary example start-up lifecycle; no vendor or collision code is changed.
