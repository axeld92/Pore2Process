# Exchange between fast and slow paths in a repeated porous cell

## Independent inputs

Use the archived, balanced N96 TRT carrier field and the same two
repeated synthetic FCC cells. Restore the paired internal seam. Use
the same no-through-flow solid walls, open ends, liquid flow, and
1e-9 m²/s liquid diffusivity. No new hydraulic solve. The previously
reported two-cell curve is a fixed benchmark and its aggregate results
are already known. No sample or statistic from that curve enters the
new fit or choice of exchange strength. Lock the predictions before
opening the reference history for scoring.

## Derive exchange from marker trajectories

Inject equal-mass markers at the inlet plane, selecting each face in
proportion to its archived inward flow and distributing points uniformly
on the selected face. Interpolate each normal velocity linearly between
opposite faces within a voxel, as in the archived tracer solver. Integrate
advection exactly within that linear field and change voxels at face
crossings. Use symmetric advection/diffusion splitting. Advance
molecular diffusion by Gaussian displacements in each coordinate with
specular reflection at closed faces and at the two exterior end faces.
Allow diffusion across the restored interior seam. Measure the first
interior crossing and the final advective exit time for each marker.
The second duration includes any subsequent excursions back across
the interior seam.

Use 65,536 markers in each case, with an 80 µs primary step, a 40 µs
check and an independent seed at 80 µs. Cap tracking at 8 s and
require at least 99.9% of markers to cross both periods. Retain and
report censoring. Run a matched zero-diffusion advection control to
distinguish redistribution by flow alone from changes caused by the
molecular-displacement term. Use independent analytical advection
and reflection checks. Compare the raw first-period pulse with the
archived one-cell response using the arrival, rise, mean, variance,
and 99% recovery limits. The raw first crossing differs slightly from
an open single-cell exit because diffusion can cross the internal seam.

## Predict the two-cell outlet without a refit

Use the existing frozen N384 one-cell delay distribution as the
marginal distribution in both periods. Transfer only the measured
pairing of first- and second-period time ranks from the trajectories.
Divide the one-cell delay distribution into equal-probability bins
with one bin per completed marker, and use each bin's conditional
mean. Map the two ranks independently to these bins. This preserves
total mass and the mean liquid transit per period and retains the
flow-derived correlation. Add the two mapped times and convolve
with the same 4 ms inlet pulse. Compare against complete persistence
with the same marginal distribution, and apply the same thresholds
as the previous [complete holdout](../residence_time_holdout/protocol.md).

Report rank correlation and the probabilities of remaining in the
fastest and slowest deciles. Check sampling and time-step sensitivity
against 0.005 maximum cumulative difference, 5% relative variance and
5% relative 99%-recovery-time difference from the primary prediction.
Report both the conservative predictor's overall qualification and
whether its change improved the late tail without degrading the
early front beyond the established limits. The zero-diffusion control
shares the same calibrated marginal to isolate the effect of changing
rank pairing. No adsorption or DWSIM calculations are included.
