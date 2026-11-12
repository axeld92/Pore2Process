# Flow- and inventory-preserving tracer translation

Declared before evaluating this stricter alternative. Use only the accepted
revision3 native data. Preserve the native run and its two-moment shape-only
candidate unchanged. This stage launches zero new native/process simulations.

For each native segment take actual fluid volume V, transport length L and flow
Q. Fix pore velocity v=Q*L/V and fluid area Afluid=V/L. Hence v*Afluid=Q and
Afluid*L=V exactly. Do not fit velocity, area, inventory, capacity or delay.
Determine only positive D from the SHORT bed32 inert-tracer mean, using the
independent absorbing-end analytical moment equation. No adsorption data or
longer-segment response is used. Reject if no representable positive D exists.

Freeze D for the longer length32 segment; use its measured Q,V,L. Test its
normalized profile (RMSE<=0.03), mean (<=10%), variance (<=20%), plus absolute
steady EXIT FLUX (<=10%) at unit inlet concentration. Also test the short
segment's variance/shape/flux; only its mean is a calibration target. Check
Julia's time/space refinement (profile change<=0.005, moment changes<=2%),
its native Q/V identities and per-area mass balance (<=1e-8).

Scale Julia's area-specific inventories and fluxes by the declared Afluid;
verify feed minus exit equals stored mass without a correction. The native
wall-population inventory remains a separately reported numerical artifact,
not physical solid uptake or added fitted liquid volume. This is a bounded
synthetic, absorbing-boundary contact closure. No intrinsic material-diffusivity,
experimental adsorption, realistic-Peclet or DWSIM service prediction is claimed.

Save a unique translations run with the immutable native-summary hash, software/
protocol snapshots, exact SI model inputs, full Julia outputs and diagnostics.
Independent --check replays Julia and rederives D, without running OpenLB.
