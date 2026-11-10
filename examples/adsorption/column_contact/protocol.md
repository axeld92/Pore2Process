# Frozen protocol: column inventory and contact audit

Recorded before this model's curve selection. Preserve the earlier finite-bath,
LDF and synthetic studies and their source hashes. No new laboratory data gate.

1. Reinspect Alhamed (2009) PDF pages 2, 5 and 6: Table 1's **bulk** density
   0.74 g/cm3 and voidage 0.427, Section 2.3 packing/distributor description,
   Figure 7 and Table 5. Density times geometry is a derived inventory, not
   the weighed mass of either Figure 7 bed. No fraction-specific mass, voidage,
   tracer residence-time distribution or pressure loss is available.
2. Integrate only the digitized Figure 7 observation window. A monotone clean-
   bed origin and linear segments are declared quadrature assumptions. Tail
   uptake remains unknown; finite-window retention is not saturation capacity.
   Calculate required missing tail if the printed inventory fully equilibrates.
3. Freeze the accepted batch Ds, particle density, source Langmuir isotherm and
   full active mass. Do not apply the previous fitted column film multiplier.
   Use the same declared superficial-velocity film correlation, multiplier 1,
   for both diameters. Table 5 DL/Pe remain **breakthrough-fit descriptors**;
   neither is an independently measured fluid-dispersion input to this model.
4. Test N equal mixed sections in series. Every section contains its share of
   the full liquid and solid inventory; only advective exchange connects them.
   There is no bypass, dead solid, additional film factor or fitted density.
   N is physical contact topology, not a mesh. It is never inferred from DL/Pe.
5. Verify analytically before selection: nonadsorbing Erlang step response,
   mean liquid residence time eps*L/u and variance tau^2/N, one-section CSTR,
   fast linear-equilibrium retarded response, zero flow, no-exchange limits,
   particle/interface and global conservation. Reuse the frozen particle code.
6. Select the one integer N in **1..64**, minimizing equal-weight concentration-
   fraction RMSE at the 24 **0.8 mm Figure 7** markers only. Search at Nr=16,
   dt=2 s; check every candidate. Reevaluate the selected N and its neighbors
   at Nr=32, dt=0.5 s; if the minimizer changes, evaluate all 64 at that level.
   Do not use 0.45 mm observations in selection. Do not add parameters after
   seeing holdout. N=1/64 boundary selection is an identifiability warning.
7. Freeze selected N for 0.45 mm. Retain the earlier engineering screens:
   RMSE <=0.10, t50 relative error <=20%, resolved t05 error <=25%.
   Image +/-2 pixel scenarios are reading sensitivity, not measurement errors.
   Plot/check 1, 4, 16 and 64 section scenarios for both sizes, not extra fits.
8. At fixed N refine radial shells and time step independently; require each
   t05/t50 change <=2% and max concentration-fraction change <=0.005. If early
   crossing needs a finer time window, save and qualify it separately. Changing
   N is model sensitivity and cannot be called numerical grid convergence.
9. At fixed frozen particle density, run porosity 0.35, 0.50 and 0.60 as chosen
   common packing scenarios for both sizes with no reselection of N. Enforce
   rho_bulk=(1-eps)*rho_particle; all solid remains active. These are not source
   uncertainty bounds. Do not reinterpret 0.74 g/cm3 as particle density in a
   forward run with the old batch fit: that would require batch requalification.
10. Save inputs, trajectories, source/artifact hashes, calibration scores,
    diagnostics and actual runtimes; independently replay. Keep numerical
    verification separate from physical comparison, even if the latter fails.
    No process-energy/material recommendation until column transfer qualifies.
