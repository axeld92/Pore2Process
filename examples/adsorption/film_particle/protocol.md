# Finite-bath and film/particle comparison protocol

Recorded 10 November 2026 before fitting the new transport model. Preserve the
previous Langmuir/LDF comparison, its failed physical screens, and all earlier
examples. This is a separate homogeneous-particle diffusion approximation,
not a reconstruction of the carbon's connected micropores.

## Evidence and partition

Use Alhamed (2009), DOI 10.1016/j.jhazmat.2009.05.002, date-stone carbon only.
Table 1 gives qmax=0.1115 kg/kg, b=139 m3/kg, bulk density 740 kg/m3 and
bed porosity 0.427. The implied particle-envelope density is
740/(1-0.427) kg/m3. Keep these source values unchanged in the primary model.

Audit all seven Table 3 fitted equilibrium loadings with the finite-bath
balance V*C + m*q(C) = V*C0. They are published fit endpoints, not independent
raw equilibrium measurements. A 5% relative discrepancy is a chosen audit
screen. The preliminary source audit is within 2.1%, so no equilibrium factor
is needed in the primary model. Also report a diagnostic capacity multiplier
fitted only to the training rows (0.8 mm/200 mg/L and 1.5 mm/200,400 mg/L),
without using it to silently replace the source isotherm.

Digitize Figure 2 experimental markers, not the pseudo-second-order fitted
lines. Conditions: 1 g dry carbon, 250 mL liquid, 200 mg/L phenol, 200 rpm,
22 +/- 1 degC. Its y axis starts at 30 mg/g; early uptake below that level is
unobserved. Omit unresolved overlapping symbols. Retain pixel centres, raster
hash, axes and a chosen +/- 2 pixel reading sensitivity. This is not an
experimental confidence interval. Liquid sampling uses 1 mL aliquots; the
baseline neglects volume/solute withdrawal. A recorded serial-withdrawal
sensitivity at the selected digitized times checks its possible influence;
it cannot reconstruct unreported sampling times.

Fit one shared effective intraparticle loading diffusivity and one batch film
coefficient at 0.8 mm using Figure 2's 0.8 and 1.5 mm curves only. Assume a
constant batch Sherwood number, kf(dp)=kf_ref*(0.8 mm/dp), as an explicit
transfer assumption. Do not fit a particle-size exponent. Withhold all 0.45
and 0.225 mm Figure 2 points, plus the other Table 3 conditions. The 1.47/1.5
mm label discrepancy is retained as a diameter-reading scenario.

Batch engineering screens: loading RMSE <= 3 mg/g on each calibration and
held-out size, and the 5% endpoint audit. Fast-particle curves already close
to equilibrium at the first visible point cannot identify their early uptake
time. Fits do not establish intrinsic diffusivities or the validity of the
chosen batch film scaling.

## Conservative alternative

Use spherical shells for an effective adsorbed loading q(r,t), with constant
Ds and centre symmetry. The surface obeys q_s=qmax*b*Cs/(1+b*Cs).
The same interface flux satisfies

    J = kf*(Cbulk-Cs) = rho_particle*Ds*dq/dr at r=R.

Couple its shell-volume-weighted mean loading to a well-mixed finite bath or
the conservative axial column balance. Neglect separate intraparticle liquid
inventory explicitly. Use implicit finite volumes; never correct mass after
solving. Python coordinates and plots; Julia owns all forward models.

Verify a linear-isotherm film limit against an analytical finite-bath solution,
a fixed-surface spherical-diffusion series, nonlinear equilibrium, zero
exchange, nonnegative states, conservation and radial/axial/time convergence.
Require normalized mass residual <= 1e-8; selected final loading/crossing
refinement changes <= 2%, with short early windows if needed. Do not relax
screens to promote a fit.

## Column prediction and qualified anchor

Freeze the batch diffusivity and source isotherm before any Figure 7 use.
Keep the reported-density-derived 7.032455 g bed inventory distinct from a
measured per-run mass. Use Figure 7's 10 cm, 11 mm, 43.8 mL/min, 150 mg/L
conditions. Set a baseline packed-bed film coefficient using Equation 19's
Wakao/Funazkri form with explicitly declared superficial-velocity Reynolds
number. An interstitial-velocity convention is a separate sensitivity,
not a silent redefinition or a fitted batch-to-column coefficient transfer.

The primary column calculation predicts both Figure 7 curves without fitting
either. Baseline axial dispersion is zero; Dax=0.5*v_interstitial*dp is a
chosen positive-dispersion scenario, not a measured closure.

A secondary, explicitly conditional route may fit **one common column film
multiplier** to the 0.8 mm Figure 7 curve only. Freeze it before predicting
0.45 mm. No per-curve column fitting, inventory fitting, diffusivity refit or
equilibrium refit is permitted in that route. Bounds: multiplier 1e-4..100.
Shared batch Ds bounds: 1e-18..1e-9 m2/s; batch kf_ref: 1e-10..1e-3 m/s.
Report bound hits and film/diffusion trade-offs rather than claiming unique
physical identification.

Report fixed active-inventory fractions 0.5, 0.75 and 1 as separate loss-of-
accessibility scenarios, preserving physical dry mass and the source isotherm.
These are not measured uncertainty or fitted bed masses. They do not select
the accepted prediction.

Column screens remain fraction RMSE <= 0.10, t50 error <= 20%, and resolved
t05 error <= 25%. Use the saved Figure 7 reading bounds. Distinguish fully
transferred predictions from the one-anchor conditional prediction, and
numerical/software passes from experimental-comparison passes. A failed
alternative remains a valid reported result.

No pressure data, new experiments, extra raw datasets or specimen reconstruction
are prerequisites. No new OpenLB/DWSIM run is needed for this transport check;
it cannot validate the previous synthetic permeability or absolute pump duty.
Reaction, heat coupling and regeneration remain deferred.
