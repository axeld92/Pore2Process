# Predict both Alhamed Figure 7 curves with the fixed path exchange

Use the qualified one-cell delay distribution and marker rank pairs from
the archived TRT exchange study. Coarsen the ranks using the first of
20, 40, 80 or 160 classes that stays within 0.005 maximum absolute
cumulative error and 5% relative error in 99% recovery time and variance
against the unbinned two-cell prediction. Use finer probability bins
in the slow tail. Each class uses the conditional mean delay and its
observed flow share. Flow between two classes is their measured joint
rank fraction times the total flow. This preserves the flow shares,
total void volume and mean transit.

Use a stationary first-order transition between successive periods.
Transfer the same normalized delays and transition probabilities to both
diameters. Scale period length by particle diameter using the synthetic
FCC reference length per diameter. Choose the nearest whole number of
periods for the reported 10 cm bed, then use its exact length. Rescale
the delay units by the void residence time per period at the reported
flow and porosity. The same rank transfer is a scenario assumption
across packing, Peclet number and much longer distances. No reconstruction
of the experimental packing or validation of that transfer is implied.

For a class of flow Q_i and delay T_i, assign liquid volume Q_i*T_i.
Allocate dry carbon proportional to that liquid volume. This assumes
the reported solid-to-liquid ratio throughout the path classes. Use
the full printed-density-derived dry mass and separate physical void
volume. Solve a conservative finite-volume liquid balance along each
class. Mix only at the period interfaces. In every liquid control
volume, couple the unchanged FilmParticle spherical diffusion operator
and Langmuir equilibrium to that volume's solid with the same implicit
film flux. Use the batch-only effective loading diffusivity and the
previously declared, unfitted column film correlation with superficial
velocity and multiplier one. Use the original Figure 7 operating
conditions and keep the batch film coefficient separate. Apply no
experimental curve-specific parameters, fitted dispersion, active
inventory multiplier, or inlet holdup.

Run through 18,000 s on both sizes. Initial constant-step checks found
that the computed 0.80 mm startup requires smaller steps. Use a first
phase of 24 s with shorter steps while retaining the same liquid and
particle state at the transition. In the 0.80 mm case compare 2/4 path
subdivisions, 96/128 radial shells, 0.04/0.02 s startup steps and 4/2 s
main steps. Start the 0.45 mm comparison with 1/2 subdivisions, 24/48
shells, 0.08/0.04 s startup steps and 4/2 s main steps. If needed use
4 subdivisions, 96 shells, 0.02 s startup steps and a 1 s main step.
Require 0.005 maximum difference in normalized outlet and 2% relative
changes in t05 and t50 under refinement. Require
normalized total mass and film/radial flux errors each below 1e-8 and
physical concentration and loading bounds. Use identity mixing on the
same path shares and volume allocation as the no-exchange control.

Freeze the inputs and simulated curves before scoring Figure 7. Use the
existing source markers, crossing interpolation and acceptance screens:
fraction RMSE at most 0.10, t50 relative error at most 20%, and resolved
t05 relative error at most 25%. Record any figure reading sensitivity
separately. Both diameters are predictions. Treat the printed bulk
density as a bulk density. Retain the distinction between source
equilibrium, derived inventory, batch calibration, and synthetic flow
transfer. No new OpenLB or DWSIM runs are part of this comparison.
