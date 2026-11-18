# Shared effective equilibrium sensitivity by bed depth

The two source equilibrium laws and their independently qualified
kinetics bracket the measured 25 cm retention at both flow rates.
Test a single effective equilibrium weight shared by every column
condition. Set the loading to the convex mixture
`q(C) = (1-w) q_L(C) + w q_F(C)` of the two printed Table 1 equations.
The candidate grid is `w = 0, 0.25, 0.50, 0.75, 1.00`. Interpolate
the already qualified primary loading diffusivities and stirred
reference film coefficients linearly in their logarithms. The
stirred film is used only for Figure 2. Keep the source terms and
weight separately traceable.

Predict all four batch sizes at each weight, with no new batch fit.
Use the existing 3 mg/g and numerical limits before allowing a
candidate into the column comparison. Freeze the input and model
settings for every eligible candidate.

Use both 10 and 15 cm, 23.3 mL/min Figure 5 conditions for one shared
calibration. Compare last-marker retained phenol per gram of total
carbon, using both curves with equal weight. Take the minimum sum
of squared errors. Use no crossing, curve-shape or other depth
information for weight selection. Freeze all candidate predictions
before reading the two training curves. Never refine the weight grid
against breakthrough results.

Freeze the selected weight before predicting the 20 cm low-flow
case and the 25 cm cases at 23.3 and 62.5 mL/min. The original
N60 flows, films, path rules, direct/slow fractions and access
resistance remain fixed. The short columns use the previously
archived depth-specific inventory and section counts. Scale the
fitted-clean-wall conductance in proportion to the number of
sections and carbon mass, keeping the film coefficient unchanged.
Use 48 shells and a 4 s main step or finer, with the established
refinement check. Seal predictions before scoring the held-out
curves. At every depth report 5% and 50% crossover, outlet/feed
RMSE and retention at the last marker and from the published
midpoint onward. The duplicate 25 cm Figure 5 markers check the
Figure 4 reading and are counted once. The weight is an effective
calibration and makes no claim about intrinsic adsorption equilibrium.

[Alhamed, Figures 2 and 4–5 and Table 1](https://doi.org/10.1016/j.jhazmat.2009.05.002)
