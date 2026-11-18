# Batch-supported particle kinetics

Reuse the conservative Julia finite-bath model with Alhamed's reported
Langmuir fit and the established Figure 2 split. Set the reference batch
film to 0.8 mm, with the same inverse-diameter scaling for all four sizes.
Run the archived shared fit and its high-film profile. On the low-film
side, test reference films of 40 and 45 µm/s and minimize the pooled
training error over log10(D) from -14 to -9 with the existing bounded
search and the original training discretization. The 40 µm/s point
checks the edge of the sampled trade-off. The 0.45 and 0.225 mm series
are withheld from every parameter search.

Freeze the resulting parameters and the four forward trajectories
before evaluating the markers. Use the established 3 mg/g center RMSE
limit for each size and retain the four ±2-pixel reading checks. Screen
mass balance and interface flux at 1e-8. Compare the loading at each
marker with a radial-shell doubling and a step halving. Require each
change to be within 2% of the reference finite-bath equilibrium
loading. Use the lowest and highest diffusivities among the qualified
profiles for the column test. This is a sampled sensitivity range.

For each transferred diffusivity, use both saved N60 25 cm conditions.
Retain their unchanged liquid and carbon inventory, isotherm, external
film, additional slow-cohort resistance, flow and geometry. Change
only D. Replay the existing time and radial refinement and require at
least the original selected resolution. Freeze and hash both flow
predictions at every tested resolution before scoring Figures 4–5.
Keep the same published 5% and 50% crossing, RMSE and last-marker
retention checks. Compare uptake added after the *published* midpoint
through the last marker using the measured outlet integral and the
simulated conservative mass inventory. Use the duplicate low-flow
curve only as a reading check.

Audit the source bulk-density-derived inventory and both published
isotherms against the seven Table 3 finite-bath endpoint fits. Compare
predictions at the 150 mg/L column feed as a separate equilibrium
sensitivity. Keep the column isotherm fixed in the kinetics transfer.

[Alhamed, Tables 1–3 and Figures 2, 4–5](https://doi.org/10.1016/j.jhazmat.2009.05.002)
