# Wall-weighted transfer into shared carbon

Use the frozen 40-bin advection-ranked wall-exposure shares and the
size-specific Péclet-matched exchange matrices from the previous
contact-column comparison. Keep the same liquid flow shares, delay
distribution, source geometry, flow, total dry mass of 7.032455155 g,
bulk void volume, equilibrium, film correlation and batch loading
diffusivity. The advection and diffusive path ranks are aligned by
their flow quantiles, as in that comparison. No Figure 7 point
selects a parameter or scenario.

Use one common radial particle loading profile per axial control
volume. Divide the dry mass evenly along the bed. Set total external
sphere area to 6M/(rho_p*d_p) from the dry inventory and diameter.
Divide that area among the flow classes in proportion to their
archived exposure shares. Keep the physical film coefficient per
unit area and the total film conductance fixed. Each path has a
separate film and surface liquid concentration. The shared radial
operator receives the area-weighted sum of surface fluxes. Uptake
and release are signed and are exactly opposite in the liquid and
solid balances. This assumes instantaneous averaging within each
radial shell, with no direct mixing of the liquid films. In the
inert limit the original carrier exchange remains the only mixing.

Use a shared-carbon control with area proportional to path flow
shares. It uses the same physical carbon and total area and
isolates the effect of weighting the conductance by wall contact.
The control and the wall-weighted case both retain the same
Péclet-matched exchange. A class with zero sampled contact has
zero direct interfacial conductance in the wall-weighted case.

Run both sizes to 18,000 s. Compare successively refined finite
volumes and radial shells, with shorter startup steps during the
first 24 s. Require normalized mass and interface flux errors
below 1e-8 and physical concentration and loading bounds. Require
2% relative convergence in times to 5% and 50% of feed, with a
maximum normalized curve change of 0.005 after 24 s. If crossings
within the first 12 s need a separate fine computation, retain the
full and short results and verify the difference at their join is
at most 0.005 of feed before using both for comparison. Flag a
crossing that remains numerically sensitive on a concentration
plateau.

Freeze all inputs and outputs before scoring the existing Figure 7
markers for both sizes. Use the earlier curve RMSE ≤0.10, midpoint
error ≤20%, and resolved early-crossing error ≤25%. Compare to the
archived separate-carbon predictions to assess the effect of sharing
solid loading. The external contact and mixing remain assumptions
for the experimental bed. Run no new OpenLB or DWSIM cases.
