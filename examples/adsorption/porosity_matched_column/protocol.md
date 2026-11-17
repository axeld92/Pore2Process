# Transfer of the 42.7%-void liquid film to the 25 cm phenol column

Use the 23.3 and 62.5 mL/min, 150 mg/L cases with 0.80 mm carbon.
Freeze the earlier 25 cm column feed, 17.581 g of dry carbon, liquid inventory,
particle size, particle diffusivity and Langmuir equilibrium. Preserve
the previous 50%-direct access scenario and its additional resistance on
the other half of the carbon. Recalculate the slow film coefficient from
the new liquid film in series with that fixed resistance. Use fully mixed
joins between the same eight-section repeats and the final four sections.

Use the carbon-area-normalized film coefficients and the corresponding
transport and wall-contact rules from the N45 and N60 TRT runs, both at
20 µs tracer step. Compare the two paired cases, with N60 as the primary
prediction. Also cross the two film coefficients with the two sets of
transport and contact rules. This measures film-only sensitivity at fixed
transport as well as the effect of changing the transport and wall-contact
rules at fixed film. The earlier correlation-based forecast serves as an
archived comparison.

Use the existing Julia column solver and its time and particle-shell
refinement checks. Verify phenol mass, interface flux, total carbon and
cohort inventory. Seal the inputs and every predicted curve before using
the digitized Figure 4 markers. For 23.3 mL/min, use the duplicate 25 cm
curve from Figure 5 as a reading check. Score time to 5% and 50% of feed,
curve RMSE at the markers, retention through the last marker and retention
at the shorter shared feed dose. Keep the previously declared ±25%
early-crossing, ±20% midpoint and 0.10 RMSE limits and show retention
against the published marker reading range. Use no breakthrough data for
parameter selection. Figure 7 remains unfitted.

[Alhamed Figures 4–5 and Table 5](https://doi.org/10.1016/j.jhazmat.2009.05.002)
