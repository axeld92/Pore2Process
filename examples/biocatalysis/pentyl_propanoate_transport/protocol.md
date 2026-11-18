# Liquid readout and bead transport protocol

Use the frozen 2:1 feed, 6.0 g catalyst inventory and 20-run
[reactor matrix](../pentyl_propanoate/protocol.json). The reactor
observations are released to the scoring functions after each runner
writes its forward predictions.

1. Apply published pure densities and binary excess-volume measurements to
   the fixed full-inventory 2025 molar conversion predictions. Evaluate
   complete liquid samples, omission of generated water and selected-phase
   readouts based on ternary measurements. Treat the ternary projection as
   a sensitivity with unresolved four-component phase behavior.
2. If the liquid measurement effects are small on the scale of the
   discrepancy, test a change of internal transport at the original
   catalyst mass. Use published pore-volume estimates, analogue optical
   particle sizes and an outer-shell distribution. Translate changes
   using first-order spherical effectiveness with the source Table B1
   reference. Check the dry pore volume on a skeletal-density basis.
   Apply relative rate changes to both frozen reactor models.
3. Compare with the uniform-sphere strong-diffusion limit and with the
   earlier catalyst exposure multiplier. Combine the two porosity
   conventions with diameter change and complete water omission.
   Seal these forward comparisons before using any measured reactor
   conversion. Retain all adverse or unavailable results.

Report the phase basis, rate factors, acid balance, grid sensitivity and
held-out errors. Preserve the earlier example files and rate coefficients.
