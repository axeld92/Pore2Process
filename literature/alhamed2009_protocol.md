# Alhamed (2009): comparison protocol

Recorded on 10 November 2026 before fitting the local column model to Figure 7.
This is a same-paper, held-out particle-size comparison, not validation against
an independently collected experimental campaign.

## Data and split

- Material: ZnCl2-activated date-stone carbon, kept separate from the
  coconut-shell carbon used in the accepted synthetic example.
- Figure 7: 10 cm bed, 11 mm column diameter, 150 mg/L phenol,
  43.8 mL/min, 22 +/- 1 degC. Open squares: 0.8 mm calibration;
  open diamonds: 0.45 mm holdout. Digitize marker centres, not fitted lines.
- Fixed source inputs: Table 1 bulk density 740 kg/m3, bed porosity 0.427,
  Langmuir capacity 0.1115 kg/kg, affinity 139 m3/kg. Interpret the explicitly
  labelled bulk density as dry mass per packed-bed volume; the implied
  particle-envelope density is 740/(1-0.427) kg/m3. This is not skeletal density.
- Fit only the 0.8 mm curve, using unweighted concentration-fraction residuals.
  Fit an effective LDF rate and axial dispersion. These are reduced-model
  parameters, not separately identified film and intraparticle diffusivities.
- Primary holdout: freeze equilibrium, inventory and dispersion; change only
  particle diameter and scale the fitted rate as (0.8/0.45)^2. A secondary
  scenario scales dispersion proportionally to diameter. Neither uses holdout
  observations in parameter estimation.

## Screens and uncertainty

The chosen engineering screens are fraction RMSE <= 0.10 and relative t50
error <= 20% on the held-out curve. Report the 5% outlet crossing separately
with a 25% relative screen when the observed crossing is resolved away from
zero. These are engineering tolerances, not statistical confidence tests.
No experimental pointwise uncertainty is reported for these plotted points.

Retain the original embedded raster, its hash, pixel axes and manually selected
marker centres. Use +/- 2 pixels in time and concentration as an explicit
figure-reading sensitivity, not experimental error. Show crossing-time ranges
from these bounds. Near-zero times on the 0.8 mm curve may be unresolved at
this resolution; do not turn them into precise relative-error checks.

Use the existing conservative Julia column solver. Check normalized mass
residual <= 1e-8. On both curves compare the baseline with doubled axial cells
and halved time step; require < 2% changes in t50 and resolved t05 before
interpreting the model comparison. Preserve failures instead of relaxing
screens. A failed diameter transfer is a result, not a validated adsorption law.

## Source audit and exclusions

Transcribe Table 3 batch fits and Table 5 column fits separately from measured
points. Do not treat the paper's fitted dispersion as an independently measured
hydrodynamic closure or its pseudo-second-order batch constants as Deff.
Check the printed Equation 10 against the plotted curves before reusing it.
Figure 6 says 20 cm while several matching Table 5 concentration rows say
10 cm; do not silently choose a height for that concentration series.
Table 5's two lowest fit R2 values contradict the blanket high-R2 claim.

No pressure-drop measurements are provided. The synthetic FCC permeability
and DWSIM pump results remain separate conditional examples. This comparison
does not calibrate their absolute hydraulic predictions or claim an industrial
particle-size optimum. Do not rerun OpenLB or DWSIM merely to relabel synthetic
results as experimental validation.
