# Pentyl propanoate transfer audit protocol

Use the same 2:1 pentan-1-ol-to-propanoic-acid feed, 6.0 g Lipozyme 435
loading, 3.3 mol/L inlet acid concentration and 20 temperature/flow
conditions as the [original benchmark](../pentyl_propanoate/protocol.json).
Predict first and open the experimental conversion CSV for scoring.
The source's concentration-based conversion is compared to molar
conversion under the same constant-liquid-volume assumption.

1. Integrate the normalized blue-dextran water outlet trace. Under
   independent convolution, subtract the nominal 7 min rectangular
   inlet moments. Include an impulse at the start and one at the end of
   that assumed inlet window as sensitivity cases. Match the relative
   variance to the existing Julia axial-dispersion model. Use the
   frozen 2026 batch kinetics and compare with the original plug-flow
   result. Treat water-to-reacting-liquid transfer as a scenario.
2. Use the 2025 printed rate coefficients and internal effectiveness
   factors. Reconstruct acid consumption with the printed fixed inlet
   concentration and zero outlet gradient. Calculate catalyst
   concentrations from the same weighed mass on bulk and liquid
   volume bases. Check the effect of using the full catalyst exposure,
   dropping the internal effectiveness factor, changing the inlet to
   a flux condition and using the low-Reynolds-number correlation.
   Check the rate time basis using the 2025 batch markers. Compare
   published Figure 6c calculated lines and experimental markers
   separately. Treat the source basis inferred from the line
   comparison as a reconstruction.
3. Test rate transfer without imposing a kinetic law when the
   measured-tracer scenario leaves the conversion gap large. For
   plug flow with proportional catalyst activity, match batch
   time to total bed catalyst mass divided by inlet acid molar flow.
   Interpolate both measured batch series at each equivalent time.
   Use Table S2 concentration loading for 2025 and the nominal 10%
   acid-mass loading for 2026. Repeat the 2026 transfer using
   conversion recalculated from the printed acid mole fractions.
   Leave forecasts outside the source time range unavailable.

Report signed errors and RMSE in percentage points. Compare the direct
transfer routes on the same 15 conditions as well as over the points
covered by each source.
