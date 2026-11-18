# Solvent-free continuous pentyl propanoate synthesis and water-tracer tests

Lorenzo, Ortega, González Garcinuño and Martín-del Valle (2025),
*Industrial & Engineering Chemistry Research*, 64, 15311–15319.
[Article](https://doi.org/10.1021/acs.iecr.5c02780) and
[supporting information](https://acs.figshare.com/articles/journal_contribution/_/29631736).
The [extraction record](extraction.json) gives the source locations, PDF hashes,
operating conditions and digitization axes.

## Packed-bed conversion

[Conversion measurements](conversion_2to1.csv) contain the 20 solvent-free
pentyl propanoate runs at a 2:1 inlet molar ratio of pentan-1-ol to propanoic
acid. The feed concentrations were 6.6 and 3.3 mol/L, respectively. The authors
dried the reagents and describe measuring initial moisture. The residual
moisture in this feed is unavailable in the article and supplement.

The column was 16 mm in diameter and 140 mm long with 6.0 g of Lipozyme435.
Jacket water circulated at 12 L/min. The paper describes the reactive runs
as quasi-isothermal within ±0.1 °C of the set temperature. Steady state was
reported at about 30 min. Effluent acid concentration was determined by
titration and checked by gas chromatography. The paper estimates concentration
error below 5%. Individual conversion uncertainties are unreported.

The conversion CSV gives the measured fractional decrease in acid
concentration, `(C_in - C_out) / C_in`. It represents the point markers in
Figure 3c. Five markers obscured by nearby series in Figure 3c were read from
the filled-circle series in Figure S1, which plots the same observations
against temperature at each flow. The other 15 markers were cross-checked
against Figure S1 within 0.6 percentage points. Readings are rounded to
0.001 in fractional conversion. `source_run_id` follows SI Table S2. Marker
coordinates use the native image in the indicated figure. The evidence record
includes the axis calibration. These measurements are reserved for reactor
comparison.

## Catalyst and water-tracer data

Section 2.2 reports a mean bead diameter of 593 µm, with size-distribution
span 0.976. The text reports a dry basis, while Figure 2 is labelled wet and
the procedure uses distilled water. We retain the reported mean as a
representative size with that basis ambiguity. The parenthetical approximation
to the diameter in Section 2.2 has an inconsistent unit. The printed 6.0 g
catalyst loading is retained independently of the accompanying approximate
volume of 2 cm³. The geometric volume of the glass column is about 28.15 mL.
Applying the reported water-bed void fraction of 0.64 to that volume would
imply about 10.1 mL of bead envelopes. The measured dry pore characterization
and the correlation-derived wet porosities are recorded separately.

The [water-tracer measurements](water_tracer.csv) are digitized
concentrations from Figure A1a. The circles correspond to blue dextran
and the diamonds to sodium azide. Water flow was 1.2 mL/min with separate
8.4 mL tracer pulses and 1 mL outlet aliquots. Blue dextran, treated as
excluded from the beads, was prepared at 1 mg/mL. Sodium azide, which can
enter the beads, was prepared at 0.1 mg/mL. The nominal pulse length
would be 7 min for injection at the stated water flow. The inlet pulse
shape, timing convention and water temperature are unreported.

The concentration and time columns use the axes as plotted. Many blue
dextran markers lie about 1 min apart on the graph. An uninterrupted
1 mL sample at the reported flow would take 0.83 min. The observed time
axis is retained. Marker centers were visually selected and checked
against the embedded raster. Values indistinguishable from the baseline
at this image scale were entered as zero. The retained resolution is
0.01 min and 0.0001 mg/mL for reproducible calculation.

Trapezoidal integration of the sampled profiles gives normalized means
of 19.09 and 24.30 min and variances of 27.78 and 55.82 min² for blue
dextran and sodium azide, respectively. Table A1 reports means of 19.07
and 24.27 min and variances of 27.83 and 56.84 min². Its moments normalize
the measured finite-pulse concentrations according to Appendix A. The
Péclet numbers, axial dispersion coefficients and wet void fractions in Table A1 are
correlation-derived values. The tracer measurements were made with water.

The digitized concentration integrals and printed water flow imply about
6.72 mg of blue dextran and 2.48 mg of sodium azide at the outlet over the
plotted intervals. The reported pulse volumes and concentrations imply
input doses of 8.4 and 0.84 mg, respectively. The concentration plot
therefore supports a normalized shape and timing comparison, with an
unresolved absolute-dose discrepancy.
