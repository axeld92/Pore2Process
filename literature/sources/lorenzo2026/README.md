# Lorenzo and Ortega (2026) batch esterification

Beatriz Lorenzo and Juan Ortega, *Industrial & Engineering Chemistry Research*
65 (32), 17092–17102. [Main article](https://doi.org/10.1021/acs.iecr.6c02375)
and [Supporting Information](https://acs.figshare.com/articles/journal_contribution/Influence_of_Temperature_Initial_Reactant_Ratio_and_Catalyst_Type_and_Concentration_on_a_More_Sustainable_Synthesis_of_Pentyl_Propanoate_in_a_BSTR/33164325).
The article is open access under CC BY 4.0. The supplement file is
`ie6c02375_si_001.pdf` from ACS Figshare item 33164325, file 67271163. Its MD5
matches the publisher's `c8e96383ec826d3710c3b507e4484e44`.

## Batch measurements

[`batch_2to1.csv`](batch_2to1.csv) transcribes the complete time series at the
nominal initial 2:1 molar feed ratio of 1-pentanol to propanoic acid, using the
Lipozyme 435 immobilised lipase. The solvent-free reaction yields pentyl
propanoate and water. Operating conditions are from main Table 1. Time, mole
fractions (`x_acid`, `x_alcohol`, `x_ester`, `x_water`) and dimensionless acid
conversion are from Supporting Information Table S2, pages 3–5, including the
printed initial composition and zero conversion at time zero for every assay.
Time is elapsed minutes, following the main article's time axis. The decimal
comma in assay 3 at 5.3 min is normalised to a decimal point in the CSV.

| Assay | Temperature (°C) | Catalyst loading (% of initial acid mass) | Rows including start |
| --- | ---: | ---: | ---: |
| 1 | 40 | 2.5 | 11 |
| 2 | 40 | 5 | 11 |
| 3 | 40 | 10 | 16 |
| 4 | 40 | 20 | 19 |
| 9 | 50 | 10 | 15 |
| 10 | 60 | 10 | 17 |
| 11 | 70 | 10 | 15 |
| 12 | 80 | 10 | 17 |

The article reports a 50 mL jacketed stirred vessel with a 40 mL charge, a
nominal 29.7 mL of 1-pentanol and 10.3 mL of acid in the selected assays, and
time zero at acid addition after temperature equilibration. Catalyst loading
is mass of immobilised preparation divided by initial mass of propanoic acid.
Supporting Information Table S1 gives measured pure-liquid densities at 25 °C
of 810.92 kg/m³ for 1-pentanol, 988.07 kg/m³ for propanoic acid and 867.93 kg/m³
for pentyl propanoate. Applying the pure-acid density to the nominal charge
gives approximately 10.177 g initial acid and 1.018 g immobilised preparation
for a 10% run. The calculation uses the density at 25 °C as a reference for
the loading and provides no operating-temperature mixture density.

The authors report duplicate runs with no significant differences, a
temperature-probe accuracy of ±0.2 °C, and estimated absolute uncertainties
of ±0.001 for mole fraction and ±0.002 for acid conversion. They estimated
water by mass balance and did not publish individual replicate time series.

### Consistency checks

All 121 rows were checked against the PDF table cells and rendered pages. The
CSV preserves the values as printed. The article defines acid conversion as
`1 - x_acid / x_acid_at_start`. Applying that formula to the printed
compositions gives differences up to 0.04661 in conversion (assay 9, 13.0 min).
The conversion reconstructed from ester produced, `x_ester / x_acid_at_start`,
agrees more closely in most rows, with a largest difference of 0.01940 (assay
2, 44.8 min). The largest deviation of the four mole fractions from a sum of
1 is 0.018 (assay 4, 168.4 min). These differences exceed the reported assay
uncertainties. Use the printed conversion as the primary kinetic observation.
Any re-derived conversion, composition closure correction or omission of
rows needs its own traceable sensitivity check.

## Equilibrium and calorimetry

[`apparent_equilibrium_2to1.csv`](apparent_equilibrium_2to1.csv) separately
transcribes the dimensionless 2:1 apparent mole-fraction equilibrium constants
from main Table 2, using the authors' assay mapping. The authors estimated
these constants from the corresponding batch runs and used a mole-fraction
model that absorbs liquid nonideality. The values share experimental evidence
with the kinetic time series. The (K_x) fit in Supporting Information Table
S3 is a separate regression to Table 2. Supporting Information Table S4
contains the authors' kinetic fits covering other assays as well as these
2:1 runs.

Main Section 5 reports a reaction enthalpy of −5.62 kJ per mole of ester. The
Calvet calorimeter run for this value was at 25 °C, a 1:1 initial molar ratio
and 10% catalyst loading, as described in Supporting Information Section
S.2.4 and Figure S5. It belongs to a separate calorimetric experiment and is
outside the selected kinetic series. Temperature-dependent apparent
mole-fraction equilibrium constants and the measured calorimetric enthalpy
therefore have distinct evidence and meanings.

## Additional equimolar loading series (40 °C)

[`batch_1to1.csv`](batch_1to1.csv) transcribes every point in SI Table S2,
assays 5–8 (printed pages 3–4), from the same verified supplement. The four
assays have respectively 11, 11, 17 and 15 rows including zero time. Main
Table 1 specifies 23.7 mL of pentanol and 16.3 mL of acid, a nominal 40 mL
charge and 2.5, 5, 10, 20% catalyst **relative to initial acid mass**.
Nominal concentrations from the pure-liquid acid density at 25 °C are
10.07, 20.13, 40.26 and 80.53 g/L. These are calculated from the charge,
not measured operating-temperature mixture concentrations. Use the same
reported conversion uncertainty and printing/closure checks as above.
Main Table 2 gives apparent Kx = 11.15 from assay 7 (40 °C, 10%); it is not
an independent equilibrium determination for prediction of that assay.
SI SHA-256: `b08bd59fb16a005584a6a48104f289b3a55f70165897a037cc5e9823d9f38285`.
