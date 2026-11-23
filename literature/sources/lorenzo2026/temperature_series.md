# Full equimolar temperature series and activity-loss challenge

[Lorenzo and Ortega (2026)](https://doi.org/10.1021/acs.iecr.6c02375),
main Tables 1–2, Sections 3–4; [Supporting Information](https://doi.org/10.1021/acs.iecr.6c02375.s001)
Table S2 pp. 5–6. Supplement SHA-256 and other batch data in
[the earlier extraction](README.md).

[`batch_1to1_temperature.csv`](batch_1to1_temperature.csv) transcribes
all 50 printed rows of assays 13–16, including starts (18, 13, 12,
7 rows). The 1:1 nominal feed is 23.7 mL pentanol and 16.3 mL acid,
40 mL nominal charge, 10% catalyst by initial acid mass. Preserve
printed compositions and conversion; no correction to the time or
chemistry columns. Combined with the preserved 40 °C equimolar and
full 2:1 series, these are all sixteen immobilised lipase assays,
225 rows including starts. Assays 17 and 18 use Amberlyst-15 and
are intentionally not enzyme-kinetic data.

[`apparent_equilibrium_1to1.csv`](apparent_equilibrium_1to1.csv) comes
from main Table 2: Kx = 11.15, 12.52, 17.46 for assays 7, 13, 14
(40, 50, 60 °C). The authors give **no equilibrium constant**
for 1:1 at 70 or 80 °C (assays 15, 16); their final conversions
are 0.435 and 0.185 at 270.1 and 270 min. The authors discuss
enzyme inactivation and exclude these two assays in their
fourteen-assay nonlinear calibration (§4.3). The apparent
equilibrium relation for 1:1 above 60 °C is an extrapolation
of active low-temperature data, not a measured plateau. The
2:1 series remains active through 80 °C and has separately
reported apparent Kx values. Table S3 estimates an empirical
temperature dependence for each feed. Table S4's author
Ping-Pong fit is a multivariable correlation, not an
independent time-resolved enzyme-activity assay. Eq. 17
relates an apparent denaturation constant to starting acid,
feed ratio and temperature; do not reinterpret its fitted
coefficients as the rate constant of irreversible loss.

The present audit uses a separately qualified reduced
reversible rate model and explicitly distinct clock-time
activity hypotheses. Do not use a low final conversion to
fit an apparent equilibrium at 70/80 °C. The kinetics and
apparent Kx from the originating curves are not independent
measurements; a withheld assay's Table 2 Kx must be withheld
too. No individual activity histories, catalyst reuse identity,
actual water activity or independent uncatalysed blank is
reported for these experiments.
