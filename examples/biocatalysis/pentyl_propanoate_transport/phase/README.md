# Does the conversion readout explain the discrepancy?

**Bottom line.** At full 6 g loading, the measured liquid-volume inputs give
at most 0.06 percentage points of full-liquid readout shift. Omitting every
mole of produced water lowers conversion by at most 1.55 points and leaves
all 20 predictions high. An ester-excluded ternary phase projection shifts the
pentanol-rich readout by at most 0.27 points. The tested readout effects are
small compared with the reaction discrepancy. Actual quaternary phase behavior
and assay-temperature properties remain unresolved.

## What was measured

The authors sampled the effluent after 10 min and then every 30 min and report
steady state after about 30 min. Section 2.3 describes NaOH titration of
remaining acid with phenolphthalein, a GC check and a check for interference
from ester hydrolysis. The stated concentration error is under 5% relative.
The reported conversion is `1 - C_sample/C_feed`. The paper leaves the sample
temperature, phase selection and outlet volumetric flow unreported. The 1 mL
aliquot size in Section 2.2 belongs to the water-tracer runs. The authors dried
the reagents without reporting a residual water concentration. The projections
assume a dry inlet and the frozen 2:1 reaction stoichiometry.
[Lorenzo et al. (2025), Sections 2.1–2.3](https://doi.org/10.1021/acs.iecr.5c02780)

For one mole of acid fed and molar conversion `X`, the outlet contains
`(1-X, 2-X, X, X)` moles of acid, pentanol, water and ester. For a sample that
retains `n_A` moles of that acid and `V_sample` mL on the same feed basis, the
concentration readout is `1 - n_A * V_feed / V_sample`. This equation keeps
liquid-volume and acid-partition effects distinct. The frozen reaction
calculation uses the earlier Eq. S1 reproduction with 6.0 g per void-liquid
volume, Eq. 4 and the reported Table B1 internal effectiveness.
[Prior full-inventory calculation](../../pentyl_propanoate_transfer/published/README.md)

## Reference-temperature volume and phase projections

The 25 °C volumetric calculation uses the [measured pure-liquid densities](properties.json).
Sum pure molar volumes for both inlet and outlet. Adding two
separate binary excess-volume terms for acid–pentanol and acid–water gives a
sensitivity to measured nonideality. The symmetric Redlich–Kister expression is
`x_A(1-x_A) sum(A_i (2x_A-1)^i)` evaluated on each binary subcomposition and
multiplied by the corresponding pair's mole total. This uses published
correlations for the two binaries. Their sum is a projection, with the
interactions involving the ester and the higher-order mixture terms unresolved.
The pure-volume estimate gives a reaction volume increment of +0.551 mL/mol of
ester at 25 °C. On the original 20 predictions, full-liquid volume ratios
change by at most 0.175% with additive volume and 0.379% with both binary terms.
[Organic densities in 2026 SI Table S1](https://api.figshare.com/v2/file/download/67271163) · [Acid–pentanol data](https://doi.org/10.33224/rrch.2020.65.3.06) · [Acid–water data](https://doi.org/10.1021/je030102f)

The dry-sample sensitivity excludes the volume of every mole of generated
water and retains every other component, including all remaining acid. Within
this additive-volume scenario, it is the largest downward readout shift from
water omission alone. The additive inlet volume is about 292.38 mL per mole of
acid, while the printed 3.3 M implies 303.03 mL per mole. The combined sensitivity
uses the latter feed basis and removes all produced water. It expresses the
effect of combining two different volume bases. The primary scenarios use a
consistent feed and product volume basis.

Amamehi's measured acid–pentanol–water tie lines at 30 and 50 °C show strongly
different water and acid concentrations in the two liquid phases. The 2026
study's Figure 9 also identifies limited miscibility in the water–ester binary
at 25 and 50 °C. The ternary sensitivity projects each outlet onto its
non-ester composition, interpolates adjacent measured tie lines by the lever
rule and locates all ester in the pentanol-rich phase. The tabulated mole
fractions are normalized for rounding. This projection omits the ester from
the equilibrium calculation and uses additive 25 °C molar volumes for both
phase readouts. The two source temperatures are evaluated separately, without
temperature or composition extrapolation. The equilibrium calculation for
the full four-component reaction mixture remains open. At 50 °C, the 40 °C,
3 mL/min model composition falls outside the measured tie-line envelope and
the predicted phase-specific readouts are unavailable.
[Measured ternary tie lines](https://doi.org/10.12962/j2964710X.v5i1.19035) · [Water–ester binary, 2026 Figure 9](https://doi.org/10.1021/acs.iecr.6c02375)

## Comparison after freezing the predictions

| Readout | Cases | RMSE (percentage points) | Largest downward shift (percentage points) | Cases high |
| --- | ---: | ---: | ---: | ---: |
| Frozen molar conversion | 20 | 10.46 | 0.00 | 20 |
| Whole liquid, additive | 20 | 10.49 | 0.00 | 20 |
| Whole liquid, binary excess projection | 20 | 10.47 | 0.06 | 20 |
| Exclude all generated water | 20 | 9.48 | 1.54 | 20 |
| Exclude water and use printed inlet basis | 20 | 8.64 | 3.44 | 20 |
| Pentanol-rich phase at 30 °C | 20 | 10.38 | 0.26 | 20 |
| Pentanol-rich phase at 50 °C | 19 | 10.16 | 0.20 | 19 |
| Water-rich phase at 30 °C | 20 | 19.36 | 2.51 | 20 |
| Water-rich phase at 50 °C | 19 | 19.99 | 0.00 | 19 |

At 60 °C and 3 mL/min the frozen prediction is 60.23% versus 49.70% measured.
The additive whole-sample readout is 60.28%. Complete omission of produced
water gives 58.74%, or 57.24% in the combined inlet-basis sensitivity. Retaining
all unreacted acid would require a sample volume equal to 79.06% of the inlet
volume on the fixed feed basis to match the observation. Across the 20 points,
the required volume ratio ranges from 39.46% to 81.70%. These inversions report
the size of the discrepancy after the comparisons and are not used by the
forward solver. Selective water-rich sampling generally raises the predicted
conversion and increases the RMSE.
[Readout comparisons](reference/comparison.csv) · [Run evidence](reference/summary.json)

**Next step.** The large residual after the tested liquid readout shifts
supports a separate check of intrabead transport at the full 6 g loading.
Phase effects on the reaction rate require independent activity or transport
evidence in the operating liquid.

## Reproduce

Run from the previously configured Python and Julia environment.

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transport.phase.run
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_transport/phase/test_phase.py
```

The runner writes the forward estimates to `predictions_before_comparison.csv`
before reading the 20 reactor observations for scoring. A phase-specific output
is blank when the ternary composition lacks an interpolable tie line.
