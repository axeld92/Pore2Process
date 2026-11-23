# Catalyst mass-proportionality audit

[Results and engineering conclusion](REPORT.md) · [Protocol](protocol.md) ·
[Analysis settings](protocol.json)

Python handles the published observations, mass-normalized progress
comparison, independent holdout splits and outputs. Julia uses the
existing reversible batch and film models and tests concentration-dependent
mass-specific rate laws. All 20 continuous-reactor observations remain
excluded from parameter estimation. The actual 6 g of catalyst is retained.

From the repository root after [setup](../../../docs/guides/setup.md):

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_loading.run
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_loading/test_loading.py
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_loading/test_loading.jl
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_loading.check
```

The runner makes a new timestamped directory under
`results/biocatalysis/pentyl_propanoate_loading/`, unless
`--output` specifies a fresh destination. It saves direct
matched-conversion observations, all loading-law fits and
whole-loading predictions, their scores, the frozen
fixed-inventory reactor predictions and the resulting comparison.
Units appear in output field names; acid conversion is a fraction
and concentration uses g/L (numerically kg/m³). A blank held-out
loading identifies an all-data fit; a blank shape-free prediction
is outside the training exposure range. The long point-by-point
`batch_predictions.csv` is generated with each run and its hash
recorded; compact summaries and predictions are retained here.

Regenerate the report figures using the retained summary:

```sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_loading.plot \
  --output results/biocatalysis/pentyl_propanoate_loading/curve_collapse.png \
  --summary examples/biocatalysis/pentyl_propanoate_loading/reference/summary.json \
  --extrapolation-output results/biocatalysis/pentyl_propanoate_loading/loading_extrapolation.png
```

Verified on macOS arm64, Python 3.12.9 and Julia 1.11.6.
Six Python and 35 Julia checks passed, including known
parameter recovery with a withheld synthetic loading and
an independent integration of the reversible rate model.
Doubling the loading-profile search resolution changes the full
primary predictions by at most 2.75e-9 in conversion.
Refining the axial reacting grids changes predictions by at
most 1.17e-5. The largest normalized acid balance error is
4.71e-10, maximum film balance error 6.65e-13 mol/L/min,
and the preceding central-film result replays within 6e-16.
Both final runs reproduced byte-identical numerical JSON and CSV
outputs. The nine repository layout tests and 31 preceding Julia
batch/film tests also passed.
The numerical screen is complete; the high-concentration rate
extrapolations remain unqualified. No new OpenLB simulation
or DWSIM process validation is claimed in this stage.
