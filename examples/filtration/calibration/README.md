# Experimental support calibration

## Crossflow source audit

The four new crossflow papers have been reviewed. [The received-paper report](../../../literature/crossflow_papers_review.md)
and [extraction](../../../literature/crossflow_extraction.json) retain six García
table rows, six Chiu (2005) and four Chiu (2006) figure readings. To reproduce
the source arithmetic and reading audit:

```sh
source scripts/env.sh
uv run python -m examples.filtration.calibration.audit_crossflow_sources --require-sources \
  --output results/crossflow_evidence/NEW_AUDIT/audit.json
```

Use an unused output path. The audit needs only standard-library Python;
it does not fit a model or call Julia, OpenLB or DWSIM. Source energy arithmetic
agrees within 0.912%, while the absolute seven-channel flow convention remains
unresolved. Pixel-reading bounds, reported repeatability and operating ranges
are separate. Existing calibration datasets and simulation results are preserved.

## Matched support and whole-membrane calibration

The matched Zhu support/whole-membrane study is also implemented. Its
[PDF-refined result and reading limits](../../../docs/zhu_pdf_refinement.md) are separate
from the earlier Biesheuvel fit below. To reproduce it:

```sh
source scripts/env.sh
uv run python -m examples.filtration.calibration.zhu --source pdf --require-sources
uv run python -m examples.filtration.calibration.zhu --check results/zhu_calibration/runs/RUN_ID
```

Replace `RUN_ID` with the printed directory. A new `--output` directory can be
specified; completed directories cannot be overwritten. Without
`--require-sources`, the retained digitized source data can be replayed when
local full text/figures are absent; any present source must match its hash.
The fit and reading-sensitivity analysis call the existing explicit Julia runtime.
`--source pdf` is now the default. Use `--source jpeg` for the preserved earlier
reading; `--check` replays the reading saved inside that run.

To audit the Figure 6 reading from the original JPEG:

```sh
uv run python -m examples.filtration.calibration.digitize_zhu --output results/zhu_digitization.json
uv run python -m examples.filtration.calibration.digitize_zhu --source pdf --output results/zhu_pdf_digitization.json
```

Each Zhu run has `inputs.json`, strict `summary.json`, `predictions.csv`,
`heldout_validation.png` and `conditional_energy.png`. Permeance is m/(Pa s),
effective impedance is Pa s/m, and Julia flux inputs/outputs are m/s. Original
figure x-values in m/h and y-values in LMH are explicitly named in the CSV;
one-based point IDs follow increasing support flux. Odd IDs train and even
IDs are held out. Reading bounds are not experimental uncertainty. The
effective fixed impedance includes possible non-support interface/access
effects and does not imply a uniquely measured selective-layer property.

The PDF reading passes both frozen nominal and conservative reading-corner
screens (4.56% and 7.88% maximum error). The earlier JPEG reading passes only
the nominal screen; its failed corner screen and 11.04% lowest-flux training
residual remain preserved. Both readings retain the low-flux model discrepancy.
The checker reports numerical integrity separately from those scientific
outcomes. No viscosity, intrinsic support k or real connected geometry is
invented. No new native OpenLB/DWSIM result or UI is claimed by this CLI.
Unavailable experimental data are not prerequisites for this computational
project; see [verification scope](../../../docs/verification_scope.md).

## Earlier Biesheuvel support closure

This separate CLI uses Biesheuvel and Verweij's measured alumina support
properties. The numerical closure, log-space fit and error metrics run in
the existing Julia package. Python selects the recorded rows, coordinates
the bridge and exports plots/data. No dependency or bridge environment was
added. The original OpenLB cases, DWSIM cases and Panel app are preserved.

The [result report](../../../docs/experimental_calibration.md) distinguishes execution
success from the **failed primary prediction screen**. This is a retrospective
same-study test of a reduced closure, not experimental validation of an OpenLB
pore geometry or of the complete membrane/process.

## Run and replay

From the project root:

```sh
source scripts/env.sh
uv sync --locked
uv run python -m examples.filtration.calibration.run --require-sources
uv run python -m examples.filtration.calibration.check results/calibration/runs/RUN_ID
```

Replace `RUN_ID` with the directory printed by the run. Optional
`--output results/calibration/my-run` selects a new directory; an existing
directory is rejected. Removing `--require-sources` permits reproduction from
the tracked, visually audited extraction if the local article files are
absent. That absence is explicitly recorded in `inputs.json`. Present source
files must always match their recorded hashes. No source is downloaded during
a simulation. The public author-thesis PDF URL and SHA-256 are retained in
[the source audit](../../../literature/biesheuvel_source_audit.json).

`calibration.check` verifies artifact hashes, replays the saved input through
Julia, compares all fitted/predicted/error/process values and checks the CSV.
Add `--require-screen-pass` to require the separate scientific outcome; it
deliberately fails for the current result. A normal run/check completing does
not imply that the prediction screen passed. Hashes establish artifact
integrity, not authenticity against deliberate rewriting of all records.

Tests, using the same Julia executable selected by the bridge:

```sh
/Users/dula/.julia/juliaup/julia-1.11.6+0.aarch64.apple.darwin14/bin/julia --startup-file=no --project=julia -e 'using Pkg; Pkg.test()'
uv run pytest -q
```

The full Python suite includes the original macOS VTK graphics checks and
requires normal host graphics access on this machine. For a calibration-only
check: `uv run pytest -q tests/test_calibration.py`.

## Fixed protocol and data boundary

[protocol.json](protocol.json) specifies the split and 25% engineering screen.
It was saved before this implementation's first fit, after the table had
already been reviewed. No blind, prospective or independent-laboratory
validation is claimed. There are three fitting points and three primary
held-out points, not replicated measurements with known noise.

The model is `k = b*d_primary^2*epsilon^3/[180*(1-epsilon)^2]`.
Only the dimensionless multiplier `b` is fitted. Minimising squared residuals
in natural-log permeability gives
`log(b) = mean(log(k_measured)-log(k_CK))`.
Equal log weights reflect an explicit modelling choice; uncertainty weights
are unavailable. Leave-one-training-case-out sensitivity is not a confidence
interval. The constant 180 and nominal particle sizes are fixed, not tuned.

| Quantity | Contract |
| --- | --- |
| Table row ID | Original one-based Table 2 row; not a Python/VTK node ID |
| Particle diameter | Approximate primary powder mean in metres; not a pore diameter or capillary radius |
| Porosity | Dimensionless volume fraction, strictly between 0 and 1; caller must assess percolation/applicability |
| Permeability | Darcy k in m²; measured, baseline CK and calibrated prediction kept separately |
| Julia vector inputs | Equal, nonzero-length `Vector{Float64}` arrays of diameter, porosity and k, one entry per selected case |
| Relative error | `k_pred/k_measured - 1`; log RMSE uses natural logs |
| Uncertainty | Unavailable measurement uncertainty is JSON null and empty CSV; no invented error bars |
| Process boundary | Existing Julia layers-in-series model, SI inputs; assumed selective resistance fixed within every design pair and across measured/predicted comparisons |

Julia functions `support_permeability`, `fit_support_multiplier` and
`permeability_errors` are usable directly with `using Pore2ProcessDemo`.
The Carman–Kozeny closure is an empirical-reference model alongside the
geometry-resolved OpenLB stage; it cannot replace that required stage.

## Outputs

Each immutable run directory contains:

- `inputs.json`: protocol, source/extraction/audit snapshot, source-file
  availability and hashes, model/bridge/lockfile hashes.
- `summary.json`: measured and predicted values, fit and held-out diagnostics,
  screening outcome, sensitivity, conditional process results, runtime versions,
  timing, timestamp and artifact hashes. JSON is strict; no NaN/Infinity.
- `predictions.csv`: nine rows with roles, SI quantities and signed errors.
- `permeability_validation.png`: all nine points and held-out errors; every
  point is shown, including out-of-domain cases.
- `conditional_energy.png`: measured-support versus fitted-support energy
  savings with the **same assumed selective resistance**.

The conditional energy comparison uses rows 3 → 2, both training rows. It
illustrates propagation of model discrepancy; it is not an independent process
validation. Nominal 2 mm supports, 100 LMH, 1 m², 1 mPa·s and 75% shaft-to-hydraulic
efficiency are explicit assumptions. Some high-selective-resistance scenarios
exceed the reported 2 MPa apparatus pressurization; JSON flags this. That
apparatus value is neither a material pressure rating nor a validated Darcy
range. No mechanical feasibility or full-scale energy performance is inferred.
