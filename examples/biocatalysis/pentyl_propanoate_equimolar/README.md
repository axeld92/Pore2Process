# Equimolar packed-bed holdout

[Engineering results](REPORT.md) · [Dual-feed protocol](protocol.json) ·
[Published-sign audit protocol](source_sign_protocol.json)

This separate example imports the previously checked batch parameters and
their equilibrium regressions. It never calls a fitting routine. Julia
uses the original axial-dispersion and liquid-film implementations with a
feed-aware rate wrapper; Python joins measured conversions **after** saving
both feeds' predictions. The actual weighed mass remains 6 g. Every case
has explicit ratio, acid concentration, residence, activity, film capacity
and internal effectiveness. No existing reference or OpenLB source changes.

From the project root after [setup](../../../docs/guides/setup.md):

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_equimolar.run --output results/biocatalysis/pentyl_propanoate_equimolar/new_run
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_equimolar.run --output results/biocatalysis/pentyl_propanoate_equimolar/new_replay
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_equimolar.source_sign_audit --previous results/biocatalysis/pentyl_propanoate_equimolar/new_run/bed_comparison.csv --output results/biocatalysis/pentyl_propanoate_equimolar/new_sign
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_equimolar.source_sign_audit --previous results/biocatalysis/pentyl_propanoate_equimolar/new_run/bed_comparison.csv --output results/biocatalysis/pentyl_propanoate_equimolar/new_sign_replay
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_equimolar/test_model.jl
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_equimolar/test_equimolar.py
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_equimolar.check
```

Use fresh output paths; the full run refuses to overwrite a prior result.
The main campaign took 149 s in one cold process and 141 s in another on
the tested Apple Silicon macOS machine (Python 3.12.9, Julia 1.11.6).
There is no runtime requirement for the private paper PDFs because the
digitized, cross-checked evidence and source hashes are retained.

Regenerate compact evidence in a fresh destination and redraw the plots:

```sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_equimolar.retain \
  --run results/biocatalysis/pentyl_propanoate_equimolar/new_run \
  --replay-run results/biocatalysis/pentyl_propanoate_equimolar/new_replay \
  --sign-run results/biocatalysis/pentyl_propanoate_equimolar/new_sign \
  --sign-replay results/biocatalysis/pentyl_propanoate_equimolar/new_sign_replay \
  --output results/biocatalysis/pentyl_propanoate_equimolar/new_compact
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_equimolar.plot \
  --source results/biocatalysis/pentyl_propanoate_equimolar/new_compact \
  --sign-source results/biocatalysis/pentyl_propanoate_equimolar/new_compact \
  --output results/biocatalysis/pentyl_propanoate_equimolar/new_figures
```

The root [reference](reference/summary.json) retains pointwise central,
age and transport comparisons, source and output hashes, the sign audit
and [cold replay record](reference/replay.json). Complete intermediate
ledgers and per-run details remain in the ignored `results/` output. Numeric
conversions are dimensionless fractions throughout the API and data files;
errors in `*_pp` fields are percentage points. Catalyst mass and dimensions
use kg and m, concentration mol/L or g catalyst/L liquid as labeled, flow
mL/min, time min and temperature °C externally. Absolute temperatures in K
enter the kinetic laws. The inferred outlet ester flow (g/h) is based on
constant total volume flow, a dry feed and the single reaction; it is not an
independent mass-flow observation.

### Verified status

- All 40 experimental cases withheld; 22 central scenarios, seven
  stationary catalyst-age hypotheses at each of two time-on-stream values,
  four transport sensitivities on five fixed model cases.
- 50 direct Julia assertions, four new Python tests and seven preceding
  kinetic tests passed. Independent fresh runs produced identical prediction
  ledgers and summaries, including the 80 source-sign bed and 450 batch-time
  data points. All previous 2:1 predictions replay exactly. The source files
  and compact result hashes pass the standalone check.
- Source sign ambiguity documented, not silently patched or fitted.
  No common validated kinetic prediction for 1:1 and 2:1; no new reacting
  OpenLB or DWSIM acceptance claim.
