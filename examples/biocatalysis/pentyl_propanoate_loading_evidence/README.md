# Independent high-loading evidence and transfer audit

[Engineering results](REPORT.md) · [Analysis settings](protocol.json)

This extension preserves the prior catalyst-mass audit. Python manages
the independently extracted progress series and holdouts; Julia fits
the shared rate shape and separate apparent kinetics/thermodynamics,
and uses the tested existing fixed-inventory reactor and film models.
The private PDFs remain in the literature inbox. The compact source
data, extraction locations and reading assumptions are committed in
the corresponding paper directories, so reproduction does not need
the private source copies. All reactor observations remain outside
calibration. Different types of loading denominators are explicit.

From the repository root after [setup](../../../docs/guides/setup.md):

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_loading_evidence.run
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_loading_evidence/test_evidence.py
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_loading_evidence/test_shared_loading.jl
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_loading_evidence.check
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_loading_evidence.plot \
  --output results/biocatalysis/pentyl_propanoate_loading_evidence/loading_evidence.png
```

The runner creates a fresh timestamped run directory under
`results/biocatalysis/pentyl_propanoate_loading_evidence/`, or a
fresh location with `--output`. The matched progress CSV and external
time predictions use minutes, fraction conversion unless identified
as percent, and g/100 g total reactants unless the source definition
is different. A blank latest time indicates the sensitivity
did not attain the target within the observed time span. A blank
total-substrate loading means the denominator is unreported.
The JSON summary records every source, fit, predictive split,
factor and reactor score. All numerical JSON/CSV outputs including
the complete point-by-point predictions are hashed; only compact
records are retained here. The frozen calibration and predictions
are written before the measured reactor conversion file is opened.

Tested with Python 3.12.9 and Julia 1.11.6 on macOS arm64. Every
one of the 54 newly transcribed pentyl batch rows was checked against
the supplement cells. Two final full runs reproduce byte-identical
numerical CSV and JSON outputs. The four new Python checks and six
prior loading checks, plus nine repository layout/integration
checks, pass. All nineteen new Julia assertions (including an
independently generated reversible progress law with a hidden
high-loading curve, cross-study progress-time laws, and the
near-irreversible limit) and the thirty-five preceding Julia
loading/film checks pass. The
shared profile's 25/49-point refinement changes conversion by
less than 6e-10. Selected near-irreversible equilibrium profiles
differ by less than 1e-6 on training when the Kx limit increases
from 1e6 to 1e9. The maximum axial refinement change is 1.17e-5,
normalized acid balance error <4.7e-10 and local film balance
error <5.2e-13 mol/L/min. The previous proportional reference
replays within 6e-16. The numerical audit is complete; the
scientific transfer remains unqualified.
