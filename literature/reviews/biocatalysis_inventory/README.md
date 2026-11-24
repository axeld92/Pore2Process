# Catalyst inventory and activity audit

[Engineering report](REPORT.md) · [Protocol](protocol.toml) ·
[Complete compact result ledger](reference.toml).

Four supplied articles and one supplement provide independent evidence
for loading, wet filling, activity history, and microreactor scale-up.
The private originals stay ignored in `literature/inbox/`; each
source-specific `metadata.json` and `evidence.toml` retain its hash,
experimental conditions, units, readings and uncertainty limitations.
Do not merge their different chemistries or copy any fitted activity
constants into the pentyl-propanoate input files.

From the repository root, using the established project environments:

```sh
source scripts/env.sh
julia --startup-file=no --project=julia literature/reviews/biocatalysis_inventory/audit.jl /tmp/inventory_audit.toml
julia --startup-file=no --project=julia literature/reviews/biocatalysis_inventory/test_audit.jl
cmp /tmp/inventory_audit.toml literature/reviews/biocatalysis_inventory/reference.toml
uv run --locked --offline python literature/reviews/biocatalysis_inventory/plot.py \
  --summary /tmp/inventory_audit.toml --output /tmp/inventory_audit.png
```

Julia performs bounded interpolation, volume/unit calculations and
profiled least-squares activity fits; Python only renders the independent
results. The activity fit uses the same declared 0–10 min markers
for both candidate forms and withholds 30 min, with the earlier 0–4 min
logarithmic slope independently checked against the paper convention.
The 600/1200 rate-grid refinement changes holdout predictions by less than
1e-5 percentage points. Reading sensitivity perturbs each training marker
by the declared allowance, not the unreported experimental variance.
The Damstrup and Strniša comparisons use one source curve each, keeping
the other loadings or geometries withheld; no pentyl-propanoate
conversion point is read by this audit.

Tested macOS arm64, Julia 1.11.6, CPython 3.12 project environment.
The 64 Julia tests cover axis/unit checks, independently generated
analytical activity predictions, and source/holdout consistency.
Two final numerical runs reproduce byte-identical TOML. Visual
source-page and supplement inspection verifies the tables and
marker conventions; the figure is rebuilt from the saved readings.
The nine existing Python repository layout/integration checks pass
(`tests/integration/test_layout.py` and
`tests/integration/test_repository_layout.py`).
