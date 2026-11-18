# Pentyl propanoate batch-to-bed transfer audit

[Results and interpretation](REPORT.md) · [Comparison protocol](protocol.md)

This follow-up uses the frozen [pentyl propanoate benchmark](../pentyl_propanoate/README.md)
and the curated [2025](../../../literature/sources/lorenzo2025/README.md)
and [2026](../../../literature/sources/lorenzo2026/README.md) source data.
The 20 packed-bed conversion measurements score the forward predictions.

Run from the project root with the existing [Python and Julia setup](../../../docs/guides/setup.md).

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transfer.mixing.run
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transfer.published.run
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transfer.run
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transfer.check
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_transfer/test_transfer.py
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_transfer/published/test_model.jl
```

Each runner writes a new run under
`results/biocatalysis/pentyl_propanoate_transfer/`. Each prediction table
is written before the packed-bed observations are read. The curated tables
and source hashes are linked from the report. The checker verifies the
retained bytes and inputs without a solver run. The direct-transfer CSV
uses empty cells for conditions outside the reported batch time range.

The [mixing methods](mixing/README.md) describe transfer of the water tracer
to the Julia axial-dispersion model. The [published model reconstruction](published/README.md)
describes the 2025 kinetic equation, units check and catalyst-volume bases.
