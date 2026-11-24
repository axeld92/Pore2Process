# Bajić plumbing-volume and dispersion audit

See the [scientific report](REPORT.md), [source evidence](../../../literature/sources/bajic2017/README.md)
and [frozen protocol](protocol.toml). The preceding Strniša TRT and Panjan
tracer examples and their numerical/measurement records are unchanged.
This audit runs an independent Julia calculation of the transit-time
moments in developed straight tubing, using measured/related-source tube
diameters and the previously reported nominal nonreactive fluid properties.
It does not run a new porous OpenLB geometry or start reaction.

From the configured project root:

```sh
source scripts/env.sh
julia --threads=1 --startup-file=no --project=julia examples/biocatalysis/bajic_plumbing_audit/test_tube_moments.jl
uv run --locked --offline python -m examples.biocatalysis.bajic_plumbing_audit.audit
uv run --locked --offline python -m examples.biocatalysis.bajic_plumbing_audit.check --require-run
uv run --locked --offline pytest -q examples/biocatalysis/bajic_plumbing_audit/test_audit.py
uv run --locked --offline python scripts/check_repository.py
```

The `audit` command writes a new compact reference and figures. Julia is
called by the existing project runtime selector; its executable and the
actual command are recorded. The original source PDFs are not required for
replaying the extracted volume and moment records. The check verifies the
hashes of the extracted sources, numerical inputs/code and outputs, and the
PDF hashes when the originals are locally available. `--require-run`
also requires the ignored full run directory
`results/biocatalysis/bajic_plumbing_audit/runs/<id>/`. Omit it for a
source-only consistency check on a fresh checkout.

The [reference summary](reference/summary.json) keeps source conventions,
the signed mean and variance residuals, all reading/diffusivity cases, the
actual radial convergence test, and the illustrative bore-volume bounds.
[Tube moments](reference/tube_moments.csv) use seconds for time, s² for
variance and m²/s for molecular diffusivity. Radii are SI in the Julia
solver; printed fixture lengths are mm and volumes µL. The inferred extra
length is a mean-only diagnostic, not a measured length or a fit of the
retained bed pressure/flow field.
