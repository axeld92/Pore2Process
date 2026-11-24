# Panjan: measurement-loop and middle-length tracer audit

Read the [scientific report](REPORT.md), [source extraction](../../../literature/sources/panjan2018/README.md),
and [protocol](protocol.toml). This is a follow-up to the numerically
qualified but experimentally failed [Strniša tracer example](../strnisa_tracer/REPORT.md).
It does **not** run a new R2 OpenLB geometry, fit R2, or start reaction.

From the existing configured project root:

```sh
source scripts/env.sh
uv run --locked --offline pytest -q examples/biocatalysis/panjan_tracer_audit/test_audit.py
uv run --locked --offline python -m examples.biocatalysis.panjan_tracer_audit.audit --derive-kernel
uv run --locked --offline python -m examples.biocatalysis.panjan_tracer_audit.check --require-run
```

`--derive-kernel` verifies the actual saved A_H16 marker file against the
preceding reference and bins its exit ages conservatively at 0.05 s in
`reference/A_impulse.csv`. It requires the previous full local run. The
resulting compact impulse and its manifest are tracked as numerical
evidence. On a checkout without the full field or markers, omit that
flag to replay the same archived impulse:

```sh
uv run --locked --offline python -m examples.biocatalysis.panjan_tracer_audit.audit
uv run --locked --offline python -m examples.biocatalysis.panjan_tracer_audit.check
uv run --locked --offline python scripts/check_repository.py
```

Running `audit` explicitly publishes a new compact reference and figures.
The actual saved original OpenLB/Julia run is unchanged. The ignored
`results/biocatalysis/panjan_tracer_audit/runs/<id>/` folder contains each
full predicted response before its scoring and a prediction manifest;
`reference/predicted_narrow_samples.csv` retains the nominal-bracket
predictions at 5 s cadence through the figure's 7 min range. The
[reference summary](reference/summary.json) retains all full-window
comparison scores and code/source/output hashes. A source-only check
does not pretend to rerun OpenLB; `--require-run` also verifies the full
forward-prediction archive if requested.

All dimensional pulse concentrations are mM, times seconds, reactor
liquid volumes µL, marker weights dimensionless exit-mass fractions,
and recovered apparent dose µmol. The original measured blank already
includes injection, tubing and detector effects; only the saved pore
age kernel and frozen diagnostic retention scenarios are convolved with
it. R2 has a fixed reported liquid inventory of 35 µL. Two independent
A-derived transit sections and a single scaled retained-memory case
are reduced-model sensitivity brackets, not geometry-specific TRT
predictions of R2. The wide R3/R4 plots are integrated only over their
visible nonzero-ended windows; no full RTD moments are supplied for them.
