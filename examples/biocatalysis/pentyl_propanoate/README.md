# Pentyl propanoate with immobilised lipase

[Scientific report](REPORT.md) · [Comparison protocol](protocol.json)

This example fits the 2:1 batch experiments of
[Lorenzo and Ortega (2026)](../../../literature/sources/lorenzo2026/README.md)
and checks the conversion prediction on the independent
[Lorenzo et al. (2025) packed-bed series](../../../literature/sources/lorenzo2025/README.md).
The water-tracer experiment provides a separate check of bed mixing.

## Reproduce the Julia reaction calculation

From the project root with the [configured environment](../../../docs/guides/setup.md),
run the following commands. The full reaction run takes the retained
[OpenLB water-flow result](mixing/reference.json) as an input.

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate/test_model.jl
uv run --locked python -m examples.biocatalysis.pentyl_propanoate.run --batch-only
uv run --locked python -m examples.biocatalysis.pentyl_propanoate.run
uv run --locked python -m examples.biocatalysis.pentyl_propanoate.check
```

A new run is written under `results/biocatalysis/pentyl_propanoate/runs/`.
`batch_fit.csv` preserves the observed series and the two batch-model fits.
`bed_predictions_before_comparison.csv` is written before the converter opens
the packed-bed observations. `bed_comparison.csv` and the corresponding PNG
compare that frozen prediction to the twenty measured cases. Fractions in
the files are dimensionless and values in the plot are percentages. The
summary includes source hashes, catalyst inventory through the protocol,
software version, error metrics in percentage points and balance diagnostics.
The `check` command replays the curated comparison without a native run.

## Reproduce the native water-flow check

Use the [OpenLB commands and conditions](mixing/README.md) to build and run
the independent tracer check. To predict using a fresh tracer summary,
pass its project-local path via `--mixing-summary`.
