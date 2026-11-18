# Liquid-film resistance and water activity

[Bottom line and results](REPORT.md) · [Protocol](protocol.md)

The film calculation adds an external liquid-film balance to the
[2025 reactor reconstruction](../pentyl_propanoate_transfer/published/README.md).
It uses the [independent liquid-property evidence](film/evidence.md) and all
6.0 g of catalyst. The [batch rate analysis](water_activity/evidence.md)
uses independently controlled initial water activities for Novozym 435.

Run from the project root after [environment setup](../../../docs/guides/setup.md).

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_film_water.film.run
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_film_water.water_activity.run
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_film_water/film/test_model.jl
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_film_water.check
```

The runners write to fresh directories under
`results/biocatalysis/pentyl_propanoate_film_water/`. Pass `--output` to
give a run a specific directory inside the checkout. The film runner
saves the predictions and liquid-film diagnostics before reading the
20 measured conversions. The batch runner reads the independent
initial-rate series. Both preserve their input hashes in the summary.

Verified on macOS arm64 with Python 3.12.9 and Julia 1.11.6. The
16 film-model checks passed. Two independent runs of each analysis
produced identical CSVs and JSON summaries.
