# Pentyl propanoate catalyst-conditioning audit

[Bottom line and results](REPORT.md) · [Published evidence](evidence.md) ·
[Protocol](protocol.md) · [Numeric inputs](inputs.json)

This independent screening step compares first-contact procedures,
published bead wetting and sorption, water-tracer exchange, reuse
and residual activity. Julia solves conservative transient bead
moisture storage and runs the previously tested film/reaction model
with independent apparent-rate sensitivities. Python keeps source
properties and scenarios separate and performs the held-out
full-6-g reactor comparison. No new OpenLB flow field, DWSIM case,
or fitted packed-bed parameters are introduced.

From the repository root after [environment setup](../../../docs/guides/setup.md):

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_conditioning.run
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_conditioning/test_conditioning.py
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_conditioning/test_moisture.jl
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_conditioning.check
```

Results are emitted to a new timestamped directory under
`results/biocatalysis/pentyl_propanoate_conditioning/`. Use
`--output` for a distinct run directory. The full archive contains
water-inventory tables, activity differences, time curves for
water transfer, the axial predictions saved before opening the
20-run comparison, and the subsequent scoring. Conversions,
activity and normalized remaining inventories are fractions.
Water storage concentration uses kg/m³ (numerically g/L).
Blank activity or equilibrium cells indicate an unavailable
calculation within the corresponding source's concentration
range. The compact [reference summary](reference/summary.json)
records hashes of the model/data files and the full run outputs.

Verified on macOS arm64 with Python 3.12.9 and Julia 1.11.6.
The five new Python checks, 19 Julia moisture checks, five earlier
water-activity checks, and 16 earlier film checks passed. The
moisture cases use 64/128 radial cells and 2/1 s time steps;
the bed cases use 512/1024 axial cells. The direct Julia diffusion
test compares against the exact spherical diffusion solution.
The central-film predictions reproduce the preceding audit, and two
complete final runs produced byte-identical summary and CSV files.
The repository catalogue and nine layout checks also passed.
This completes the conditioning *screen*; transfer of the actual
conditioning history and rate into the process case is unresolved.
