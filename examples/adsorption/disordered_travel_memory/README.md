# Within-path travel-time memory

[Protocol](protocol.md) · [Results](REPORT.md) · [Compact record](reference_summary.json)

Run from the checkout root with the [N160 reference trajectories](../disordered_n160/REPORT.md)
and [independent operator probes](../disordered_section_exchange/REPORT.md)
available locally. The runner uses the Julia 1.11.6 executable recorded
in the N160 tracer plan.

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/adsorption/disordered_travel_memory/test_model.jl
uv run python -m examples.adsorption.disordered_travel_memory.run
uv run python -m examples.adsorption.disordered_travel_memory.run --check "$memory_summary"
uv run python -m examples.adsorption.disordered_travel_memory.plot
```

Set `memory_summary` to the printed summary path. The check verifies
hashes, repeats all four Julia calculations in a temporary directory and
checks their products byte for byte.

The Julia entry point is `solve.jl training.f64 probe_passages.f64 output_directory`.
The input files contain row-major little-endian Float64 first-passage
times in seconds. The training input includes archived sections 1–2.
The probe input includes independent sections 1–8. Each prediction CSV
contains recovery, mean time in seconds and variance in seconds squared
at sections 3–8. `rules.f64` stores 18 column-major 40 × 40 tables. For
each of six seams, it stores the transition matrix and the weighted
first and second moments of the local delay, in seconds and seconds
squared. Rows are donors and columns are receivers, ordered by the ten
local transit-time paths and then the four elapsed-time ranks.
