# Section-specific exchange and transit times

[Protocol](protocol.md) · [Results](REPORT.md) · [Compact record](reference_summary.json)

Run from the checkout root with the saved [N160 cohort](../disordered_n160/REPORT.md)
available locally. The fresh probes reuse the Julia 1.11.6 executable and the
saved balanced TRT face flows. Tracking 16,384 markers for each of the four
cases is the main cost. The model and scoring replay are much faster.

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/adsorption/disordered_section_exchange/test_model.jl
uv run python -m examples.adsorption.disordered_section_exchange.run
uv run python -m examples.adsorption.disordered_section_exchange.run --check "$section_exchange_summary"
uv run python -m examples.adsorption.disordered_section_exchange.plot
```

Set `section_exchange_summary` to the printed summary path. The read-only
check verifies the saved flow, tracking code, local duration extraction,
predictions, archived arrivals and scoring. `julia` in the test command
should resolve to the same executable recorded in the probe plan.

The Julia solver accepts two row-major little-endian Float64 files.
`training.f64` contains cumulative first-passage times at sections 1 and 2.
`probe_durations.f64` contains local transit times through sections 2–8.
All times are in seconds. The resulting six-row prediction CSVs give tracer
recovery, mean arrival time in seconds, and variance in seconds squared.
The per-section rule tables have donor rows and receiver columns.
`counts` and `transition` give probe counts and transition probabilities.
`delay` and `delay2` store the transit-time moments multiplied by the
transition probability, in seconds and seconds squared respectively.
