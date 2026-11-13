# Disordered packing path exchange and wall-contact holdout

Use the saved [N160 cohort](../disordered_n160/REPORT.md). The
[protocol](protocol.md) specifies the split and the acceptance limits.
[Results](REPORT.md)

The [variance audit](VARIANCE_AUDIT.md) uses the original frozen prediction
to distinguish changes in path exchange and local transit times from
missing travel-time memory.

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=julia examples/adsorption/disordered_exchange/test_model.jl
uv run python -m examples.adsorption.disordered_exchange.run
uv run python -m examples.adsorption.disordered_exchange.run --check "$holdout_summary"
uv run python -m examples.adsorption.disordered_exchange.plot "$holdout_summary"
```

Set `holdout_summary` to the path printed by the run. The Python
runner extracts the two-section training files and verifies the
archived data. Julia estimates and propagates the model. From the
example directory the Julia entry point is
`julia --startup-file=no solve.jl training.f64 output_directory`.
The input file contains little-endian Float64 rows with four columns
in seconds: first and second cumulative passage times, followed by
the two absolute first-wall-hit clocks. NaN means no observed contact.
The output directory must exist. Transition and wall-free tables use
donor rows and receiver columns. The `pair.csv` and `receiver.csv`
rows cover sections 3–8 with columns for total survival, fast survival,
the two conditional capture fractions, inert mean arrival in seconds,
inert variance in seconds squared, and fast-path persistence.
