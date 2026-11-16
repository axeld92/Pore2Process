# Repeat-join sensitivity at fixed carbon inventory

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the checkout root. The [fixed-carbon N160 column archive](../disordered_memory_column/REPORT.md)
provides the prior curves and uses the held-out tracer rules and the
batch-qualified particle model. The current run reuses its Julia solver
and all archived inputs. The archived mixed join is the numerical control.

```sh
source scripts/env.sh
julia --startup-file=no --project=julia examples/adsorption/disordered_repeat_join/test_model.jl
uv run python -m examples.adsorption.disordered_repeat_join.run
uv run python -m examples.adsorption.disordered_repeat_join.run --check "$join_summary"
uv run python -m examples.adsorption.disordered_repeat_join.plot "$join_summary"
```

Set `join_summary` to the path printed by the simulation. The read-only
check verifies source and curve hashes, phenol conservation, refinement
and Figure 7 scores. The full run evaluates the unreactive and reactive
column with both joins. The original delay-network solver remains the
implementation for both cases.
