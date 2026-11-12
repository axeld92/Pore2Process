# Predicting Alhamed Figure 7 with path exchange and particle uptake

The [protocol](protocol.md) specifies the archived tracer closure,
the column transfer assumptions, the conserved inventories and the
scores. The liquid balances and particle uptake are solved in
[Julia](../../../julia/adsorption/ExchangeColumn.jl). The coordinator
writes both predictions and hashes before it scores Figure 7.
It reads the qualified tracer and particle archives and uses the
same local project Julia runtime as the earlier adsorption cases.
[Results and verification](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=julia julia/adsorption/test_exchange_column.jl
uv run python -m examples.adsorption.exchange_column.run
uv run python -m examples.adsorption.exchange_column.run --check <summary-path>
uv run python -m examples.adsorption.exchange_column.plot --summary <summary-path>
```

Replace `<summary-path>` with the printed completed path. The run saves
the path classes, stationary flows, source assumptions, physical units,
outputs and scores under `results/adsorption/exchange_column/`.
A marker transition count is indexed by upstream class and receiving
class, in order of increasing time rank. The Julia mixing matrix
maps upstream concentrations to receiving concentrations weighted
by their contributions to the receiving flow. The column history
uses seconds, kg/m³ dissolved phenol, and kg solute inventory or
cumulative transfer. Read-only replay reconstructs the tracer
closure and recomputes the experimental comparisons.
