# Reproduce the particle-kinetics and equilibrium checks

[Results and source audit](REPORT.md) · [Kinetics protocol](protocol.md) · [Freundlich protocol](freundlich_protocol.md)

Run from the project root with the existing Python–Julia environment,
the Alhamed PDF and the saved N60 data from the
[25 cm packing-film comparison](../../porosity_matched_column/README.md).
The column stage checks the archived input and transport-file paths.
Each batch stage prints its new run directory. Pass that directory to
its matching transfer command.

```sh
source scripts/env.sh

uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.run
BATCH_RUN=results/adsorption/breakthrough_hypotheses/batch/runs/20261117T234340Z-2095512d
uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.transfer "$BATCH_RUN"

uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.freundlich
FREUNDLICH_BATCH=results/adsorption/breakthrough_hypotheses/batch/freundlich_batch/20261117T235228Z-c703c8e5
uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.freundlich_transfer "$FREUNDLICH_BATCH"
```

Replay the tracked checks against the local runs that produced the report.
The result directories contain the original input plans, numerical
refinements, trajectory digests and scores.

```sh
source scripts/env.sh
BATCH_RUN=results/adsorption/breakthrough_hypotheses/batch/runs/20261117T234340Z-2095512d
TRANSFER_RUN=results/adsorption/breakthrough_hypotheses/batch/transfer/20261117T234456Z-82ae7115
FREUNDLICH_BATCH=results/adsorption/breakthrough_hypotheses/batch/freundlich_batch/20261117T235228Z-c703c8e5

uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.run --check "$BATCH_RUN"
uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.transfer "$TRANSFER_RUN" --check
uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.freundlich --check "$FREUNDLICH_BATCH"
FREUNDLICH_TRANSFER=results/adsorption/breakthrough_hypotheses/batch/freundlich_transfer/20261117T235324Z-835f21b3
uv run --locked python -m examples.adsorption.breakthrough_hypotheses.batch.freundlich_transfer "$FREUNDLICH_TRANSFER" --check
```
