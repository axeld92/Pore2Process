# Phenol breakthrough hypothesis tests

The study varies wall flow, carbon uptake histories and batch-qualified
particle kinetics. It preserves the total dry carbon calculated from
Alhamed's reported bulk density and holds the published breakthrough curves
out of parameter selection. See the [protocol](protocol.md), [results](REPORT.md)
and [compact run record](reference_summary.json).

Run from the repository root with the configured Julia environment, the
Alhamed source in the literature inbox and the archived
[N60 porosity-matched film](../porosity_matched_film/README.md) and
[25 cm column](../porosity_matched_column/README.md) runs. The runners check
the upstream inputs against their saved digests. The archived trajectories
and new bulk results are kept under `results/`.

```sh
source scripts/env.sh
uv sync --locked
uv run python -m examples.adsorption.breakthrough_hypotheses.screen
uv run python -m examples.adsorption.breakthrough_hypotheses.run
```

The first command after setup tests parallel annulus/core flow without exchange.
The second tests the periodic baseline and split carbon histories using
the saved N60 delay and contact rules. It records refinement results and
hashes the predicted curves before reading Figures 4 and 5. Pass the printed
run directory with `--check` to replay the score without solving.

The [hydraulic benchmark](benchmark/README.md) audits published annular water
flux and shows how the no-slip wall alters estimates of wall flow. The
[confined stock-TRT study](hydrodynamics/README.md) supplies water flow,
void and solid fractions and a forward-tracer radial transfer matrix. The
column runner uses a half-particle-diameter annulus. Its JSON input has
`cases.q233` and `cases.q625` objects. Each requires these fields.

| Field | Meaning |
| --- | --- |
| `flow` | Fraction of total signed water flow in the wall annulus |
| `liquid` | Fraction of total void volume in the annulus |
| `carbon` | Fraction of total solid volume in the annulus, at uniform particle density |
| `length_m` | Tracer transfer length in metres |
| `transition` | Forward radial tracer probabilities, with wall/core inlet rows and wall/core outlet columns |

Pass the saved handoff JSON to evaluate confined flow with shared and
separate carbon histories. Set `wall_base` to the directory printed
by the [confined-flow commands](hydrodynamics/README.md#commands).

```sh
uv run python -m examples.adsorption.breakthrough_hypotheses.run --wall "$wall_base/column_handoff.json"
```

The [particle kinetics study](batch/README.md) freezes the stirred-bath
alternatives and transfers those supported by withheld particle sizes.
Each column comparison uses the full published measurement window for
that condition. The late uptake metric integrates from the **published**
50% breakthrough time to the last published marker.

The [conditional equilibrium test](batch/conditional/README.md) selects
one equilibrium weight on the two shortest beds and freezes it for three
other conditions. The [joint protocol](coupled_protocol.md) reuses it
with the confined handoff and separate radial carbon histories.
Set `calibration_base` to the printed conditional run directory.

```sh
uv run python -m examples.adsorption.breakthrough_hypotheses.run \
  --wall "$wall_base/column_handoff.json" --calibration "$calibration_base"
```

Use the `runs` entries in the compact record to locate the exact saved
outputs. Set `forecast_summary` to a periodic, confined or joint column
`summary.json` to replay its scores. The plot uses the saved runs named
in the compact record.

```sh
uv run python -m examples.adsorption.breakthrough_hypotheses.run --check "$(dirname "$forecast_summary")"
uv run python -m examples.adsorption.breakthrough_hypotheses.plot
uv run pytest -q tests/adsorption/test_breakthrough_hypotheses.py
```
