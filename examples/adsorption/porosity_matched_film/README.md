# Clean-wall phenol transfer at 42.7% void fraction

[Protocol](protocol.md) · [Report](REPORT.md) · [Compact record](reference_summary.json)

Run from the project root using the installed OpenLB 1.9r0 archive and
the [previous column and tracer archives](../disordered_film_conductance/README.md).
The native command prints the case directory. Set `dense_run` to that
directory for the subsequent commands. The copied rectangular OpenLB
application and its build live under `.local/adsorption/porosity_matched_film/`.

```sh
source scripts/env.sh
uv run python -m examples.adsorption.porosity_matched_film.native --grid 45 --pilot
uv run python -m examples.adsorption.porosity_matched_film.native --base "$dense_run" --grid 45 --warm-grid 45 --warm-trial
uv run python -m examples.adsorption.porosity_matched_film.native --base "$dense_run" --grid 60 --warm-grid 45 --warm-trial --threads 4
uv run python -m examples.adsorption.porosity_matched_film.handoff "$dense_run" --project 45
uv run python -m examples.adsorption.porosity_matched_film.handoff "$dense_run" --project 60
uv run python -m examples.adsorption.porosity_matched_film.handoff "$dense_run" --compact 45 60
uv run python -m examples.adsorption.porosity_matched_film.run "$dense_run" --trace 45 20 6144
uv run python -m examples.adsorption.porosity_matched_film.run "$dense_run" --trace 60 40 8190 --threads 3
uv run python -m examples.adsorption.porosity_matched_film.run "$dense_run" --trace 60 20 8190
uv run python -m examples.adsorption.porosity_matched_film.run "$dense_run" --seal 45 60
uv run python -m examples.adsorption.porosity_matched_film.run "$dense_run" --check
uv run pytest -q tests/adsorption/test_porosity_matched_film.py
```

The warm starts use the converged fields from `_trial` directories.
The pilot fields supply initial moments even if their Mach or density
screen misses. Every accepted native output still passes the full
validation. The native command accepts a separate `--flow q233` or
`--flow q625` when repeating one solve. Each case records the executable
command and initial-field digest.

The geometry contains sphere centers in metres as x, y, z columns with
periodic coordinates. The binary velocity and balanced faces use the
[original packing conventions](../disordered_packing/README.md), with
rectangular grid extents recorded in each case. Each trajectory row
contains eight successive first-passage times followed by eight
section-specific wall clocks, in seconds. NaN marks a section with no
wall contact. The replay verifies the full local archive and recalculates
the conductance and escape estimates.
