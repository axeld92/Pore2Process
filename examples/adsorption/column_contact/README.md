# Column contact and inventory audit

[Scientific report: purpose, theory, results and conclusions](REPORT.md).

This separate extension preserves the batch-qualified particle model and all
earlier examples. [The report](../../../docs/verification/adsorption/column_contact.md) explains the
source-density audit, common contact model, held-out comparison and limits.

From the existing macOS ARM64 project root:

```sh
source scripts/env.sh
uv sync --locked
julia --project=julia julia/adsorption/test_film_particle.jl
julia --project=julia julia/adsorption/test_column_contact.jl
uv run python -m examples.adsorption.column_contact.study
uv run python -m examples.adsorption.column_contact.check
uv run pytest -q
```

The study selects N using Figure 7's 0.8 mm curve only. It does not refit the
batch model, isotherm, film coefficient, active mass or bed density. Selection
is deliberately repeated on a fresh run; no optimization cache is hidden.
Unique output folders live under `results/adsorption/column_contact/runs/`.
The checker replays every saved case and calibration evaluation, verifies
hashes and SI CSV exports, and writes the compact documentation reference.

Direct Julia use (same `julia/` project):

```julia
include("julia/adsorption/FilmParticle.jl")
include("julia/adsorption/ColumnContact.jl")
using .ColumnContact
r = solve_contact(sections=1, diffusivity_m2_s=1.846035489738998e-12,
                  film_m_s=4.540819586217862e-5)
```

`sections` is a **physical mixing model choice**, not a spatial resolution.
Only time and radial discretization are numerical refinements. The particle
implementation is reused unchanged; the new wrapper supplies equal section
inventories and directed advective coupling to its conservative integrator.
It exposes no dispersion, bypass, inactive-solid or inventory-fit parameter.
