# Finite-bath and film/particle adsorption check

Separate date-stone-carbon model and source audit. It preserves the original
synthetic OpenLB–Julia–DWSIM example and the earlier failed Langmuir/LDF study.
The shared fit uses only the 0.8 and 1.5 mm Figure 2 curves; smaller batch
particles are withheld. Column diffusivity, source equilibrium and physical
inventory remain frozen. An optional single column film multiplier uses only
0.8 mm Figure 7; 0.45 mm remains withheld.

The batch comparison passes the chosen 3 mg/g screen. Both unfitted and
one-anchor column routes fail the experimental comparison. Do not treat a
successful numerical replay as column validation.

From the project root:

```sh
source scripts/env.sh
uv sync --locked
uv run python -m adsorption.film_particle.data
julia --project=julia julia/adsorption/test_film_particle.jl
uv run python -m adsorption.film_particle.run
uv run python -m adsorption.film_particle.check
uv run pytest -q
```

No new dependencies or environment changes are required. The existing bridge
sets the selected Julia executable/project before JuliaCall imports. Calls
stay on the main thread. The numerical package is usable directly from Julia:

```julia
include("julia/adsorption/FilmParticle.jl")
using .FilmParticle
result = solve_bath()
endpoint = bath_equilibrium()
```

Results are saved in unique `results/adsorption/film_particle/runs/` folders,
with a separate `latest.json`. Input snapshots, SI trajectories, separate
physical/active masses, sample losses, conservation diagnostics, fit history,
numerical refinements, source hashes and comparison figures are retained.
See [the equations, data contract and verified findings](../../docs/film_particle_adsorption.md).
