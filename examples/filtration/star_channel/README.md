# Explicit star-channel and biological-feed showcase

This extends the existing examples without replacing their saved results.
The channel has six explicit points (2.3 mm outer and 1.4 mm inner radii), a
0.3 m active length and an unambiguous area/perimeter/flow convention. It is a
synthetic geometry, not a reconstruction of a paper's specimen.

OpenLB uses its stock D3Q19 descriptor, TRT collision, exact-difference forcing,
halfway bounce-back, STL geometry facilities, plane integrals and ValueTracer.
The dynamics alias assembles existing components through OpenLB's `Tuple`;
it adds no numerical implementation. Pressure-driven supports use stock TRT
and Bouzidi walls. All source, build and output additions are outside OpenLB.

Julia independently verifies the duct conductance using a polygon Poisson
reference, and translates native OpenLB measurements into a distributed
permeating channel and conservative dilute COD balance. DWSIM actually solves,
saves and reloads the stock equipment/recycle flowsheet. COD is an external
oxygen-demand tracer; DWSIM retains the water-carrier property model.

From the project root:

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/bootstrap_julia.py
uv run python -m examples.filtration.graded_support.run --protocol star_channel/support_protocol.json --output results/graded_support_trt --workers 2
uv run python -m examples.filtration.star_channel.native_run
uv run python -m examples.filtration.star_channel.run
uv run python -m examples.filtration.star_channel.check
uv run python -m examples.filtration.star_channel.reference
uv run pytest -q
uv run panel serve app.py --address 127.0.0.1 --port 5006 --show
```

The first two simulation commands can take substantial time. Their outputs
have unique run directories. `latest.json` is published only after acceptance.
Use `--require-sources` with `star_channel.run` to require the original supplied
article Markdown/images; otherwise retained source extractions remain usable.
Existing source files are hash-checked when present.
After an interrupted support run, `--reuse-native /absolute/path/to/prior/run`
reuses only measurements with matching binary, inputs and artifact hashes;
analysis, connectivity checks and actual DWSIM cases are repeated. The support
launcher permits 3,600 seconds per native solve; this wall-time bound is
separate from its unchanged steady-state and accuracy criteria.

The interface has separate capillary-network and filtration tabs. A filtration
run launches an actual DWSIM child process, keeping Julia calls on the server
event loop. Inputs and downloads belong to each browser session. Editing inputs
marks a result stale; validation or simulator failure preserves the last success.

The new biological case is **clean initial hydraulics with illustrative COD
retention**. The measured source retention is transferred with a 0.5–0.8
sensitivity. The paper's critical-flux measurements occur at much higher Re;
critical/sustainable flux for this new geometry and speed remains unavailable.
No laboratory data collection, inaccessible raw data or specimen reconstruction
is required to reproduce this computational example.

See [the scientific report](../../../docs/star_channel_showcase.md) for the verified
range, results, assumptions, failed configurations and exact checks.
