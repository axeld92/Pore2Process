# Spatially varying periodic shear control

This diagnostic removes solid surfaces while retaining a prescribed velocity
that changes across the box. It asks whether the audited OpenLB bulk scalar
model can transport a smooth concentration profile at the application flow
and molecular diffusivity without the FCC particle boundaries.

The velocity is a smooth, divergence-free cosine shear in an all-fluid
periodic cube. It is not a wall-bounded channel validation, a resolved FCC
simulation, or an accepted column dispersion coefficient. See the frozen
[protocol](protocol.md) for physical scales, the comparison tolerances and
the independent continuum reduction.

From the repository root:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.channel.run
uv run python -m examples.adsorption.transport_root_cause.channel.run --check
uv run python -m examples.adsorption.transport_root_cause.channel.plot
uv run pytest tests/adsorption/test_transport_channel.py -q
```

The five native cases use stock regularized TRT/D3Q27 at a matched physical
time step: uniform sinusoidal transport at N32, constant preservation under
shear at N32, and sinusoidal transport under shear at N32/48/64. Python's
Fourier calculation is solely an independent reference for this separable
verification problem; it is not a second production pore solver.

Run folders and their exact native commands are saved below
`results/adsorption/transport_root_cause/channel/`. The `latest.json` pointer
identifies the completed guarded summary. Replay verifies current sources,
source snapshots, native outputs and numerical analysis without executing a
new simulation. Build products stay in `.local/adsorption/transport_root_cause_channel`.
Original OpenLB source is checked against the pinned distribution archive.

All concentration fields are kg/m3, coordinates are metres, velocity is m/s
and indices are zero-based. `inventory.csv` records every native time step.
`initial.csv`, `peak.csv` and `final.csv` hold complete x-y sections at k=0
(z=dx/2); each row
also reports the maximum concentration variation over the homogeneous z
direction. `peak.csv` corresponds to the greatest observed departure from
the prescribed continuum interval (or the initial state when no larger
departure occurs); `raw.json` supplies its step and physical time.

See [REPORT.md](REPORT.md) for the actual executed results and limits. A
small [captured reference summary](reference.json) retains the native summary
hash and checked results alongside the source; full outputs remain in the
run directory.
