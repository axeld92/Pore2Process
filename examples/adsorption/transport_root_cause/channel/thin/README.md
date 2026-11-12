# Thin periodic shear follow-up

The full-cube shear fields are independent of z. This separately guarded
follow-up uses four periodic z layers, first verifies that reduction against
the full N32 cube, and then checks N96, N128 and an N128 half-time-step case.
It preserves the stock three-dimensional OpenLB operator and the physical
x-y problem; it is a diagnostic cost reduction, not a representation of a
three-dimensional packed bed.

See the frozen [protocol](protocol.md) for the equivalence gate and comparison
criteria. The parent full-cube study must complete first.

From the repository root:

```sh
source scripts/env.sh
uv run python -m examples.adsorption.transport_root_cause.channel.thin.run
uv run python -m examples.adsorption.transport_root_cause.channel.thin.run --check
uv run python -m examples.adsorption.transport_root_cause.channel.thin.plot
```

Run outputs, source snapshots, executable hashes and exact commands are below
`results/adsorption/transport_root_cause/channel/thin/`. Current source and
native artifact hashes are checked on replay. No new DWSIM or adsorption
calibration is performed, and no successful all-fluid control by itself
qualifies the FCC tracer closure.

The [report](REPORT.md) and [captured reference summary](reference.json)
record the completed equivalence, spatial and time-step checks.
