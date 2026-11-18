# Liquid readout and internal bead transport

[Bottom line and results](REPORT.md) · [Protocol](protocol.md)

The [phase calculation](phase/README.md) checks how the reported
concentration conversion could change with sample volume or phase
selection. The [particle calculation](particle/README.md) tests the
full catalyst mass under independently sourced transport scenarios.

From the project root with the [configured environment](../../../docs/guides/setup.md):

```sh
source scripts/env.sh
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transport.phase.run
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transport.particle.run
uv run --locked --offline python -m examples.biocatalysis.pentyl_propanoate_transport.check
uv run --locked --offline pytest -q examples/biocatalysis/pentyl_propanoate_transport/phase/test_phase.py
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate_transport/particle/test_model.jl
```

Runners create new folders under `results/biocatalysis/pentyl_propanoate_transport/`.
Both write `predictions_before_comparison.csv` before opening the
reactor conversion measurements. The checker verifies curated outputs and
input hashes from the two run records. Fractions in the CSVs range from
zero to one for the reactor conversions. The phase projections are only
available where the measured tie lines support interpolation.
