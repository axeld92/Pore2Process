# Python–Julia capillary-network reference

[Scientific report: purpose, theory, results and conclusions](REPORT.md).

This launcher uses the existing sparse Julia solver and Python export pipeline.
No second implementation is introduced.

```sh
uv run python -m examples.hydraulics.capillary_network.run --output results/default
uv run panel serve app.py --address 127.0.0.1 --port 5006 --show
```

See [setup and model derivation](../../../docs/guides/setup.md).
