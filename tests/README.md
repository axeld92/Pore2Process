# Tests by application

| Folder | Coverage |
| --- | --- |
| `hydraulics/` | Original network, Python–Julia contract and numerical/export reference |
| `filtration/` | Calibration, supports, circuits, capacity policies and saved process records |
| `adsorption/` | Conservative uptake/transport and literature comparisons |
| `integration/` | Shared interface, path discovery and provenance integrity |

From the configured root, use for example:

```sh
uv run pytest --collect-only -q
uv run pytest tests/adsorption -q
uv run pytest tests/integration/test_layout.py -q
```

Some tests call actual Julia; the UI regression needs normal host graphics.
They are distinct from the expensive external OpenLB/DWSIM study commands in
each example README. Moving tests does not change their scientific thresholds.
