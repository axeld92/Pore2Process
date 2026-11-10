# Documentation

The maintained scientific explanation of each example is its `REPORT.md`.
Start with the [example catalogue](../examples/README.md). Use the corresponding
`README.md` for commands and `example.toml` for machine-readable identity,
prerequisites, reference files and validation status.

| Documentation | Contents |
| --- | --- |
| [Guides](guides/README.md) | Setup, reproduction, research directions and adding examples |
| [Architecture](architecture/README.md) | Layout, software/data boundaries and simulator policies |
| [Verification records](verification/README.md) | Dated execution, acceptance and source-audit records |
| [History](history/README.md) | Original handoffs, planning and relocation registries |
| [Literature](../literature/README.md) | Paper metadata, extractions and shared evidence reviews |

Example figures and compact numerical references live in each example's
`figures/` and `reference/` folders. Full runs remain under their existing
ignored `results/` locations. Installed tools and build caches remain in `.local/`.

Run `uv run python scripts/check_repository.py` for documentation, metadata,
dependency and resource-path checks. It launches no numerical or process solver.
