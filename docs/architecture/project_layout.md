# Project layout

Runnable examples are grouped by application. Their scientific report,
commands, protocols, figures and compact reference now live together.

```text
examples/<area>/<case>/
  README.md                 Setup, run/check commands and contracts
  REPORT.md                 Maintained scientific explanation
  example.toml              Identity, prerequisites and validation status
  protocol.json / .md        Existing declared study protocol
  figures/                  Report figures and screenshots
  reference/                Compact accepted observations
  evidence/                 Optional case-specific transferred evidence
  openlb/                   Project native case, where applicable
src/pore2process_demo/       Shared Python bridge, interface and runtime helpers
julia/                      Independent numerical models and explicit environment
tests/<area>/               Application tests; shared checks under integration/
literature/sources/<paper>/  Paper metadata, single-paper extractions and points
literature/reviews/         Intact multi-paper comparisons
literature/inbox/           Private original articles and supplements
docs/guides/                Setup, reproduction and contribution guidance
docs/architecture/          Shared boundaries, policies and this layout
docs/verification/          Dated execution and audit records
docs/history/               Original handoffs, planning and path registries
templates/example/          Starting files for a new example
```

The original synthetic adsorption case lives directly in `examples/adsorption/`;
its later comparisons are children. Existing entry points retain their module
names. The capillary-network CLI and root `app.py` remain available.

## Navigation and status

Start at the [example catalogue](../../examples/README.md) or
[documentation index](../README.md). Each scientific `REPORT.md` is the maintained
explanation; dated records in `verification/` preserve execution details and
earlier findings. A software/numerical pass does not change a failed physical
comparison. Per-example `example.toml` records that distinction and its actual
prerequisites. The [template guide](../guides/adding_examples.md) defines the schema.

## Locations deliberately retained

Original OpenLB sources, Julia modules, dependency locks, installed tools and
native build caches retain their locations. Accepted run folders, their source
snapshots and artifacts stay in their original ignored `results/` directories.
Only new examples adopt the suggested `results/<area>/<case>/runs/<id>/` convention.
The article inbox and unrelated branding work are untouched.

Every pre-existing scientific JSON/CSV file moved by this organization retains
its exact bytes. Old paths inside frozen protocols, extractions and references
are provenance, not broken data that should be rewritten. `project_path()`
resolves those explicit aliases for current code. Human navigation uses the new
locations and per-paper metadata; shared extractions are not duplicated.

## Historical source integrity

The ordered [first example relocation](../history/reorganization.json) and
[resource relocation](../history/resource_organization.json) registries retain
the pre-move Git revisions, path maps and reviewed source guards.
`verified_source_path()` verifies current guards through each applicable move
and retrieves exact archived before bytes when a historical hash needs them.
It rejects unrelated implementation edits and corrupted archive bytes.
New run fingerprints include both registries and the resource resolver.

The earlier layout record is preserved verbatim in
[history](../history/layout_before_resource_organization.md). Full Git history
is needed to materialize archived sources in a new clone. No symlink forest or
duplicate legacy resource files is required at the old locations.

## Organization checks

```sh
source scripts/env.sh
uv run python scripts/check_repository.py
uv run pytest --collect-only -q
uv run pytest tests/integration/test_layout.py tests/integration/test_repository_layout.py -q
```

These check navigation, metadata, paths, test discovery and provenance rather
than regenerate the scientific studies. Numerical changes remain subject to
their own model-specific verification.
