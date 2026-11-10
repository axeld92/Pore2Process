# Adding an example

Copy [the example template](../../templates/example/README.md) into a suitably
named `examples/hydraulics/`, `examples/filtration/` or `examples/adsorption/`
subfolder. Another application area can be added when the research question
requires it. Do not overwrite an existing case or restructure completed runs.

## Files and responsibilities

| File/folder | Purpose |
| --- | --- |
| `README.md` | Prerequisites, exact setup/run/check commands and file contracts |
| `REPORT.md` | Purpose, theory/equations, assumptions, results, discussion, conclusions, references |
| `example.toml` | Stable ID, area, status, prerequisites, entry modules and evidence links |
| `protocol.json` or `protocol.md` | Declared physical/numerical inputs, split, screens and comparison scope |
| `run.py` | Coordinates actual execution using independent numerical models |
| `check.py` | Verifies saved evidence without silently replacing it with new simulations |
| `figures/` | Figures referenced by the report |
| `reference/` | Compact retained observations with explicit run/provenance links |
| `evidence/` | Optional example-specific transferred inputs or multi-source scenarios |

The template's run/check stubs deliberately fail until implemented. Its result
and validation states start as unavailable/not assessed. An empty result is not
a scientific pass. Existing examples retain their established entry names,
including `study.py` where already used; metadata records the actual entry point.

## Metadata schema, version 1

Paths in `reference_files` and `protocol_files` are relative to the example.
`id` is its path below `examples/`, and `area` is the first path component.
`prerequisites` lists other example IDs needed for the documented workflow;
these are not claims of independent experimental evidence. The graph must be
acyclic. `run_module`/`check_module` identify existing Python modules, but this
metadata is not an automatic job runner. An empty checker module means the
run guide must identify the appropriate check route.

Allowed lifecycle states are `planned`, `in_progress`, `complete`, `closed`.
Validation states include `not_assessed`, `not_applicable`, `passed`, `failed`,
`conditional`, `scenario_only`, `mixed`. Record batch/column or other subchecks
separately when their outcomes differ. These labels summarize the report's
evidence and never replace its scope or acceptance thresholds.

## Implementation and evidence

Keep custom numerical models usable independently of the interface. Reuse stock
OpenLB facilities, at least TRT for porous media, and the installed examples'
lifecycle. Keep all project sources/builds outside the vendor tree. Select the
physical application before choosing which effective quantities reach DWSIM.

Declare units, coordinates, indexing, geometry/area conventions and unavailable
inputs. Seek published experimental inputs with specimen, source location and
conditions retained; no new laboratory measurements or unavailable raw dataset
is a prerequisite. Do not subtract unrelated specimens to infer a layer property.

Validate physical equations against suitable analytical/numerical/literature
references, checking conservation and convergence as appropriate. Preserve
failed screens. Record actual simulator execution separately from saved-result
replay and browser inspection. Numerical changes require new verification;
documentation-only changes require appropriate document checks.

## Checks and publication

1. Add the report/run guide to the [catalogue](../../examples/README.md).
2. Run `uv run python scripts/check_repository.py` for links, metadata and paths.
3. Collect the grouped tests with `uv run pytest --collect-only -q` and run the
   checks relevant to the change. External native/process studies are explicit
   commands in the run guide, never an accidental consequence of metadata checking.
4. Keep full runs under ignored `results/<area>/<example>/runs/<run-id>/` for
   new cases. Existing accepted result paths stay where they are.
5. Publish source, locks, curated evidence and a compact reference; keep full
   article files, runtime installations, binaries and routine runs out of Git.
