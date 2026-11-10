# Project layout

Organized on 10 November 2026. Runnable cases are grouped by application;
shared numerical and runtime code has a stable location.

```text
examples/
  hydraulics/
    capillary_network/       Original Python–Julia reference
    capillary_benchmark/     Native capillary benchmark
    darcy_pump/              Measured property → actual DWSIM pump
    tool_checks/             Independent native/runtime probes
  filtration/
    calibration/
    layered_support/
    graded_support/
    once_through/
    crossflow_recycle/
    literature_benchmark/
    star_channel/
    operating_map/
  adsorption/
    alhamed/                 Initial reduced-model comparison
    film_particle/           Batch/column film–particle comparison
    column_contact/          Inventory/contact interpretation
src/pore2process_demo/        Shared interface, bridge and utilities
  runtime/                   Shared external OpenLB builder and DWSIM pump runtime
julia/                       Existing numerical modules, tests and environment
scripts/setup/               Standard-library native-tool setup and inventory
scripts/                     Python–Julia bootstrap and preflight
```

The [example catalogue](../examples/README.md) links to each case. Imports use
the installed `examples` package, so launch commands do not depend on adding the
current directory to `PYTHONPATH`. Run `uv sync --locked` once after the move to
refresh the editable package. `app.py` and the original `pore2process` CLI remain
available. The shared Julia project and dependency locks are retained.

For example:

```sh
source scripts/env.sh
uv sync --locked
uv run python -m examples.hydraulics.capillary_network.run --output results/default
uv run python -m examples.filtration.operating_map.check
uv run panel serve app.py --address 127.0.0.1 --port 5006 --show
```

A saved-study checker requires that study's completed result files. To create
new results, follow the [ordered reproduction guide](reproduction.md) and each
case README. The [filtration overview](../examples/filtration/README.md) preserves
its closed status; adsorption retains its own work-in-progress evidence and
failed physical comparisons.

## Shared responsibilities

- `paths.py` finds the checkout independently of example nesting and working
  directory. `PORE2PROCESS_ROOT` can explicitly select a checkout.
- `io.py` contains the existing strict JSON and SHA-256 helpers.
- `runtime/openlb.py` contains the existing external builder and vendor audit.
- `runtime/hydraulics.py` contains the existing native capillary measurement and
  hydraulic process checks, extracted without changing their numerical bodies.
- `runtime/dwsim.py` contains the existing shared water-pump adapter runtime;
  process-specific adapters remain with their cases.
- `provenance.py` resolves legacy source paths for historical verification and
  includes shared helpers in new run fingerprints.

Original OpenLB files, installed tool locations and scientific models are not
restructured. C++/Julia numerical files are preserved byte-for-byte. The native
include path and .NET project root lookup are adjusted for the new source depth.

## Existing evidence and source history

Existing `results/`, `.local/`, `literature/inbox/`, curated JSON/CSV extractions,
Julia modules and manifests remain at their existing locations. Old run IDs,
input snapshots, source snapshots, hashes and result values are not rewritten.

[history/reorganization.json](history/reorganization.json) records the old/new
source paths, the immutable pre-move Git revision and reviewed post-move hashes.
When a historical source hash differs solely because of the recorded relocation,
the checker verifies the current guarded implementation and exact archived
pre-move bytes. It rejects subsequent source edits or corrupted archived bytes.
This preserves strict historical replay without presenting old hashes as hashes
of the reorganized source. New runs fingerprint the current paths and extracted
shared modules.

The pre-move source is kept in Git history; local verification materializes only
requested source files in ignored `.local/reorganization/source_cache/`. A
shallow clone must fetch the relevant history before historical replay. Existing
run artifacts and native/simulator evidence continue through their original
hash, conservation and numerical validators.

The original local development history is preserved on
`codex/checkpoint-before-reorganization-20261110`. The organization branch has
GitHub's published history as its ancestry, avoiding replacement of the public
overview or loss of the local checkpoint. The original handoff and pre-move local
README are preserved verbatim in `docs/history/`.

## Checks for this move

The post-move Python suite passed **133 tests**, including installed launchers,
root discovery, historical source integrity and existing export checks. Direct
Julia checks passed **902 package assertions** and **147 adsorption assertions**.
Julia numerical files, C++ source/headers and dependency locks remain unchanged.

Saved filtration evidence replay passed 91 DWSIM cases and 682 artifacts; saved
synthetic and literature adsorption comparisons also replayed with their prior
physical pass/fail conclusions preserved. These replays reuse saved native and
simulator evidence rather than regenerate the scientific studies.

The relocated native builder and common DWSIM adapter built successfully. A
small fresh capillary case completed, and the Panel filtration page completed
an actual DWSIM solve/save/reload with its download enabled. Automated browser
download capture timed out; a newly downloaded ZIP was not inspected during
this move. ZIP generation/roundtrip checks passed in the Python suite.

An additional multi-resolution capillary rerun was unnecessary for the move.
It reached its timeout at resolution 40 and was not continued or counted as a
completed benchmark. Its partial outputs are retained under
`results/reorganization/`; the accepted studies and their result records remain
in their original locations. No further full scientific reruns are needed to
finish this folder organization.
