# Resource organization — 10 November 2026

This maintenance change applies the shared documentation roles, example-local
figures/references, per-paper navigation, grouped tests and future-example
template. It does not change a numerical model or repeat a scientific study.

## Structure and preserved evidence

- Shared documentation is divided into guides, architecture, dated verification
  records and history. Example `REPORT.md` files are the maintained scientific
  explanations; the corresponding `README.md` files retain run commands.
- 124 existing files move to their new locations. Original scientific JSON/CSV,
  image files, protocols and numerical implementations retain their bytes.
- The 16 examples have metadata with prerequisites and distinct validation
  status. Thirteen paper folders include the eight membrane/crossflow sources,
  adsorption sources/candidates and the existing film-correlation citation.
  Multi-paper records remain shared, and no new paper extraction is claimed.
- Tests are grouped under hydraulics, filtration, adsorption and integration.
  A copyable template starts as planned/unverified, with explicit run/check stubs.
- Accepted run folders, source snapshots, installed tools, article inbox and
  unrelated branding work retain their original locations.

## Resource and source integrity

`project_path()` resolves explicit old resource names. The second ordered
relocation registry records the pre-move revision and current source guards.
Historical source verification checks both sides of each applicable move,
requiring exact archived bytes and rejecting unrelated implementation edits.
New run fingerprints include case metadata and both registries. Reference
outputs are excluded from incidental source-directory scans, while explicit
scientific input references remain fingerprinted.

## Verification scope

The maintenance checks cover links, metadata/prerequisites, test collection,
resource digests and guarded historical source resolution. Targeted provenance
tests include deliberate unreviewed-source and corrupted-archive failures.
No new OpenLB/DWSIM solve, Julia scientific study or browser-physics acceptance
is required by this folder change. The original physical comparison outcomes,
including failed adsorption columns, remain unchanged.

## Completed maintenance checks

| Check | Outcome |
| --- | --- |
| Repository navigation/metadata | 16 examples, 13 paper folders, 136 Markdown files; no errors |
| Grouped test discovery | 137 tests collected; the original 133 remain discoverable |
| Targeted path/provenance tests | 9 passed; includes four new resource/transition checks |
| Original numerical/data/image files | 134 checked; all retain their exact bytes |
| Pre-move source/resource resolution | 379 before digests resolved with exact-byte checks |
| Existing Python example/test numeric literals | 126 files checked; no changes |
| Real data-reader smoke checks | Protocol alias, 21-point Zhu PDF reading and frozen adsorption evidence pass |
| Standard-library preflight | Isolated `--help` invocation passes without project imports |

Source-only inspection of eight retained result records checks 93 unique
recorded source pairs. Eighty-eight resolve; five older whole-source hashes
already differed from the pre-move current implementation. Comparing against
the pre-move digests/registry finds **zero organization regressions**. The
pre-existing differences concern older versions of the shared DWSIM runner,
Julia package entry module, filtration interface, visualization and channel
exporter. They are not accepted by weakening guards, and this maintenance
record does not claim a fresh scientific revalidation of those old versions.

The unchanged numerical/data hashes and source comparisons are maintenance
evidence. No native/process study or Julia scientific suite was rerun. The
existing local Panel service was restarted to load current resource paths;
startup/HTTP checks are distinct from its previously recorded browser acceptance.

[Layout guide](../architecture/project_layout.md) ·
[Example catalogue](../../examples/README.md) ·
[Relocation registry](../history/resource_organization.json).
