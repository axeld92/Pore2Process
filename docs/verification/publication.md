# Initial source publication

> **Dated verification record.** This retains the execution/audit history.
> Use the [maintained scientific report](../../examples/README.md) for the current explanation
> and the example README for current commands.

Prepared 9 November 2026 after the completed explicit-channel, biological-feed
sensitivity and local-interface handoff.

## Snapshot and provenance

The public source is a snapshot of the completed local implementation, including
its earlier examples, curated literature extractions and compact verification
references. The original development checkout and commit history are preserved
locally. The publication is based on GitHub's existing documentation branch and
retains its initial README/roadmap commits.

Numerical source, simulator adapters, protocols, tests, locks and existing
scientific records were copied without modification. The public landing page
adds navigation and the current result; the previous detailed README is retained
as [setup_and_examples.md](../guides/setup.md), with links adjusted for its
new location. [reproduction.md](../guides/reproduction.md) lists the ancestor stages that
are required before launching the connected filtration showcase.

The snapshot excludes full article files, installed Python/Julia/OpenLB/DWSIM
runtimes, local caches, native binaries and routine generated runs. Only the
specifically selected project logo is included from the separate branding work.
Pinned third-party distributions are obtained through the existing setup tools.

## Verification scope

The completed development handoff records:

- 49 native TRT cases: 17 channel, 24 support and eight capillary regression cases.
- 28 actual DWSIM process cases solved, saved, reloaded and checked.
- 878 Julia assertions and 61 Python tests.
- Browser interaction, geometry rotation, separate session state, invalid-input
  recovery and downloaded ZIP content checks.
- 2,070 original OpenLB files checked against the pinned distribution archive.

Publication checks compare final UI and screenshot hashes with the recorded
browser evidence, inspect the upload for accidental credentials and excluded
artifacts, validate documentation links and exercise the numerical reference
from a separate checkout. These checks do not constitute a repeat of all native
simulations or new physical validation.

The [implementation status](first_simulation.md) and
[scientific report](filtration/star_channel.md) distinguish historical BGK evidence,
current stock TRT results, published measurements and scenario assumptions.
Neither sustainable flux nor physical fouling time is claimed for the new
geometry and low-Re operating range.

## Checks repeated for publication

The separate publication checkout completed `uv sync --locked --offline` using
cached packages and bootstrapped Julia 1.11.6 with PythonCall/JuliaCall 0.9.37.
Direct Julia tests passed **878 assertions**; `pytest -q` passed **61 tests**
in 129.15 seconds. This timing includes a fresh local environment and is not a
solver-performance claim.

The default CLI generated all five required JSON/CSV/PNG/VTP exports and returned
`Q_in = 7.85398163397448e-10 m3/s`,
`k = 2.454369260617024e-12 m2`, with normalized conservation error
`1.97e-15`. The independent reference-check command also passed. Existing
scientific references were retained rather than replaced with publication-run
timing observations. All relative Markdown file targets resolved.

The requested logo was copied without image edits. The upload was checked for
credential patterns and excluded environment/article/simulator artifacts. No
new native OpenLB or DWSIM study was run for this publication; their completed
handoff evidence remains the source of those verification claims.

## Filtration closure update — 10 November 2026

The closed higher-Re filtration continuation is now included: native application,
Julia capacity/momentum-screen functions, source audit, conditional operating map,
compact measured/simulator reference, scientific figure and tests. The stock-DWSIM
adapter's input ceiling increases from 5,000 to 20,000 iterations for this study;
existing defaults and scientific acceptance thresholds remain unchanged.

The publication was prepared in a separate checkout. Final native/application
and process implementation hashes match the accepted filtration runs. The
isolated environment installs with `uv sync --locked --offline` and bootstraps
the existing Julia runtime. Repeated checks pass **902 Julia assertions** and
**74 Python tests** (112.43 seconds for that Python suite). Full saved-evidence
replay passes **91 DWSIM cases**, **682 artifacts**, **1,620 capacity decisions**
and **1,458 discrete policy selections**, including the upstream native/support
checks. Replay launches no new external simulations and does not rewrite the
accepted run folders. Existing tool distributions were read from the local
cache; no simulator distribution is published.

The README and roadmap mark filtration hydraulics closed for now, preserve the
selected logo and eight article references, and identify isothermal adsorption
as the next case in progress. Active adsorption code, tests, source extractions
and draft case reports are excluded. The original development checkout and its
concurrent work remain untouched.

## Example organization — 10 November 2026

The development workspace now groups cases under `examples/hydraulics/`,
`examples/filtration/` and `examples/adsorption/`. This update includes the
existing adsorption implementation and comparison reports as work in progress;
it does not change their physical validation conclusions. Shared runtime and
setup helpers have stable locations, and documentation uses the new commands.

The pre-move source snapshot is preserved in the published Git ancestry.
Accepted run folders and their source snapshots are retained, with guarded
legacy-source replay documented in [project_layout.md](../architecture/project_layout.md).
That guide also records the targeted relocation checks and the incomplete
extra benchmark. The logo, eight membrane/crossflow references, dependency locks
and original numerical implementations are preserved.
