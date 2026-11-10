# Higher-Re hydraulics and a qualified operating map

This extension preserves the completed capillary, support, crossflow and
low-Re star-channel examples. It uses stock OpenLB D3Q27/TRT in a separate
application folder. The numerical model remains usable without Panel.

Verified on 9–10 November 2026: 16 native cases, 91 actual DWSIM
solve/save/reload cases, 902 Julia assertions and 74 Python tests. Independent
replay covers 682 local artifacts and rechecks upstream evidence. Native
verification is for a fully developed laminar branch, approximately Re=1000–2000;
physical transition and sustainable filtration remain outside that claim.

From the project root:

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/preflight.py --output results/operating_map/preflight.json
uv run python -m operating_map.evidence
uv run python -m operating_map.native
uv run python -m operating_map.run
uv run python -m operating_map.check
uv run python -m operating_map.reference
julia --project=julia -e 'using Pkg; Pkg.test()'
uv run pytest -q
```

Native verification can take tens of minutes or longer. Four independent
single-core native processes may run; Julia stays on the calling thread.
Only after every numerical screen passes is `native/latest.json` published.
The first native pass used a 3,600-second ceiling; the longer-domain case
was restarted from its own measured field with a 7,200-second allowance.
This is a new initial estimate, not an exact population checkpoint, and all
stationarity/refinement criteria remain the same.
The process matrix uses two independent DWSIM child processes; its own pointer
is published after replay checks. Per-solve native and DWSIM limits are 7,200
and 180 seconds, separate from numerical acceptance criteria.

Interrupted, unaccepted runs can resume in place with:

```sh
uv run python -m operating_map.native --resume /absolute/path/to/native/run
uv run python -m operating_map.run --resume /absolute/path/to/process/run
```

Requests, binary and initial-field hashes must match before reusing native
results; geometry is regenerated and byte-checked. New native warm files
include measured non-equilibrium stress, converted to the target lattice units
and supplied to stock `iniRegularized`. Every warm run still evolves and must
pass both stationarity tracers and the minimum observation interval. Native interrupted histories
are retained under `attempts/`. Accepted runs are immutable. Process resumption
requires unchanged source/input snapshots and revalidates actual saved outputs.

The primary biological-feed source is Chiu et al. (2006), Fig. 6 and §3.2.
Garcia (2008), Chiu (2005) and Frederic (2018) have separate roles; no feeds or
specimens are pooled. Existing supplied source files are hash-checked; when
absent, retained extractions can be replayed with the absence recorded.

The map compares three supports, three mean fluxes, three circulation levels
and three fixed-impedance scenarios. Source-capacity fractions, transfer factors,
Reynolds/velocity mappings and interpolation/monotonic policies are evaluated
separately from hydraulics. Nine additional 1.2-times-hydraulic-viscosity cases at 12 LMH/fixed factor 1
keep inlet Re fixed, so flow increases with viscosity. DWSIM retains its water
carrier thermodynamics. These are paired sensitivities, not a variable over
which the main selection is optimized. At most 90 distinct process cases are
solved and saved; a policy selection is not another DWSIM solve.

Outputs include full input/source snapshots, comparison/selection CSVs, SI axial
profiles, actual `.dwxmz` archives and flowsheets, a measured native-field VTP/STL,
verification hashes and a scientific figure. Native CSV `j,k` are local block
indices; physical y/z=(index-0.5)*dx. VTP IDs are zero-based. Unavailable values
are JSON null and blank CSV cells.

The existing Panel showcase keeps its independently verified low-Re controls.
The higher-Re map is inspected through its generated report, figure and CSVs;
no changed process velocity is presented as a new native field.

See [the scientific report](../docs/operating_window.md) for actual acceptance
results, assumptions and limits.
