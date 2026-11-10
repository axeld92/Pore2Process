# OpenLB and DWSIM compatibility — 8 November 2026

> **Dated verification record.** This retains the execution/audit history.
> Use the [maintained scientific report](../../../examples/hydraulics/tool_checks/REPORT.md) for the current explanation
> and the example README for current commands.

Both engines ran natively on this Mac and can be coordinated by Python. This
completes the tool-compatibility task. It does not complete their coupling or
select a real material/application. The existing capillary-network solver,
Python environment/lock, Julia environment/manifest, and Panel application were
preserved.

## Machine and setup

Project: `/Users/dula/Documents/ChatGPT/Pore2Process`. macOS 26.6.2, arm64,
Python 3.12.9. No applicable AGENTS.md was found in the project or its ancestors.
The initial Git working tree was clean at `5199f33` on `codex/first-simulation`.
An existing OpenLB source tree in another research project was inspected only
for configuration clues; this task did not build or edit that tree.

| Component | Observed / tested |
| --- | --- |
| GCC C/C++ | Homebrew 15.2.0; used for the OpenLB build |
| Apple clang | 21.0.0; discovered, not used for this build |
| GNU make | 3.81; used |
| tinyxml2 / zlib | 11.0.0 / 1.2.12; linked Homebrew tinyxml2 and macOS zlib |
| Open MPI / CMake | 5.0.8 / 4.0.0; discovered, unused |
| Git | 2.49.0 |
| OpenLB | Official 1.9r0 source archive; native ARM64 executable |
| DWSIM | Official 10.2.15 Apple Silicon portable distribution; API reports 10.2.15.0 |
| .NET | Private SDK 10.0.403 and runtime 10.0.14, ARM64 |

DWSIM was not found in the standard application locations; dotnet/Mono were
not on PATH. Official archives and a private SDK were installed under ignored
`.local/compatibility/`. No system Python or existing Julia installation was
replaced. The archived source/package URLs and pinned SHA256/SHA512 values are
in the bootstrap scripts and reference JSON. Repeating both bootstraps passed
hash checks and reused their installations. Initial downloads required external
network access, which was approved. No outstanding download/install blocker
remains for these probes.

Exact setup, launch, and validation commands are in
[compatibility/README.md](../../../examples/hydraulics/tool_checks/README.md). Generated
`results/compatibility/environment.json` records the bounded target-machine
inventory, source records, and OpenLB linked library paths. Dependency updates
on this Mac can affect the OpenLB binary; rerun the inventory and benchmarks
after compiler/library changes. These source/SDK pins do not replace a future
portable build environment for other operating systems.

## OpenLB numerical check

The unmodified official `examples/laminar/poiseuille2d` ran in serial double
precision using BGK / D2Q9, Re=10, relaxation time 0.8, domain 2 m × 1 m,
characteristic maximum speed 1 m/s, and density 1 in the example's 2D convention.
It imposes a parabolic velocity inlet, zero-gauge-pressure outlet, and
interpolated no-slip walls. It is a channel-flow compatibility benchmark rather
than the later porous-material model. OpenLB's analytical reference is
uₓ(y)=4Umax(y/H)(1−y/H), with Δp=8ρνUmaxL/H² and a linear pressure profile.
The solver performs lattice-Boltzmann collision/streaming; it does not produce
the solution by evaluating those analytical expressions.

| Quantity | N=32 | N=64 |
| --- | --- | --- |
| Inlet / outlet boundary nodes | 31 / 31 | 63 / 63 |
| Reported convergence iteration | 19,373 | 75,582 |
| Relative velocity L2 error | 0.00279288 | 0.000699683 |
| Relative pressure L2 error | 0.00737092 | 0.00184006 |

Refinement reduced velocity error by 3.9916 and pressure error by about 4.006.
Both cases passed the explicitly implemented acceptance checks (velocity <1%,
pressure <5%, finite metrics, correct boundaries, reported convergence, and
decreasing errors). Two levels are useful evidence; they do not establish
convergence for an arbitrary porous geometry. OpenLB's stop test monitors maximum
velocity. A separate integrated inlet/outlet flux/conservation diagnostic and
permeability export are not implemented in this probe.

An exploratory N=40 run revealed a real limitation: the stock geometry reported
no material-4 outlet nodes. Its velocity converged but density drifted and
relative pressure L2 error reached 4.96933. That run is **rejected** as a physical
benchmark. A resolution-sensitive boundary-classification issue is observed;
floating-point geometry selection is a suspected cause, not a confirmed root
cause. N=32/64 use binary-exact spacing and have the expected outlet counts.
The runner now checks both boundary counts and pressure error, preventing a low
velocity error from hiding this failure. Robust boundary selection across
requested resolutions must be addressed before using an adapted case for
permeability. Exploratory logs remain under `results/compatibility/openlb/n40/`.

Two upstream/tool quirks were also recorded:

- The archive's version rules detect the enclosing Pore2Process Git repository
  and would label OpenLB with its commit. The build explicitly supplies
  `GIT_VERSION=1.9r0`; the archive hash establishes provenance.
- In this example `GNUPLOT_ENABLED=0` enters a profile-PNG branch. The installed
  gnuplot 6.0.2 fails because it expects `libcerf.2.dylib`, while the installed
  libcerf is 3.3. The final numerical runs set `GNUPLOT_ENABLED=1`, which records
  EOC data without invoking PNG rendering, and `VTK_ENABLED=0`. This task did not
  repair the unrelated gnuplot installation or verify OpenLB graphics. Console
  error norms, rather than the example's EOC data file, supply the checks.

## DWSIM automation and numerical check

The supported [.NET automation route](https://dwsim.org/help/main/57-automation.html)
was used through `DWSIM.Automation.FluentAPI.dll`. Python launches a small C#
adapter as a separate process; it does not load CLR into the Python/JuliaCall
process. The macOS portable archive contains the assemblies but no standalone
MCP-server executable was found in it. The compatibility route needs no COM
registration, Windows installation, or desktop automation.

The flowsheet uses pure water and Steam Tables (IAPWS-IF97): feed → pump →
product, plus a connected energy stream. Inputs are 298.15 K, 100,000 Pa inlet,
10 kg/s, and 75% pump efficiency. The baseline outlet is 1,000,000 Pa, matching
the operating point in DWSIM's [official pump example](https://dwsim.org/help/api/pressure-changers/pump.html).
The second outlet is 1,900,000 Pa, doubling pressure rise rather than absolute
pressure. The actual DWSIM solver ran in each case; the adapter then saved,
loaded, and solved each `.dwxmz` file, checking preserved pressure and mass flow.

| Quantity | Baseline | Double pressure rise |
| --- | --- | --- |
| Shaft power | 12,035.5356965 W | 24,071.0713931 W |
| Mass imbalance | 0 kg/s | 0 kg/s |
| Energy-balance residual | 0.000783894 W | 0.000753268 W |
| Pump formula relative error | 0 | 0 |

Feed density was 997.0474354 kg/m³. Baseline product temperature was
298.2386130 K. Both streams remained liquid. The separate physical check uses
P=ṁΔp/(ρη) with density returned by the property package. This checks pump
behavior and units, not an independent validation of the water property package.
Energy residuals are below 1e-6 of shaft power; shaft power scales by exactly
two within floating-point precision. The boundary uses W and J/kg after explicit
conversion from DWSIM's kW and kJ/kg; pressures are absolute Pa and efficiency
is a fraction in JSON. Inputs, results, runtime version, timings, logs, and
saved flowsheets are preserved under `results/compatibility/dwsim/`.

The DWSIM GUI was not launched or visually checked. SDK builds emitted a
macOS CSSM diagnostic in the restricted terminal, but compilation, solves, and
save/reload checks all returned success. Developer-certificate generation is
disabled in the runner. No GUI or rendering pass is inferred from these checks.

## Remaining work

This task establishes two executable engines and a Python-to-DWSIM process
boundary. It does not compute a porous-material permeability with OpenLB,
transfer OpenLB output through Julia, or drive DWSIM with that property. Those
are the next milestones in [integration_roadmap.md](../../history/planning/integration_roadmap.md).
The N=40 boundary issue must be resolved or guarded in the later benchmark;
the optional gnuplot installation remains broken. Neither prevents rerunning
the verified compatibility cases above.
