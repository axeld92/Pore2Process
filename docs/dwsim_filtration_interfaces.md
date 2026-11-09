# DWSIM hydraulic filtration interface decision

Inspected the installed target assemblies on macOS 26.6.2 ARM64, 8–9 November
2026. DWSIM FluentAPI reports **10.2.15.0**, with the private .NET **10.0.403**
SDK (10.0.14 runtime). CPython remains **3.12.9** and Julia **1.11.6**. No
additional runtime or package installation was needed. Build outputs live in
`.local/process_filtration/`; the installed DWSIM and original OpenLB remain
untouched. This is an actual target-machine check, not the hosted inventory.

## Interfaces examined

Reflection of the installed FluentAPI establishes these available builders:

| Stock component | Relevant supported interface | Use here |
| --- | --- | --- |
| Pump | Outlet-pressure mode, efficiency, material and energy connections | Pressurize the full fresh feed |
| Splitter | `WithSplitRatios(double[])`, numbered product connections | Permeate/feed mass recovery; remaining mass goes to retentate |
| Valve | `WithPressureDrop(Quantity)`, outlet-pressure and Kv alternatives | Separate channel loss, permeate TMP and retentate letdown |
| Mixer | Pressure behavior and numbered feed connections | Available for a later recycle loop |
| CrossflowUF | Concentration/diafiltration modes, VCF, sieving, specified flux, TMP | Examined; not selected for this hydraulic closure |

The reproduction runner supports `--inspect` and saves
`adapter_build/installed_api.json` from the actual installed assembly on each
run. `Inner` exposes the supported underlying flowsheet API if a later recycle
requires a generic logical operation. A recycle was not implemented or tested
in this milestone; presence of mixer/generic interfaces does not establish
recycle convergence.

The subsequent [crossflow extension](crossflow_recycle.md) now verifies that
interface in actual solved/saved/reloaded loops, including seed and tolerance
sensitivity. This page retains the once-through interface decision above.

## Why the splitter and valves

The official [CrossflowUF reference](https://dwsim.org/help/api/bioprocess/crossflow-uf.html)
describes prescribed flux/VCF/sieving and identifies TMP as a reported input,
unused in the flux calculation. It therefore does not supply our required
OpenLB-derived pressure–flux law. Its chemical separation capabilities are
also unnecessary for the current pure-water question. We preserve that stock
model as a candidate for a later bioprocess case; no custom membrane unit or
replacement DWSIM dynamics have been created.

The stock [splitter](https://dwsim.org/help/api/mixing-separation/splitter.html)
distributes material at inlet conditions. Stock
[valves](https://dwsim.org/help/api/pressure-changers/valve.html) in pressure-drop
mode then apply explicit hydraulic loads with isenthalpic flashes. Julia sets
the TMP and feed flow before the flowsheet solve. The pressure drop is a
linear hydraulic resistance supplied externally, not a valve Kv correlation
or a DWSIM prediction of filtration selectivity. The inlet channel-loss valve
is a separate declared lump, not the membrane's axial permeability used twice.

The [official unit-operation guide](https://dwsim.org/help/main/24-unit-operations.html)
documents splitter port order and mass/energy propagation. Actual installed
builder methods, solver outputs and the saved XML connections are checked as
the target-version evidence; online documentation alone is not an execution
test. These primary references were consulted during implementation.

## Data boundary

`process_filtration/model.py` calls the same Julia package through the existing
single bridge. JSON passes absolute pressure in Pa, temperature in K, mass
flow in kg/s, TMP/channel loss in Pa, and dimensionless mass recovery and pump
efficiency. All inputs must be finite; recovery is strictly between zero and
one and efficiency in (0,1]. The two nonempty product branches are intentional.

The adapter converts DWSIM's kW to W and kJ/kg to J/kg exactly once. Each stream
returns its pressure, temperature, mass/actual volume flow, density, liquid
viscosity, enthalpy, vapor fraction and water fraction. Recovery is a mass ratio;
reference volume = mass / common fresh-feed density. Actual flashed product
volume is also retained and need not equal that reference volume.

Seven streams, five stock units and one energy stream form the flowsheet.
The saved `.dwxmz` XML is inspected for all twelve directed connections,
including splitter outlet 0 = permeate branch and outlet 1 = retentate branch.
This is separate from the pore/network zero-based node convention, which is
unchanged. No new array-index conversion or mutable simulation state is added.
