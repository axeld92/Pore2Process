# Connected hydraulic design comparison

[Scientific report: purpose, theory, results and conclusions](REPORT.md).

This small technical demonstration runs OpenLB, Python, Julia and DWSIM on the
same Mac. It compares two synthetic straight-pore designs, derives permeability
from the actual OpenLB results, calculates a fixed-flow porous pressure load in
Julia, and supplies that pressure requirement to an actual DWSIM water pump.
It does not select or validate a real material.

From `/Users/dula/Documents/ChatGPT/Pore2Process`, after the Python/Julia setup
and the pinned tool downloads in [compatibility/README.md](../tool_checks/README.md):

```sh
source scripts/env.sh
uv run python -m examples.hydraulics.darcy_pump.run
uv run python -m examples.hydraulics.darcy_pump.check
```

The first command builds the external case and DWSIM adapter, obtains the water
properties from DWSIM, runs both pore geometries at N=24 and 40, solves both
DWSIM pump cases, checks acceptance and renders the comparison figure. It makes
fresh solves; it has no analytic-permeability shortcut or solver cache. It needs
several minutes on this Mac. No packages are installed by this command.
`--output /absolute/directory` selects a different output root.

Every attempt creates a separate `results/coupling/runs/<UTC>-<id>/` directory.
`latest.json` points to the most recent completed, passing run. A failed attempt
keeps its logs and partial artifacts and does not replace that pointer. Check
the command's exit status; an older passing pointer does not make a new attempt
successful. All generated results remain ignored in Git.

The checker verifies recorded artifact hashes, recomputes permeability and the
pressure load using real Julia, and rechecks the saved DWSIM results and input
mapping. It does not re-solve DWSIM. The original adapter creates, solves, saves,
reloads and solves each `.dwxmz` flowsheet. To re-run the entire calculation,
use the first command again. Both commands accept an explicit output/summary
location; `uv run python -m examples.hydraulics.darcy_pump.check /absolute/path/summary.json` selects
an immutable snapshot.

## Model and data contract

- Pore radius changes from 50 to 40 µm. Physical pore length stays 400 µm.
  The representative cell is a fixed 200 µm square: area 4e-8 m². Its porosity
  changes with radius. Its transverse center is at y=z=R; the cell extends
  from R−100 to R+100 µm. Translation of the cell origin has no physical effect.
- Current OpenLB code uses the unmodified D3Q19 descriptor, stock TRT dynamics
  (magic parameter 3/16), cylinder/cuboid
  indicators, Bouzidi wall and local pressure boundary implementations.
  The case follows the official example lifecycle, in separate project files;
  see [the code boundary](../../../docs/openlb_code_policy.md).
  Previously accepted BGK comparison runs remain historical; the new
  [showcase](../../../docs/star_channel_showcase.md) uses freshly verified TRT evidence.
- The imposed OpenLB gauge pressure is 0.1→0 Pa. Initial velocity is zero;
  no analytical velocity or permeability is prescribed. Julia infers k from
  integrated midplane flow and interior pressure taps at L/4 and 3L/4.
- SI coordinates are x=iΔx, y=(j+0.5)Δx, z=(k+0.5)Δx. CSV indices are zero-based
  int64; physical fields are float64. Positive flux is +x on every plane.
  The number of axial intervals is explicit; it changes when radius changes.
- DWSIM's Steam Tables (IAPWS-IF97) supplies liquid-water density and dynamic
  viscosity at 298.15 K, 100000 Pa absolute. Those actual values configure
  OpenLB, Julia's Darcy model and the mass-to-volume mapping. There is no
  implicit substitution of 1000 kg/m³ or 0.001 Pa s for those water properties.
- The idealized repeated medium is a 1 m slab with 1 m² area and Q=1e-5 m³/s
  (0.6 L/min), at a pump efficiency of 0.75. The large area represents many
  pores in parallel; the per-pore Reynolds number remains small. All geometries
  are straight, aligned, saturated, and free of junction or entrance losses.
- Julia calculates Δp=μLQ/(kA). Python sets DWSIM mass flow to ρQ and **absolute**
  discharge pressure to 100000 Pa + Δp. DWSIM calculates shaft power, outlet
  state and the energy stream. The separately computed porous load is an
  external requirement on the pump; a porous unit operation is not implemented
  in the flowsheet. No other equipment losses or elevation changes are added.
- Inlet properties are held constant for the Darcy translation. Actual DWSIM
  outlet density and viscosity changes must each remain below 0.01%; process
  maximum pore Reynolds number must remain below 0.05. The assumption and its
  checked range are recorded for each design.
- DWSIM kW and kJ/kg are converted to W and J/kg at the adapter boundary;
  efficiency is a fraction in JSON and percent in the DWSIM API.

Independent Hagen–Poiseuille theory predicts k=πR⁴/(8A_cell). At fixed cell
area, k40/k50=0.4096; at fixed flow and efficiency the pressure and power ratios
should be 1/0.4096=2.44140625. These expressions provide test oracles, and never
supply the permeability used by the process run.

## Acceptance and artifacts

Both designs must pass the hydraulic field, coordinate, pressure, integrated
mass-conservation and finite-value checks. N=40 permeability and profile errors
must be below 1%, and permeability error must decrease from N=24 to 40. The fine
all-plane mass spread must be below 0.1%. The permeability ratio must agree with
R⁴ within 1%; the actual DWSIM power ratio must agree with the inverse measured
permeability ratio within 1e-9 relative.

Each DWSIM solve checks liquid state, pressure, shared inlet properties, mass
conservation, the energy stream, and power against QΔp/η (1e-9 relative).
The energy-balance tolerance is 0.1% of shaft power: the small pressure rises
produce tiny enthalpy differences, so the water flash's numerical residual is
reported separately. It is not claimed to be exact energy closure.

| File | Content |
| --- | --- |
| `fluid_anchor/` | DWSIM property-state probe, inputs, output, log and saved flowsheet |
| `<design>/n24/`, `n40/` | Actual OpenLB raw JSON, SI plane/midplane/history CSVs, log and Julia-derived analysis |
| `<design>/translation.json` | Measured k, property values, Darcy inputs/drop, originating pore run and assumptions |
| `<design>/dwsim/` | Exact DWSIM input, output, log and saved/reloaded `water_pump.dwxmz` |
| `summary.json` | Full snapshots, acceptance results, ratios, source/binary/runtime hashes, versions and timings |
| `verification.json` | Separate artifact-hash and Julia translation verification |
| `design_comparison.png` | Permeability → required pressure → actual DWSIM power |

The original Julia/Panel network app remains a separate, independent example.
This milestone adds a command-line workflow and scientific figure; it does not
add an end-to-end UI or launch the DWSIM desktop GUI.
