# Does a salt pulse sink in a horizontal glass-bead bed?

The paired simulations send the same 0.75 mg NaCl dose through the same synthetic 4 mm bead packing at 0.3 mL/min. The variable under test is the excess weight of the dissolved salt. The outlet comparison uses the homogeneous coarse bed in section 2 of Berkowitz et al. The [scientific report](REPORT.md) gives the outcome. [Experimental inputs](../../../literature/sources/berkowitz2009/evidence.toml) · [Berkowitz et al. (2009)](https://doi.org/10.1029/2008WR007342)

## Reproduce

Run from the project root on the macOS ARM64 environment used for the previous [homogeneous-bed case](../berkowitz_homogeneous/README.md). This uses its pinned OpenLB 1.9r0 installation and the same Julia and Python environments. The new external native example consumes the earlier field and a voxel-for-voxel copy of its pore mask. The original OpenLB source stays in the vendor directory.

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/bootstrap_julia.py
uv run python -m examples.adsorption.breakthrough_hypotheses.hydrodynamics.native --build
uv run python -m examples.hydraulics.berkowitz_buoyancy.build
```

On a fresh checkout, generate the carrier and conservative faces for these runs. These commands reuse seed 31 and the 28 mm packing period. Existing verified fields in the appropriate results directories can be reused.

```sh
passive=examples.hydraulics.berkowitz_homogeneous.run
for n in 6 8 12; do
    uv run python -m "$passive" flow coarse --grid "$n"
    uv run python -m "$passive" faces coarse --grid "$n"
    if [ "$n" != 12 ]; then uv run python -m "$passive" tracer coarse --grid "$n"; fi
done
for n in 8 12; do
    uv run python -m "$passive" flow coarse --grid "$n" --low-compressibility
    uv run python -m "$passive" faces coarse --grid "$n" --low-compressibility
    uv run python -m "$passive" tracer coarse --grid "$n" --low-compressibility
done
```

Run the primary pair, the time-step check and the coarse-grid controls. The number after N is the number of water-grid cells per bead diameter. The same stock TRT model is used throughout. Lowering its relaxation time to 0.60 reduces artificial compressibility. The unbounded cases reproduce the initial particle-density interpolation. The bounded cases enforce the physical upper inlet concentration while conserving mass.

```sh
run=examples.hydraulics.berkowitz_buoyancy.run
uv run python -m "$run" N12_bounded --grid 12 --low-mach --sweeps 96 --threads 5
uv run python -m "$run" N12_off --grid 12 --low-mach --sweeps 96 --gravity 0 --threads 3
uv run python -m "$run" N8_bounded --grid 8 --low-mach --sweeps 96 --threads 4
uv run python -m "$run" N8_kernel --grid 8 --low-mach --sweeps 96 --kernel 0.0015 --threads 3
uv run python -m "$run" N8_half_off --grid 8 --low-mach --interval 30 --sweeps 48 --gravity 0 --threads 2
uv run python -m "$run" N8_half --grid 8 --low-mach --interval 30 --sweeps 48 --threads 3
uv run python -m "$run" N8_low_off --grid 8 --low-mach --sweeps 96 --gravity 0 --threads 2
uv run python -m "$run" N8_low_down --grid 8 --low-mach --sweeps 96 --unbounded --threads 4
uv run python -m "$run" N8_down --grid 8 --no-low-mach --sweeps 24 --unbounded --threads 4
uv run python -m "$run" N8_off --grid 8 --no-low-mach --sweeps 24 --gravity 0 --threads 2
uv run python -m "$run" N6_down --grid 6 --no-low-mach --sweeps 24 --unbounded --threads 4
uv run python -m "$run" N6_off --grid 6 --no-low-mach --sweeps 24 --gravity 0 --threads 2
uv run python -m "$run" N6_up --grid 6 --no-low-mach --sweeps 24 --gravity -1 --unbounded --threads 2
```

Check the water response separately by applying the saved salt distribution from an early stage of injection/flush and again at one hour. Each independent native solve starts from the archived fresh-water carrier. Collect the results after the cases complete.

```sh
audit=examples.hydraulics.berkowitz_buoyancy.audit
uv run python -m examples.hydraulics.berkowitz_buoyancy.hydrostatic
uv run python -m "$audit" N8_low_down early --threads 3
uv run python -m "$audit" N12_bounded early
uv run python -m "$audit" N12_bounded hour
julia --startup-file=no --project=julia examples/hydraulics/berkowitz_buoyancy/test_transport.jl
uv run pytest -q examples/hydraulics/berkowitz_buoyancy examples/hydraulics/berkowitz_homogeneous/test_benchmark.py
uv run python -m examples.hydraulics.berkowitz_buoyancy.report
```

The results are under `results/hydraulics/berkowitz_buoyancy/`. The checked-in [compact numerical record](reference_summary.json) and [figure](pulse.png) preserve the comparison. To regenerate them explicitly, pass `--output examples/hydraulics/berkowitz_buoyancy/reference_summary.json --plot examples/hydraulics/berkowitz_buoyancy/pulse.png` to the report command.

## Data and numerical method

- **Geometry.** The tube has an inside diameter of 27 mm. Its axis is x, positive z is upward, and z = 0 bisects the tube in the reported heights. Nine identical 28 mm packing periods form a continuous 252 mm hydraulic domain, periodic along x. The pulse enters at x = 0 and its first passage is recorded at 50, 100 and 200 mm. Markers leave at 200 mm. A 52 mm fresh-water guard separates the outlet from the next periodic inlet image. The program checks the furthest upstream marker excursion against that guard.
- **Fluid motion.** The external case uses stock OpenLB D3Q27, TRT, a magic parameter of 3/16, Kupershtokh forcing and bounce-back on the stationary glass and tube. The physical axial pressure gradient comes from the saved pure-water solution. The salt contributes only excess weight relative to fresh water. For local feed fraction `f`, the additional vertical acceleration is `−g f Δρ_feed / ρ_water`. A pressure gradient can balance a uniform density increment in a closed vertical direction. Spatial differences in density give rise to the circulation tested here. Uptake, particle settling and reaction do not enter these salt equations.
- **Two-way coupling.** Julia advances a single finite injection of markers over five minutes, followed by a fresh-water flush. Each marker carries the same salt mass. The water advects the markers on balanced shared pore faces. Molecular diffusion uses independent Gaussian steps reflected at solid faces. At each coupling midpoint, a trilinear kernel of nominal width 2 mm reconstructs the local salt concentration from the marker mass and connected liquid volume. In the bounded runs, any oversupply is shared among nearby liquid bins with room below feed concentration. The native solver receives this spatial field and returns a new three-dimensional water velocity. The lattice-density diagnostic tracks the pressure variable of the water solver. The independently specified salt-density increment enters its body force.
- **Time and throughput.** Water inertia relaxes on a faster scale than the salt pulse. The calculation takes a finite number of native sweeps per frozen salt field, then projects the velocity change onto shared faces that conserve liquid volume. A multiple of the original pressure-driven flow supplies the small remaining pump correction at fixed discharge. The live water field is compared with longer frozen-plume solves. The salt clock advances in physical seconds. Native sweeps are a multirate approximation to the evolving water response. The coupled first-passage result includes the actual injection clock.
- **Density law.** The feed is 0.500 kg NaCl/m³ and water density is 997.54 kg/m³. At this dilution, the salt's partial molar volume gives `Δρ = c (1 − ρ_water V̄ / M)`. The 25 °C limiting molar volume of 16.40 cm³/mol is used at 23 °C, giving 0.3600 kg/m³ excess density at the inlet. The molecular salt diffusivity is 1.4653 × 10⁻⁹ m²/s. [Harned and Hildreth, Table I](https://doi.org/10.1021/ja01146a043) · [Water and transport calculation](../berkowitz_homogeneous/protocol.toml)

### Outputs

- `inputs.json` stores physical and computational settings and input hashes. Each case has its own native and Julia logs.
- `entries.f64` and `outlet_z.f64` contain little-endian Float64 seconds and metres, one per marker. Each row in `arrivals.f64` contains three little-endian Float64 first-arrival times in seconds at the three downstream planes, measured from the start of injection.
- `history.csv` gives the midpoint salt height and split between halves, mass accounting, native water flow and pressure-density diagnostics, and conservation of the face velocities. The active marker count excludes salt already collected at the outlet.
- `salt_{early,hour}.bin` holds local feed fractions. The corresponding water files hold physical velocity components in m/s, ordered by repeat and by pore row in the archived field. The sample nearest the target coupling boundary is saved at its preceding midpoint (e.g. 570 s and 3570 s for 60 s updates). The large transient fields stay under the ignored results directory.
