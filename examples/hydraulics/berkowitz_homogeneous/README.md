# Salt spreading through homogeneous glass beads

This benchmark predicts the width of an inert salt pulse after 200 mm of fine or coarse glass beads at 0.3 mL/min. The independent comparison uses the two dispersivities reported by Berkowitz and colleagues for homogeneous packings. The measured curves in their Figure 2 belong to their layered column. The [results](REPORT.md) and [experimental source record](../../../literature/sources/berkowitz2009/README.md) describe the distinction.

## Reproduce the runs

Run from the repository root. The tested environment was macOS 26.6.2 ARM64, CPython 3.12.9, NumPy 2.5.4, Julia 1.11.6 and OpenLB 1.9r0, compiled with Homebrew GCC 15.2.0 and OpenMP. The [native setup guide](../tool_checks/README.md) describes the required archive and compiler. The project builds the sphere case outside the vendor tree and verifies the source and executable hashes before use.

```sh
source scripts/env.sh
uv sync --locked
uv run python scripts/bootstrap_julia.py
uv run python -m examples.adsorption.breakthrough_hypotheses.hydrodynamics.native --build

case_module=examples.hydraulics.berkowitz_homogeneous.run
for grid in 6 8 10; do
  uv run python -m "$case_module" all fine --grid "$grid" --threads 4
done
for grid in 6 8 12 16 20; do
  uv run python -m "$case_module" all coarse --grid "$grid" --threads 4
done
for beads in fine coarse; do
  uv run python -m "$case_module" all "$beads" --seed 47 --threads 4
  uv run python -m "$case_module" tracer "$beads" --dt-factor 0.5 --threads 4
  uv run python -m "$case_module" tracer "$beads" --diffusion-factor 1.05 --threads 4
  uv run python -m "$case_module" all "$beads" --low-compressibility --threads 4
done
uv run python -m "$case_module" all fine --period 10.5 --threads 4
uv run python -m "$case_module" all coarse --period 14 --threads 4
uv run python -m examples.hydraulics.berkowitz_homogeneous.report
julia --threads=3 --startup-file=no --project=julia examples/hydraulics/berkowitz_homogeneous/test_transport.jl
uv run pytest -q examples/hydraulics/berkowitz_homogeneous/test_benchmark.py
```

The stage commands `flow`, `faces` and `tracer` can also be run individually. A completed carrier can supply several tracer checks. Each realization has a separate directory under `results/hydraulics/berkowitz_homogeneous/`. The report command writes a fresh `comparison.json` and `outlet_pulses.png` there. The compact [reference calculation](reference_summary.json) and [outlet pulse figure](outlet_pulses.png) are versioned. Stock OpenLB emits one construction frame per material update. The runner keeps these in temporary storage and retains the converged field, flow history and log. Allow several gigabytes of scratch space per active fine-bed run.

## Packed tube and flowing water

The internal diameter is 27 mm, with 1 mm or 4 mm glass beads. Water fills approximately 39% of the volume. The experimental superficial speed is 0.524 mm/min and the mean pore speed is 1.344 mm/min. The 200 mm section therefore holds 44.66 mL of water and has a nominal travel time of 148.86 min. Each synthetic geometry includes the entire tube cross-section. The baseline flow cell spans 5.25 fine-bead diameters or seven coarse-bead diameters along the axis. The longer realizations double these periods. Each bed has a second packing at the original period.

The generator seeds staggered cylindrical shells for the fine spheres. In the confined coarse tube it seeds helical shells and a central row. It perturbs the positions and relaxes overlapping spheres in stages. The independent geometry audit checks every pair with the axial periodic image and tests every sphere against the tube wall. Packing order and correlations are model assumptions. The paper reports bead size, packing procedure and approximate porosity. Individual bead positions and single-bed outlet traces are unavailable. [Berkowitz et al., section 2](https://doi.org/10.1029/2008WR007342)

The existing external sphere case uses stock sphere and cylinder indicators, D3Q27, TRT with relaxation time 0.85 and magic parameter 3/16, Kupershtokh body forcing and bounce-back on glass and tube walls. Flow is periodic only along the tube axis. It adjusts the axial force until the mean flow is within 0.2% of the prescribed rate and checks for steady velocity and energy. The records retain the actual flow and the saved velocity field.

The compiled sphere case uses a reference diameter of 0.8 mm. To reuse that executable at the experimental dimensions, all lengths are scaled by $s=d_*/d$. With the starred values in native units, velocities scale by $b=(\nu_*/\nu)/s$ and throughput by $s^2 b$. Both cases therefore have the experimental particle Reynolds number. Python divides the native velocity by $b$ and uses the physical spacing for transport. The reported tracer Péclet numbers use the physical bead diameter, pore speed and salt diffusivity.

The liquid properties are rounded pure-water values at 23 °C and approximately atmospheric pressure. The salt input is 500 mg/L or 8.56 mmol/L. Harned and Hildreth's Table II gives NaCl mutual diffusivities of 1.544 and 1.555 × 10⁻⁹ m²/s near this concentration at 25 °C. We used 1.545 × 10⁻⁹ m²/s and scaled it with absolute temperature divided by water viscosity to obtain 1.4653 × 10⁻⁹ m²/s at 23 °C. A separate run increases that diffusivity by 5% to cover dilution. Constant fluid density is an assumption of the passive-tracer calculation. [IAPWS liquid-water properties](https://iapws.org/technical-guidance/release/LiquidWater) · [Harned and Hildreth, Tables I–II](https://doi.org/10.1021/ja01146a043)

At eight cells per bead, a second TRT calculation lowers the relaxation time from 0.85 to 0.60. Each control starts from its own converged baseline velocity and the corresponding rescaled lattice-density perturbation. The native solver reconverges at the same physical flow. The tracer is then rerun using the new field. This checks sensitivity to the weak compressibility of the numerical water solver. Salt-induced density changes in the physical experiment are a separate question.

## Tracer and outlet measurement

The OpenLB output contains seven little-endian doubles per fluid voxel, ordered as zero-based axial and transverse indices, three native physical velocity components in m/s and native density in kg/m³. Python forms the shared face fluxes. The existing Julia pressure projection restores local conservation and preserves the achieved axial throughput. The record gives both the local correction and the fraction of pore space reachable only through diagonals. The particle tracker uses the connected six-face graph and records how much native axial flow belongs to excluded cells.

Julia reads three one-based neighbor IDs or signed face rates per voxel and side. Zero denotes a solid neighbor and rates are in s⁻¹ after division by the physical voxel volume. A marker is injected at the axial cut with probability proportional to forward liquid throughput. It follows the locally divergence-free velocity between faces and diffuses using independent Gaussian displacements with reflection at solid faces. The analytic subcell drift and face-crossing logic are shared with the existing Julia transport code. Axial paths retain their positions through periodic joins and can backtrack upstream of the injection cut. Travel is recorded at the first downstream crossing of the faces nearest 50, 100 and 200 mm. The files retain the exact distances and all arrival times in seconds.

For an ideal advection–dispersion pulse reaching a downstream plane, $\bar t=L/u$ and $\sigma_t^2=2D_L L/u^3$. The equivalent dispersivity is therefore $\alpha=L\sigma_t^2/(2\bar t^2)$. The calculation uses the marker mean and unbiased sample variance. Sixteen batches give a Monte Carlo standard error for the estimate. This quantifies sampling variability within a synthetic run. The authors give no uncertainty for their fitted single-packing dispersivities.

The physical injection lasted five minutes and delivered 0.75 mg of NaCl. The impulse calculation supplies bed dispersivity directly. For the collection plot, each marker arrival is convolved with the five-minute rectangular injection and integrated in five-minute bins. Integration over all bins must recover the injected mass. The dashed comparison is a one-dimensional prediction reconstructed from a published dispersivity and the nominal water inventory. Source measurements enter through the two independent scalar comparisons.
