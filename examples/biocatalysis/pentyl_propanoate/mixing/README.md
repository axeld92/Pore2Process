# Packed-bed water-tracer check

The experiment was a pulse of blue dextran through the same tube used
for pentyl propanoate production, with water flowing at 1.2 mL/min.
Blue dextran is too large to sample the bead interiors in the
published interpretation. Sodium azide is a separate bead-access
experiment and has a different exit-time distribution. The comparison
uses the blue-dextran measurements in Figure A1a and the printed
moments in Table A1. [Published measurements](../../../../literature/sources/lorenzo2025/README.md)

## Flow and tracer models

- [flow.py](flow.py) builds a project-owned application against the
  repository's pinned OpenLB 1.9r0 archive. It reuses the existing
  periodic packing application with the present bead diameter and
  water properties. The OpenLB flow uses stock D3Q27 TRT, bounce-back
  bead walls and Kupershtokh forcing. No vendor files are modified.
- The seeded geometry is a periodic block of 128 nonoverlapping
  593 µm spheres. Its continuous void fraction is 0.64, taken from
  the paper's aqueous tracer correlation. The block is 8.53 mm long
  and about 2.13 mm wide. It represents interbead flow through a
  synthetic packing. The narrow periodic cross-section omits tube
  wall layers, end distributors and external plumbing. The reported
  6 g catalyst inventory is used separately for the reaction model.
- The superficial velocity is the reported water flow divided by
  the area of the 16 mm ID tube. The reference does not specify a
  tracer temperature. The simulation assumes water at 25 °C with
  viscosity 0.890 mPa s and density 997.05 kg/m³. [IAPWS liquid water properties](https://iapws.org/relguide/LiquidWater.pdf)
- Native axial throughput is matched on each mesh. The existing
  Julia projection closes the net flow through every block face at
  the native throughput and closes every connected control volume.
  Tiny diagonal-only pockets are counted in OpenLB and removed from
  the six-face tracer graph. Their volume and native axial flow are
  recorded.
- [tracer.jl](tracer.jl) releases markers in proportion to the
  incoming water flux and reuses the existing exact within-voxel
  advection and reflecting Brownian moves. Marker positions unwrap
  through successive copies of the same periodic packing until they
  first cross a plane 140 mm downstream. The repeated packing and
  its interbead flow are the assumptions being tested. Each reported
  case tracks 8,192 markers for up to 90 min.
- The tracer diffusivity is a Stokes–Einstein estimate of
  \(6.13\times10^{-12}\ \mathrm{m^2/s}\) using the paper's approximate 80 nm
  dextran diameter and the assumed water temperature. A second
  calculation doubles it. The pulse model assumes the 8.4 mL dose
  enters at the water flow rate over 7 min with uniform concentration.
  Both the ideal instantaneous result and the rectangular-pulse
  result are recorded. The experiment's actual injection time
  profile was unreported.

## Reproduce

Run from the project root with the environment established in the
project setup. The external build requires GCC 15, OpenMP and the
locally installed OpenLB 1.9r0 archive. The validated run used Julia
1.11.6, NumPy 2.5.4 and six threads per run.

```sh
source scripts/env.sh
uv run python -m examples.biocatalysis.pentyl_propanoate.mixing.flow
uv run python -m examples.biocatalysis.pentyl_propanoate.mixing.benchmark
julia --startup-file=no --project=julia examples/biocatalysis/pentyl_propanoate/mixing/test_tracer.jl
uv run python -m pytest -q examples/biocatalysis/pentyl_propanoate/mixing/test_benchmark.py
```

The compact [reference](reference.json) records the tested inputs and
results. Flow fields, temporary executables and marker samples remain
under the ignored project-local build and results directories. The
benchmark exports `tracer_curve.csv` with the prediction and Figure
A1a traces on the same measurement time axis.

The comparison uses a 5% mean and 10% variance screen as an
engineering check against printed moments. Mesh and time-step
comparisons apply those limits to pulse-free moments. These screens
are numerical criteria. Figure A1 and Table A1 report experimental
results without uncertainties on the tracer moments.

## Verified result

| Blue-dextran exit response | Mean (min) | Variance (min²) |
| --- | ---: | ---: |
| Synthetic 14 cm bed, 64-cell transverse mesh, 0.10 s marker step | 15.021 | 1.847 |
| Same simulation with the assumed 7 min injection pulse | 18.521 | 5.930 |
| Reported Table A1 moments | 19.070 | 27.830 |
| Normalized Figure A1a digitization | 19.091 | 27.780 |

The pulse simulation gives a mean within 2.9% of Table A1. Its exit
times have a standard deviation of 2.44 min compared with 5.28 min
in the experiment. The normalized cumulative tracer recovery
curves differ by up to 18.3 percentage points. The overall
water-tracer benchmark fails the variance check. At 1.2 mL/min the
digitized blue-dextran concentration area implies 6.72 mg out
against the nominal injected dose of 8.4 mg. The raw concentration
comparison therefore also includes a dose mismatch. [Digitization and experimental conditions](../../../../literature/sources/lorenzo2025/README.md)

The OpenLB flow converged on all three meshes and the native
throughput was within 0.11% of the target. The 64-cell transverse
mesh used 33.33 µm voxels, with a maximum Mach number of 0.00375.
Its face-rate projection changed the field by 1.38% RMS and kept
the throughput on every plane within \(2.0\times10^{-14}\) relative
error. All markers exited in each simulation. The variance changed
by 11.2% between the 32- and 48-cell meshes at matched time step.
Refining the 48-cell mesh to 64 cells changed it by 1.07% at 0.10 s.
Halving the fine-grid marker time step changed it by 0.82%.
Doubling the assumed molecular diffusivity on the 48-cell mesh
changed it by 2.31%.

The 0.64 void fraction was inferred from the published water trace.
Agreement in the mean provides a consistency check on liquid volume.
The large discrepancy in width remains after the tested numerical
changes. If the input pulse really lasted 7 min, any pulse shape
confined to that interval could contribute at most 12.25 min² of
additional variance to the calculated bed response. Even that bound
would give a variance of only 14.10 min². The experiment leaves the
actual pulse profile and the contributions of tube-scale packing
and plumbing unresolved.

The reaction example uses the pulse-free variance divided by squared
mean, 0.0081843, from `water_tracer.pulse_free_relative_variance` as
a mixing sensitivity of the reacting liquid. This transfer assumes
the same dimensionless spread across the liquid and flow rates.
