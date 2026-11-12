# Stationary-wall stock TRT comparison

Test the archived stationary disk at N=128 and four axial layers using its
original 11,293 steps, time step, 0.4 mm radius, 12.9373 mm/s peak speed and
molecular diffusivity of 1e-9 m²/s. Each setup gets a uniform 0.5 kg/m³ run and
the original angular mode with initial amplitude 0.25 kg/m³. Compare all
native fields at the same coordinates and times to the saved, independently
converged Julia radial solution.

| Setup | Fluid collision | Wall reflection |
| --- | --- | --- |
| Archived comparison | Combined regularized TRT | Full-way bounce-back |
| Direct, full-way | Direct TRT | Full-way bounce-back |
| Regularized, half-way | Combined regularized TRT | Stock Bouzidi operator at q=0.5 |
| Direct, half-way | Direct TRT | Stock Bouzidi operator at q=0.5 |

Keep stock D3Q27, first-order equilibrium, odd relaxation for physical
diffusivity and TRT magic 0.3 in every comparison. The q=0.5 route uses the
existing project helper to set stock registered wall-link distances. This
moves the timing of reflection to the same step and retains the voxel
mid-link approximation to the wall. It does not use interpolated cut
distances. Keep the original OpenLB source intact. Use no clipping,
concentration correction, transport parameter fit or bed run.

Require the original screens: uniform maximum departure at most 1e-9 kg/m³,
profile RMS divided by initial amplitude at most 1% both globally and in the
outer 50 µm at mid and final time, engineering excursion beyond initial
bounds at most 0.02 kg/m³, and relative closed inventory drift at most 1e-8.
Also report the 1e-8 kg/m³ physical bound screen. Full-way inventory
includes fluid and reflected populations temporarily stored in active solid
links. Half-way inventory is fluid only; the reflected copies left in solid
cells are outside that balance. Read back native geometry, velocity,
initialization, link coverage, history extrema and axial invariance.

If an alternative passes these screens at N=128, repeat that alternative
at N=64 for a matched two-case resolution check. If none passes, archive
the comparison as a negative qualification. Hold the full bed for separate
authorization regardless of outcome.
