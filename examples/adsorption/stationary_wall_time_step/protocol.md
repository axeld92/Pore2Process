# Stationary circular wall at half the time step

Use the archived regularized TRT half-way cases in
`20261111T133159Z-7f69a322` as the controls and run their exact OpenLB 1.9r0
binary again. Retain the D3Q27 lattice, wall links at q=0.5, 128 cells across
the enclosing square, four axial layers, prescribed physical velocity, and
initial concentrations. Keep molecular diffusivity at 1e-9 m²/s and TRT
magic at 0.3. Halve dt and double the step count. The stock application's
relaxation parameters and lattice velocity follow the new time scaling.

The archived binary saves a snapshot at floor(steps/2). The original 11,293
steps are odd, so the halved-step run's automatic midpoint is one refined
step later than the archived Julia midpoint. A truncated repeat of the
varying case ends at exactly the Julia midpoint. Check that its full
history agrees with the full-duration run through that step and use its
final field for the midpoint comparison.

Compare the uniform run to 0.5 kg/m³ throughout and the varying runs to the
same archived Julia radial mode at 0.10737255542652152 s and
0.21476412830883945 s. Keep the established uniform departure (1e-9 kg/m³),
global and outer-50-µm RMS (1% of the 0.25 kg/m³ initial amplitude at both
times), physical bounds (1e-8 kg/m³, reported separately), engineering
bounds (0.02 kg/m³ excursion) and closed liquid inventory (1e-8 relative)
checks. Audit native inputs, histories, wall links and reconstructed
diffusivity. This is a two-level time-step comparison on a single disk.
