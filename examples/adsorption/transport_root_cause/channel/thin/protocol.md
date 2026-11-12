# Thin periodic slab: frozen follow-up protocol

This follows the full-cube periodic-shear study. Its fields and velocity are
exactly independent of z, and the native full-cube output checks that
invariance. Replacing N periodic z layers with four periodic z layers
therefore changes computational depth/inventory, but does not change the
continuum x-y solution. This reduction is only an all-fluid diagnostic;
it does not represent or reduce a three-dimensional FCC bed.

Keep the parent H, mean U, D=1e-9 m2/s, shear fraction 0.75 and sinusoidal
initial condition. The stock regularized TRT/D3Q27 collision, streaming,
unit conversion, geometry coordinates and parameters stay identical.
Four layers exceed the stencil radius and retain ordinary stock periodicity.
Peak snapshots are buffered and written once; this affects diagnostic I/O
only. No production numerical operator is replaced or changed.

Run N32 at the parent physical dt=1.90174557964e-5 s and actual parent duration
0.21478314576454161 s (11294 steps). Before interpreting higher resolution,
require thin32 versus full32 agreement <=1e-11 kg/m3 in initial/final/peak
concentration cross-sections and every-step minimum/maximum, equal peak step,
matching coordinates and velocity, z variation <=1e-12 kg/m3, and normalized
inventory histories agreeing <=1e-11. Their absolute inventories scale with
the depth ratio 4/32; no fluid volume is silently reused.

Then run N96 and N128 at the same physical dt/duration and N128 at half dt,
with22588 steps so all actual observation times match exactly. The native
runner rounds end/dt only after requiring the ratio within1e-6 of an integer;
this avoids a one-step float-rounding difference. Record actual sampled
maximum lattice speed, physical depth, both relaxation times and native time.
The continuous maximum lattice-speed bound must remain <=0.05. N128 has
approximately0.0397 at the parent dt and0.0198 at half dt.

Compare each against the same independently spatially/time-converged Fourier
reference, with the parent strict/engineering bounds, mass and1% RMS checks.
Compare N128 full/half-dt final fields directly; call time sensitivity small
for this profile only if their RMS difference/initial amplitude <=1%. Do
not infer a formal convergence order from nonmonotone or unresolved results.
Even a successful all-fluid case does not qualify the FCC tracer or a
measured column dispersion coefficient.

Keep the parent native summary hash, source snapshots, executable hash,
commands, all-step extrema/inventories, initial/peak/final cross-sections,
independent reference-convergence results, build metadata and artifact hashes.
No original OpenLB source is modified and no new experimental data are used.
