# N96 disordered-packing tracer refinement

Refine the exact archived N64 packing to 96 voxels per reference
length in a 768 × 96 × 96 periodic OpenLB domain. Use the saved
32 sphere centres, 0.8 mm sphere diameter, carbon inventory and
original forward-flowing tracer cut. Reuse the stock OpenLB
1.9r0 D3Q27 TRT, Kupershtokh force, BounceBack and sphere
indicators from the disordered N64 study. Interpolation of the
converged N64 field supplies only the initial condition. The
N96 lattice must pass the original convergence, flow-match,
Mach, density and independently rasterized geometry checks.

Preserve the existing conservative least-squares face projection
and its local balance and throughput gates. Match the two Péclet
numbers to the archived measured interstitial-velocity basis,
using the converged N96 velocity. Keep the same perfect-sink
eight-section Julia tracker and the same 8,192 inlet seeds for
the direct comparison. Run both particle-matched conditions at
80 and 40 µs. Require at least 99.9% inert recovery.

Evaluate the maximum absolute N64-to-N96 difference over
sections 3–8 in total wall-free survival, conditional capture
among incoming survivors, and wall-free survival among initially
fast tracer. The grid screen is 0.02 of the relevant population
for each metric. Check both matched time steps and retain the
80 µs comparison as primary. Also check 80 versus 40 µs on N96
against the previous 0.01 time-step screen. If the sample-based
comparison is borderline, retain it and record any supplemental
marker ensemble separately.

The 7.032455155 g dry carbon inventory remains the common
column-scale basis. The tracer checks wall contact only. No
adsorption calculation or Figure 7 evaluation is part of this
study.
