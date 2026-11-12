# Prospective flat-wall Poiseuille scalar control

10 November 2026. Test the existing stock scalar configuration at straight
impermeable walls, separately from the three-dimensional FCC geometry and the
previous no-solid periodic shear. No adsorption, new fitted coefficient or
process prediction follows from this diagnostic.

Keep Lx=H=0.0013891725543141453 m, mean speed U=0.012937295182166382 m/s,
Dm=1e-9 m2/s, dt=1.90174557964e-5 s, 11294 complete time steps (the actual
duration of the full-cube shear parent), and initial
C=.5+.25*sin(2*pi*x/Lx) kg/m3. Prescribe the analytical divergence-free velocity
u=(6*U*(y/H)*(1-y/H),0,0) at fluid-cell centres. This is prescribed transport
in a known velocity field, not a new Navier–Stokes validation or a fitted flow.
Keep stock D3Q27, CombinedAdvectionDiffusionRLBdynamics wrapping stock TRT,
first-order equilibrium and full-way BounceBack. No vendor changes, clipping,
artificial diffusion, custom collision, descriptor or transport update.

Use N fluid cells in x/y, dx=H/N, periodic x/z, and physical walls at y=0,H.
Fluid centres are y=(j+.5)*dx, j=0..N-1. The two solid node rows are at
y=-dx/2,H+dx/2, so the walls are at the fluid/solid midlinks. OpenLB's maintainer
explicitly identifies both its full-way and halfway bounce-back wall positions
as midlink: https://www.openlb.net/forum/topic/bounceback-boundary/ . Full-way
reflection adds numerical storage and a two-step return; do not silently
count the complete solid rows as physical liquid or fit an effective H.

For nonperiodic y, record physical liquid mass plus the actual wall populations
that arrived from fluid: at the lower wall only directions with c_y=-1, at the
upper wall c_y=+1. Restore the OpenLB quadrature weights before summing these
stored populations. These are exactly the populations that reflect into fluid
on the following step. Other solid populations face the unused exterior and
are excluded from the closed liquid-plus-active-wall inventory. Check this
kinetic inventory <=1e-8 relative; report liquid inventory and active wall
storage separately. This is diagnostic accounting, not a replacement flux.

First compare N32 with 32 versus 4 periodic z layers. Require initial/final
and selected peak x-y fields and every-step extrema agreement <=1e-11 kg/m3;
equal peak step; physical volume and inventory scaling by 4/32; final z
variation <=1e-12 kg/m3. Thin z changes depth, not the invariant x-y solution.
Only after that equivalence passes run N64 and N96 with four z layers.
Also run an initially constant .5 case at N32 and a diffusion-only N32 control
whose independent solution is the decaying x sine wave with no y dependence.

The independent continuum oracle uses one Fourier mode in x and Neumann
diffusion in y. Even extension over a period 2H enforces zero normal derivative
at y=0,H. NumPy FFT diffusion and exact pointwise Poiseuille phase evolution
with Strang splitting provide a verification oracle, never a second pore
production solver. Independently refine transverse reference points and time
steps; require each sampled maximum change/amplitude <=1e-4. Verify reflection
symmetry and the zero-flow analytical solution before assessing native fields.
If the oracle does not qualify, label native profile accuracy unavailable.

Retain every-step extrema, liquid/active-wall inventories, x-y initial/final/
largest-excursion fields, physical coordinates/velocity, native time, sampled
maximum lattice speed, units, input snapshots, source/build/artifact hashes and
the pinned OpenLB distribution integrity check. Screen sinusoidal ranges
[.25,.75] strictly and [.23,.77] as the engineering allowance; profile RMS/
initial amplitude <=1%. Constant and diffusion checks remain separately
reported. Passing a flat-wall diagnostic cannot qualify curved FCC transport,
an open-boundary tracer, long-time dispersion or adsorption service.
