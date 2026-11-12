# Periodic shear-flow isolation: frozen protocol

This control isolates bulk scalar transport in a spatially varying velocity
from solid-wall effects. It is an all-fluid periodic box, not a physical
wall-bounded Poiseuille-channel benchmark or a substitute adsorption model.
No result qualifies the FCC tracer closure or a process prediction.

The cube has side H=0.0013891725543141453 m, matching the existing d800 FCC
repeat. The fixed mean velocity U=0.012937295182166382 m/s is its coarse-grid
interstitial velocity. Prescribe the exactly divergence-free field
u=(U[1+a cos(2 pi y/H)],0,0). Set a=0.75 for shear and zero for the uniform
control. All three directions are periodic. No solids, inlet or outlet are
present. Diffusivity D=1e-9 m2/s retains the application's physical regime.
This cosine shear has smooth gradients and cannot reproduce every gradient
or constriction in a real FCC velocity field.

Set C(x,y,z,0)=0.5+0.25 sin(2 pi x/H) kg/m3. The continuum maximum principle
keeps C in [0.25,0.75] for all time; its periodic integral remains constant.
A separate constant C=0.5 shear control checks exact preservation.

Use dt=1.90174557964e-5 s and requested end=0.214764128309 s, with ceil(end/dt)
steps and actual time explicitly recorded. Compare N=32,48,64 shear boxes
and N32 uniform/constant controls. The OpenLB descriptor is stock D3Q27 and
the collision is stock CombinedAdvectionDiffusionRLBdynamics composed with
AdvectionDiffusionTRTdynamics, exactly the audited FCC bulk model. Set the
stock TRT magic parameter to 0.3 and odd relaxation time to
0.5+3 D dt/dx2; derive even relaxation consistently. No operator, boundary,
limiter, artificial diffusion or dynamics is implemented in this application.

The independent continuum oracle retains the initial single x Fourier mode:
C=0.5+0.25 Im[exp(i k x) f(y,t)], k=2 pi/H,
f_t=D(f_yy-k2 f)-i k u(y) f, f(y,0)=1. A NumPy Fourier Strang calculation
uses exact diffusion/pointwise-advection substeps. Check Fourier resolution
128→256 at 1024 substeps, then time 1024→2048 at 256 points; both changes
must be <1e-4 times initial amplitude at native sampled coordinates. The
uniform case additionally has the exact translating, exponentially decaying
sine solution, independent of the split calculation.

Native checks: finite values, maximum global inventory drift <=1e-8,
constant error <=1e-9 kg/m3, z invariance <=1e-12 kg/m3, and continuum bounds
with a separately reported 1e-8 strict tolerance and existing 0.02 kg/m3
engineering allowance. Compare RMS and maximum native error with the converged
continuum oracle and inspect the trend; do not claim formal order from an
unresolved or nonmonotone trend. RMS/initial-amplitude <=1% is the predeclared
profile comparison. Peak snapshots contain a complete x-y cross-section
and explicit z-invariance diagnostics; all-step extrema retain timing.

Preserve source snapshots, executable hash, exact command, inventory history,
initial/peak/final cross-sections, version/build metadata and artifact hashes.
No new experimental input is needed. Original OpenLB files stay unchanged.
