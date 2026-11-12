# Initial normal-compatibility control

Declared before execution on 10 November 2026. The original sine initial
concentration has a nonzero normal derivative at most curved impermeable
particle surfaces. Its weak continuum solution must still respect the maximum
principle, but an unresolved initial boundary adjustment could excite numerical
errors. This separate diagnostic changes only the initial concentration. It
does not qualify or repair the original case, change molecular diffusivity,
change velocity, fit parameters or replace a boundary/collision operator.

Use the same qualified d800 N32 native flow (optionally N64 afterward), physical
D=1e-9 m²/s, dt=1.9017455796408345e-5 s, end request
0.2147547125969688 s, periodic domain, stock D3Q27 CombinedRLB/TRT and stock
full-way BounceBack. The application is copied from the frozen instrumented
scalar case; only initial-data construction and corresponding metadata change.

For the nearest existing stock FCC sphere, let d be its stock signed distance
and s the nearest radial surface point. Reuse `IndicatorSphere3D` signedDistance,
getCenter/getRadius and stock `olb::normalize`; do not implement a new indicator,
lattice, collision or boundary. Define C0(x)=.5+.25 sin(2 pi x/L). Set

    C(x) = (1-w) C0(s_x) + w C0(x)
    w = stock PolynomialStartScale<T,T>(1,1)(t)
      = 6 t^5 - 15 t^4 + 10 t^3,
    t = clamp(d/delta, 0, 1).

The nearest sphere image is found using the same stock indicators as geometry.
Coordinates are periodically wrapped before evaluation. The fixed physical
buffer is delta=0.4*(L/sqrt(2)-diameter), approximately 72.9 micrometres. It is
less than half the minimum sphere gap, so buffers do not overlap. It is held
fixed in metres if N64 is tested. Inside the sphere use the same normal-constant
extension to initialize active wall populations. The radial projection at an
exact sphere centre is undefined; only inactive deep-solid data would need the
original C0 value there (even cell-centred N32/N64 never samples a centre).

The profile preserves the original wall trace and tangential variation, has
zero wall-normal derivative at the surface, joins C0 with two continuous
derivatives at delta, and is a convex combination within [.25,.75] kg/m³.
Independent Python verification uses four FCC basis centres and minimum-image
coordinates, rather than the native loop of sphere images. Before interpreting
the simulation it checks native sampled values, periodicity, projection/range,
normal finite differences and the transition derivative. Geometry IDs and
velocities must still match the real native source exactly.

Report initial physical-liquid range/mean/SD and full-run extrema, conservation,
peak location and final moments. Compare with the unchanged sine case at the
same grid and dt. A smaller excursion is evidence about sensitivity to initial
compatibility, not evidence that a different initial condition fixes the
original profile or qualifies open-boundary breakthrough/column dispersion.
