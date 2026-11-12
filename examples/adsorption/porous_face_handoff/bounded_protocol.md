# Bounded transport follow-up on the same balanced faces

The centered cases on N32 and N64 produced negative tracer concentration.
Keep the archived centered cases and the identical projected carrier faces.
Test one conservative flux-corrected alternative on the N64 and N96 saved
geometries, keeping physical diffusivity, initial states and times fixed.
Limit antidiffusive face increments using each cell's current local range
and its allowable positive and negative increments. Use a donor-cell
low-order update with physical diffusion, a fourth-order centered face
value where the full axial stencil is liquid, and two-point centering at
other faces. Advance with three strong-stability-preserving Runge–Kutta
stages. Do not clip concentrations or increase molecular diffusivity.

Check boundedness and mass at each step and run the uniform control on
N96. Check the implementation against the exact straight-flow sine and
repeat the N64 varying case with half the time-step cap. Retain the same
0.1% time-step RMS and 5% spatial variance/amplitude plus 0.05 rad phase
screens. Record the physical diffusion dissipation and the net difference
between it and observed variance loss. Flag any contribution above 15%
of observed variance loss as unresolved numerical spreading. The variance
difference is a discretization diagnostic and does not give a continuum
dispersion coefficient. Do not increase the geometry size beyond one
periodic cell.
