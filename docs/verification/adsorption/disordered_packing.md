# Matched ordered and disordered tracer

Cohort `20261112T161615Z-52964a01` used stock OpenLB 1.9r0
D3Q27 TRT and a full eight-section disordered sphere packing.
The matched N64 comparison reduced wall-free passage of
initially fast tracer from 93.27% to 76.40% and from 88.65%
to 69.59% at the 0.80 and 0.45 mm Péclet numbers. All base
cases recovered the full inert dose. The disordered N48 to
N64 changes exceeded the 2-point grid screen. This
establishes the direction in the tested packing and leaves
the precise magnitude grid-sensitive.

The N64 temporal screens passed after doubling the 0.45 mm
marker ensembles. Read-only replay verified 34 sources and
942 artifacts.
[Report and reproduction](../../../examples/adsorption/disordered_packing/REPORT.md)
