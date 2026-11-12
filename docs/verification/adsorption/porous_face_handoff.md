# OpenLB-to-Julia face-flow transfer, 11 November 2026

The [report](../../../examples/adsorption/porous_face_handoff/REPORT.md)
records a verified throughput-preserving reconstruction of saved TRT
fields on one periodic FCC cell. All three meshes preserved uniform
tracer. The N96 face corrections were 0.362% RMS and the axial flow
matched the native volume throughput within 7.34e-14 relative.

The centered transport produced negative concentrations. The conservative
bounded alternative respected concentration bounds, with 17.1% final
variance sensitivity between N64 and N96 and an estimated 31.8% net
discretization contribution to N96 variance loss. Tracer spreading remains
unqualified. The study generated zero full-bed and process runs.
[Run and replay guide](../../../examples/adsorption/porous_face_handoff/README.md)
