# Short inlet–outlet inert tracer benchmark

Open the axial seam of one 1.389172554 mm N96-derived periodic FCC
cell. Retain the same saved TRT-derived balanced N192 and N288 face
flows and solid masks. Keep transverse periodicity and impermeable
solid faces. Replace each positive seam flow with equal positive
carrier flow entering the inlet cell and leaving the outlet cell.
Require all cut faces to be forward-flowing and preserve the archived
throughput, local continuity and liquid volume.

Initialize the pore liquid with zero tracer. Prescribe incoming total
tracer flux equal to Q_face times 1 kg/m³ for a 4 ms rectangular pulse,
then zero. At the outlet use Q_face times the cell concentration with
zero diffusive flux. There is no adsorption or reaction. Use molecular
diffusivity 1e-9 m²/s. Use the existing bounded interior flux corrections
and three-stage SSP time integration with explicit conservative boundary
increments. Check the open update on a constant field and an independent
discrete donor-cell pulse solution before applying it to the porous case.

Use a common 80 µs time step on N192 and N288, below both tested
stability caps. An N96 scouting calculation at 200 µs is used solely
to choose an adequate observation window. Define one transit as liquid
volume divided by archived throughput. Use one-transit checkpoints,
each rounded up to a whole common time step. Start with at least four
transits, allow up to ten in the scout, and use one extra transit beyond
the first scout checkpoint that passes the pulse recovery and tail
screens (with a minimum of five). If no checkpoint passes, use the
ten-transit window and report any remaining tail limitation. Both
comparison meshes use exactly the same physical end time. Record
outlet mass and raw first/second exit moments as they are integrated.
No extrapolated mass or moments are used.

Define early arrival as the first time when 5% of the injected mass has
exited. Also report 10% arrival and mass recovery. Compare exit time
measured from pulse onset and variance of the finite pulse, with mean
difference at most 5% and variance difference at most 10% relative to
N288. Require whole-history mass balance within 1e-8 of injected tracer,
concentrations within 0–1 kg/m³ to 1e-8, at least 99.9% recovered, and
at most 1% change in each recovered first and second raw moment over the
last nominal transit. Report incomplete moments and bound the claim if
these pulse-tail checks fail. This is a one-cell transport test with
periodic developed carrier flow. Do not fit dispersion or run a full bed.
