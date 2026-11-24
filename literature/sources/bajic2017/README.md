# Bajić et al. (2017): microreactor assembly and nominal liquid inventory

[Process Biochemistry 52, 63–72](https://doi.org/10.1016/j.procbio.2016.09.021)
· [metadata](metadata.json) · [evidence](evidence.toml).
The user supplied the original PDF in the private inbox. The recorded
SHA-256 identifies it; the original article is not distributed in Git.

## Relevant experimental details

- Section 2.2 and Fig. 1 (pp. 64–65) show 15 mm thick PMMA plates,
  drilled ports in the upper plate, PEEK high-pressure fittings and PFA
  tubing of 0.50 mm internal and 1.59 mm external diameter.
- Sections 2.2 and 3.2 distinguish narrow rectangular channels from wider
  channels with triangular inlet and outlet regions. The latter have 4 mm
  pillars. The dye visualization in Fig. 3 uses the wider empty geometry
  at 1,000 µL/min; it does not supply a quantitative narrow-bed pulse RTD.
- Section 2.2 describes estimating channel volume with water filling;
  Section 2.5 and Table 1 report channel and inferred void volumes. All
  twelve Table 1 channel volumes agree with nominal rectangular
  length × width × depth within 0.041 µL. The wide-channel width refers
  to its rectangular section. The tabulated values therefore cannot
  establish a separate wetted fitting or distributor volume.

## Transfer boundary

This is a reactor assembly precedent, not a new specimen of Strniša's
or Panjan's narrow bead packing. Its immobilized preparation and particle
shape differ. The bore IDs and effective wet lengths of the ports, fitting
cavities and tubing lengths relative to the bypass are not reported.
[Strniša's own assembly description](../strnisa2018/plumbing_geometry.toml)
independently confirms the tubing ID. [Panjan's thesis](../panjan2018/plumbing_geometry.toml)
reports thicker PMMA plates (20 mm) but not a tubing ID, so using the
same 0.50 mm ID there remains a related-source scenario.

See the [plumbing audit](../../../examples/biocatalysis/bajic_plumbing_audit/REPORT.md)
for the independent circular-tube dispersion calculation and the
comparison against the retained tracer discrepancies. It does not
assume the actual wet bore diameter equals the tubing diameter or
transfer the wider channel's distributors to the narrow designs.
