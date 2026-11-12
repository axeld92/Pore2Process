# Balanced Julia face flows on the disk, 11 November 2026

The [report](../../../examples/adsorption/disk_face_flow/REPORT.md) records
seven Julia disk executions using the original prescribed shear, physical
diffusivity and archived radial reference. The uniform tracer stayed at
0.500 kg/m³ exactly. The 160×1024 working mesh passed all 1% varying-profile
limits, with a worst RMS error of 0.235% of the initial amplitude. The
320×2048 confirmation changed the sampled profile by at most 0.139%.

The initial cohort `20261111T153035Z-2fda74bc` and spatial confirmation
`20261111T153229Z-96c67149` passed read-only replay. The confirmation
qualifies the disk after the initial cohort's spatial-change screen missed
near the wall. The study uses prescribed face flows and generated zero
full-bed or process runs. [Run guide](../../../examples/adsorption/disk_face_flow/README.md)
