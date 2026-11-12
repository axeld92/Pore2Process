# Qualified-reference continuation

The original native cohort `20261110T212641Z-c0b585ff` completed full32,
thin32, constant32, diffusion32 and thin64. Its native solver and data remain
unchanged. It stopped while assessing thin64 because the independent continuum
reference's spatial refinement change was 1.1232e-4 of amplitude, exceeding the
unchanged 1e-4 tolerance. That failure concerns oracle resolution, not a failed
native execution. Do not declare a native profile pass from that assessment.

Freeze the original runner/reference. In separate reference/assessment files,
increase the independent transverse reference refinement from 1024/2048 to
2048/4096 points on physical H (4096/8192 on the even period); retain separate
2048/4096 time refinements and the same Neumann symmetry/derivative checks.
Reassess the preserved five native cases without rerunning them. Hash every
parent artifact and original source snapshot before this continuation.

Use the unchanged, hash-verified original native executable for only the missing
thin96 case with the original protocol's parameters. Retain actual commands,
new source snapshots, the parent artifact manifest and distinct native counts:
five completed parent cases, one new completed case, and two earlier interrupted
native attempts whose exit status is unavailable. A successful analysis must
state six completed cases and two interrupted attempts, not eight successful
cases. The two interrupted folders are preserved separately and have no live
orphan process according to the recorded narrow ps/lsof check.

The original physical wall convention, active-wall inventory, acceptance
screens and prohibition of application/dispersion qualification remain intact.
