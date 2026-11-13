# Project instructions

Preserve completed examples and unrelated files.

## Computational scope

- This is a computational-tool project. Evidence comes from published papers
  and their available supporting information; no new laboratory measurements
  or unavailable raw datasets are prerequisites.
- Use synthetic or literature-informed geometries. Verification against
  analytical solutions, numerical convergence and literature-based comparisons
  are required as appropriate; reproducing a specific experimental specimen's
  connected 3D pore geometry is outside the current acceptance scope.
- Keep unreported experimental inputs unavailable or make clearly labelled
  scenario assumptions with sensitivity checks. Do not invent uncertainty or
  treat obtaining raw data, extra repositories or author contact as a gate.

## Scientific evidence

- Seek experimental inputs in scientific articles and retain source, specimen,
  figure/table, conditions, units and uncertainty with each extraction.
- The user can provide article Markdown and supplements in `literature/inbox/`
  using ZHAW/DTU access. Identify the exact paper and missing data when needed.
- Keep illustrative assumptions distinct from extracted measurements. Do not infer
  separate layer resistances by subtracting measurements of unrelated specimens.

## OpenLB

- Porous-media simulations must use at least TRT collision dynamics. Reuse
  stock OpenLB TRT or a suitable more advanced model; retain older BGK runs
  as historical evidence and revalidate new results after a dynamics change.
- Reuse OpenLB's dynamics, descriptors/velocity sets, geometry indicators,
  functors, data structures and other existing facilities. Do not recreate them.
- Follow the structure and lifecycle of the official examples: mesh creation,
  geometry preparation, lattice preparation, initial/boundary values, simulation
  and results. Use the installed version's examples as the API reference.
- Never modify original OpenLB files. Keep project simulation code, headers,
  configuration, build products and results in separate project folders.
- Revalidate numerical results after changing a simulation case or discretization.

## Version control

- Use `main` for local work and GitHub pushes unless the user asks for a separate branch.
- Commit substantial project code changes to Git and push the commits to the
  configured GitHub repository.
- Keep source, tests, protocols and compact verification records traceable.
  Leave unrelated workspace files and bulk generated run data out of commits.
