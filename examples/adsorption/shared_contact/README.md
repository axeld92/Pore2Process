# Wall-weighted transfer into shared carbon

The [protocol](protocol.md) fixes the Péclet-matched exchange matrices,
particle properties, dry carbon and film conductance before scoring.
The [Julia solver](../../../julia/adsorption/SharedContactColumn.jl)
uses the original radial diffusion operator and allows each flow
class to exchange solute with the shared inventory in its axial slice.
[Results and verification](REPORT.md)

From the project root:

```sh
source scripts/env.sh
export OPENBLAS_NUM_THREADS=1
julia --startup-file=no --project=julia julia/adsorption/test_shared_contact_column.jl
uv run python -m examples.adsorption.shared_contact.run
uv run python -m examples.adsorption.shared_contact.run --check <summary-path>
uv run python -m examples.adsorption.shared_contact.plot --summary <summary-path>
```

Use the summary path printed by the forward run. Each result archive
retains full and any separately resolved short histories. The
replay reconstructs the archived weights and matrices, validates
the histories' mass and interfacial checks, and reproduces the
numerical and Figure 7 comparisons.

The model stores one liquid concentration per path and axial slice.
Paths are ordered by transit rank, and the transition matrix rows
index receiving paths. Particle loading is kg phenol per kg dry
carbon. All local quantities in the solver use SI units. Gross film
uptake and release are integrated separately and close against the
total adsorbed inventory.

The separate plateau check uses the full-run summary and the same
frozen physical parameters:

```sh
uv run python -m examples.adsorption.shared_contact.plateau <summary-path>
uv run python -m examples.adsorption.shared_contact.plateau <summary-path> --check
```
