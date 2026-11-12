# Curated engineering accuracy evidence

The JSON records retain numerical outcomes and exact source/native/assessment
summary hashes. Full native fields, histories, input snapshots and builds stay
under ignored `results/adsorption/transport_accuracy/`.

- `screen_reference.json`: four actual native baseline cases; uniform and sharp
  front results are separate.
- `smooth_input_reference.json`: four actual native runs with two known input
  durations; input mean and variance are audited and removed from the response.
- `regularized_inlet_reference.json`: two actual short stock-inlet diagnoses;
  these short windows cannot qualify transport tails.
- `profile_reference.json`: four actual periodic smooth-profile controls;
  free-flow analytical comparisons and porous concentration bounds are separate.
- `diagnosis_reference.json`: four subsequent bounded controls using the exact
  preserved profile binary; stationary water, time-step and matched-grid results,
  native/scalar field audits and the unchanged unqualified service state.

Producing a reference by `check` only verifies saved evidence and replays Julia
analysis. It does not execute OpenLB or DWSIM and does not promote a failed model.
