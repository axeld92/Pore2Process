import numpy as np

from examples.adsorption.disordered_section_exchange.run import errors


def test_variance_gate_covers_all_held_out_sections_and_both_error_signs():
    observed = np.column_stack((np.arange(3.,9.),np.full(6,10.)))
    predicted = np.column_stack((np.ones(6),observed))
    predicted[2,2] = 8.9
    predicted[5,2] = 10.8
    result = errors(observed,predicted)
    assert result["mean_passed"]
    assert not result["variance_passed"]
    assert result["variance_peak_section"] == 5
    assert result["variance_relative_error"][-1] > 0.
