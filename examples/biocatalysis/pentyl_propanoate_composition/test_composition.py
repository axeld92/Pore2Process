import numpy as np
import pytest

from examples.biocatalysis.pentyl_propanoate_composition import run as r


def test_complete_source_matrix_and_independent_equilibrium():
    plan=r.read_json(r.PLAN)
    groups=r.read_batch(plan)
    assert len(groups)==16 and sum(map(len,groups.values()))==225
    assert [len(groups[a]) for a in (13,14,15,16)]==[18,13,12,7]
    assert groups[15][-1]["printed"]==.435
    assert groups[16][-1]["printed"]==.185
    assert groups[12][-1]["printed"]==.902
    eq=r.equilibrium_rows()
    assert {int(row["assay_id"]) for row in eq}=={3,7,9,10,11,12,13,14}


def test_paired_groups_partition_the_data_without_equilibrium_leakage():
    plan=r.read_json(r.PLAN)
    groups=r.read_batch(plan)
    test=[a for fold in plan["holdouts"]["paired_folds"].values() for a in fold]
    assert sorted(test)==sorted(groups)
    seen=[]
    class Dummy:
        @staticmethod
        def fit_equilibrium(temps,values):
            seen.append(list(zip(temps,values)))
            return type("fit",(),{"log_K_reference":1.5,"inverse_temperature_slope":.4,"log_K_rmse":.02})
    eq=r.equilibrium_training(Dummy(),r.equilibrium_rows(),[a for a in groups if a not in (3,7)])
    assert eq["1"]["assays"]==[13,14]
    assert len(seen[0])==2
    assert 3 not in eq["2"]["assays"]


def test_no_per_curve_activity_parameter():
    assert len(r.LAWS)==13 and len(r.PARAMETERS)==16
    # The published batch dose is not a reference concentration factor.
    groups=r.read_batch(r.read_json(r.PLAN))
    assert groups[8][0]["loading_pct"]==20
    assert groups[16][0]["loading_pct"]==10
    assert np.isclose(groups[13][0]["printed_ratio"],.497/.503)
