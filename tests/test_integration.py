from dataclasses import asdict, replace
import io
import json
from pathlib import Path
import subprocess
from zipfile import ZipFile

import numpy as np
import pytest
import pyvista as pv

from pore2process_demo.bridge import ROOT, ColumnInputs, NetworkInputs, predict_column_flow, run_case
from pore2process_demo.export import summary, write_exports, zip_export
from pore2process_demo.units import from_display, microlitres_per_minute, millilitres_per_minute
from pore2process_demo.visualization import pressure_profile


@pytest.fixture(scope="module")
def default():
    return run_case()


def test_direct_julia_matches_bridge(default):
    config = json.loads((ROOT / ".local/julia_runtime.json").read_text())
    code = 'using Pore2ProcessDemo; s=run_case(); println(join(vcat(s.Q_in_m3_s,s.permeability_m2,s.pressures_Pa),","))'
    output = subprocess.check_output([config["executable"],"--startup-file=no",f"--project={ROOT / 'julia'}","-e",code],text=True)
    direct = np.fromstring(output.strip(),sep=",")
    np.testing.assert_allclose(direct,np.r_[default.Q_in_m3_s,default.permeability_m2,default.pressures_Pa],rtol=1e-12,atol=0)


def test_contract_and_analytical_reference(default):
    s = default
    assert s.coordinates_m.shape == (96,3)
    assert s.edges.shape == (224,2)
    for array in (s.coordinates_m,s.pressures_Pa,s.flows_m3_s,s.radii_m,s.lengths_m,s.net_outflows_m3_s,s.conductances_m3_Pa_s):
        assert array.dtype == np.float64
        assert not array.flags.writeable
    assert s.edges.dtype == np.int64
    np.testing.assert_array_equal(s.edges[:3],[[0,1],[0,6],[0,24]])
    assert s.edges.min() == 0 and s.edges.max() == 95
    np.testing.assert_array_equal(s.inlet_nodes,np.arange(0,96,6))
    np.testing.assert_array_equal(s.outlet_nodes,np.arange(5,96,6))
    assert len(s.interior_nodes) == 64
    for k in range(4):
        for j in range(4):
            for i in range(6):
                np.testing.assert_allclose(s.coordinates_m[i+6*(j+4*k)],np.array([i,j+.5,k+.5])*1e-3)
    np.testing.assert_allclose(s.flows_m3_s,s.conductances_m3_Pa_s*(s.pressures_Pa[s.edges[:,0]]-s.pressures_Pa[s.edges[:,1]]),atol=0)
    x,mean = pressure_profile(s)
    np.testing.assert_allclose(mean,100-100*x/.005,atol=1e-7)
    assert s.Q_in_m3_s == pytest.approx(7.853981633974485e-10,rel=1e-9)
    assert s.permeability_m2 == pytest.approx(np.pi*(50e-6)**4/(8*(1e-3)**2),rel=1e-9)
    assert s.max_interior_imbalance_m3_s <= 1e-9*7.853981633974485e-10
    assert s.flow_mismatch_m3_s <= 1e-9*7.853981633974485e-10
    assert predict_column_flow(s) == pytest.approx(1.227184630308513e-8,rel=1e-9)


@pytest.mark.parametrize("change,qfactor,kfactor",[({"p_in_Pa":200.},2,1),({"radius_m":100e-6},16,16),({"viscosity_Pa_s":2e-3},.5,1),({"p_in_Pa":0.,"p_out_Pa":100.},-1,1)])
def test_changed_inputs(default,change,qfactor,kfactor):
    other = run_case(replace(default.inputs,**change))
    assert other.Q_in_m3_s == pytest.approx(default.Q_in_m3_s*qfactor,rel=1e-9)
    assert other.permeability_m2 == pytest.approx(default.permeability_m2*kfactor,rel=1e-9)
    assert not np.shares_memory(default.pressures_Pa,other.pressures_Pa)


def test_units():
    s = from_display(radius_um=50,spacing_mm=1,viscosity_mPa_s=1)
    assert s.radius_m == pytest.approx(50e-6)
    assert s.spacing_m == .001 and s.viscosity_Pa_s == .001
    assert microlitres_per_minute(1e-9) == pytest.approx(60.)
    assert millilitres_per_minute(1e-6) == pytest.approx(60.)
    assert millilitres_per_minute(None) is None


def test_worker_calls_rejected_before_julia(default):
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=1) as pool:
        with pytest.raises(RuntimeError,match="main server thread"):
            pool.submit(run_case).result(timeout=5)
        with pytest.raises(RuntimeError,match="main server thread"):
            pool.submit(predict_column_flow,default).result(timeout=5)


def test_export_snapshot_and_roundtrip(default,tmp_path):
    column = ColumnInputs(length_m=.04)
    write_exports(default,tmp_path,column)
    assert {p.name for p in tmp_path.iterdir()} == {"summary.json","nodes.csv","throats.csv","pressure_profile.png","network.vtp"}
    saved = json.loads((tmp_path / "summary.json").read_text())
    assert saved["inputs"] == asdict(default.inputs)
    assert saved["units"]["flow"] == "m^3/s"
    assert saved["Q_in_m3_s"] == default.Q_in_m3_s
    assert saved["column"]["network_run_id"] == default.run_id
    assert saved["column"]["inputs"] == asdict(column)
    assert saved["column"]["flow_m3_s"] == pytest.approx(predict_column_flow(default)/2)
    assert saved["versions"]["PythonCall"] == saved["versions"]["juliacall"]
    nodes = np.genfromtxt(tmp_path / "nodes.csv",delimiter=",",names=True)
    throats = np.genfromtxt(tmp_path / "throats.csv",delimiter=",",names=True)
    np.testing.assert_allclose(nodes["pressure_Pa"],default.pressures_Pa)
    np.testing.assert_array_equal(throats["from_node_id_zero_based"],default.edges[:,0])
    np.testing.assert_allclose(throats["flow_from_to_m3_s"],default.flows_m3_s)
    mesh = pv.read(tmp_path / "network.vtp")
    assert mesh.n_points == 96 and mesh.n_lines == 224
    np.testing.assert_allclose(mesh.points,default.coordinates_m)
    np.testing.assert_allclose(mesh.point_data["pressure_Pa"],default.pressures_Pa)
    np.testing.assert_allclose(mesh.cell_data["flow_from_first_to_second_m3_s"],default.flows_m3_s)
    assert (tmp_path / "pressure_profile.png").read_bytes().startswith(b"\x89PNG")
    archive = ZipFile(zip_export(default,column))
    assert json.loads(archive.read("summary.json"))["run_id"] == default.run_id


def test_zero_and_validation(tmp_path):
    zero = run_case(NetworkInputs(p_in_Pa=73.,p_out_Pa=73.))
    assert zero.permeability_m2 is None and predict_column_flow(zero) is None
    assert np.max(np.abs(zero.flows_m3_s)) < 1e-20
    write_exports(zero,tmp_path)
    content = (tmp_path / "summary.json").read_text()
    assert "NaN" not in content and "Infinity" not in content
    assert json.loads(content)["column"]["flow_m3_s"] is None
    for change in ({"nx":1},{"nx":2.5},{"ny":True},{"radius_m":.0005},{"viscosity_Pa_s":0},{"p_in_Pa":float("nan")}):
        with pytest.raises(ValueError) as error:
            run_case(replace(NetworkInputs(),**change))
        assert "Stacktrace:" not in str(error.value)
