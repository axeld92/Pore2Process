"""One native field export for the study; no rescaled velocity claimed native."""
from pathlib import Path
import shutil
import numpy as np
import pyvista as pv
from pore2process_demo.io import read_json as read

def native_field(native_path,directory):
    native_path=Path(native_path);directory=Path(directory);directory.mkdir()
    n=read(native_path);c=n['cases'][n['accepted_case']];base=native_path.parent/n['accepted_case']
    values=np.atleast_1d(np.genfromtxt(base/'section.csv',delimiter=',',names=True));x=c['raw']['dx_m']
    mesh=pv.PolyData(np.column_stack((np.full(len(values),x),values['y_m'],values['z_m'])))
    for name in ('ux_m_s','uy_m_s','uz_m_s','density_kg_m3'):mesh.point_data[name]=values[name]
    mesh.point_data['point_id_zero_based']=np.arange(len(values),dtype=np.int64)
    mesh.field_data['coordinate_units']=['m'];mesh.field_data['native_measured_Re']=[c['Re_mean']]
    mesh.field_data['scope']=['Actual native high-Re developed laminar branch; not a rescaled process field.']
    mesh.save(directory/'native_section.vtp')
    for name in ('section.csv','planes.csv','history.csv','raw.json','channel.stl','request.json'):shutil.copy2(base/name,directory/name)
