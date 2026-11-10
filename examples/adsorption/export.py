"""Export the actual measured native field, never a rescaled process velocity."""
import numpy as np
import pyvista as pv
import json

def export_native(directory,path):
    a=np.genfromtxt(directory/'field.csv',delimiter=',',names=True);r=json.loads((directory/'raw.json').read_text())
    mesh=pv.PolyData(np.column_stack([a[k] for k in ('x_m','y_m','z_m')]))
    mesh.point_data['point_id_zero_based']=np.arange(len(a),dtype=np.int64)
    mesh.point_data['lattice_index_zero_based']=np.column_stack([a[k] for k in ('i','j','k')]).astype(np.int64)
    mesh.point_data['native_velocity_m_s']=np.column_stack([a[k] for k in ('ux_m_s','uy_m_s','uz_m_s')])
    mesh.point_data['native_density_kg_m3']=a['density_kg_m3']
    mesh.field_data['coordinate_units']=['m']
    mesh.field_data['scope']=['Native fluid cell centres at reference diameter/properties; no rescaled process field.']
    mesh.field_data['native_diameter_m']=[r['particle_diameter_m']]
    mesh.field_data['native_viscosity_Pa_s']=[r['viscosity_Pa_s']]
    mesh.save(path)
