"""Immutable, session-specific exports with SI contracts and native provenance."""
import csv
from datetime import datetime, timezone
from io import BytesIO
import json
from pathlib import Path
import shutil
from zipfile import ZipFile, ZIP_DEFLATED

import numpy as np
import pyvista as pv

from pore2process_demo.io import sha
from pore2process_demo.io import write_json
from .geometry import polygon
from .plot import profile_figure

def export_completed(directory,case,native_path,source_path,protocol,versions):
    directory=Path(directory);directory.mkdir(parents=True)
    native_path=Path(native_path);native=json.loads(native_path.read_text())
    source_path=Path(source_path)
    versions={**versions,'DWSIM':case['dwsim']['version'],
              'OpenLB':native['cases'][native.get('accepted_case','grid48')]['raw']['version'],
              'DWSIM_property_package':case['dwsim']['property_package']}
    a=case['translation']['axial'];f=case['translation']['feed']
    summary={'status':'passed','run_id':directory.parent.name,'timestamp_utc':datetime.now(timezone.utc).isoformat(),
        'input_snapshot':case['request'],'case':case,'protocol':protocol,'versions':versions,
        'native_summary_sha256':sha(native_path),'support_summary_sha256':sha(source_path),
        'units':{'pressure':'Pa gauge TMP; DWSIM stream pressures absolute Pa','flow':'m3/s, positive +x',
            'length':'m','COD':'kg/m3 (COD equivalent mass, not a chemical species)','shaft_power':'W','specific_shaft_energy':'kWh/m3'},
        'indexing':'Native local block j,k are integers with y=(j-0.5)*dx, z=(k-0.5)*dx. Exported VTP point IDs start at zero.',
        'assumptions':protocol['scope'],'feed_evidence':protocol['feed'],
        'disclaimer':'Synthetic computational demonstration. COD retention is a transferred scenario assumption; critical flux and sustainable operating time are unavailable for the new geometry and speed.'}
    write_json(directory/'summary.json',summary)
    for source,target in [(native_path,'native_summary.json'),(source_path,'support_summary.json')]:
        shutil.copy2(source,directory/target)
    with (directory/'axial_profile.csv').open('w',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['x_m','TMP_Pa','local_flux_m_s','channel_flow_positive_x_m3_s','COD_kg_m3'])
        writer.writerows(zip(a['x_m'],a['pressure_Pa'],a['local_flux_m_s'],a['channel_flow_m3_s'],f['axial_COD_kg_m3'],strict=True))
    profile_figure(case['translation']).savefig(directory/'pressure_profile.png',dpi=150)
    accepted=native.get('accepted_case','grid48')
    section=native_path.parent/accepted/'section.csv'
    shutil.copy2(section,directory/'native_section.csv')
    shutil.copy2(native_path.parent/accepted/'channel.stl',directory/'native_periodic_channel_m.stl')
    values=np.genfromtxt(section,delimiter=',',names=True)
    sample_x=native['cases'][accepted]['raw']['dx_m']
    mesh=pv.PolyData(np.column_stack((np.full(len(values),sample_x),values['y_m'],values['z_m'])))
    for name in ('ux_m_s','uy_m_s','uz_m_s','density_kg_m3'):mesh.point_data['native_'+name]=values[name]
    mesh.point_data['point_id_zero_based']=np.arange(len(values),dtype=np.int64)
    mesh.field_data['native_section_x_m']=[sample_x]
    mesh.field_data['coordinate_units']=['m'];mesh.field_data['velocity_scope']=['Native low-Re verification field; no rescaled process velocity is represented as an OpenLB solve.']
    mesh.save(directory/'native_section.vtp')
    # Keep the actual simulator archive, not an image or a surrogate.
    shutil.copy2(directory.parent/'crossflow.dwxmz',directory/'process.dwxmz')
    shutil.copy2(directory.parent/'flowsheet.png',directory/'flowsheet.png')
    provenance={p.name:sha(p) for p in directory.iterdir() if p.is_file()}
    write_json(directory/'manifest.json',provenance)
    result=BytesIO()
    with ZipFile(result,'w',ZIP_DEFLATED) as archive:
        for p in sorted(directory.iterdir()):
            if p.is_file():archive.write(p,p.name)
    result.seek(0)
    return summary,result

def native_pane(native_path,existing=None):
    import panel as pn
    native_path=Path(native_path);s=json.loads(native_path.read_text())
    accepted=s.get('accepted_case','grid48')
    c=s['cases'][accepted];section=np.genfromtxt(native_path.parent/accepted/'section.csv',delimiter=',',names=True)
    p=pv.Plotter(off_screen=True,window_size=(700,430));p.set_background('#f6f9fb')
    wall=pv.read(native_path.parent/accepted/'channel.stl');wall.points*=1e3
    # Clip the extended STL to the actual periodic segment for display.
    length=c['raw']['actual_axial_planes']*c['raw']['dx_m']*1e3
    wall=wall.clip_box((0,length,0,4.6,0,4.6),invert=False)
    p.add_mesh(wall,color='#879aa6',opacity=.25,show_edges=True)
    points=pv.PolyData(np.column_stack((np.full(len(section),c['raw']['dx_m']*1e3),section['y_m']*1e3,section['z_m']*1e3)))
    points['Native axial velocity (m/s)']=section['ux_m_s']
    p.add_mesh(points,scalars='Native axial velocity (m/s)',render_points_as_spheres=True,
        point_size=7,cmap='viridis',show_scalar_bar=False)
    p.view_isometric();p.reset_camera()
    if existing is None:
        axes={'fontsize':12,'digits':2,'xticker':{'ticks':[length/2],'labels':[f'x {length/2:.2f}']},
              'yticker':{'ticks':[1.15,2.3,3.45],'labels':['y 1.15','y 2.30','y 3.45']},
              'zticker':{'ticks':[1.15,2.3,3.45],'labels':['z 1.15','z 2.30','z 3.45']}}
        pane=pn.pane.VTK(p.ren_win,height=420,sizing_mode='stretch_width',orientation_widget=False,axes=axes)
        view=pn.Column(pane,pane.construct_colorbars(),
            pn.pane.Markdown(f'Periodic window: x = 0–{length:.4g} mm; sampled field at x = {c["raw"]["dx_m"]*1e3:.4g} mm.'),sizing_mode='stretch_width')
    else:
        view=existing;pane=view[0];pane.object=p.ren_win
        pane.color_mappers=list(pane.get_color_mappers(infer=True))
        view[1].object=pane.construct_colorbars().object;pane.reset_camera()
    return view,p
