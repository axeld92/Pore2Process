"""Independent initial-data oracle using minimum-image FCC basis coordinates."""
import json
from pathlib import Path

import numpy as np


def evaluate(points,side,mean=.5,amplitude=.25):
    p=np.atleast_2d(np.asarray(points,dtype=float))
    x=np.mod(p,side)
    basis=side*np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
    delta=np.mod(x[:,None,:]-basis[None,:,:]+side/2,side)-side/2
    radius=side/np.cbrt(2*np.pi/(3*.4))/2
    norm=np.sqrt(np.sum(delta**2,axis=2));signed=norm-radius
    pick=np.argmin(signed,axis=1);rows=np.arange(len(x))
    d=signed[rows,pick];n=norm[rows,pick]
    direction=np.divide(delta[rows,pick],n[:,None],out=np.zeros_like(x),where=n[:,None]>0)
    surface=x-d[:,None]*direction
    buffer=.4*(side/np.sqrt(2)-2*radius)
    t=np.clip(d/buffer,0,1)
    w=10*t**3-15*t**4+6*t**5
    base=mean+amplitude*np.sin(2*np.pi*x[:,0]/side)
    trace=mean+amplitude*np.sin(2*np.pi*surface[:,0]/side)
    value=(1-w)*trace+w*base
    value[n==0]=base[n==0]
    return value,dict(signed_distance=d,surface=surface,normal=direction,
                     radius=radius,buffer=buffer,original=base)


def validate_initial(directory):
    directory=Path(directory)
    meta=json.loads((directory/'initial_profile_metadata.json').read_text())
    side=meta['side_m'];buffer=meta['buffer_m'];amplitude=meta['amplitude_kg_m3']
    a=np.genfromtxt(directory/'initial.csv',delimiter=',',names=True)
    coords=np.column_stack([a[k+'_m'] for k in 'xyz'])
    expected,data=evaluate(coords,side)
    error=float(np.max(abs(a['C_kg_m3']-expected)))
    fluid=a['material']==1
    p=np.genfromtxt(directory/'initial_profile_probes.csv',delimiter=',',names=True)
    pcoords=np.column_stack([p[k+'_m'] for k in 'xyz'])
    pexpected,_=evaluate(pcoords,side)
    probe_error=float(np.max(abs(p['C_kg_m3']-pexpected)))
    periodic_error=max(float(np.max(abs(p[k]-p['C_kg_m3']))) for k in ('periodic_x_C','periodic_y_C','periodic_z_C'))
    wall_trace_error=0.;normal_error=0.;transition_error=0.;transition_derivative_error=0.
    h=meta['probe_normal_step_m']
    for group in np.unique(p['group_id']):
        q=p[p['group_id']==group]
        if len(q)!=6:raise ValueError('Incomplete native compatibility probes.')
        q=np.sort(q,order='normal_distance_m')
        wall_trace_error=max(wall_trace_error,abs(q['C_kg_m3'][1]-(.5+.25*np.sin(2*np.pi*q['x_m'][1]/side))))
        normal_error=max(normal_error,abs((q['C_kg_m3'][2]-q['C_kg_m3'][0])/(2*h))*buffer/amplitude)
        transition_error=max(transition_error,abs(q['C_kg_m3'][4]-(.5+.25*np.sin(2*np.pi*q['x_m'][4]/side))))
        original=.5+.25*np.sin(2*np.pi*q['x_m'][3:6]/side)
        transition_derivative_error=max(transition_derivative_error,
            abs((q['C_kg_m3'][5]-q['C_kg_m3'][3]-original[2]+original[0])/(2*h))*buffer/amplitude)
    checks=dict(buffer_nonoverlap=0<buffer<meta['minimum_sphere_gap_m']/2,
        fixed_buffer=abs(buffer-data['buffer'])<=1e-15,
        native_formula=error<=1e-12,native_probe_formula=probe_error<=1e-12,
        periodicity=periodic_error<=1e-12,wall_trace=wall_trace_error<=1e-12,
        zero_normal_derivative=normal_error<=1e-7,transition_value=transition_error<=1e-12,
        transition_derivative=transition_derivative_error<=1e-7,
        range=bool(np.min(a['C_kg_m3'])>=.25-1e-13 and np.max(a['C_kg_m3'])<=.75+1e-13),
        stock_geometry=bool(np.all(data['signed_distance'][fluid]>-1e-14) and np.all(data['signed_distance'][~fluid]<=1e-14)))
    return dict(checks=checks,passed=all(checks.values()),metadata=meta,
        max_native_initial_formula_error_kg_m3=error,max_native_probe_formula_error_kg_m3=probe_error,
        max_periodic_difference_kg_m3=periodic_error,max_wall_trace_error_kg_m3=wall_trace_error,
        max_normal_derivative_times_buffer_over_amplitude=normal_error,
        max_transition_error_kg_m3=transition_error,
        max_transition_derivative_error_times_buffer_over_amplitude=transition_derivative_error,
        fluid_initial_minimum_kg_m3=float(np.min(a['C_kg_m3'][fluid])),
        fluid_initial_maximum_kg_m3=float(np.max(a['C_kg_m3'][fluid])),
        fluid_initial_mean_kg_m3=float(np.mean(a['C_kg_m3'][fluid])),
        fluid_initial_sd_kg_m3=float(np.std(a['C_kg_m3'][fluid])),
        max_change_from_original_sine_kg_m3=float(np.max(abs(a['C_kg_m3'][fluid]-data['original'][fluid]))))
