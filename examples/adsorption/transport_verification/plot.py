"""Publication-style static diagnostics from checked, saved fields."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

from pore2process_demo.io import sha,write_json
from pore2process_demo.paths import ROOT
from examples.adsorption.transport_accuracy.run import validate_record
from . import reference as ref
from .analysis import load,profile,SHELL_M,RMS_LIMIT,EXCURSION_LIMIT
from .run import HERE,OUTPUT


def latest(stage):
    return ROOT/json.loads((OUTPUT/f'{stage}_latest.json').read_text())['summary']


def plot_fields(stretch,disk,output):
    r,_,_,f,_,q=load(stretch/'native/stretch96');x,y,z=q;i=r['n']//2
    mode=np.load(stretch/'reference_modes.npy',allow_pickle=False)[-1]
    fine=(np.arange(256)+.5)*ref.L/256
    continuum=.5+ref.A*np.imag(np.exp(1j*ref.K*x[i])*ref.sample(mode,fine,fine))
    at_nodes=.5+ref.A*np.imag(np.exp(1j*ref.K*x[i])*ref.sample(mode,y,z))
    r2,mask,_,f2,_,q2=load(disk/'native/rotation128');q2x,q2y,_=q2
    X,Y=np.meshgrid(fine,fine,indexing='ij')
    exact=np.where(np.hypot(X-ref.L/2,Y-ref.L/2)<=ref.R,ref.disk_profile(X,Y,ref.TIMES[-1]),np.nan)
    exact_nodes=profile('disk',q2,ref.TIMES[-1],ref.U,ref.A)[:,:,0]
    material=mask[:,:,0]==1;native_disk=np.where(material,f2['final'][:,:,0],np.nan)
    data=[(continuum,f['final'][i],f['final'][i]-at_nodes,'Stretching, N96',('y','z')),
          (exact,native_disk,np.where(material,native_disk-exact_nodes,np.nan),'Circular wall, N128',('x','y'))]
    fig,axs=plt.subplots(2,3,figsize=(11.5,7.1),layout='constrained');Lmm=ref.L*1000
    cmap=plt.get_cmap('coolwarm').copy();cmap.set_bad('#555b62');extent=(0,Lmm,0,Lmm)
    span=max(.25,max(float(np.nanmax(abs(field-.5))) for row in data for field in row[:2]))
    span=np.ceil(span*20)/20
    diffspan=max(float(np.nanmax(abs(row[2]))) for row in data)
    for row,(exact,native,err,name,labels) in enumerate(data):
        for col,(field,title) in enumerate(zip((exact,native,err),('Continuum reference','OpenLB','OpenLB − reference'))):
            ax=axs[row,col];limit=diffspan if col==2 else span;center=0 if col==2 else .5
            image=ax.imshow(field.T,extent=extent,origin='lower',cmap=cmap,vmin=center-limit,vmax=center+limit,interpolation='nearest')
            ax.set(title=f'{name}\n{title}',xlabel=f'{labels[0]} (mm)',ylabel=f'{labels[1]} (mm)')
            if row==1:
                ax.add_patch(Circle((Lmm/2,Lmm/2),ref.R*1000,fill=False,color='black',linewidth=.7))
            if col==2:delta=image
            else:concentration=image
    fig.colorbar(concentration,ax=axs[:,:2].ravel().tolist(),shrink=.82,label='Concentration (kg/m³)')
    fig.colorbar(delta,ax=axs[:,2].tolist(),shrink=.82,label='Difference (kg/m³)')
    fig.savefig(output,dpi=175);plt.close(fig)


def plot_errors(stretch,disk,output):
    fig,axs=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    ns=[32,64,96];nc=[32,64,96,128]
    for key,label in [('middle','Halfway'),('final','Final')]:
        axs[0,0].plot(ns,[100*stretch['cases'][f'stretch{n}']['fields'][key]['RMS_over_amplitude'] for n in ns],'o-',label=label)
    for key,label in [('RMS_over_amplitude','Rotation, all liquid'),('wall_shell_RMS_over_amplitude','Rotation, outer 50 µm')]:
        axs[0,1].plot(nc,[100*disk['cases'][f'rotation{n}']['fields']['final'][key] for n in nc],'o-',label=label)
    axs[0,1].plot([32,128],[100*disk['cases'][f'diffusion{n}']['fields']['final']['RMS_over_amplitude'] for n in (32,128)],'s--',label='Diffusion only')
    for ax,title in zip(axs[0],('Stretching field error','Circular-wall field error')):
        ax.axhline(RMS_LIMIT*100,ls='--',color='black',lw=1,label='1% screen')
        ax.set(yscale='log',xlabel='Cells across the square or cube',ylabel='RMS (% of initial amplitude)',title=title)
        ax.legend(fontsize=8)
    for ax,record,grids,prefix,title in [(axs[1,0],stretch,ns,'stretch','Wall-free excursions'),(axs[1,1],disk,nc,'rotation','Circular-wall excursions')]:
        base=ROOT/record['path'];base=base.parent
        for n in grids:
            _,_,_,_,h,_=load(base/'native'/f'{prefix}{n}')
            e=np.maximum(0,np.maximum(.25-h['min_concentration_kg_m3'],h['max_concentration_kg_m3']-.75))
            ax.plot(h['time_s'],e,label=f'N{n}')
        ax.axhline(EXCURSION_LIMIT,ls='--',color='black',lw=1,label='Allowance')
        ax.set(xlabel='Time (s)',ylabel='Excursion (kg/m³)',title=title);ax.legend(fontsize=8)
    fig.savefig(output,dpi=175);plt.close(fig)


def stretching_gradient_regions(final,mode,dx):
    """Compare final field error in the sharpest reference gradients."""
    q=(np.arange(final.shape[0])+.5)*dx
    k=2*np.pi*np.fft.fftfreq(len(mode),d=ref.L/len(mode));F=np.fft.fft2(mode)
    f=ref.sample(mode,q,q);phase=np.exp(1j*ref.K*q[:,None,None]);exact=.5+ref.A*np.imag(phase*f)
    derivatives=(1j*ref.K*f,ref.sample(np.fft.ifft2(1j*k[:,None]*F),q,q),
                 ref.sample(np.fft.ifft2(1j*k[None,:]*F),q,q))
    norm=np.sqrt(sum(np.imag(phase*a)**2 for a in derivatives))*ref.A
    cutoff=np.quantile(norm,.9);high=norm>=cutoff;error=(final-exact)/ref.A
    return dict(gradient_decile_cutoff_kg_m4=float(cutoff),
        high_gradient_RMS_over_amplitude=float(np.sqrt(np.mean(error[high]**2))),
        remaining_RMS_over_amplitude=float(np.sqrt(np.mean(error[~high]**2))))


def disk_background(base):
    """Diagnostic superposition at identical times, not a corrected solution."""
    d=base/'native';r,mask,u,wave,_,q=load(d/'rotation128')
    c,other,other_u,constant,_,_=load(d/'constant128')
    for k in ('n','nz','dx_m','dt_s','steps','speed_m_s','diffusivity_m2_s'):
        if r[k]!=c[k]:raise ValueError('Disk background case is not matched')
    np.testing.assert_array_equal(mask,other);np.testing.assert_array_equal(u,other_u)
    radius=np.hypot(q[0][:,None,None]-ref.L/2,q[1][None,:,None]-ref.L/2)
    fluid=mask==1;wall=fluid&(radius>=ref.R-SHELL_M);result={}
    for key,t in zip(('middle','final'),ref.TIMES):
        err=wave[key]-profile('disk',q,t,ref.U,ref.A);background=constant[key]-.5
        def measures(sel):
            rms=lambda x:float(np.sqrt(np.mean(x[sel]**2))/ref.A)
            return dict(original_RMS_over_amplitude=rms(err),background_RMS_over_amplitude=rms(background),
                        diagnostic_wave_only_RMS_over_amplitude=rms(err-background))
        result[key]=dict(all_liquid=measures(fluid),outer_shell=measures(wall))
    return dict(fields=result,production_correction_applied=False)


def compact_summary(sp,cp,st,cu):
    modes=np.load(sp.parent/'reference_modes.npy',allow_pickle=False)
    continuum=[]
    for t,f in zip(ref.TIMES,modes):
        k=2*np.pi*np.fft.fftfreq(len(f),d=ref.L/len(f));F=np.fft.fft2(f)
        a=(1j*ref.K*f,np.fft.ifft2(1j*k[:,None]*F),np.fft.ifft2(1j*k[None,:]*F))
        # Exact maximum over the x-phase of |Im(exp(ikx) a)|^2.
        g2=(sum(abs(v)**2 for v in a)+abs(sum(v*v for v in a)))/2
        if np.max(abs(f))>1+1e-4:raise ValueError('Continuum reference violates its range')
        continuum.append(dict(time_s=t,minimum_inverse_gradient_length_m=float(1/np.sqrt(np.max(g2))),
             concentration_sd_kg_m3=float(ref.A*np.sqrt(np.mean(abs(f)**2)/2))))
    native={}
    for path,study in ((sp,st),(cp,cu)):
        entries={}
        for name,c in study['cases'].items():
            raw=c['raw'];shape=(raw['n'],raw['n'],raw['nz']);d=path.parent/'native'/name
            fluid=np.fromfile(d/'material.u8',dtype='u1').reshape(shape)==1
            final=np.fromfile(d/'final.f64',dtype='<f8').reshape(shape)
            entries[name]=dict(grid=list(shape),range_kg_m3=c['range_kg_m3'],maximum_excursion_kg_m3=c['maximum_excursion_kg_m3'],
                closed_inventory_relative_drift=c['closed_inventory_relative_drift'],physical_liquid_relative_drift=c['physical_liquid_relative_drift'],
                fields=c.get('fields'),checks=c['checks'],final_liquid_sd_kg_m3=float(np.std(final[fluid])),
                dx_m=raw['dx_m'],peak_time_s=c['peak_time_s'],native_wall_s=raw['wall_s'],
                fluid_volume_m3=c['initial_fluid_volume_m3'],analytic_volume_m3=c.get('analytic_fluid_volume_m3'))
            if name.startswith('stretch'):
                entries[name]['final_reference_gradient_regions']=stretching_gradient_regions(final,modes[-1],raw['dx_m'])
            if raw['kind']=='disk':
                q=[(np.arange(n)+.5)*raw['dx_m'] for n in shape]
                peak=np.fromfile(d/'peak.f64',dtype='<f8').reshape(shape)
                def maximum_point(field):
                    idx=np.unravel_index(np.argmax(np.where(fluid,field,-np.inf)),shape)
                    gap=ref.R-np.hypot(q[0][idx[0]]-ref.L/2,q[1][idx[1]]-ref.L/2)
                    return dict(ijk=[int(i) for i in idx],smooth_wall_gap_m=float(gap),value=float(field[idx]))
                exact=profile('disk',q,ref.TIMES[-1],raw['speed_m_s'],raw['amplitude_kg_m3'])
                entries[name]['peak_high_concentration']=maximum_point(peak)
                entries[name]['final_maximum_absolute_error']=maximum_point(abs(final-exact))
        native[study['stage']]=dict(summary=str(path.relative_to(ROOT)),summary_sha256=sha(path),cases=entries)
    return dict(schema='pore2process.transport-verification-reference.v1',stages=native,
        reference_convergence=st['reference'],continuum_stretching=continuum,depth_comparison=cu['depth_comparison'],
        disk_matched_background=disk_background(cp.parent),
        native_runs=st['new_native_runs']+cu['new_native_runs'],new_bed_runs=0,new_DWSIM_runs=0,
        application_transport_qualified=False,source_sha256=sha(Path(__file__)))


def main():
    p=argparse.ArgumentParser();p.add_argument('--stretching',type=Path);p.add_argument('--curved',type=Path);a=p.parse_args()
    sp=(a.stretching or latest('stretching')).resolve();cp=(a.curved or latest('curved')).resolve()
    st,cu=validate_record(sp),validate_record(cp)
    if cu['parent_summary_sha256']!=sha(sp):raise ValueError('Curved stage belongs to another stretching record')
    out=HERE/'figures';out.mkdir(exist_ok=True)
    plot_fields(sp.parent,cp.parent,out/'field_comparison.png')
    plot_errors(dict(st,path=str(sp.relative_to(ROOT))),dict(cu,path=str(cp.relative_to(ROOT))),out/'mesh_comparison.png')
    write_json(HERE/'reference_summary.json',compact_summary(sp,cp,st,cu))
    write_json(OUTPUT/'figures.json',dict(inputs={str(p.relative_to(ROOT)):sha(p) for p in (sp,cp)},
        analysis_source_sha256=sha(Path(__file__)),outputs={str(p.relative_to(ROOT)):sha(p) for p in out.glob('*.png')},visual_inspection='pending'))
    print('Saved',out)


if __name__=='__main__':main()
