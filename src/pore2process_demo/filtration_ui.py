"""Per-session explicit-channel study; Julia on event loop, DWSIM subprocess."""
import asyncio
import copy
from datetime import datetime, timezone
import html
from io import BytesIO
from pathlib import Path
from uuid import uuid4

import panel as pn

from .bridge import runtime, runtime_versions
from examples.filtration.star_channel import model
from examples.filtration.star_channel.export import export_completed, native_pane
from examples.filtration.star_channel.plot import profile_figure
from examples.filtration.star_channel.run import evidence, support_evidence
from pore2process_demo.paths import ROOT; from pore2process_demo.io import read_json as read
from pore2process_demo.io import write_json
from pore2process_demo.provenance import verified_source_path
from examples.filtration.crossflow_recycle import dwsim
from examples.filtration.crossflow_recycle.model import validate, topology


class FiltrationApp:
    def __init__(self):
        self.completed=None;self.archive=None;self.plotter=None;self.vtk_container=None
        self.session_id=uuid4().hex
        self.controls={
            'support':pn.widgets.Select(label='Support interior radius',options={'4 µm · straight':'straight','5 µm · widened':'widened_5um','6 µm · widened':'widened_6um'},value='straight'),
            'flux_LMH':pn.widgets.FloatInput(label='Mean filtration flux (LMH)',value=12.,step=1.,start=8.,end=12.),
            'inlet_Re':pn.widgets.Select(label='Inlet Reynolds number',options=[10.,20.,30.],value=20.),
            'rejection':pn.widgets.FloatSlider(label='Assumed local COD rejection',start=.5,end=.8,step=.05,value=.65),
            'viscosity_multiplier':pn.widgets.Select(label='Assumed viscosity multiplier',options=[1.,1.2],value=1.),
            'fixed_multiplier':pn.widgets.Select(label='Fixed non-support impedance multiplier',options=[.5,1.,2.],value=1.),
        }
        self.run_button=pn.widgets.Button(label='Run filtration process',color='primary',sizing_mode='stretch_width')
        self.download=pn.widgets.FileDownload(label='Download completed process (.zip)',filename='filtration.zip',disabled=True,sizing_mode='stretch_width')
        self.download_slot=pn.Column(self.download,margin=0,sizing_mode='stretch_width')
        self.status=pn.pane.Alert('Loading verified study evidence…',alert_type='info',sizing_mode='stretch_width')
        self.metrics=pn.pane.Markdown('Run a case to calculate the actual DWSIM recycle and pump duties.',sizing_mode='stretch_width')
        self.profile=pn.Column(sizing_mode='stretch_width');self.network=pn.Column(sizing_mode='stretch_width')
        self.source_view=pn.Column(sizing_mode='stretch_width')
        try:
            self.native_path,self.parent_path=evidence();self.source_path=support_evidence()
            self.native=read(self.native_path);self.source=read(self.source_path)
            self.protocol=read(ROOT/'examples/filtration/star_channel/protocol.json');self.h=model.settings(read(self.parent_path))
            # A source hash mismatch cannot silently reuse stale geometry results.
            from pore2process_demo.io import sha
            for name,digest in self.native['source_sha256'].items():
                if sha(verified_source_path(name,digest))!=digest:raise ValueError('Channel implementation changed. Rerun examples.filtration.star_channel.native_run.')
            self.ready=True
            self.status.object='Ready. This runs Julia scale translation and an actual DWSIM solve, save and reload.'
            self._sources()
            try:
                view,self.plotter=native_pane(self.native_path);self.network.objects=[view];self.vtk_container=view
            except Exception as e:
                self.network.objects=[pn.pane.Alert('3D view unavailable: '+html.escape(str(e)),alert_type='warning')]
        except Exception as e:
            self.ready=False;self.run_button.disabled=True
            self.status.object='Study evidence unavailable: '+html.escape(str(e))+'. Run the setup commands in examples/filtration/star_channel/README.md.'
            self.status.alert_type='warning'
        self.run_button.on_click(self._on_run)
        for widget in self.controls.values():widget.param.watch(self._stale,'value')
        self.view=self._layout()
        if pn.state.curdoc is not None:pn.state.on_session_destroyed(lambda context:self.close())

    def _layout(self):
        for widget in self.controls.values():widget.sizing_mode='stretch_width'
        controls=pn.Column('### Filtration inputs',*self.controls.values(),self.run_button,self.download_slot,
            'All supports retain 4 µm openings. The fixed effective non-support impedance remains unchanged within each comparison.',
            'Geometry controls select verified OpenLB measurements. Run recalculates Julia and solves DWSIM; new CFD geometries require the documented native verification commands.',
            width=300,styles={'background':'#f1f6f7','padding':'16px','border-radius':'12px'})
        results=pn.Column(self.status,self.metrics,self.profile,
            pn.Row(pn.Column('### Verified OpenLB geometry and velocity',self.network,
                'Coordinates are in mm; +x is axial flow. Rotate by dragging. The short periodic segment is shown at physical scale. Colours are the actual low-Re OpenLB verification field, not a rescaled field claimed as another simulation.',sizing_mode='stretch_width'),sizing_mode='stretch_width'),
            sizing_mode='stretch_width',styles={'flex':'1 1 650px','min-width':'0'})
        return pn.Column('# Filtration: explicit pores, channel and process',
            'A six-point synthetic star channel → Julia distributed pressure/flux and dilute COD balance → an actual DWSIM water-carrier recycle. COD and a retention sensitivity come from Chiu et al. (2006). This is a clean initial hydraulic case; sustainable flux at this geometry and speed is unavailable.',
            pn.FlexBox(controls,results,sizing_mode='stretch_width'),
            pn.Accordion(('Evidence, comparisons and assumptions',self.source_view),active=[0],sizing_mode='stretch_width'),
            max_width=1600,sizing_mode='stretch_width',margin=20)

    def _sources(self):
        g=self.native['geometry'];f=self.protocol['feed']
        notes=(f"**Channel:** six points, outer radius 2.3 mm, inner radius 1.4 mm; open area {g['open_area_m2']*1e6:.4g} mm², "
            f"wetted perimeter {g['wetted_perimeter_m']*1e3:.5g} mm, hydraulic diameter {g['hydraulic_diameter_m']*1e3:.5g} mm; active length 0.3 m. "
            'Parallel replication supplies 1 m² membrane area. OpenLB uses stock TRT collision and stock exact-difference forcing.\n\n'
            '**Feed evidence:** biological-treatment effluent, module-inlet COD 68.7 mg/L (reported range 65.3–72.1); '
            '0.2 µm membrane, pH 7.5, 20–25 °C. Lowest measured point: Re 1388, critical flux 21.4 LMH. '
            'That critical-flux point is outside our verified Reynolds-number range and is retained as source context. No local or sustainable flux limit is inferred for the synthetic channel.\n\n'
            '**Retention:** 0.65 was measured below critical flux at 1.6 m/s and 0.1 bar. Its use at our different speed/TMP is a transfer assumption; 0.5–0.8 is a scenario sensitivity. '
            'COD is balanced in Julia as an external lumped tracer. DWSIM uses pure-water properties. Fresh concentration is calculated to retain the stated module-inlet concentration.\n\n'
            '**Limits:** no validated fouling time, particle transport, turbulence, osmosis, reactions or cleaning model. The common non-support impedance is effective and cannot be identified uniquely as selective-layer resistance.')
        self.source_view.objects=[pn.pane.Markdown(notes)]
        support_view=ROOT/'examples/filtration/graded_support/figures/graded_support_geometries.png'
        if support_view.exists():
            self.source_view.append(pn.Accordion(('Support geometries · unchanged 4 µm openings',
                pn.pane.PNG(support_view,sizing_mode='stretch_width')),active=[],sizing_mode='stretch_width'))
        latest=ROOT/'results/star_channel/process/latest.json'
        if latest.exists():
            path=Path(read(latest)['summary'])
            self.source_view.append(pn.pane.PNG(path.parent/'study.png',sizing_mode='stretch_width'))
        old=ROOT/'examples/filtration/literature_benchmark/figures/crossflow_benchmark.png'
        if old.exists():
            self.source_view.append(pn.Accordion(('Earlier TiO₂ comparison · unresolved source flow convention',pn.Column(
                'Historical literature-flow comparison with earlier BGK support measurements. It uses a different feed and operating basis. New process runs use the separately revalidated TRT supports.',
                pn.pane.PNG(old,sizing_mode='stretch_width'))),active=[],sizing_mode='stretch_width'))

    def _stale(self,event=None):
        if self.completed is not None:
            self.status.object='Inputs changed. Display and download belong to the previous completed process '+self.completed['run_id']+'. Press Run filtration process to update.'
            self.status.alert_type='warning'

    def _request(self):
        return {name:widget.value for name,widget in self.controls.items()}

    async def _on_run(self,event=None):
        if self.run_button.disabled:return
        self.run_button.disabled=True;self.run_button.loading=True
        for w in self.controls.values():w.disabled=True
        self.status.object='Running Julia translation and DWSIM. The simulator runs in a separate process; inputs are locked for this run.'
        self.status.alert_type='info'
        try:
            # Async Bokeh callback executes on the event loop. Julia remains a
            # synchronous main-thread call; only the DWSIM child is awaited.
            await asyncio.sleep(.05)
            request=copy.deepcopy(self._request());jl=runtime()
            translated=model.translate(jl,self.native,self.source,self.protocol,request)
            a=model.inputs(self.h,self.source,translated)
            run_id=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'-'+uuid4().hex[:8]
            directory=ROOT/'results/star_channel/sessions'/self.session_id/run_id
            directory.mkdir(parents=True);write_json(directory/'inputs.json',a)
            if not dwsim.RUNNER.is_file():raise ValueError('DWSIM adapter is not built. Run python -m examples.filtration.star_channel.run first.')
            command=[str(dwsim.DOTNET),str(dwsim.RUNNER),str(dwsim.DISTRIBUTION),str(directory/'inputs.json'),str(directory)]
            with (directory/'stdout.log').open('w') as log:
                process=await asyncio.create_subprocess_exec(*command,cwd=dwsim.DISTRIBUTION,env=dwsim.environment(),stdout=log,stderr=asyncio.subprocess.STDOUT)
                try:code=await asyncio.wait_for(process.wait(),180)
                except asyncio.TimeoutError:
                    process.kill();await process.wait();raise ValueError('DWSIM exceeded the 180-second limit; the failed process log is retained.')
            if code:raise ValueError('DWSIM solve failed; diagnostic log: '+str(directory/'stdout.log'))
            actual=read(directory/'summary.json')
            checks=validate(actual,a,translated['process'],self.h,self.source)
            case={'id':model.key(request),'request':request,'translation':translated,'dwsim':actual,
                'validation':checks,'saved_connections':topology(directory/'crossflow.dwxmz')}
            case['verified_baseline_comparison']=self._baseline_comparison(case)
            completed,archive=export_completed(directory/'export',case,self.native_path,self.source_path,self.protocol,runtime_versions())
            completed['run_id']=run_id
            fig=pn.pane.Matplotlib(profile_figure(translated),tight=True,sizing_mode='stretch_width',height=370)
            self.completed=completed;self.archive=archive;self.profile.objects=[fig]
            self.metrics.object=self._metrics(case,run_id)
            # Build a complete embedded link before mounting it. Updating an
            # initially disabled FileDownload can retain its old blob URL.
            self.download=pn.widgets.FileDownload(label='Download completed process (.zip)',
                filename='pore2process-filtration-'+run_id+'.zip',file=BytesIO(archive.getvalue()),
                embed=True,sizing_mode='stretch_width')
            self.download_slot.objects=[self.download]
            self.status.object='Completed process '+run_id+'. Actual DWSIM solve, save/reload and balance checks passed. COD retention is illustrative; sustainable flux remains unavailable.'
            self.status.alert_type='success'
        except Exception as e:
            prior=' The previous successful result and download remain available.' if self.completed is not None else ' No successful process result yet.'
            self.status.object='Run failed: '+html.escape(str(e))+prior;self.status.alert_type='danger'
        finally:
            self.run_button.disabled=not self.ready;self.run_button.loading=False
            for w in self.controls.values():w.disabled=False

    def _metrics(self,case,run_id):
        a=case['translation']['axial'];f=case['translation']['feed'];v=case['validation'];p=case['dwsim']['powers']
        comparison=case.get('verified_baseline_comparison')
        saving='' if comparison is None else f"**Saving versus the verified straight-support case with the same inputs:** {comparison['saving_fraction']:.2%}  \n"
        return (f"**Completed `{run_id}` · {case['request']['flux_LMH']:g} LMH · {case['request']['support']}**\n\n"+saving+
            f"**Pump shaft energy:** {v['specific_shaft_energy_kWh_m3_reference_volume']:.6g} kWh/m³ · "
            f"Fresh-feed pump {p['feed_W']:.5g} W · circulation pump {p['circulation_W']:.5g} W  \n"
            f"**Mean TMP:** {a['mean_TMP_Pa']/1000:.4g} kPa · axial channel drop {a['channel_drop_Pa']:.4g} Pa · "
            f"fixed effective non-support resistance share {a['selective_and_other_fixed_fraction']:.1%}  \n"
            f"**Inlet Re:** {a['inlet_Re']:g} · mean inlet speed {a['inlet_velocity_m_s']:.5g} m/s · assumed viscosity {a['viscosity_Pa_s']*1000:.5g} mPa s  \n"
            f"**Maximum local flux:** {case['translation']['maximum_local_flux_LMH']:.4g} LMH · **Sustainable flux: unavailable** at this geometry/speed  \n"
            f"**COD (external Julia tracer):** fresh {f['fresh_COD_kg_m3']*1000:.4g}, module inlet {f['module_inlet_COD_kg_m3']*1000:.4g}, "
            f"permeate {f['permeate_COD_kg_m3']*1000:.4g}, purge {f['retentate_COD_kg_m3']*1000:.4g} mg/L  \n"
            f"Overall water recovery 50% · recycle iterations {case['dwsim']['recycle']['iterations_taken']} · "
            f"maximum relative mass residual {v['max_mass_relative_residual']:.2e} · enthalpy residual / shaft power {v['max_enthalpy_relative_residual']:.2e}. "
            'Cooling heat is exported separately; electricity for cooling is not included.')

    def _baseline_comparison(self,case):
        from pore2process_demo.io import sha
        latest=ROOT/'results/star_channel/process/latest.json'
        if not latest.is_file():return None
        path=Path(read(latest)['summary']);study=read(path)
        if study.get('native_sha256')!=sha(self.native_path) or study.get('source_sha256')!=sha(self.source_path):return None
        request={**case['request'],'support':'straight'}
        baseline=next((c for c in study['cases'] if c['request']==request),None)
        if baseline is None:return None
        power=baseline['validation']['total_shaft_power_W']
        return {'matrix_summary_path':str(path),'matrix_summary_sha256':sha(path),
            'baseline_case_id':baseline['id'],'baseline_total_shaft_power_W':power,
            'saving_fraction':1-case['validation']['total_shaft_power_W']/power}

    def close(self):
        if self.plotter is not None:self.plotter.close();self.plotter=None
