"""
    crossflow_recycle_load(tmp_Pa, permeate_flow_m3_s, recovery,
        module_to_permeate_ratio, channel_resistance_Pa_s_m3,
        feed_pump_efficiency, circulation_pump_efficiency)

Reference-density hydraulic closure for a pressurized clean-water recycle loop.
The feed pump raises only fresh feed through TMP. A second pump raises the
full module flow through a separate linear upstream channel loss R*Q_module.
Retentate splits into a depressurized purge and a high-pressure return.
Recovery = permeate/fresh feed; module/permeate ratio is independent of it.
All Q values use the same fresh-feed density; DWSIM reports actual flashed
volumes and property-corrected pump powers separately. This does not predict
rejection, axial membrane pressure variation or a pump-curve operating point.
"""
function crossflow_recycle_load(tmp_Pa, permeate_flow_m3_s, recovery,
        module_to_permeate_ratio, channel_resistance_Pa_s_m3,
        feed_pump_efficiency, circulation_pump_efficiency)
    tmp=positive("transmembrane pressure_Pa",tmp_Pa)
    Qp=positive("permeate reference flow_m3_s",permeate_flow_m3_s)
    y=positive("overall permeate/fresh-feed recovery",recovery)
    y<1 || throw(ArgumentError("Recovery must be strictly between zero and one."))
    C=positive("module/permeate flow ratio",module_to_permeate_ratio)
    C>1/y || throw(ArgumentError("Module/permeate ratio must exceed 1/recovery for a positive recycle return."))
    R=finite("channel resistance_Pa_s_m3",channel_resistance_Pa_s_m3)
    R>=0 || throw(ArgumentError("Channel resistance must be nonnegative."))
    ef=positive("feed-pump efficiency",feed_pump_efficiency)
    ec=positive("circulation-pump efficiency",circulation_pump_efficiency)
    max(ef,ec)<=1 || throw(ArgumentError("Pump efficiencies must be <= 1."))
    Qf=positive("fresh-feed reference flow_m3_s",Qp/y)
    Qm=positive("module reference flow_m3_s",C*Qp)
    Qret=positive("module retentate reference flow_m3_s",Qm-Qp)
    Qreturn=positive("recycle return reference flow_m3_s",Qm-Qf)
    Qpurge=positive("purge reference flow_m3_s",Qf-Qp)
    dp=finite("channel pressure drop_Pa",R*Qm)
    Wf=positive("feed shaft power_W",Qf*tmp/ef)
    Wc=finite("circulation shaft power_W",Qm*dp/ec)
    W=positive("total shaft power_W",Wf+Wc)
    e=positive("specific shaft energy_kWh_m3",W/Qp/3.6e6)
    (; transmembrane_pressure_Pa=tmp,permeate_reference_flow_m3_s=Qp,
       fresh_feed_reference_flow_m3_s=Qf,module_reference_flow_m3_s=Qm,
       module_retentate_reference_flow_m3_s=Qret,recycle_reference_flow_m3_s=Qreturn,
       purge_reference_flow_m3_s=Qpurge,overall_recovery_fraction=y,
       module_to_permeate_ratio=C,channel_resistance_Pa_s_m3=R,channel_pressure_drop_Pa=dp,
       membrane_permeate_split_fraction=Qp/Qm,purge_retentate_split_fraction=Qpurge/Qret,
       loop_return_fraction=Qreturn/Qm,feed_shaft_power_W=Wf,circulation_shaft_power_W=Wc,
       total_shaft_power_W=W,specific_shaft_energy_kWh_m3=e,
       membrane_hydraulic_dissipation_W=Qp*tmp,purge_throttle_dissipation_W=Qpurge*tmp,
       channel_hydraulic_dissipation_W=Qm*dp,pump_loss_W=Wf*(1-ef)+Wc*(1-ec))
end
