"""
    once_through_filtration(tmp_Pa, permeate_flow_m3_s, recovery, efficiency,
                            channel_drop_Pa=0.0)

Incompressible clean-water hydraulic closure. Recovery = permeate/feed at one
common reference density (also the mass recovery). Both external products and
the fresh feed have the same pressure. A separate upstream lump represents
channel loss; the membrane sees uniform feed-side pressure = product P + TMP.
All retentate pressure is throttled away; no recycle or pressure recovery.
These are declared equipment assumptions, not additional OpenLB predictions.
"""
function once_through_filtration(tmp_Pa, permeate_flow_m3_s, recovery, efficiency,
                                 channel_drop_Pa=0.0)
    tmp=positive("transmembrane pressure_Pa",tmp_Pa)
    Qp=positive("permeate reference flow_m3_s",permeate_flow_m3_s)
    y=positive("permeate/feed recovery",recovery)
    y<1 || throw(ArgumentError("Recovery must be strictly between zero and one for two nonempty product streams."))
    eta=positive("pump efficiency",efficiency)
    eta<=1 || throw(ArgumentError("Pump efficiency must be <= 1."))
    channel=finite("separate channel pressure drop_Pa",channel_drop_Pa)
    channel>=0 || throw(ArgumentError("Channel pressure drop must be nonnegative."))
    Qf=positive("feed reference flow_m3_s",Qp/y)
    Qr=positive("retentate reference flow_m3_s",Qf-Qp)
    rise=positive("pump pressure rise_Pa",tmp+channel)
    hydraulic=positive("hydraulic power_W",Qf*rise)
    shaft=positive("shaft power_W",hydraulic/eta)
    specific=positive("specific shaft energy_kWh_m3",shaft/Qp/3.6e6)
    (; permeate_reference_flow_m3_s=Qp,feed_reference_flow_m3_s=Qf,
       retentate_reference_flow_m3_s=Qr,recovery_fraction=y,
       transmembrane_pressure_Pa=tmp,channel_pressure_drop_Pa=channel,
       pump_pressure_rise_Pa=rise,hydraulic_power_W=hydraulic,shaft_power_W=shaft,
       membrane_dissipation_W=Qp*tmp,retentate_throttle_dissipation_W=Qr*tmp,
       channel_dissipation_W=Qf*channel,pump_loss_W=shaft-hydraulic,
       specific_shaft_energy_kWh_m3=specific)
end
