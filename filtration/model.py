"""SI boundary to the real Julia layer model; no Python numerical solver."""
from pore2process_demo.bridge import runtime

FIELDS = ("support_resistance_m_inv", "selective_resistance_m_inv", "fixed_resistance_m_inv",
          "total_resistance_m_inv", "support_pressure_drop_Pa", "selective_pressure_drop_Pa",
          "fixed_pressure_drop_Pa", "total_pressure_drop_Pa", "flow_m3_s", "shaft_power_W",
          "specific_shaft_energy_kWh_m3", "selective_resistance_fraction", "support_resistance_fraction")


def layer_load(k, settings, selective_resistance, viscosity):
    result = runtime().filtration_load(k, settings["support_thickness_m"], selective_resistance,
        viscosity, settings["flux_m_s"], settings["membrane_area_m2"],
        settings["pump_efficiency_fraction"],
        fixed_resistance_m_inv=settings["other_fixed_resistance_m_inv"])
    return {**{name:float(getattr(result,name)) for name in FIELDS},
            "selective_dominates":bool(result.selective_dominates),
            "selective_at_least_90_percent":bool(result.selective_at_least_90_percent)}
