"""Independent actual DWSIM water-pump compatibility probe."""
from datetime import datetime, timezone
import json
import math
from pore2process_demo.runtime.dwsim import OUTPUT, build_adapter, solve_case

def main():
    build_record=build_adapter(OUTPUT)
    cases = []
    for name, pressure in (("baseline", 1e6), ("double_pressure_rise", 1.9e6)):
        directory = OUTPUT / name
        directory.mkdir(parents=True, exist_ok=True)
        inputs = {"Temperature_K": 298.15, "Inlet_pressure_Pa": 1e5,
                  "Outlet_pressure_Pa": pressure, "Mass_flow_kg_s": 10.0,
                  "Pump_efficiency_fraction": 0.75,
                  "Scope": "Independent water-pump compatibility probe; no OpenLB-derived input yet."}
        print(f"Solving real DWSIM water-pump case: {name}...", flush=True)
        result=solve_case(directory,inputs)
        numeric = [value for value in result.values() if isinstance(value, (int, float))]
        if not all(math.isfinite(value) for value in numeric):
            raise RuntimeError("Non-finite DWSIM result.")
        checks = [abs(result["outlet_pressure_Pa"] - pressure) <= 1e-6,
                  abs(result["mass_balance_error_kg_s"]) <= 1e-9 * inputs["Mass_flow_kg_s"],
                  abs(result["energy_balance_error_W"]) <= 1e-6 * result["shaft_power_W"],
                  result["pump_power_relative_error"] <= 1e-9,
                  abs(result["energy_stream_power_W"] - result["shaft_power_W"]) <= 1e-6,
                  result["inlet_vapor_fraction"] == 0, result["outlet_vapor_fraction"] == 0,
                  result["save_reload_solve_checked"] is True,
                  (directory / "water_pump.dwxmz").is_file()]
        if not all(checks):
            raise RuntimeError(f"DWSIM acceptance failed for {name}: {checks}; inspect summary/log.")
        cases.append(result)
        print(f"{name}: {result['shaft_power_W']:.6f} W; checks passed", flush=True)
    ratio = cases[1]["shaft_power_W"] / cases[0]["shaft_power_W"]
    if not math.isclose(ratio, 2.0, rel_tol=1e-9):
        raise RuntimeError("Doubling pressure rise did not double pump power.")
    summary = {"timestamp_utc": datetime.now(timezone.utc).isoformat(), "cases": cases,
               "power_ratio": ratio, **build_record,
               "contract": "JSON field suffixes give SI units; efficiency is a fraction; internal kW/kJ converted at adapter boundary"}
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(f"Passed: {OUTPUT / 'summary.json'}", flush=True)


if __name__ == "__main__":
    main()
