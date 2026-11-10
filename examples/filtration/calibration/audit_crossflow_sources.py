"""Audit published crossflow data and figure readings, without fitting a model.

This checks source integrity and arithmetic. Physical inconsistencies are
reported as unresolved evidence limits, rather than corrected or hidden.
Only the Python standard library is required.
"""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path


from pore2process_demo.paths import ROOT
EXTRACTION = ROOT / "literature/crossflow_extraction.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def affine(pixel, anchors):
    p0, p1, v0, v1 = anchors
    if p0 == p1:
        raise ValueError("Axis anchors must have distinct pixel coordinates.")
    return v0 + (pixel - p0) * (v1 - v0) / (p1 - p0)


def reading(pixel, anchors, bound):
    # Independently perturb marker and both tick centres; these are reading
    # sensitivity bounds, not experimental confidence intervals.
    values = [
        affine(pixel + dp, [anchors[0] + d0, anchors[1] + d1, *anchors[2:]])
        for dp, d0, d1 in itertools.product((-bound, bound), repeat=3)
    ]
    return {"nominal": affine(pixel, anchors), "reading_interval": [min(values), max(values)]}


def audit(data, require_sources=False):
    source_checks = []
    for source in data["sources"]:
        for item in source["files"]:
            path = ROOT / item["path"]
            if path.exists():
                if sha256(path) != item["sha256"]:
                    raise ValueError(f"Source changed: {item['path']}")
                status = "hash_matches"
            elif require_sources:
                raise FileNotFoundError(path)
            else:
                status = "unavailable_locally; replaying retained extraction"
            source_checks.append({"path": item["path"], "status": status})

    table = data["garcia2008_table1"]
    area = table["filtration_area_m2"]
    target = data["illustrative_target_flux_LMH"]
    rows = []
    for row in table["rows"]:
        u, dp, flux, q = (row[k] for k in ("velocity_m_s", "pressure_loss_Pa", "critical_flux_LMH", "flow_L_h"))
        if not all(math.isfinite(v) and v > 0 for v in (u, dp, flux, q, area, target)):
            raise ValueError("Table values and scenario flux must be finite and positive.")
        q_si = q / 3_600_000
        flux_si = flux / 3_600_000
        energy = dp * q_si / (flux_si * area)
        rows.append({
            "row_id": row["row_id"],
            "velocity_m_s": u,
            "critical_flux_m_s": flux_si,
            "flow_m3_s_as_tabulated": q_si,
            "paper_energy_recomputed_J_m3": energy,
            "paper_energy_relative_difference": energy / row["reported_energy_J_m3"] - 1,
            "implied_total_open_area_m2_if_flow_is_bundle_total": q_si / u,
            "target_flux_LMH": target,
            "nominal_target_meets_tabulated_critical_flux": target <= flux,
            "fixed_target_circulation_energy_J_m3_using_table_flow": dp * q / (target * area),
        })

    figures = {}
    for key, figure in data["figure_readings"].items():
        figures[key] = [
            {"point_id": p["point_id"],
             "x": reading(p["pixel_xy"][0], figure["x_anchors"], figure["reading_bound_px"]),
             "critical_flux_LMH": reading(p["pixel_xy"][1], figure["y_anchors"], figure["reading_bound_px"])}
            for p in figure["markers"]
        ]

    # This combines dimensions reported for the same nominal module family
    # in different papers; it is an explicit consistency check, not a measured
    # replacement for Garcia's Q column.
    family = data["chiu2006_module"]
    family_area = family["hydraulic_diameter_m"] * family["filtration_area_m2"] / (4 * family["active_length_m"])
    for row in rows:
        row["family_area_to_table_implied_area_ratio"] = family_area / row["implied_total_open_area_m2_if_flow_is_bundle_total"]

    admissible = [r for r in rows if r["nominal_target_meets_tabulated_critical_flux"]]
    return {
        "scope": "Source arithmetic and reading audit only; no Julia/OpenLB/DWSIM model execution or validation.",
        "source_checks": source_checks,
        "garcia_rows": rows,
        "maximum_absolute_paper_energy_relative_difference": max(abs(r["paper_energy_relative_difference"]) for r in rows),
        "lowest_tabulated_admissible_velocity_m_s": min((r["velocity_m_s"] for r in admissible), default=None),
        "family_total_open_area_m2_from_Dh_S_and_length": family_area,
        "absolute_flow_basis": "unresolved; no automatic per-channel or factor-seven correction",
        "figure_readings": figures,
        "biological_maximum_digitized_critical_flux_LMH": max(p["critical_flux_LMH"]["nominal"] for p in figures["chiu2006_fig6"]),
        "biological_target_within_measured_critical_flux_range": any(target <= p["critical_flux_LMH"]["nominal"] for p in figures["chiu2006_fig6"]),
        "independent_validation": "None: Garcia Table 1 credits Chiu 2005 for critical flux. Biological feed is a different system.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--require-sources", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    data = json.loads(EXTRACTION.read_text())
    result = audit(data, args.require_sources)
    result["extraction_sha256"] = sha256(EXTRACTION)
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as handle:
            handle.write(text)
        print(f"Saved source audit to {args.output}")
        print(f"Paper energy arithmetic: maximum relative difference {result['maximum_absolute_paper_energy_relative_difference']:.3%}")
        print("Absolute bundle-flow basis remains unresolved; this is not a simulation validation.")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()
