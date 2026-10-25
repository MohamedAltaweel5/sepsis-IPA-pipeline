"""Emit the thesis placeholder map (spec section 13.3).

Maps each thesis [XX%]/[XX] placeholder to the computed value in results.json, so
placeholders can be filled MECHANICALLY - never by hand, never guessed. Writes a
human-readable Markdown table and a machine-readable JSON.
"""
from __future__ import annotations
from ..prep.constants import rel
import json
from pathlib import Path

from ..prep.constants import OUTPUTS_DIR


def build_placeholder_map(results: dict) -> list[dict]:
    b = results["prep"]["baseline_benchmark"]
    f = results["phase3"]["stage1_fidelity"]["fidelity_mean"]
    gate = results["phase3"]["stage1_fidelity"]["gate"]["gate_passed"]
    ev = results["phase3"]["stage2_evaluation"]["candidates"]
    top = next((c for c in ev if c["passed_stage2"]), None)

    rows = [
        ("§4.5.2 / abstract", "Baseline 1-hour compliance (real)",
         f"{b['compliance_1h_pct']}%", "prep.baseline_benchmark.compliance_1h_pct"),
        ("§4.5.2", "Baseline median time-to-antibiotics (real)",
         f"{b['tta_median_h']} h", "prep.baseline_benchmark.tta_median_h"),
        ("§4.5.2", "Baseline emergency-phase duration (real, median)",
         f"{b['ephase_median_h']} h", "prep.baseline_benchmark.ephase_median_h"),
        ("§2.10 / §4.5.1", "Fidelity - control-flow distance (CFLD)",
         f"{f['control_flow_tv_distance']}", "fidelity_cfld"),
        ("§2.10 / §4.5.1", "Fidelity - cycle-time relative EMD",
         f"{f['cycle_time_relative_emd']}", "fidelity_cycletime_distance"),
        ("§2.10 / §4.5.1", "Fidelity gate passed",
         f"{gate}", "fidelity_gate_passed"),
    ]
    if top:
        rows += [
            ("§4.5.2 / abstract", "Predicted 1-hour compliance (top candidate)",
             f"{top['predicted_compliance_1h']}%", f"phase3...{top['candidate_id']}.predicted_compliance_1h"),
            ("§4.5.2", "Predicted median TTA (top candidate)",
             f"{top['predicted_median_tta_h']} h", f"phase3...{top['candidate_id']}.predicted_median_tta_h"),
            ("§4.5.2", "Predicted scoped cycle time (top candidate) [NOT the §7.1 emergency phase]",
             f"{top['predicted_scoped_cycle_time_h']} h", f"phase3...{top['candidate_id']}.predicted_scoped_cycle_time_h"),
            ("§4.5.2 / abstract", "Composite risk reduction (top candidate)",
             f"{top['predicted_risk_reduction']}%", f"phase3...{top['candidate_id']}.predicted_risk_reduction"),
        ]
    return [{"thesis_location": a, "quantity": b_, "value": c, "results_key": d}
            for (a, b_, c, d) in rows]


def emit(results: dict, out_dir: Path = OUTPUTS_DIR) -> dict:
    rows = build_placeholder_map(results)
    (out_dir / "placeholder_map.json").write_text(json.dumps(rows, indent=2))

    lines = ["# Thesis placeholder map",
             "",
             "> Auto-generated from `results.json`. Every value is computed by the "
             "pipeline. **All Phase III (predicted/*fidelity*) values are "
             "MODEL-PREDICTED, not empirical (C6).**", "",
             "| Thesis location | Quantity | Value | results.json key |",
             "|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['thesis_location']} | {r['quantity']} | "
                     f"**{r['value']}** | `{r['results_key']}` |")
    (out_dir / "placeholder_map.md").write_text("\n".join(lines) + "\n")
    return {"placeholder_map_md": rel(out_dir / "placeholder_map.md"),
            "placeholder_map_json": rel(out_dir / "placeholder_map.json"),
            "n_placeholders": len(rows)}
