"""Pipeline orchestrator: prep -> Phase I -> Phase II -> Phase III -> report.

Currently implemented: PREP (load, quality, scope, benchmark) with all section-12
checkpoints. Later phases are appended here as they are built. Every computed
value is written to outputs/results.json (spec section 13.1) so the thesis
[XX%] placeholders can be filled mechanically - never by hand.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.prep.load import load_event_log, lift_case_attributes
from src.prep.quality import quality_report
from src.prep.scope import scope_er_core, variant_stats, benchmark_summary
from src.phase1_discovery.discover import run_discovery
from src.phase1_discovery.views import bottleneck_summary, save_views
from src.phase2_diagnosis.detect import detect_all
from src.phase2_diagnosis.heuristics import generate_candidates
from src.phase2_diagnosis.filter import filter_candidates
from src.phase2_diagnosis.score import score_and_rank
from src.phase3_simulation.fidelity import run_repeated_fidelity
from src.phase3_simulation.evaluate import evaluate_all
from src.report import figures as report_figures
from src.report import numbers as report_numbers


def run_prep(df=None) -> dict:
    if df is None:
        df = load_event_log()
    ca = lift_case_attributes(df)
    q = quality_report(df)
    scoped = scope_er_core(df)
    vs = variant_stats(scoped)
    bench = benchmark_summary(df)

    return {
        "log_totals": {
            "n_cases": int(df["case:concept:name"].nunique()),
            "n_events": int(len(df)),
            "n_activities": int(df["concept:name"].nunique()),
            "different_resources_total": 1,   # no org:resource column (spec 5, 12)
            "different_transitions_total": int(
                df["lifecycle:transition"].nunique()
            ) if "lifecycle:transition" in df.columns else 1,
        },
        "cohort": {
            "Hypotensie_true": int(ca["Hypotensie"].sum()),
            "InfectionSuspected_true": int(ca["InfectionSuspected"].sum()),
            "SIRSCriteria2OrMore_true": int(ca["SIRSCriteria2OrMore"].sum()),
        },
        "quality": q,
        "scope_er_core": vs,
        "baseline_benchmark": bench,
    }


def run_phase1(scoped_df) -> dict:
    """Phase I - discovery + human-interpretable views (spec section 8)."""
    disc = run_discovery(scoped_df)
    bs = bottleneck_summary(scoped_df)
    figs = save_views(scoped_df)
    return {
        "discovery": disc,
        "bottlenecks": bs["slowest_edges"],
        "view_figures": figs,
    }


def run_phase2(df, scoped_df) -> dict:
    """Phase II - benchmark diagnosis + candidate generation (spec section 9)."""
    detection = detect_all(df)
    generation = generate_candidates(scoped_df, detection)
    filtering = filter_candidates(generation["candidates"])
    scoring = score_and_rank(filtering["permitted"])
    return {
        "detection": detection,
        "generation": generation,
        "filtering": {
            "rules_applied": filtering["rules_applied"],
            "n_permitted": filtering["n_permitted"],
            "n_rejected": filtering["n_rejected"],
            "rejected": filtering["rejected"],
        },
        "permitted_candidates": filtering["permitted"],
        "scoring": scoring,
    }


def run_phase3(df, scoped_df, permitted_candidates, detection=None,
               k_splits=8, n_replications=3, eval_n=1000, eval_reps=20) -> dict:
    """Phase III - baseline fidelity gate (Stage 1) then candidate eval (Stage 2)."""
    from src.phase3_simulation.evaluate import build_eval_model
    from src.phase3_simulation import bias_check

    fidelity = run_repeated_fidelity(scoped_df, k_splits=k_splits,
                                     n_replications=n_replications)
    gate_passed = fidelity["gate"]["gate_passed"]
    # Reproducibility of the timing-model comparison (paper Table 11): the alternative
    # independent per-edge sampler under the SAME splits and seeds. Not used for gating.
    fidelity_edge = run_repeated_fidelity(scoped_df, k_splits=k_splits,
                                          n_replications=n_replications, mode="edge")
    timing_model_comparison = {
        "block": {"fidelity_mean": fidelity["fidelity_mean"],
                  "fidelity_std": fidelity["fidelity_std"],
                  "per_split": fidelity["per_split"],
                  "pass_fraction": fidelity["per_split_gate_pass_fraction"]},
        "edge": {"fidelity_mean": fidelity_edge["fidelity_mean"],
                 "fidelity_std": fidelity_edge["fidelity_std"],
                 "per_split": fidelity_edge["per_split"],
                 "pass_fraction": fidelity_edge["per_split_gate_pass_fraction"]},
        "note": "Gate uses the block model only; edge is reported for comparison.",
    }
    evaluation = evaluate_all(scoped_df, permitted_candidates, gate_passed,
                              n=eval_n, replications=eval_reps, detection=detection)
    # Validity control: selection-bias check on the empirical antibiotics-delay pool.
    eval_model = build_eval_model(scoped_df)
    bias_followup = bias_check.compute(df, scoped_df, eval_model, permitted_candidates)
    return {
        "stage1_fidelity": fidelity,
        "timing_model_comparison": timing_model_comparison,
        "stage2_evaluation": evaluation,
        "candidate_timing_sample_bias": {"followup": bias_followup},
    }


def tool_versions() -> dict:
    """Record exact tool versions for reproducibility (spec section 4)."""
    import platform
    v = {"python": platform.python_version()}
    for mod in ["pm4py", "scipy", "numpy", "pandas", "matplotlib", "yaml"]:
        try:
            m = __import__(mod)
            v[("PyYAML" if mod == "yaml" else mod)] = getattr(m, "__version__", "unknown")
        except Exception:
            v[mod] = "not installed"
    v["prosimos"] = "not used (control-flow-first DES built in-repo; see Phase III)"
    return v


def main():
    df = load_event_log()
    scoped = scope_er_core(df)
    results = {"prep": run_prep(df)}
    results["phase1"] = run_phase1(scoped)
    results["phase2"] = run_phase2(df, scoped)
    results["phase3"] = run_phase3(df, scoped, results["phase2"]["permitted_candidates"],
                                   detection=results["phase2"]["detection"])
    # Convenience keys the thesis placeholder map references directly (spec 13.1).
    b = results["prep"]["baseline_benchmark"]
    results["baseline_compliance_1h"] = b["compliance_1h_pct"]
    results["baseline_compliance_3h"] = b["compliance_3h_pct"]
    results["baseline_median_tta_h"] = b["tta_median_h"]
    f = results["phase3"]["stage1_fidelity"]
    results["fidelity_cfld"] = f["fidelity_mean"]["control_flow_tv_distance"]
    results["fidelity_cycletime_distance"] = f["fidelity_mean"]["cycle_time_relative_emd"]
    results["fidelity_gate_passed"] = f["gate"]["gate_passed"]

    # Reproducibility block (spec sections 4, 13.1): versions, seeds, replication counts.
    results["tool_versions"] = tool_versions()
    results["run_config"] = {
        "scope": "ER-core (6 activities)",
        "train_test_split_frac": 0.30,
        "fidelity_split_seeds": f["split_seeds"],
        "fidelity_replications_per_split": f["n_replications_per_split"],
        "candidate_eval_n_per_rep": results["phase3"]["stage2_evaluation"]
            .get("conditions", {}).get("n_per_rep"),
        "candidate_eval_replications": results["phase3"]["stage2_evaluation"]
            .get("conditions", {}).get("replications"),
        "candidate_eval_seeds": (results["phase3"]["stage2_evaluation"]["candidates"][0]
            ["seeds"] if results["phase3"]["stage2_evaluation"]["candidates"] else []),
    }

    # Report: monochrome IEEE figures (§13.2) + placeholder map (§13.3).
    results["report"] = {
        "figures": report_figures.generate_all(df, scoped, results),
    }
    results["report"].update(report_numbers.emit(results))

    out = ROOT / "outputs" / "results.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2, default=str))
    print(f"Wrote {out}")
    # Compact console summary (full detail in results.json).
    print(json.dumps({k: results[k] for k in results if not isinstance(results[k], dict)},
                     indent=2, default=str))


if __name__ == "__main__":
    main()
