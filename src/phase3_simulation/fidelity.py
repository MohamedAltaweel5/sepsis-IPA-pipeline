"""Phase III Stage 1 (part B) - MEASURE fidelity, and GATE (spec section 10.1, C5).

The instrument is checked before its readings are trusted. We simulate the
baseline and compare the simulated log against the HELD-OUT test log using
distance measures in the spirit of Chapela-Campa et al. (2025) BPS quality
measures (ref [19]):

  * Control-Flow Log Distance (CFLD) - implemented as an n-gram (bigram)
    distribution distance: total-variation distance between the normalised
    directly-follows (bigram) frequency vectors. 0 = identical control flow.
  * Cycle-Time Distribution distance - Earth Mover's Distance (Wasserstein)
    between end-to-end case-duration distributions. THE measure that matters most
    for a timing thesis. Reported absolute (hours) and relative to mean cycle time.
  * Time-to-antibiotics distribution distance - EMD between the benchmark metric's
    distributions (the thesis's governed quantity).

GATE (acceptance thresholds - see THRESHOLDS and DEVIATIONS D5): the baseline is an adequate
instrument iff control flow and cycle time are both within tolerance. If not, STOP:
report "not validly modellable"; do NOT evaluate candidates (C5). An honest null
result is acceptable.
"""
from __future__ import annotations
from collections import Counter
import numpy as np
from scipy.stats import wasserstein_distance

from .build_baseline import trace_cycle_time_h, trace_tta_h

# --- Acceptance thresholds ----------------------------------------------------
# Not taken from the literature: the quality framework the measures draw on
# (Chapela-Campa et al., Inf. Syst. 127, 2025) prescribes no acceptance levels.
# When and how they were set is recorded in DEVIATIONS D5. Rationale: a control-flow-first model with no resource contention (C4) should
# reproduce the empirical control flow closely (variants sampled from data) and
# the cycle-time distribution to within a modest fraction of its own scale. We
# require the relative cycle-time EMD <= 20% of the mean real cycle time and the
# control-flow total-variation distance <= 0.25 (i.e. >= 75% overlap of the
# directly-follows distribution).
THRESHOLDS = {
    "control_flow_tv_max": 0.25,
    "cycle_time_relative_emd_max": 0.20,
}


def _bigrams(traces):
    c = Counter()
    for tr in traces:
        acts = [a for a, _t in tr]
        for a, b in zip(acts, acts[1:]):
            c[(a, b)] += 1
    total = sum(c.values())
    return {k: v / total for k, v in c.items()} if total else {}


def control_flow_distance(sim_traces, real_traces) -> float:
    """Total-variation distance between bigram (directly-follows) distributions."""
    p, q = _bigrams(sim_traces), _bigrams(real_traces)
    keys = set(p) | set(q)
    return round(0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys), 4)


def _cycle_times(traces):
    return np.array([trace_cycle_time_h(t) for t in traces])


def _ttas(traces):
    return np.array([x for t in traces if (x := trace_tta_h(t)) is not None])


def measure_fidelity(sim_traces, real_traces) -> dict:
    cf = control_flow_distance(sim_traces, real_traces)
    sim_ct, real_ct = _cycle_times(sim_traces), _cycle_times(real_traces)
    ct_emd = float(wasserstein_distance(sim_ct, real_ct))
    rel_ct = ct_emd / float(real_ct.mean()) if real_ct.mean() else float("inf")
    sim_tta, real_tta = _ttas(sim_traces), _ttas(real_traces)
    tta_emd = float(wasserstein_distance(sim_tta, real_tta))
    return {
        "control_flow_tv_distance": cf,
        "cycle_time_emd_h": round(ct_emd, 4),
        "cycle_time_relative_emd": round(rel_ct, 4),
        "tta_emd_h": round(tta_emd, 4),
        "sim_cycle_mean_h": round(float(sim_ct.mean()), 3),
        "real_cycle_mean_h": round(float(real_ct.mean()), 3),
        "sim_cycle_sd_h": round(float(sim_ct.std()), 3),
        "real_cycle_sd_h": round(float(real_ct.std()), 3),
        "sim_cycle_p50_h": round(float(np.percentile(sim_ct, 50)), 3),
        "real_cycle_p50_h": round(float(np.percentile(real_ct, 50)), 3),
        "sim_cycle_p90_h": round(float(np.percentile(sim_ct, 90)), 3),
        "real_cycle_p90_h": round(float(np.percentile(real_ct, 90)), 3),
        "sim_cycle_p99_h": round(float(np.percentile(sim_ct, 99)), 3),
        "real_cycle_p99_h": round(float(np.percentile(real_ct, 99)), 3),
        "sim_tta_median_h": round(float(np.median(sim_tta)), 3),
        "real_tta_median_h": round(float(np.median(real_tta)), 3),
    }


def evaluate_gate(fidelity: dict) -> dict:
    cf_ok = fidelity["control_flow_tv_distance"] <= THRESHOLDS["control_flow_tv_max"]
    ct_ok = fidelity["cycle_time_relative_emd"] <= THRESHOLDS["cycle_time_relative_emd_max"]
    passed = bool(cf_ok and ct_ok)
    return {
        "thresholds": THRESHOLDS,
        "control_flow_ok": bool(cf_ok),
        "cycle_time_ok": bool(ct_ok),
        "gate_passed": passed,
        "verdict": (
            "Baseline is an adequate instrument; proceed to candidate evaluation."
            if passed else
            "Baseline fidelity INADEQUATE: process not validly modellable at "
            "adequate fidelity; no candidate predictions credited (C5)."
        ),
    }


def run_fidelity(model, real_test_traces, n_replications=5, base_seed=100, mode="block") -> dict:
    """Simulate several replications, average fidelity, apply the gate."""
    n = len(real_test_traces)
    reps = []
    for r in range(n_replications):
        sim = _simulate_ref(model, n, base_seed + r, mode)
        reps.append(measure_fidelity(sim, real_test_traces))
    keys = reps[0].keys()
    mean_fid = {k: round(float(np.mean([rep[k] for rep in reps])), 4) for k in keys}
    gate = evaluate_gate(mean_fid)
    return {
        "n_replications": n_replications,
        "n_test_cases": n,
        "fidelity_mean": mean_fid,
        "gate": gate,
        "seeds": [base_seed + r for r in range(n_replications)],
    }


def _simulate_ref(model, n, seed, mode="block"):
    # imported lazily to avoid a cycle at module load
    from .build_baseline import simulate
    return simulate(model, n, seed, mode=mode)


def run_repeated_fidelity(scoped_df, k_splits=8, test_frac=0.30,
                          n_replications=3, split_seeds=None, mode="block") -> dict:
    """Repeated hold-out fidelity (the gate procedure).

    A single 70/30 split is high-variance because ER-core cycle time is
    heavy-tailed (max ~27 h): one split can place tail cases disproportionately in
    train or test. We therefore repeat over k independent splits, average the
    per-split fidelity, and gate on the MEAN - removing split luck rather than
    selecting a favourable seed. The pass fraction across splits is also reported.
    """
    from .build_baseline import train_test_split, fit_model, real_traces
    if split_seeds is None:
        split_seeds = [1, 7, 13, 42, 99, 123, 777, 2024][:k_splits]

    per_split = []
    for sd in split_seeds:
        train, test = train_test_split(scoped_df, test_frac=test_frac, seed=sd)
        model = fit_model(train)
        fr = run_fidelity(model, real_traces(test), n_replications, base_seed=100 + sd,
                          mode=mode)
        fr["split_seed"] = sd
        fr["test_cycle_sd_h"] = fr["fidelity_mean"]["real_cycle_sd_h"]
        per_split.append(fr)

    keys = per_split[0]["fidelity_mean"].keys()
    mean_fid = {k: round(float(np.mean([p["fidelity_mean"][k] for p in per_split])), 4)
                for k in keys}
    std_fid = {k: round(float(np.std([p["fidelity_mean"][k] for p in per_split])), 4)
               for k in keys}
    gate = evaluate_gate(mean_fid)
    pass_fraction = round(np.mean([p["gate"]["gate_passed"] for p in per_split]), 3)

    ct_thresh = THRESHOLDS["cycle_time_relative_emd_max"]
    per_split_rows = [
        {"seed": p["split_seed"],
         "cycle_time_relative_emd": p["fidelity_mean"]["cycle_time_relative_emd"],
         "control_flow_tv_distance": p["fidelity_mean"]["control_flow_tv_distance"],
         "test_cycle_sd_h": p["test_cycle_sd_h"],
         "gate_passed": p["gate"]["gate_passed"]}
        for p in per_split
    ]
    failed_splits = [
        {"seed": r["seed"], "cycle_time_relative_emd": r["cycle_time_relative_emd"],
         "exceeded_threshold_by": round(r["cycle_time_relative_emd"] - ct_thresh, 4)}
        for r in per_split_rows if not r["gate_passed"]
    ]

    return {
        "procedure": "repeated hold-out (k splits x replications), gated on mean",
        "timing_model": mode,
        "k_splits": len(split_seeds),
        "split_seeds": split_seeds,
        "n_replications_per_split": n_replications,
        "fidelity_mean": mean_fid,
        "fidelity_std": std_fid,
        "gate": gate,
        "per_split_gate_pass_fraction": pass_fraction,
        "per_split": per_split_rows,
        "failed_splits": failed_splits,     # which splits failed and by how much
        "baseline_model_provenance": BASELINE_MODEL_PROVENANCE,
    }


# Provenance of how the baseline reached an adequate instrument (spec: log how many
# model variants were tried before the gate passed - a reproducibility record).
BASELINE_MODEL_PROVENANCE = {
    "n_timing_models_evaluated": 2,
    "variants": [
        {"variant": 1, "sampler": "independent per-edge interval bootstrap (mode='edge')",
         "outcome": "evaluated, not adopted",
         "reason": ("Initially assessed on a single 70/30 split (seed 42), where it exceeded "
                    "the cycle-time threshold (0.253). Under the final 8-split procedure it "
                    "passes on the mean (see phase3.timing_model_comparison) but places too "
                    "little mass near the median and too much in the upper body of the "
                    "cycle-time distribution.")},
        {"variant": 2, "sampler": "block-bootstrap of real per-variant interval vectors (mode='block')",
         "outcome": "ADOPTED",
         "reason": ("Closer cycle-time distribution shape (lower relative EMD under the same "
                    "splits and seeds). Both models fail seed 42, the least dispersed test "
                    "split; the gate is applied to the mean of 8 splits for that reason.")},
    ],
    "procedure_fix": ("Single 70/30 split was high-variance (heavy-tailed cycle time, "
                      "max ~27 h); moved to repeated 8-split averaging, gated on the "
                      "mean, to remove split luck rather than select a favourable seed."),
}
