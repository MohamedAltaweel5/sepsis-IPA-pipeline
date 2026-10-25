"""Phase III Stage 1 (part A) - Build the baseline BPS model (spec section 10.1).

Control-flow-first discrete-event model (C4), parameterised entirely from the log:

  * Control flow  = empirical VARIANT distribution (branching probabilities from
                    variant frequencies, spec 10.1).
  * Timing        = per-EDGE empirical inter-event interval distributions
                    (bootstrap resampling). We model INTER-EVENT INTERVALS, never
                    activity durations (C7: single timestamp per event).
  * Resource model= group-level only, no individual actors (C4). The log has
                    different_resources_total = 1, so there is no contention to
                    model; intra-trace timing is independent across cases. Case
                    arrival is therefore immaterial to the intra-trace metrics
                    (time-to-antibiotics, cycle time) and is held fixed - stated,
                    not hidden (spec 10.1 arrival caveat).

The model is fitted on a TRAIN split and later compared against a held-out TEST
split, so fidelity is not measured on the fitting data (spec 10.1).
"""
from __future__ import annotations
import numpy as np

from ..prep.constants import CASE, ACT, TS


def case_sequences(scoped_df):
    """dict case_id -> list[(activity, timestamp)] in temporal order (ER-core)."""
    seqs = {}
    for cid, g in scoped_df.groupby(CASE, sort=False):
        g = g.sort_values(TS)
        seqs[cid] = list(zip(g[ACT].tolist(), g[TS].tolist()))
    return seqs


def train_test_split(scoped_df, test_frac=0.30, seed=42):
    cases = np.array(sorted(scoped_df[CASE].unique()))
    rng = np.random.default_rng(seed)
    rng.shuffle(cases)
    n_test = int(round(len(cases) * test_frac))
    test = set(cases[:n_test])
    train = set(cases[n_test:])
    return (scoped_df[scoped_df[CASE].isin(train)].copy(),
            scoped_df[scoped_df[CASE].isin(test)].copy())


def fit_model(train_scoped) -> dict:
    """Fit the generative model: variant probabilities + per-edge interval pools."""
    seqs = case_sequences(train_scoped)

    variant_counts: dict[tuple, int] = {}
    variant_vectors: dict[tuple, list] = {}   # variant -> list of real interval-vectors
    edge_intervals: dict[tuple, list] = {}
    for _cid, ev in seqs.items():
        acts = tuple(a for a, _t in ev)
        variant_counts[acts] = variant_counts.get(acts, 0) + 1
        vec = [(tb - ta).total_seconds() for (a, ta), (b, tb) in zip(ev, ev[1:])]
        variant_vectors.setdefault(acts, []).append(np.asarray(vec, dtype=float))
        for (a, ta), (b, tb) in zip(ev, ev[1:]):
            edge_intervals.setdefault((a, b), []).append((tb - ta).total_seconds())

    total = sum(variant_counts.values())
    variants = [(list(v), n / total) for v, n in variant_counts.items()]
    # Global interval pool for any (rare) edge unseen after a transform.
    global_pool = [s for pool in edge_intervals.values() for s in pool]

    return {
        "variants": variants,                                  # (sequence, prob)
        "variant_vectors": variant_vectors,                    # joint timing per variant
        "edge_intervals": {k: np.asarray(v) for k, v in edge_intervals.items()},
        "global_interval_pool": np.asarray(global_pool),
        "n_train_cases": total,
    }


def _sample_interval(model, edge, rng) -> float:
    """Bootstrap-sample one inter-event interval (seconds) for an edge."""
    pool = model["edge_intervals"].get(edge)
    if pool is None or len(pool) == 0:
        pool = model["global_interval_pool"]
    return float(pool[rng.integers(len(pool))])


def simulate(model, n_cases, seed=0, mode="block") -> list:
    """Generate n_cases simulated ER-core traces.

    Returns list of traces; each trace is list[(activity, rel_seconds)] with the
    first (ER Registration) event at t=0 (absolute arrival immaterial - see module
    docstring).

    mode="block" (default, most faithful): sample a variant by its frequency, then
    reuse a real train-case interval-vector OF THAT VARIANT, preserving the joint
    within-case timing structure. mode="edge": independent per-edge bootstrap
    (kept for transformed adjacencies where a joint vector no longer applies).
    """
    rng = np.random.default_rng(seed)
    seqs = [tuple(v) for v, _p in model["variants"]]
    probs = np.array([p for _v, p in model["variants"]])
    idx = rng.choice(len(seqs), size=n_cases, p=probs)
    traces = []
    for i in idx:
        acts = seqs[i]
        if mode == "block":
            vecs = model["variant_vectors"][acts]
            vec = vecs[rng.integers(len(vecs))]
            times = np.concatenate([[0.0], np.cumsum(vec)]) if len(vec) else np.array([0.0])
            trace = [(acts[k], float(times[k])) for k in range(len(acts))]
        else:
            t = 0.0
            trace = [(acts[0], 0.0)]
            for a, b in zip(acts, acts[1:]):
                t += _sample_interval(model, (a, b), rng)
                trace.append((b, t))
        traces.append(trace)
    return traces


# --- Metric extractors shared by fidelity.py and evaluate.py ----------------

def trace_cycle_time_h(trace) -> float:
    return (trace[-1][1] - trace[0][1]) / 3600.0


def trace_tta_h(trace):
    """Hours from ER Registration (t=0) to first IV Antibiotics, or None."""
    for act, t in trace:
        if act == "IV Antibiotics":
            return t / 3600.0
    return None


def real_traces(scoped_df) -> list:
    """Real scoped cases as relative-second traces (for fidelity comparison)."""
    out = []
    for _cid, ev in case_sequences(scoped_df).items():
        t0 = ev[0][1]
        out.append([(a, (t - t0).total_seconds()) for a, t in ev])
    return out
