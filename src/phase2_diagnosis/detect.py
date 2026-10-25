"""Phase II Step 1 - Benchmark deviation detection (spec section 9.1).

Compare observed timing/occurrence against the EXTERNAL SSC 2021 standard
(`config/ssc_benchmark.yaml`), NEVER against the process's own average (C2).

This is PURE MEASUREMENT. No training, no learned threshold, no accuracy to
report (nothing is predicted) - only deviation counts, magnitudes, compliance
rates, each weighted by SSC recommendation strength (C1).

Two metric types:
  * interval_from_reference: hours from ER Registration to first occurrence;
    deviation = case whose interval exceeds target_hours; magnitude = hours over.
  * occurrence: deviation = case in which the activity never occurs.

Denominator convention (matches the design checkpoint 284/823 = 34.5%): for
interval metrics the population is cases that HAVE the activity (so timing is
defined); negative-interval cases (GW, LZ) are retained and, being < target,
count as compliant. Untreated cases are reported separately, not as timing
deviations (their interval is undefined - spec 6.4).
"""
from __future__ import annotations
from pathlib import Path
import yaml

from ..prep.constants import CASE, ACT, TS, CONFIG_DIR, REFERENCE_ACTIVITY


def load_benchmark(path: Path = CONFIG_DIR / "ssc_benchmark.yaml") -> dict:
    return yaml.safe_load(open(path))


def _first_ts(g, act):
    s = g.loc[g[ACT] == act, TS]
    return s.iloc[0] if len(s) else None


def _reference_ts(g):
    return _first_ts(g, REFERENCE_ACTIVITY)


def detect_provision(df, provision: dict) -> dict:
    """Measure one SSC provision against the log. Returns a deviation record."""
    activity = provision["activity"]
    metric = provision["metric"]
    strength = provision["strength"]
    weight = float(provision["weight"])

    n_cases = df[CASE].nunique()

    if metric == "occurrence":
        cases_with = set(df.loc[df[ACT] == activity, CASE].unique())
        n_present = len(cases_with)
        n_missing = n_cases - n_present
        return {
            "provision_id": provision["id"],
            "activity": activity,
            "metric": metric,
            "strength": strength,
            "weight": weight,
            "population": n_cases,
            "compliant_n": n_present,
            "compliance_pct": round(100 * n_present / n_cases, 1),
            "deviation_n": n_missing,
            "weighted_deviation": round(weight * n_missing, 2),
            "target": "performed",
        }

    if metric == "interval_from_reference":
        target = float(provision["target_hours"])
        intervals = {}     # case -> hours (population: cases that have the activity)
        for cid, g in df.groupby(CASE, sort=False):
            ref = _reference_ts(g)
            act_ts = _first_ts(g, activity)
            if ref is not None and act_ts is not None:
                intervals[cid] = (act_ts - ref).total_seconds() / 3600.0
        pop = len(intervals)
        deviations = {c: h - target for c, h in intervals.items() if h > target}
        mags = sorted(deviations.values(), reverse=True)
        weighted = weight * sum(mags)  # strength-weighted total hours-over
        return {
            "provision_id": provision["id"],
            "activity": activity,
            "metric": metric,
            "strength": strength,
            "weight": weight,
            "target_hours": target,
            "population": pop,                       # e.g. antibiotics -> 823
            "compliant_n": pop - len(deviations),    # e.g. 284
            "compliance_pct": round(100 * (pop - len(deviations)) / pop, 1),  # 34.5
            "deviation_n": len(deviations),          # 539
            "deviation_median_h_over": round(sorted(mags)[len(mags) // 2], 3) if mags else 0.0,
            "deviation_max_h_over": round(max(mags), 3) if mags else 0.0,
            "weighted_deviation_h": round(weighted, 2),
            "untreated_note": (
                f"{n_cases - pop} cases never receive {activity}; interval undefined, "
                "excluded from timing (not counted as deviations)."
            ),
        }

    raise ValueError(f"Unknown metric {metric!r} for provision {provision['id']!r}")


def detect_all(df, benchmark: dict | None = None) -> dict:
    """Run detection for every governed provision; rank by weighted deviation."""
    if benchmark is None:
        benchmark = load_benchmark()
    records = [detect_provision(df, p) for p in benchmark["provisions"]]
    # Strong-recommendation breaches outweigh weak ones (spec 9.1).
    ranked = sorted(
        records,
        key=lambda r: -(r.get("weighted_deviation_h", 0) or r.get("weighted_deviation", 0)),
    )
    return {
        "standard": benchmark["meta"]["standard"],
        "provisions": records,
        "ranked_by_weighted_deviation": [r["provision_id"] for r in ranked],
        "undiagnosable": benchmark.get("undiagnosable", []),
    }
