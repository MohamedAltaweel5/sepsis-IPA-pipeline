"""Scope definition and benchmark timing (spec section 7, 11, 12).

The pipeline operates on the ER-CORE scope, justified by benchmark coverage: the
SSC governs only the ER window. This module produces:

  * the ER-window bound (first Admission NC/IC, exclusive) per case;
  * the ER-core scoped log (6 activities within the window);
  * variant statistics on the scoped log;
  * the benchmark timing table (time-to-antibiotics, compliance, emergency phase).

All timing uses INTER-EVENT intervals from ER Registration (spec C7: activity
durations cannot be measured - single timestamp per event).
"""
from __future__ import annotations
import pandas as pd
import numpy as np

from .constants import (
    CASE, ACT, TS, ER_CORE_ACTIVITIES, ADMISSION_ACTIVITIES, REFERENCE_ACTIVITY,
)


# --- Temporal + activity scoping (section 7.1, 7.2) -------------------------

def er_window_bound(g: pd.DataFrame):
    """First Admission (NC or IC) timestamp for a case, or NaT if never admitted."""
    adm = g.loc[g[ACT].isin(ADMISSION_ACTIVITIES), TS]
    return adm.iloc[0] if len(adm) else pd.NaT


def scope_er_core(df: pd.DataFrame) -> pd.DataFrame:
    """Return the ER-core scoped log: events before first admission (exclusive)
    whose activity is one of the 6 ER-core activities."""
    # Per-case ER-window bound = first admission timestamp (NaT if never admitted).
    adm = df[df[ACT].isin(ADMISSION_ACTIVITIES)]
    bound = adm.groupby(CASE, sort=False)[TS].min()
    b = df[CASE].map(bound)
    within_window = b.isna() | (df[TS] < b)
    return df[within_window & df[ACT].isin(ER_CORE_ACTIVITIES)].copy()


def variant_stats(scoped: pd.DataFrame) -> dict:
    """Variant count and dominant variant on the scoped log."""
    counts: dict[tuple, int] = {}
    for _cid, g in scoped.groupby(CASE, sort=False):
        v = tuple(g.sort_values(TS)[ACT])
        counts[v] = counts.get(v, 0) + 1
    ordered = sorted(counts.items(), key=lambda kv: -kv[1])
    n_cases = scoped[CASE].nunique()
    top_v, top_n = ordered[0]
    return {
        "n_events": int(len(scoped)),
        "n_variants": len(counts),
        "n_cases": int(n_cases),
        "dominant_variant": list(top_v),
        "dominant_count": int(top_n),
        "dominant_pct": round(100 * top_n / n_cases, 1),
        "top_variants": [
            {"variant": list(v), "count": int(n), "pct": round(100 * n / n_cases, 1)}
            for v, n in ordered[:10]
        ],
    }


# --- Benchmark timing (section 11, 12) --------------------------------------

def _first_ts(g, act):
    s = g.loc[g[ACT] == act, TS]
    return s.iloc[0] if len(s) else pd.NaT


def time_to_antibiotics(df: pd.DataFrame) -> pd.Series:
    """Hours from first ER Registration to first IV Antibiotics, per treated case.

    Computed on the FULL log (not the ER-window-bounded log): 'treated' means the
    case received antibiotics at all (823 cases). Negative-interval cases (GW, LZ)
    are retained here so the compliance denominator is 823, matching the design
    checkpoint; detect.py flags/handles them per the strength weighting.
    """
    out = {}
    for cid, g in df.groupby(CASE, sort=False):
        reg = _first_ts(g, REFERENCE_ACTIVITY)
        abx = _first_ts(g, "IV Antibiotics")
        if pd.notna(reg) and pd.notna(abx):
            out[cid] = (abx - reg).total_seconds() / 3600.0
    return pd.Series(out, name="tta_hours")


def emergency_phase_duration(df: pd.DataFrame) -> pd.Series:
    """Hours from first ER Registration to first Admission (NC or IC), per case.

    The ER-window performance metric (section 7.1). Full case duration is NOT used
    (dominated by ward stay the benchmark does not govern). Non-negative only.
    """
    out = {}
    for cid, g in df.groupby(CASE, sort=False):
        reg = _first_ts(g, REFERENCE_ACTIVITY)
        adm = er_window_bound(g)
        if pd.notna(reg) and pd.notna(adm):
            h = (adm - reg).total_seconds() / 3600.0
            if h >= 0:
                out[cid] = h
    return pd.Series(out, name="ephase_hours")


def benchmark_summary(df: pd.DataFrame) -> dict:
    """The baseline benchmark table (spec section 12 headline numbers)."""
    tta = time_to_antibiotics(df)
    ep = emergency_phase_duration(df)
    n_treated = len(tta)
    return {
        "cases_with_antibiotics": int(n_treated),
        "untreated_cases": int(df[CASE].nunique() - n_treated),
        "tta_median_h": round(float(tta.median()), 3),
        "tta_mean_h": round(float(tta.mean()), 3),
        "tta_p25_h": round(float(tta.quantile(0.25)), 3),
        "tta_p75_h": round(float(tta.quantile(0.75)), 3),
        "tta_p90_h": round(float(tta.quantile(0.90)), 3),
        "tta_max_h": round(float(tta.max()), 3),
        "compliance_1h_n": int((tta <= 1).sum()),
        "compliance_1h_pct": round(100 * float((tta <= 1).mean()), 1),
        "compliance_3h_n": int((tta <= 3).sum()),
        "compliance_3h_pct": round(100 * float((tta <= 3).mean()), 1),
        "ephase_n": int(len(ep)),
        "ephase_median_h": round(float(ep.median()), 3),
        "ephase_mean_h": round(float(ep.mean()), 3),
    }
