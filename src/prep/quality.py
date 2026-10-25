"""Data-quality fixes (spec section 6). Every action is counted and reported;
nothing is silently dropped.

Design rules honoured:
  1. Negative treatment intervals -> EXCLUDE 2 cases (GW, LZ) from benchmark timing
     computation. IV Antibiotics recorded BEFORE ER Registration = timestamp error.
     They may remain in discovery. (They still count toward the treated denominator
     and, being <1 h, as compliant - see detect.py; that reproduces the design's
     284/823 = 34.5% checkpoint.)
  2. org:group == "?" -> 294 events, all Return ER. Systematic, out of scope. Noted only.
  3. Extreme case-duration outliers -> irrelevant once ER-scoped (section 7.1). Noted only.
  4. Untreated cases -> 227 (21.6%) have no IV Antibiotics. Clinical reality, NOT defects.
     Kept in discovery; excluded from time-to-antibiotics stats (undefined for them).
"""
from __future__ import annotations
import pandas as pd

from .constants import CASE, ACT, TS, GROUP, REFERENCE_ACTIVITY


def _first_ts(g: pd.DataFrame, act: str):
    s = g.loc[g[ACT] == act, TS]
    return s.iloc[0] if len(s) else pd.NaT


def find_negative_interval_cases(df: pd.DataFrame) -> list[str]:
    """Cases where first IV Antibiotics precedes first ER Registration (timestamp error)."""
    bad = []
    for cid, g in df.groupby(CASE, sort=False):
        reg = _first_ts(g, REFERENCE_ACTIVITY)
        abx = _first_ts(g, "IV Antibiotics")
        if pd.notna(reg) and pd.notna(abx) and (abx - reg).total_seconds() < 0:
            bad.append(cid)
    return sorted(bad)


def quality_report(df: pd.DataFrame) -> dict:
    """Compute the quality actions and return a JSON-serialisable report."""
    neg = find_negative_interval_cases(df)

    unknown_group = df[df[GROUP] == "?"] if GROUP in df.columns else df.iloc[0:0]
    unknown_group_acts = sorted(unknown_group[ACT].unique().tolist())

    cases_with_abx = set(df.loc[df[ACT] == "IV Antibiotics", CASE].unique())
    all_cases = set(df[CASE].unique())
    untreated = sorted(all_cases - cases_with_abx)

    return {
        "negative_interval_cases": neg,                       # -> ['GW', 'LZ']
        "negative_interval_excluded_from_timing": len(neg),   # -> 2
        "unknown_group_events": int(len(unknown_group)),      # -> 294
        "unknown_group_activities": unknown_group_acts,       # -> ['Return ER']
        "cases_with_antibiotics": len(cases_with_abx),        # -> 823
        "untreated_cases": len(untreated),                    # -> 227
        "total_cases": len(all_cases),                        # -> 1050
    }
