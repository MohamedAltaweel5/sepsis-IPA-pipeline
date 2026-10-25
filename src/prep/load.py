"""XES loading and normalisation (spec section 5, 6.5).

Loads the gzipped/plain XES event log into a canonical pandas DataFrame and lifts
the clinical attributes recorded only on each case's ER Registration event up to
case level. Uses PM4Py (spec ref [5]).

Anonymisation note (spec section 5): absolute timestamps are randomised but
intra-trace intervals are preserved. Every benchmark metric is an intra-trace
interval, so this does not affect the analysis. We never rely on absolute dates.
"""
from __future__ import annotations
import warnings
import pandas as pd

from .constants import DATA_XES, CASE, ACT, TS, COHORT_ATTRS


def load_event_log(path=DATA_XES) -> pd.DataFrame:
    """Read the XES log into a canonical, timestamp-parsed, sorted DataFrame.

    Timestamps are parsed to UTC (a fixed offset) so that interval arithmetic
    measures true elapsed time and is invariant to the source timezone.
    """
    import pm4py
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df = pm4py.read_xes(str(path))
    if not isinstance(df, pd.DataFrame):
        df = pm4py.convert_to_dataframe(df)
    df[TS] = pd.to_datetime(df[TS], utc=True)
    df = df.sort_values([CASE, TS], kind="stable").reset_index(drop=True)
    return df


def lift_case_attributes(df: pd.DataFrame, attrs=COHORT_ATTRS) -> pd.DataFrame:
    """Lift attributes from each case's ER Registration event to case level.

    Returns one row per case indexed by case id. Clinical attributes (SIRS
    criteria, checkboxes, triage fields) are recorded at event level but only on
    the ER Registration event of each case (exactly 1,050 occurrences), so they
    are genuinely case-level and are read from that event (spec section 5, 6.5).
    """
    from .constants import REFERENCE_ACTIVITY
    reg = df[df[ACT] == REFERENCE_ACTIVITY].drop_duplicates(subset=[CASE], keep="first")
    keep = [c for c in attrs if c in reg.columns]
    out = reg.set_index(CASE)[keep].copy()
    # Normalise the boolean-ish attributes to real booleans (XES yields strings/objects).
    for c in keep:
        out[c] = out[c].map(_as_bool)
    return out


def _as_bool(v) -> bool:
    if isinstance(v, bool):
        return v
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return False
    return str(v).strip().lower() in {"true", "1", "yes"}
