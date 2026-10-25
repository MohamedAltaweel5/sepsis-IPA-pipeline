"""Validation checkpoints (spec section 12) - the correctness test.

Every number here was computed from the real log during design. If any assertion
fails, data handling is wrong and must be fixed before trusting anything
downstream. Run:  python -m pytest tests/test_checkpoints.py -v
(or plain:  python tests/test_checkpoints.py)
"""
from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.prep.load import load_event_log, lift_case_attributes
from src.prep.quality import quality_report
from src.prep.scope import scope_er_core, variant_stats, benchmark_summary
from src.prep.constants import CASE, ACT

import pytest


@pytest.fixture(scope="module")
def df():
    return load_event_log()


# --- Log level --------------------------------------------------------------

def test_log_totals(df):
    assert df[CASE].nunique() == 1050
    assert len(df) == 15214
    assert df[ACT].nunique() == 16


def test_single_resource_and_transition(df):
    # different_resources_total = 1: no individual org:resource column exists.
    assert "org:resource" not in df.columns
    # different_transitions_total = 1: a single lifecycle transition value.
    if "lifecycle:transition" in df.columns:
        assert df["lifecycle:transition"].nunique() == 1


# --- Cohort attributes (at ER Registration) ---------------------------------

def test_cohort_attributes(df):
    ca = lift_case_attributes(df)
    assert int(ca["Hypotensie"].sum()) == 54
    assert int(ca["InfectionSuspected"].sum()) == 848
    assert int(ca["SIRSCriteria2OrMore"].sum()) == 853


# --- Scoping ----------------------------------------------------------------

def test_er_core_scope(df):
    scoped = scope_er_core(df)
    vs = variant_stats(scoped)
    assert scoped[ACT].nunique() == 6
    assert vs["n_events"] == 5542
    assert vs["n_variants"] == 69
    assert vs["dominant_count"] == 289
    assert vs["dominant_pct"] == 27.5
    assert vs["dominant_variant"] == [
        "ER Registration", "ER Triage", "ER Sepsis Triage",
        "LacticAcid", "IV Liquid", "IV Antibiotics",
    ]


def test_quality(df):
    q = quality_report(df)
    assert q["negative_interval_cases"] == ["GW", "LZ"]
    assert q["unknown_group_events"] == 294
    assert q["unknown_group_activities"] == ["Return ER"]


# --- Baseline benchmark (the headline numbers) ------------------------------

def test_benchmark_headline(df):
    b = benchmark_summary(df)
    assert b["cases_with_antibiotics"] == 823
    assert b["untreated_cases"] == 227
    # median 1.88 h, mean 2.19 h
    assert b["tta_median_h"] == pytest.approx(1.88, abs=0.01)
    assert b["tta_mean_h"] == pytest.approx(2.19, abs=0.01)
    # THE single most important number: 284/823 = 34.5%
    assert b["compliance_1h_n"] == 284
    assert b["compliance_1h_pct"] == pytest.approx(34.5, abs=0.1)
    assert b["compliance_3h_n"] == 595
    assert b["compliance_3h_pct"] == pytest.approx(72.3, abs=0.1)
    # Emergency phase n = 810 (median/mean reported; see prep notes for the
    # ~0.15 h gap vs the design's 3.26/3.54 in this secondary stat).
    assert b["ephase_n"] == 810


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
