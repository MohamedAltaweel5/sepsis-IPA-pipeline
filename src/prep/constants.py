"""Shared constants for the sepsis IPA pipeline.

Single source of truth for column names, the ER-core activity set, and the
scoping decisions from the implementation spec (sections 7 and 11). No logic
here - just declarations - so every module scopes identically.
"""
from __future__ import annotations
from pathlib import Path

# --- Repository paths -------------------------------------------------------
ROOT = Path(__file__).resolve().parents[2]
DATA_XES = ROOT / "data" / "Sepsis_Cases_-_Event_Log.xes"
CONFIG_DIR = ROOT / "config"
OUTPUTS_DIR = ROOT / "outputs"

# --- XES column names -------------------------------------------------------
CASE = "case:concept:name"
ACT = "concept:name"
TS = "time:timestamp"
GROUP = "org:group"

# --- Activities -------------------------------------------------------------
# The full 16-activity vocabulary (spec section 5).
ALL_ACTIVITIES = [
    "ER Registration", "ER Triage", "ER Sepsis Triage", "CRP", "Leucocytes",
    "LacticAcid", "IV Antibiotics", "IV Liquid", "Admission NC", "Admission IC",
    "Return ER", "Release A", "Release B", "Release C", "Release D", "Release E",
]

# ER-CORE scope (spec section 7.2): exactly these 6 activities.
# CRP and Leucocytes are EXCLUDED - repeated monitoring measurements, not SSC
# bundle elements; they cause the variant explosion. Declared scoping decision.
ER_CORE_ACTIVITIES = [
    "ER Registration", "ER Triage", "ER Sepsis Triage",
    "LacticAcid", "IV Liquid", "IV Antibiotics",
]

# Activities that close the ER window (spec section 7.1): the ER window is bounded
# at the FIRST of these (exclusive).
ADMISSION_ACTIVITIES = ["Admission NC", "Admission IC"]

# The reference event from which all benchmark intervals are measured.
REFERENCE_ACTIVITY = "ER Registration"

# Case-level clinical attributes lifted from the ER Registration event (spec 6.5).
COHORT_ATTRS = ["Hypotensie", "InfectionSuspected", "SIRSCriteria2OrMore"]
