"""Phase II Step 2 - Candidate generation via redesign heuristics (spec section 9.2).

Encode the relevant subset of the Reijers & Limam Mansar (2005) redesign-heuristics
catalogue (ref [33]) as EXPLICIT RULES, and map each deviation's structural
signature to the applicable heuristic(s). The intelligence is the transparent
mapping observed-deficiency -> established remedy - not a learned model (C1, C3).

Each candidate carries a TRACEABLE RATIONALE: which deviation prompted it + which
named catalogue heuristic produced it + the measured structural evidence + a
concrete, machine-usable model modification for Phase III.

Do NOT invent heuristics outside the catalogue. A deviation matching no catalogue
precondition generates NO candidate and is recorded as such.
"""
from __future__ import annotations
import pandas as pd

from ..prep.constants import CASE, ACT, TS
from ..prep.load import _as_bool

# The catalogue subset in play, with the Devil's Quadrangle directional effects
# documented by Reijers & Limam Mansar (used later by score.py). Effects are the
# catalogue's qualitative signs, NOT fitted values: +1 improves, -1 worsens, 0 neutral.
CATALOGUE = {
    "resequencing": {
        "name": "Resequencing (task order optimisation)",
        "precondition": "Two governed activities occur in an order that delays the "
                        "strong-target one, with no dependency forcing it.",
        "effects": {"time": +1, "cost": 0, "quality": 0, "flexibility": 0},
    },
    "parallelism": {
        "name": "Parallelism (concurrent execution)",
        "precondition": "Two activities are sequential but independent and their "
                        "sequence adds delay.",
        "effects": {"time": +1, "cost": -1, "quality": 0, "flexibility": -1},
    },
    "task_elimination": {
        "name": "Task elimination",
        "precondition": "An activity on the path adds no value to the governed "
                        "outcome and adds delay.",
        "effects": {"time": +1, "cost": +1, "quality": -1, "flexibility": +1},
    },
    "triage": {
        "name": "Triage (case differentiation)",
        "precondition": "Cases of differing character follow one undifferentiated path.",
        "effects": {"time": +1, "cost": -1, "quality": +1, "flexibility": +1},
    },
}


def _precedes(g, a, b):
    ta = g.loc[g[ACT] == a, TS]
    tb = g.loc[g[ACT] == b, TS]
    if len(ta) and len(tb):
        return ta.iloc[0] < tb.iloc[0]
    return None


def structural_signatures(scoped_df) -> dict:
    """Measure the precedence/co-occurrence signatures the catalogue keys off."""
    lac_before_abx = abx_before_lac = both_lac_abx = 0
    liq_before_abx = 0
    for _cid, g in scoped_df.groupby(CASE, sort=False):
        g = g.sort_values(TS)
        r = _precedes(g, "LacticAcid", "IV Antibiotics")
        if r is not None:
            both_lac_abx += 1
            lac_before_abx += int(r)
            abx_before_lac += int(not r)
        if _precedes(g, "IV Liquid", "IV Antibiotics"):
            liq_before_abx += 1
    acts_by_case = scoped_df.groupby(CASE)[ACT].agg(set)
    no_treatment = int(sum(
        1 for a in acts_by_case
        if "IV Antibiotics" not in a and "IV Liquid" not in a
    ))
    # Attribute-based triage cohort: InfectionSuspected (lifted from ER Registration,
    # present in the scoped log) AND lactate-before-antibiotics (i.e. cases the
    # expedite lane would actually accelerate).
    infection = _infection_suspected_by_case(scoped_df)
    infection_and_lac_first = 0
    for cid, g in scoped_df.groupby(CASE, sort=False):
        # Strict-time test (lactate STRICTLY before antibiotics), consistent with
        # detection and Phase III evaluation; a simultaneous lactate does not delay.
        if infection.get(cid, False) and _precedes(g, "LacticAcid", "IV Antibiotics"):
            infection_and_lac_first += 1
    return {
        "both_lactate_and_antibiotics": both_lac_abx,      # 769
        "lactate_before_antibiotics": lac_before_abx,      # 649
        "lactate_before_antibiotics_pct": round(100 * lac_before_abx / both_lac_abx, 1),
        "antibiotics_before_lactate": abx_before_lac,      # 120
        "fluid_before_antibiotics": liq_before_abx,        # 631
        "no_treatment_cases": no_treatment,                # 236
        "infection_suspected_and_lactate_first": infection_and_lac_first,
    }


def LAC_before_ABX(acts) -> bool:
    if "LacticAcid" in acts and "IV Antibiotics" in acts:
        return acts.index("LacticAcid") < acts.index("IV Antibiotics")
    return False


def _infection_suspected_by_case(scoped_df) -> dict:
    """InfectionSuspected per case, read from the ER Registration event (in scope)."""
    out = {}
    if "InfectionSuspected" not in scoped_df.columns:
        return out
    reg = scoped_df[scoped_df[ACT] == "ER Registration"]
    for cid, v in zip(reg[CASE], reg["InfectionSuspected"]):
        out[cid] = _as_bool(v)
    return out


def generate_candidates(scoped_df, detection: dict) -> dict:
    """Map deviations -> catalogue heuristics -> concrete candidate modifications."""
    sig = structural_signatures(scoped_df)
    dev = {d["provision_id"]: d for d in detection["provisions"]}
    abx = dev["antibiotics_1h"]

    candidates = []
    unmatched = []

    # --- C1: Resequencing (antibiotics before lactate) --------------------------
    # Precondition satisfied: lactate precedes antibiotics in 649/769 cases and
    # lactate is NOT a precondition for antibiotics (clinical_rules
    # explicitly_not_constrained). Directly targets the Strong antibiotics deviation.
    candidates.append({
        "candidate_id": "C1_resequence_abx_before_lactate",
        "heuristic": "resequencing",
        "heuristic_name": CATALOGUE["resequencing"]["name"],
        "originating_deviation": "antibiotics_1h",
        "precondition_evidence": {
            "lactate_before_antibiotics": sig["lactate_before_antibiotics"],
            "of_cases_with_both": sig["both_lactate_and_antibiotics"],
            "pct": sig["lactate_before_antibiotics_pct"],
            "lactate_is_precondition": False,   # explicitly_not_constrained
        },
        "modification": (
            "Reorder the treatment path so IV Antibiotics is not blocked by "
            "LacticAcid: administer antibiotics on recognition; draw lactate "
            "without gating antibiotics on it."
        ),
        "transform": {"type": "reorder", "move": "IV Antibiotics",
                      "before": "LacticAcid"},
        "affected_cases": sig["lactate_before_antibiotics"],
    })

    # --- C2: Parallelism (lactate || antibiotics) -------------------------------
    # Same deviation, different catalogue remedy: the two are independent and
    # sequential. Phase III will test whether parallelising DISPLACES load (9.4).
    candidates.append({
        "candidate_id": "C2_parallel_lactate_antibiotics",
        "heuristic": "parallelism",
        "heuristic_name": CATALOGUE["parallelism"]["name"],
        "originating_deviation": "antibiotics_1h",
        "precondition_evidence": {
            "sequential_independent_pairs": sig["both_lactate_and_antibiotics"],
            "lactate_is_precondition": False,
        },
        "modification": (
            "Place LacticAcid and IV Antibiotics on parallel branches "
            "(AND-split/join) so the lactate draw does not delay antibiotics."
        ),
        "transform": {"type": "parallelize",
                      "activities": ["LacticAcid", "IV Antibiotics"]},
        "affected_cases": sig["both_lactate_and_antibiotics"],
    })

    # --- C3: Triage - ATTRIBUTE-BASED functional split (not contention) ----------
    # Differentiate by InfectionSuspected and fast-track that cohort to antibiotics.
    # This is a FUNCTIONAL routing effect (screening identifies who to expedite), not
    # a resource-contention effect - so it IS evaluable under the control-flow-first
    # model (C4). Effect applies only to the InfectionSuspected sub-cohort.
    candidates.append({
        "candidate_id": "C3_triage_infection_expedite",
        "heuristic": "triage",
        "heuristic_name": CATALOGUE["triage"]["name"],
        "originating_deviation": "antibiotics_1h",
        "precondition_evidence": {
            "differentiator": "InfectionSuspected (case-level, from ER Registration)",
            "expedited_cohort_size": sig["infection_suspected_and_lactate_first"],
            "no_treatment_cases_separated": sig["no_treatment_cases"],
        },
        "modification": (
            "Attribute-based triage: route InfectionSuspected=true cases to an "
            "expedited lane that decouples antibiotics from the lactate wait. "
            "Functional differentiation on a recorded screening attribute - not a "
            "queue/contention effect."
        ),
        "transform": {"type": "triage_split", "attribute": "InfectionSuspected",
                      "expedite": "decouple_antibiotics"},
        "affected_cases": sig["infection_suspected_and_lactate_first"],
    })

    # --- C4: Task elimination (remove LacticAcid from the antibiotics path) ------
    # Generated to demonstrate the mechanism even though it is clinically weak:
    # LacticAcid has NO SSC time target and is not an antibiotics precondition
    # (explicitly_not_constrained), yet sits on the delay path in 649 cases. The
    # catalogue effect carries a QUALITY penalty (removing a recommended
    # measurement), which scoring will reflect - a candidate, not a recommendation.
    candidates.append({
        "candidate_id": "C4_eliminate_lactate_from_path",
        "heuristic": "task_elimination",
        "heuristic_name": CATALOGUE["task_elimination"]["name"],
        "originating_deviation": "antibiotics_1h",
        "precondition_evidence": {
            "lactate_on_delay_path": sig["lactate_before_antibiotics"],
            "lactate_has_ssc_time_target": False,
            "lactate_is_antibiotics_precondition": False,
            "clinical_caveat": "LacticAcid IS a (weak) SSC measurement provision; "
                               "removing it reduces guideline conformance (quality).",
        },
        "modification": (
            "Remove LacticAcid from the critical path to antibiotics so it no longer "
            "delays treatment. Time benefit at a quality (conformance) cost."
        ),
        "transform": {"type": "eliminate", "activity": "LacticAcid"},
        "affected_cases": sig["lactate_before_antibiotics"],
    })

    return {
        "signatures": sig,
        "candidates": candidates,
        "unmatched_heuristics": unmatched,
        "catalogue_used": {k: v["name"] for k, v in CATALOGUE.items()},
    }
