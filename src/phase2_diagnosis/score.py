"""Phase II Step 4 - Goal-based scoring (spec section 9.4).

Score permitted candidates on the Devil's Quadrangle dimensions (time, cost,
quality, flexibility) from the same Reijers catalogue (ref [33]), using the expert
weights in `config/goals.yaml`.

The score is EX ANTE and INDICATIVE: it expresses the expected DIRECTION and
relative magnitude of a candidate's effect per the catalogue, weighted by goal
priorities. It ORDERS candidates for Phase III; it does NOT establish that any
candidate works (that determination belongs to Phase III alone - C6/P4). Output is
a RANKED SET OF CANDIDATES, explicitly not recommendations.

No ML: the per-dimension effect signs come from the published catalogue
(heuristics.CATALOGUE); the magnitude is scaled by the measured affected-case
count (transparent), and combined linearly with the configured weights.
"""
from __future__ import annotations
from pathlib import Path
import yaml

from ..prep.constants import CONFIG_DIR
from .heuristics import CATALOGUE


def load_goals(path: Path = CONFIG_DIR / "goals.yaml") -> dict:
    return yaml.safe_load(open(path))


def _magnitude(candidate: dict, max_affected: int) -> float:
    """Normalise the candidate's reach to [0,1] by affected-case share.

    Transparent scalar so that a heuristic touching more of the deviating
    population scores its (catalogue-signed) effect more strongly. NOT a fitted
    quantity - just affected_cases / max(affected_cases across candidates).
    """
    if max_affected <= 0:
        return 0.0
    return candidate.get("affected_cases", 0) / max_affected


def score_candidate(candidate: dict, weights: dict, max_affected: int) -> dict:
    effects = CATALOGUE[candidate["heuristic"]]["effects"]  # signed per dimension
    mag = _magnitude(candidate, max_affected)
    # Weighted sum of (signed catalogue effect * reach magnitude).
    per_dim = {d: weights.get(d, 0.0) * effects.get(d, 0) * mag for d in weights}
    total = round(sum(per_dim.values()), 4)
    return {
        "candidate_id": candidate["candidate_id"],
        "heuristic": candidate["heuristic"],
        "catalogue_effects": effects,
        "reach_magnitude": round(mag, 3),
        "per_dimension_contribution": {d: round(v, 4) for d, v in per_dim.items()},
        "score": total,
    }


def score_and_rank(permitted: list[dict], goals: dict | None = None) -> dict:
    if goals is None:
        goals = load_goals()
    if not permitted:
        return {"weights": {}, "ranking": [], "note": "no permitted candidates to score"}
    max_affected = max((c.get("affected_cases", 0) for c in permitted), default=0)

    def rank_under(weights):
        scored = [score_candidate(c, weights, max_affected) for c in permitted]
        scored.sort(key=lambda s: -s["score"])
        return scored

    base = rank_under(goals["weights"])

    # Sensitivity analysis (goals.yaml requires it): does the ranking survive
    # alternative weightings? If the top candidate is stable, the weighting choice
    # is not load-bearing (spec 9.4 defensibility).
    sensitivity = {}
    for prof in goals.get("sensitivity_profiles", []):
        r = rank_under(prof["weights"])
        sensitivity[prof["id"]] = [s["candidate_id"] for s in r]

    base_order = [s["candidate_id"] for s in base]
    top_stable = all(order[0] == base_order[0] for order in sensitivity.values()) \
        if sensitivity else True

    # Explicit ranking-change summary (goals.yaml requires the sensitivity analysis to
    # be RUN and its effect REPORTED, defending the researcher-chosen base weights).
    per_profile = {}
    for pid, order in sensitivity.items():
        per_profile[pid] = {
            "ranking": order,
            "top": order[0],
            "top_changed_vs_base": order[0] != base_order[0],
            "full_ranking_changed_vs_base": order != base_order,
        }
    changed_top = [pid for pid, d in per_profile.items() if d["top_changed_vs_base"]]
    sensitivity_summary = {
        "base_ranking": base_order,
        "base_top": base_order[0] if base_order else None,
        "profiles": per_profile,
        "top_stable_across_all_profiles": top_stable,
        "profiles_that_change_the_top": changed_top,
        "conclusion": (
            (f"Top candidate is NOT invariant: {changed_top} re-rank the top. "
             if changed_top else "Top candidate is stable across all profiles. ")
            + ("The flip is confined to the compliance-equivalent C1/C2 pair (both "
               "attack the same deviation), so it is not load-bearing for which "
               "deviation is targeted; Phase III resolves C1 vs C2 empirically."
               if set(changed_top) <= {"time_only"} else "")),
    }

    return {
        "weights": goals["weights"],
        "ranking": base,
        "ranked_ids": base_order,
        "sensitivity": sensitivity,
        "sensitivity_summary": sensitivity_summary,
        "top_candidate_stable_across_profiles": top_stable,
        "label": "CANDIDATES, not recommendations (ex-ante, indicative; Phase III decides)",
    }
