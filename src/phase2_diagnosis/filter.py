"""Phase II Step 3 - Constraint filtering (spec section 9.3).

Check each candidate against the clinical/domain rules in
`config/clinical_rules.yaml`, encoded as verifiable precedence properties. Discard
any candidate that would violate a mandatory ordering; RECORD the violated rule
for each rejection. In a clinical domain this safety gate is NOT optional (C8).

A candidate is expressed as an ordering transform. We evaluate the ordering the
transform would impose and test it against each precedence rule:
  rule(before -> after) is violated iff the candidate forces `after` to precede
  `before`. Reordering/parallelising lactate vs antibiotics touches no rule
  (lactate order is explicitly_not_constrained); a hypothetical
  antibiotics-before-screening candidate WOULD be caught by
  `screening_before_antibiotics`.
"""
from __future__ import annotations
from pathlib import Path
import yaml

from ..prep.constants import CONFIG_DIR


def load_rules(path: Path = CONFIG_DIR / "clinical_rules.yaml") -> dict:
    return yaml.safe_load(open(path))


def _imposed_orderings(transform: dict) -> list[tuple[str, str]]:
    """Return the (X before Y) orderings the transform explicitly forces.

    Only 'reorder' imposes a new precedence. 'parallelize' REMOVES an ordering
    (imposes none). 'triage_split' partitions cases and imposes no intra-path
    reordering of governed activities.
    """
    t = transform.get("type")
    if t == "reorder":
        return [(transform["move"], transform["before"])]  # move BEFORE before
    return []


def check_candidate(candidate: dict, rules: dict) -> dict:
    """Return the candidate annotated with permit/reject + any violated rule."""
    imposed = _imposed_orderings(candidate["transform"])
    violations = []
    for rule in rules["rules"]:
        if rule.get("configurable") and rule.get("_disabled"):
            continue
        if rule["type"] != "precedence":
            continue
        before = rule["before"]
        afters = rule["after"] if isinstance(rule["after"], list) else [rule["after"]]
        # Violated iff the candidate forces one of `afters` to precede `before`.
        for (x, y) in imposed:  # x forced before y
            if x in afters and y == before:
                violations.append({
                    "rule_id": rule["id"],
                    "detail": f"candidate forces {x} before {y}, but rule requires "
                              f"{before} before {x}",
                })
    permitted = len(violations) == 0
    out = dict(candidate)
    out["permitted"] = permitted
    out["violated_rules"] = violations
    return out


def filter_candidates(candidates: list[dict], rules: dict | None = None) -> dict:
    if rules is None:
        rules = load_rules()
    checked = [check_candidate(c, rules) for c in candidates]
    permitted = [c for c in checked if c["permitted"]]
    rejected = [c for c in checked if not c["permitted"]]
    return {
        "rules_applied": [r["id"] for r in rules["rules"]],
        "permitted": permitted,
        "rejected": rejected,
        "n_permitted": len(permitted),
        "n_rejected": len(rejected),
    }
