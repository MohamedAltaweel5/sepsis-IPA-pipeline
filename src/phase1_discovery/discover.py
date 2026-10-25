"""Phase I - Discovery (spec section 8).

Discover the real ER-core process from the scoped log using the PM4Py Inductive
Miner (spec ref [5]), producing a SOUND Petri net. Soundness is required by
Phase III: a net with deadlocks or improper termination cannot be simulated (C5
depends on a simulable baseline).

Outputs (the 'structured model' half of the dual representation):
  * the sound Petri net (+ initial/final markings), persisted as PNML;
  * conformance (token-based fitness, precision) against the scoped log;
  * a structural + soundness summary.

No ML, no learned thresholds (C1). The Inductive Miner is a deterministic,
algorithmic discovery method.
"""
from __future__ import annotations
from ..prep.constants import rel
import warnings
from pathlib import Path

import pm4py

from ..prep.constants import CASE, ACT, TS, OUTPUTS_DIR


def discover_sound_net(scoped_df, noise_threshold: float = 0.0):
    """Inductive Miner -> (net, im, fm). noise_threshold=0.0 => IM (not IMf):
    guarantees a sound, perfectly-fitting model of the scoped behaviour."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        net, im, fm = pm4py.discover_petri_net_inductive(
            scoped_df, noise_threshold=noise_threshold,
            activity_key=ACT, timestamp_key=TS, case_id_key=CASE,
        )
    return net, im, fm


def conformance(scoped_df, net, im, fm) -> dict:
    """Token-based replay fitness + precision against the scoped log."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = pm4py.fitness_token_based_replay(
            scoped_df, net, im, fm,
            activity_key=ACT, timestamp_key=TS, case_id_key=CASE,
        )
        prec = pm4py.precision_token_based_replay(
            scoped_df, net, im, fm,
            activity_key=ACT, timestamp_key=TS, case_id_key=CASE,
        )
    return {
        "log_fitness": round(float(fit["log_fitness"]), 4),
        "percentage_fit_traces": round(float(fit["percentage_of_fitting_traces"]), 2),
        "precision": round(float(prec), 4),
    }


def net_structure(net, im, fm) -> dict:
    is_sound, _diag = pm4py.check_soundness(net, im, fm)
    return {
        "n_places": len(net.places),
        "n_transitions": len(net.transitions),
        "n_arcs": len(net.arcs),
        "n_silent_transitions": sum(1 for t in net.transitions if t.label is None),
        "sound": bool(is_sound),
    }


def run_discovery(scoped_df, out_dir: Path = OUTPUTS_DIR / "models") -> dict:
    """Discover, verify soundness, measure conformance, persist the model."""
    out_dir.mkdir(parents=True, exist_ok=True)
    net, im, fm = discover_sound_net(scoped_df)
    struct = net_structure(net, im, fm)
    conf = conformance(scoped_df, net, im, fm)

    pnml_path = out_dir / "er_core_petri_net.pnml"
    pm4py.write_pnml(net, im, fm, str(pnml_path))

    fig_path = OUTPUTS_DIR / "figures" / "fig6_1_discovered_petri_net.png"
    fig_path.parent.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            # Rendered at 600 dpi (MDPI figure requirement) via the graphviz visualizer.
            from pm4py.visualization.petri_net import visualizer as pn_vis
            gviz = pn_vis.apply(net, im, fm, parameters={"format": "png"})
            gviz.graph_attr["dpi"] = "600"
            pn_vis.save(gviz, str(fig_path))
        except Exception:
            fig_path = None

    return {
        "model_file": rel(pnml_path),
        "model_figure": rel(fig_path) if fig_path else None,
        "structure": struct,
        "conformance": conf,
        "discovery_algorithm": "Inductive Miner (IM, noise_threshold=0.0)",
    }
