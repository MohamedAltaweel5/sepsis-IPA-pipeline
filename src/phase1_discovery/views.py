"""Phase I - Human-interpretable views (spec section 8, dual output #2).

The interpretable half of the dual representation: a directly-follows graph with a
performance overlay, the variant list, and a performance/bottleneck summary. These
are a SUPPORTING capability, not a novelty claim. They also carry the domain
expert's entry point for supplying clinical rules and goal weights in deployment.

Performance = INTER-EVENT intervals (mean seconds on each directly-follows edge).
Activity durations are NOT computed - single timestamp per event (C7).
"""
from __future__ import annotations
import warnings
from pathlib import Path

import pm4py

from ..prep.constants import CASE, ACT, TS, OUTPUTS_DIR


def performance_dfg(scoped_df):
    """Directly-follows graph annotated with mean inter-event time per edge (seconds)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dfg, start, end = pm4py.discover_performance_dfg(
            scoped_df, activity_key=ACT, timestamp_key=TS, case_id_key=CASE,
        )
        freq_dfg, _s, _e = pm4py.discover_dfg(
            scoped_df, activity_key=ACT, timestamp_key=TS, case_id_key=CASE,
        )
    return dfg, freq_dfg, start, end


def bottleneck_summary(scoped_df, top_k: int = 8) -> dict:
    """Rank directly-follows edges by mean inter-event delay (the ER-window
    bottlenecks). Pure measurement; localises where time concentrates."""
    perf_dfg, freq_dfg, _s, _e = performance_dfg(scoped_df)
    edges = []
    for edge, stats in perf_dfg.items():
        # pm4py returns a stats dict per edge: mean/median/max/min/sum/stdev (seconds).
        edges.append({
            "from": edge[0],
            "to": edge[1],
            "mean_delay_h": round(float(stats["mean"]) / 3600.0, 3),
            "median_delay_h": round(float(stats["median"]) / 3600.0, 3),
            "frequency": int(freq_dfg.get(edge, 0)),
        })
    edges.sort(key=lambda e: -e["mean_delay_h"])
    return {"slowest_edges": edges[:top_k], "all_edges": edges}


def save_views(scoped_df, out_dir: Path = OUTPUTS_DIR / "figures") -> dict:
    """Persist the interpretable views (Petri net map + performance DFG)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    perf_dfg, freq_dfg, start, end = performance_dfg(scoped_df)
    paths = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        freq_path = out_dir / "view_dfg_frequency.png"
        pm4py.save_vis_dfg(freq_dfg, start, end, str(freq_path))
        paths["dfg_frequency"] = str(freq_path)
        try:
            perf_path = out_dir / "fig6_2_dfg_performance.png"
            pm4py.save_vis_performance_dfg(perf_dfg, start, end, str(perf_path))
            paths["dfg_performance"] = str(perf_path)
        except Exception as e:  # graphviz optional; don't fail the pipeline on a figure
            paths["dfg_performance_error"] = str(e)
    return paths
