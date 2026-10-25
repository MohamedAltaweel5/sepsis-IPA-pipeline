"""Result figures (spec section 13.2) - STRICT MONOCHROME, IEEE style.

Thin black rules, grayscale fills + hatching to distinguish series, serif type, no
decorative colour - to match the existing conceptual figures. Each function
returns the saved path. Numbers come from the log / results.json; nothing here is
invented.
"""
from __future__ import annotations
from ..prep.constants import rel
from pathlib import Path
import json
import warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..prep.constants import OUTPUTS_DIR, CASE, ACT, TS

FIGDIR = OUTPUTS_DIR / "figures"

# --- IEEE monochrome style ---------------------------------------------------
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "black",
    "axes.grid": False,
    "figure.dpi": 150,
    "savefig.dpi": 600,   # MDPI: at least 600 dpi
    "savefig.bbox": "tight",
    "patch.edgecolor": "black",
    "hatch.linewidth": 0.6,
})
GRAYS = ["0.15", "0.45", "0.70", "0.88"]


def _tta_hours_full(df):
    """Registration -> first IV Antibiotics (hours), all treated cases (as prep)."""
    out = []
    for _cid, g in df.groupby(CASE, sort=False):
        g = g.sort_values(TS)
        reg = g.loc[g[ACT] == "ER Registration", TS]
        abx = g.loc[g[ACT] == "IV Antibiotics", TS]
        if len(reg) and len(abx):
            out.append((abx.iloc[0] - reg.iloc[0]).total_seconds() / 3600.0)
    return np.array(out)


def fig_tta_distribution(df, path=FIGDIR / "fig5_tta_distribution.png"):
    """Fig 5.x - baseline time-to-antibiotics ECDF with 1 h and 3 h lines."""
    tta_all = _tta_hours_full(df)                 # includes GW/LZ negatives
    x = np.sort(tta_all[tta_all >= 0])            # curve on non-negative support
    y = np.arange(1, len(x) + 1) / len(x)
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    ax.step(x, y, where="post", color="black", linewidth=1.2)
    for h, ls in [(1.0, "--"), (3.0, ":")]:
        ax.axvline(h, color="black", linewidth=0.8, linestyle=ls)
        frac = (tta_all <= h).mean()              # canonical: negatives count as <=1h
        ax.annotate(f"{h:.0f} h\n{frac*100:.1f}%", xy=(h, 0.05),
                    xytext=(h + 0.6, 0.12), fontsize=7)
    ax.set_xlim(0, 8)
    ax.set_xlabel("Time to antibiotics (h)")
    ax.set_ylabel("Cumulative fraction of treated cases")
    ax.set_title("Baseline time-to-antibiotics (ECDF)", fontsize=9)
    fig.savefig(path); plt.close(fig)
    return rel(path)


def fig_bottlenecks(bottlenecks, path=FIGDIR / "fig6_2_bottlenecks.png"):
    """Fig 6.2 - delay localisation: top directly-follows edges by mean delay."""
    edges = [b for b in bottlenecks if b["frequency"] >= 20][:7]  # drop rare/noise edges
    labels = [f"{b['from']}→{b['to']}" for b in edges]
    means = [b["mean_delay_h"] for b in edges]
    y = np.arange(len(edges))
    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    ax.barh(y, means, color=GRAYS[1], edgecolor="black", linewidth=0.8)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlabel("Mean inter-event delay (h)")
    ax.set_title("ER-window delay localisation", fontsize=9)
    for yi, m in zip(y, means):
        ax.text(m + 0.02, yi, f"{m:.2f}", va="center", fontsize=6.5)
    fig.savefig(path); plt.close(fig)
    return rel(path)


def fig_candidate_structure(path=FIGDIR / "fig6_3_candidate_structure.png"):
    """Fig 6.3 - before/after control-flow schematic for the top candidate (C1)."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.5, 2.8))
    seq_before = ["Sepsis\nTriage", "LacticAcid", "IV Liquid", "IV\nAntibiotics"]
    seq_after = ["Sepsis\nTriage", "IV\nAntibiotics", "LacticAcid", "IV Liquid"]

    def draw(ax, seq, title, highlight):
        ax.set_title(title, fontsize=8, loc="left")
        n = len(seq)
        for i, s in enumerate(seq):
            fc = "0.80" if s.replace("\n", " ") == highlight else "white"
            ax.add_patch(plt.Rectangle((i * 2, 0), 1.6, 1, fill=True,
                                       facecolor=fc, edgecolor="black", linewidth=0.9))
            ax.text(i * 2 + 0.8, 0.5, s, ha="center", va="center", fontsize=6.5)
            if i < n - 1:
                ax.annotate("", xy=(i * 2 + 2, 0.5), xytext=(i * 2 + 1.6, 0.5),
                            arrowprops=dict(arrowstyle="->", lw=0.8))
        ax.set_xlim(-0.3, n * 2); ax.set_ylim(-0.3, 1.3); ax.axis("off")

    draw(a1, seq_before, "(a) Baseline: antibiotics wait behind lactate", "IV Antibiotics")
    draw(a2, seq_after, "(b) C1 resequencing: antibiotics on recognition", "IV Antibiotics")
    fig.savefig(path); plt.close(fig)
    return rel(path)


def fig_simulation_results(evaluation, path=FIGDIR / "fig6_4_simulation_results.png"):
    """Fig 6.4 - (a) baseline vs candidates 1-h compliance; (b) the TRANSFER GAP:
    time-to-antibiotics reduction vs emergency-phase reduction per candidate."""
    cands = sorted(evaluation["candidates"], key=lambda c: c["candidate_id"])  # C1..C4
    ids = [c["candidate_id"].split("_")[0] for c in cands]
    base_comp = cands[0]["baseline_compliance_1h"]
    comp = [c["predicted_compliance_1h"] for c in cands]
    # Matched statistics (medians); means are in results.json and Table 13.
    tta_red = [c["transfer_gap"]["tta_median_reduction_pct"] for c in cands]
    ep_red = [c["transfer_gap"]["ct_median_reduction_pct"] for c in cands]
    x = np.arange(len(cands))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.8), constrained_layout=True)
    # (a) compliance
    ax1.axhline(base_comp, color="black", linewidth=0.9, linestyle="--",
                label=f"Baseline {base_comp:.1f}%")
    ax1.bar(x, comp, width=0.6, color=GRAYS[2], edgecolor="black",
            linewidth=0.8, hatch="///", label="Predicted")
    for xi, v in zip(x, comp):
        ax1.text(xi, v + 1, f"{v:.1f}", ha="center", fontsize=6.5)
    ax1.set_xticks(x); ax1.set_xticklabels(ids, fontsize=7)
    ax1.set_ylabel("Compliance ≤1 h (%)"); ax1.set_ylim(0, 100)
    ax1.set_title("(a) Predicted 1-hour compliance", fontsize=8.5)
    ax1.legend(fontsize=6.5, frameon=False, loc="upper right")
    # (b) transfer gap
    w = 0.38
    ax2.bar(x - w/2, tta_red, w, color=GRAYS[0], edgecolor="black",
            linewidth=0.8, label="Time-to-antibiotics")
    ax2.bar(x + w/2, ep_red, w, color=GRAYS[2], edgecolor="black",
            linewidth=0.8, hatch="\\\\", label="Scoped cycle time")
    ax2.set_xticks(x); ax2.set_xticklabels(ids, fontsize=7)
    ax2.set_ylabel("Reduction vs baseline (%)"); ax2.set_ylim(0, 100)
    ax2.set_title("(b) Transfer gap (medians)", fontsize=8.5)
    ax2.legend(fontsize=6.5, frameon=False, loc="upper right")
    fig.savefig(path); plt.close(fig)
    return rel(path)


def fig_fidelity_overlay(real_tta, sim_tta, path=FIGDIR / "fig6_5_fidelity.png"):
    """Fig 6.5 - simulated vs real time-to-antibiotics distribution overlay (ECDF)."""
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    for data, lbl, style in [(real_tta, "Real (held-out)", dict(color="black", lw=1.3, ls="-")),
                             (sim_tta, "Simulated", dict(color="0.45", lw=1.3, ls="--"))]:
        x = np.sort(data); y = np.arange(1, len(x) + 1) / len(x)
        ax.step(x, y, where="post", **style, label=lbl)
    ax.set_xlim(0, 8)
    ax.set_xlabel("Time to antibiotics (h)")
    ax.set_ylabel("Cumulative fraction")
    ax.set_title("Fidelity: simulated vs real TTA", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="lower right")
    fig.savefig(path); plt.close(fig)
    return rel(path)


def generate_all(df, scoped_df, results: dict) -> dict:
    """Produce every result figure. Recomputes distributions where needed."""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    from ..phase3_simulation.build_baseline import (
        fit_model, simulate, real_traces, trace_tta_h)

    paths = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        paths["fig5_tta"] = fig_tta_distribution(df)
        paths["fig6_2_bottlenecks"] = fig_bottlenecks(results["phase1"]["bottlenecks"])
        paths["fig6_3_candidate"] = fig_candidate_structure()
        paths["fig6_4_simulation"] = fig_simulation_results(
            results["phase3"]["stage2_evaluation"])
        # Fidelity overlay from a fresh baseline sim vs the real scoped log.
        model = fit_model(scoped_df)
        sim = simulate(model, len(scoped_df[CASE].unique()), seed=100)
        real = real_traces(scoped_df)
        sim_tta = np.array([x for t in sim if (x := trace_tta_h(t)) is not None])
        real_tta = np.array([x for t in real if (x := trace_tta_h(t)) is not None])
        paths["fig6_5_fidelity"] = fig_fidelity_overlay(real_tta, sim_tta)
    return paths
