# Thesis-writing handoff

Companion to `DEVIATIONS.md` and `outputs/results.json` (the complete verbatim
results). Everything below is computed by `run_pipeline.py`; nothing is hand-entered.

---

## 1. Figures produced

All in `outputs/figures/`. Strict monochrome, IEEE style (serif, thin black rules,
grayscale + hatching). Dimensions in pixels (300 dpi where matplotlib-drawn).

| # | File | Dimensions | Shows |
|---|---|---|---|
| Fig 5.x | `fig5_tta_distribution.png` | 1011×829 | Baseline time-to-antibiotics ECDF over treated cases, with 1 h (34.5%) and 3 h (72.3%) reference lines. |
| Fig 6.1 | `fig6_1_discovered_petri_net.png` | 1388×671 | Discovered sound Petri net (Inductive Miner) of the ER-core process. **pm4py default styling, not monochrome** — regenerate if strict IEEE needed. |
| Fig 6.2 | `fig6_2_bottlenecks.png` | 1287×875 | ER-window delay localisation: top directly-follows edges by mean inter-event delay (freq ≥ 20). |
| Fig 6.3 | `fig6_3_candidate_structure.png` | 873×754 | Before/after control-flow schematic for the top candidate (C1 resequencing): antibiotics moved ahead of lactate. |
| Fig 6.4 | `fig6_4_simulation_results.png` | 1802×793 | 2-panel: (a) predicted 1 h compliance C1–C4 vs baseline; (b) **transfer gap** — TTA reduction vs **scoped-cycle-time** reduction per candidate. |
| Fig 6.5 | `fig6_5_fidelity.png` | 1011×828 | Fidelity overlay: simulated vs real time-to-antibiotics ECDF. |
| (support) | `fig6_2_dfg_performance.png` | 1286×1011 | pm4py performance DFG (working view, not a thesis figure). |
| (support) | `view_dfg_frequency.png` | 1563×1046 | pm4py frequency DFG (working view, not a thesis figure). |

Note: `fig6_1` and the two `(support)` views are pm4py-rendered (colour). All
matplotlib figures (5.x, 6.2, 6.3, 6.4, 6.5) are strict monochrome.

---

## 2. Phase I — discovered model statistics

| Property | Value |
|---|---|
| Discovery algorithm | Inductive Miner (IM, noise_threshold = 0.0) — deterministic, no ML |
| Soundness | **true** |
| Token-replay fitness | **1.000** (100% fitting traces) |
| Precision | **0.8452** |
| Places / Transitions / Arcs | 15 / 16 / 40 |
| Silent transitions | 10 |
| Variant coverage of scoped log | 69 variants; dominant 289 (27.5%); top-10 variants cover ~85% |
| Model file | `outputs/models/er_core_petri_net.pnml` |

Interpretation: IM reads the six ER-core activities as a largely concurrent block
(big AND-split), which is why precision is 0.845 (it admits reorderings) while
fitness is perfect. The 69 variants come from genuine treatment-ordering freedom.

---

## 3. Phase II — full candidate output (successes AND the rest)

**Deviation detection (external SSC 2021):**

| Provision | Activity | Target | Strength | Compliance | Deviations | Weighted |
|---|---|---|---|---|---|---|
| antibiotics_1h | IV Antibiotics | ≤1 h | Strong | 34.5% (284/823) | 539 | 1128.5 h |
| fluids_3h | IV Liquid | ≤3 h | Weak | 78.1% (588/753) | 165 | 412.3 h |
| lactate_measured | LacticAcid | occurrence | Weak | 81.9% (860/1050) | 190 | 95.0 |
| sepsis_screening | ER Sepsis Triage | occurrence | Strong | 99.9% (1049/1050) | 1 | 1.0 |

Undiagnosable (no activity / no data): blood cultures, vasopressors, 30 mL/kg volume.

**Candidates generated (all originate from `antibiotics_1h`):**

| ID | Heuristic | Transform | Affected | Filter result |
|---|---|---|---|---|
| C1_resequence_abx_before_lactate | Resequencing | reorder: IV Antibiotics before LacticAcid | 649 | **permitted** (no rule violated) |
| C2_parallel_lactate_antibiotics | Parallelism | parallelize {LacticAcid, IV Antibiotics} | 769 | **permitted** |
| C3_triage_infection_expedite | Triage (attribute) | triage_split on InfectionSuspected → decouple antibiotics | 642 | **permitted** |
| C4_eliminate_lactate_from_path | Task elimination | eliminate LacticAcid | 649 | **permitted** |

**Rejections:** none (0 rejected). Rules applied: `registration_first`,
`triage_before_treatment`, `screening_before_antibiotics`. C1/C2/C4 touch lactate
ordering/presence, which the clinical rules leave `explicitly_not_constrained`
(lactate is not an antibiotics precondition), so none are rejected. `unmatched_heuristics`
is empty — all four catalogue heuristics matched a structural signature.

**Goal-based scores under ALL profiles** (score = Σ weight·effect·reach; reach =
affected/max_affected; effects are the Reijers Devil's-Quadrangle signs ±1):

| Cand | reach | base | time_only | balanced | quality_led | cost_constrained |
|---|---|---|---|---|---|---|
| C1 | 0.844 | 0.4220 | 0.8440 | 0.2110 | 0.2532 | 0.3376 |
| C2 | 1.000 | 0.3000 | **1.0000** | −0.2500 | 0.1000 | 0.0500 |
| C3 | 0.835 | **0.6679** | 0.8349 | 0.4174 | **0.6679** | 0.3339 |
| C4 | 0.844 | 0.3376 | 0.8440 | **0.4220** | 0.0000 | **0.4220** |

(C3 reach = 642/769 = 0.835 after the strict-time reconciliation; base score 0.6679. Regenerated from `results.json`; ordering unchanged — C3 is still base-top, the S1 story holds.)

**Ranking per profile (top → bottom):**
- base: **C3** > C1 > C4 > C2
- time_only: **C2** > C1 > C4 > C3
- balanced: **C4** > C3 > C1 > C2
- quality_led: **C3** > C1 > C2 > C4
- cost_constrained: **C4** > C1 > C3 > C2

Top candidate changes in 3 of 4 sensitivity profiles ⇒ the base weighting is **not
load-bearing**. (See surprise S1: the ex-ante top is C3, which is *not* the Phase III winner.)

---

## 4. Phase III — fidelity detail

**Pre-registered thresholds** (`fidelity.THRESHOLDS`, set before running):
`control_flow_tv_max = 0.25`, `cycle_time_relative_emd_max = 0.20`.

**Mean over 8 splits × 3 reps:** control-flow TV **0.0521**, cycle-time relative EMD
**0.1356**, TTA EMD 0.2287 h. sim/real cycle mean 1.900/1.980 h; sim/real TTA median
1.856/2.006 h. **Gate: PASSED.** Pass fraction 0.875 (7/8).

**Per-split (all 8):**

| Split seed | cycle-time rel-EMD | control-flow TV | gate |
|---|---|---|---|
| 1 | 0.1242 | 0.0580 | pass |
| 7 | 0.1269 | 0.0544 | pass |
| 13 | 0.1041 | 0.0428 | pass |
| **42** | **0.2016** | 0.0479 | **FAIL** (over by 0.0016) |
| 99 | 0.1668 | 0.0494 | pass |
| 123 | 0.1099 | 0.0608 | pass |
| 777 | 0.1163 | 0.0538 | pass |
| 2024 | 0.1353 | 0.0500 | pass |

Only seed 42 fails, by 0.0016 on cycle-time — the mean passes comfortably at 0.136.

**Baseline model variants tried before the gate passed: 2.**
1. Independent per-edge interval bootstrap → **REJECTED** (overestimated cycle time; single-split rel-EMD ~0.22–0.30).
2. Block-bootstrap of real per-variant interval vectors → **ACCEPTED** (rel-EMD 0.136). *More faithful instrument, not a threshold change.*
Procedure fix: single 70/30 split → repeated 8-split averaging (heavy-tailed cycle time, max ~27 h).

---

## 5. Risk-reduction formula (expression, not a number)

Composite, weight-weighted across **all governed provisions**, with a **FIXED
denominator** (so candidates are comparable) and **path-local** damage:

```
risk              = Σ_all-governed-provisions  weight · noncompliant_fraction
base_risk (FIXED) = 1.0·(1−0.345) + 0.5·(1−0.781) + 0.5·(1−0.819) + 1.0·(1−0.999) = 0.856
risk_reduction(%) = 100 · (base_risk − candidate_risk) / base_risk
```

Weights: `antibiotics_1h = sepsis_screening = 1.0`, `fluids_3h = lactate_measured =
0.5`. The candidate replaces its touched provisions' non-compliance; everything else
stays at baseline. **Path-local damage:** eliminating an activity removes it only on
the affected path, and residual compliance is recomputed from cases that still
perform it — not set to 0 globally.

| Cand | Provisions changed | base_risk | cand_risk | **risk ↓** |
|---|---|---|---|---|
| C1 / C2 | antibiotics 34.3%→80.3% | 0.856 | 0.397 | **53.59%** |
| C3 | antibiotics 34.3%→79.8% | 0.856 | 0.403 | 52.97% |
| C4 | antibiotics 34.3%→80.3% **+ lactate 81.9%→20.1%** (649 lost, 211 remain) | 0.856 | 0.706 | **17.49%** |

C4 eliminates LacticAcid on the antibiotics path (649 cases); 211 cases still
measure it, so lactate compliance falls to 20.1% (not 0%). The lactate loss
(weight 0.5) partly offsets the antibiotics gain, dropping C4 to **17.49%** — it now
agrees with Devil's Quadrangle `quality: −1`. Exposed as `risk_reduction_formula` +
`risk_components` (with `damaged_detail`) on every candidate. *History: antibiotics-only
gave C4=70%; per-candidate-denominator + global-0% gave 6.75%; final fixed-denominator
+ path-local gives 17.49%.*

---

## 6. What surprised me / what is fragile

**S1 — Ex-ante scoring and simulation disagree on the winner (the most important surprise).**
Phase II base ranking puts **C3 (triage) top** (0.6679), because the triage catalogue
signs reward quality+flexibility. But Phase III shows C3 **marginally worst** on the
governed metric (79.8% vs 80.3%). The ex-ante "best" is not the simulated best. This
is a live vindication of C5/the fidelity gate — scoring is a hypothesis, simulation
is the test — but it also means **do not report the Phase II ranking as a result**;
report it as an ordering that Phase III overturns.

**S2 — The four "distinct" candidates are one intervention in disguise — but the
composite risk separates them.** C1/C2/C4 are compliance-equivalent (all 80.316%,
identical to 3 d.p. — verified genuine via `saturation_check`: shared effective
affected set of 649; C2's extra 120 cases are antibiotics-already-first and
unchanged). Mechanically they all "stop antibiotics waiting behind the lactate
draw." **But the composite risk metric does separate C4** (17.5% vs 53.6%):
eliminating lactate damages a provision. So on 1-h compliance the diversity is
illusory; on the fixed-denominator composite it is not. Report both.

**S3 — The transfer gap is large and robust (your best finding).**
TTA improves ~79% but the ER-core **scoped cycle time** only ~19% (gap ~60 pp) across all candidates.
Fixing the governed step does **not** fix the pathway; the binding constraint
transfers to IV Liquid (~45%) / LacticAcid (~30%). This *confirms* the spec's §10.2
hypothesis that relieving one bottleneck displaces the constraint — but note it is a
**control-flow critical-path** shift, not a queue displacement (no capacity model; C4).

**S4 — Naming collision RESOLVED (was a real thesis-error risk).** (see DEVIATIONS D10)
`prep` emergency-phase (Registration→Admission = **3.22 h**) and the Phase III
end-to-end span (ER-core scoped cycle time = **1.92 h**) were both near-"emergency
phase". The Phase III quantity is now renamed `scoped_cycle_time_h` everywhere
(code, results.json, Fig 6.4 legend); no `ephase` token remains under `phase3`. The
transfer gap's 19% is a reduction in **scoped cycle time**, not the 3.22 h metric.

**S5 — The fidelity gate is borderline on cycle time.** Mean rel-EMD 0.136 vs 0.20
threshold looks comfortable, but one of eight splits fails and the earlier per-edge
model failed outright. Control flow is excellent (0.05); the *timing* fidelity is the
weaker leg. If a reviewer demanded a 0.10 cycle-time threshold, the gate would fail.
State the threshold and its rationale explicitly.

**S6 — Predicted 80.3% rests on a redesign premise.** It assumes antibiotics can be
administered essentially immediately after sepsis screening, grounded in the 83 real
cases where that happens (median 0.09 h). The bias check (`candidate_timing_sample_bias`)
shows that pool supplies 65% of affected draws, but excluding it entirely moves the
result only 2 pp (80.3→78.3), and the baseline cap bounds optimism. Still, this is a
**model-predicted** number contingent on the premise (C6) — never report it as achieved.

**S7 — C3's attribute gate barely discriminates.** 98.8% of *treated* cases are
InfectionSuspected, so "expedite the infection-suspected cohort" ≈ "expedite
everyone." The attribute-based triage is a valid mechanism demonstration but is
near-degenerate on this cohort — its 0.5 pp shortfall vs C1 is the only trace of the
6-case difference. Do not over-claim triage as a differentiated intervention here.
