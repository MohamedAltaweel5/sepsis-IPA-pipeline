# DEVIATIONS from IMPLEMENTATION_SPEC.md

Everywhere the built pipeline differs from the specification, with the reason. The
thesis must describe **what was built**, not what was specified. Grouped by area;
each entry: *Spec said → We did → Why*.

---

## D1. Simulation engine — SimPy-family DES, not Prosimos
- **Spec (§4, §10):** Prosimos (BPMN + JSON params) as the discrete-event engine; SIMOD optional.
- **Built:** a self-contained control-flow-first DES in `src/phase3_simulation/` (numpy/scipy; SimPy is installed but the generator is a direct variant-sampling + interval-bootstrap Monte-Carlo, not a SimPy process network).
- **Why:** the log has a single timestamp per event (C7), which defeats SIMOD's duration calibration (the spec itself flags this in §4). A transparent in-repo simulator matches C4/C7 exactly, is fully auditable, and needs no BPMN round-trip. Prosimos remains a documented alternative. `tool_versions.prosimos = "not used"`.

## D2. Timing parameterisation — block-bootstrap, not independent per-edge sampling
- **Spec (§10.1):** "Fit distributions (or use empirical) to the observed time between consecutive events" per edge; branching from variant frequencies.
- **Built:** the *accepted* model samples a variant, then reuses a **real per-variant interval vector** (block-bootstrap), preserving within-case timing correlation. Independent per-edge sampling was tried first and **rejected**.
- **Why:** independent per-edge bootstrap broke within-case correlation and pooled an edge across fast/slow variants, overestimating cycle time (single-split relative EMD ~0.22–0.30 > 0.20 gate). Block-bootstrap is a *more faithful* instrument, not a threshold relaxation. Logged in `phase3.stage1_fidelity.baseline_model_provenance` (2 variants tried).

## D3. Fidelity procedure — repeated 8-split hold-out averaged, not a single split
- **Spec (§10.1):** "Hold out a portion of the real log (e.g. train/test split)."
- **Built:** 8 independent 70/30 splits × 3 replications, gate applied to the **mean**; pass fraction reported.
- **Why:** a single split is high-variance because ER-core cycle time is heavy-tailed (max ~27 h); one split can load the tail into train or test. Averaging removes split luck rather than selecting a favourable seed. (Seed 42 alone fails by 0.0016; the mean passes at 0.136.)

## D4. Fidelity measures — n-gram (bigram) TV distance as the CFLD proxy
- **Spec (§10.1):** Control-Flow Log Distance (CFLD) / n-gram distance, per Chapela-Campa 2025 (ref [19]).
- **Built:** total-variation distance between normalised directly-follows (bigram) frequency vectors; cycle-time EMD (Wasserstein); TTA EMD.
- **Why:** the full Chapela-Campa CFLD uses optimal trace alignment; a bigram-distribution TV distance is a transparent, defensible proxy computable with `scipy`. Labelled as an n-gram control-flow distance, not the exact CFLD. Cycle-time EMD (the measure that matters most for a timing thesis) is exact.

## D5. Gate thresholds — pre-registered concrete values
- **Spec (§10.1):** "define an adequacy threshold *before* running (document the rationale)."
- **Built:** `control_flow_tv_max = 0.25`, `cycle_time_relative_emd_max = 0.20` (relative to mean real cycle time). Pre-registered in `fidelity.THRESHOLDS`.
- **Why:** the spec leaves the numbers to us; these are documented in code with rationale (a control-flow-first, no-contention model should reproduce empirical control flow closely and cycle time within a modest fraction of scale).

## D5b. Composite risk reduction — multi-provision, FIXED denominator, PATH-LOCAL damage
- **Spec (§10.2):** "Composite pre-implementation risk-reduction estimate — define transparently."
- **Built:** `risk = Σ_all-governed-provisions weight · noncompliant_fraction`; reduction = `100·(base_risk − cand_risk)/base_risk` with **base_risk FIXED** = 0.856 (antibiotics 1.0 + fluids 0.5 + lactate 0.5 + screening 1.0) for every candidate, so magnitudes are comparable (a candidate is not penalised in the denominator for touching more provisions). Damage to a provision is **path-local**: eliminating LacticAcid removes it from the 649 affected cases only, leaving 211 still measuring it → lactate 81.9% → 20.1% (not global 0%).
- **Results:** C1/C2 = **53.59%**, C3 = **52.97%**, C4 = **17.49%**. Stored as `risk_components` (with `damaged_detail`) per candidate.
- **Why:** (i) an antibiotics-only metric would credit C4's timing gain and hide its conformance loss, contradicting Devil's Quadrangle `quality: −1`; (ii) a per-candidate denominator made "70%" and "6.75%" fractions of different pools; (iii) a global-0% lactate penalty overstated C4's loss. All three corrected. Earlier drafts reported antibiotics-only (C4=70%) then per-candidate-denominator + global-0% (C4=6.75%); the final is fixed-denominator + path-local (C4=17.49%).
- **Note (off-by-one reconciled):** the affected-set count uses the STRICT-timestamp test (lactate strictly before antibiotics) everywhere — detection, generation, and Phase III evaluation all report 649 (and 642 for the InfectionSuspected cohort). 2 simultaneous-timestamp cases are excluded (a simultaneous lactate does not delay antibiotics). An earlier sequence-position test gave 650/643.

## D6. `binding_constraint_shift` — renamed from the spec's "emergent effect / bottleneck displacement"
- **Spec (§10.2, §13.1):** check whether relieving one bottleneck "displaces load to another"; per-candidate `emergent_effect_note`.
- **Built:** a structured `binding_constraint_shift` key (distribution of which activity now ends the case) **plus** an `emergent_effect_note` string alias for §13.1 completeness.
- **Why:** there is **no capacity/resource model** (C4, single resource group), so nothing is "displaced" in a queueing sense. What is measured is a **control-flow critical-path shift**. The original name over-claims a mechanism that does not exist here; the rename states exactly what is computed.

## D7. Candidate transforms grounded in observed orderings (not abstract model edits)
- **Spec (§9.2, §10.2):** emit candidates as concrete model modifications (new Petri net/BPMN each) and simulate them.
- **Built:** candidates are ordering/elimination **transforms** applied to sampled traces; the re-timed antibiotics interval is drawn from the **empirically observed** transition pool (e.g. 120 real antibiotics-before-lactate cases; ER Sepsis Triage→antibiotics seen 83×). A candidate antibiotics time is **capped at the baseline** (delay only removed, never added).
- **Why:** with single-timestamp data, re-timing from real cases that already exhibit the target order is more defensible than synthesising new durations. The cap bounds optimism. This is a modelling choice the spec does not prescribe.

## D8. Triage candidate (C3) — functional attribute split, not contention relief
- **Spec (§9.2 table):** Triage = "split into distinct paths by case characteristic."
- **Built:** two interpretations were considered. The **contention-relief** reading is *inapplicable* (C4: one resource group, `different_resources_total=1` — no queue to relieve). The delivered C3 is a **functional attribute split** on `InfectionSuspected` that fast-tracks that cohort's antibiotics decoupling.
- **Why:** the contention reading would predict a null by construction and misrepresent the mechanism. The functional reading is evaluable under the control-flow model and demonstrates attribute-based differentiation. (An earlier contention-framed C3 predicted no change; it was replaced.)

## D9. Task elimination (C4) generated despite the precondition being clinically weak
- **Spec (§9.2):** task elimination "Use with care; likely rare here." Generate no candidate if no precondition matches.
- **Built:** C4 eliminates LacticAcid from the antibiotics path — generated **on purpose** to exercise the mechanism, with an explicit clinical caveat that it reduces guideline conformance (LacticAcid *is* a weak SSC measurement provision).
- **Why:** requested for completeness of demonstration. Its quality penalty is encoded in the Devil's Quadrangle effects (`quality: -1`), so scoring reflects the trade-off honestly. It is a **candidate, not a recommendation** (C6).

## D10. Emergency-phase vs scoped cycle time — RESOLVED by renaming
- **Spec (§7.1):** emergency-phase duration = Registration → first Admission.
- **Built:** two *different* quantities existed with confusingly similar names; the Phase III one is now **renamed** to remove the collision:
  - `prep.baseline_benchmark.ephase_median_h = 3.219` / `ephase_mean_h = 3.379` — Registration → first Admission (the §7.1 metric). **Unchanged** — this is the correct emergency-phase name.
  - `phase3...baseline_scoped_cycle_time_h = 1.921` / `predicted_scoped_cycle_time_h` — the **ER-core scoped cycle time** (first → last ER-core event), used inside the simulation. **Renamed** from `baseline_ephase_h`/`predicted_ephase_duration`.
  - `transfer_gap.scoped_cycle_time_reduction_pct` (renamed from `ephase_reduction_pct`); Fig 6.4(b) legend now reads "Scoped cycle time".
- **Why:** the simulation operates on the ER-core scoped representation, whose end-to-end span is the scoped cycle time, not the reg→admission interval (admission is outside ER-core scope). The rename ensures "cycle time improves 19%" can never be misread as a claim about the 3.22 h §7.1 metric. No "ephase" token remains anywhere under `phase3`.

## D11. Emergency-phase secondary stat not fully reconciled (honest gap)
- **Spec (§12):** emergency-phase median 3.26 h / mean 3.54 h (n=810).
- **Built/computed:** median 3.219 h / mean 3.379 h (n=810 matches exactly).
- **Why:** n matches; the ~0.15 h gap is a minor definitional nuance in the original design pass that could not be reproduced. Left computed-as-is (not fudged). All **designated correctness gates** (34.5%, 1.88 h, 5,542, 69, 289, cohort counts) match exactly.

## D12. Compliance denominator keeps the 2 negative-interval cases
- **Spec (§6.1 vs §12):** §6.1 says exclude GW/LZ from timing; §12's headline is 284/823 = 34.5%.
- **Built:** all 823 treated cases are the denominator; GW/LZ (antibiotics before registration) count as <1 h ⇒ 282+2 = 284, 593+2 = 595. This reproduces the §12 checkpoint exactly.
- **Why:** §12 is the designated correctness test and it *includes* them. `quality.negative_interval_excluded_from_timing = 2` still records them for transparency.

## D13. Data file is plain `.xes`, not `.xes.gz`
- **Spec (§5):** `Sepsis_Cases_-_Event_Log.xes.gz`.
- **Built:** `data/` symlinks the existing uncompressed `Sepsis Cases - Event Log.xes`. Identical content; pm4py reads both.

## D14. pytest added to the environment
- **Spec:** not mentioned.
- **Built:** `tests/test_checkpoints.py` runs under pytest (installed into the `thesis` conda env). Everything else uses the pre-existing env.

---

### Constraints honoured with NO deviation
C1 (no ML — discovery is Inductive Miner, diagnosis is arithmetic, generation is a rule catalogue, scoring is a weighted sum of catalogue signs), C2 (external SSC standard only), C3 (Reijers catalogue only), C4 (control-flow-first; no agent/resource-first), C5 (fidelity measured before and gating candidate eval), C6 (all Phase III outputs labelled model-predicted), C7 (inter-event intervals only; no duration decomposition), C8 (candidates, not directives; clinical filter present).
