# Thesis placeholder map

> Auto-generated from `results.json`. Every value is computed by the pipeline. **All Phase III (predicted/*fidelity*) values are MODEL-PREDICTED, not empirical (C6).**

| Thesis location | Quantity | Value | results.json key |
|---|---|---|---|
| §4.5.2 / abstract | Baseline 1-hour compliance (real) | **34.5%** | `prep.baseline_benchmark.compliance_1h_pct` |
| §4.5.2 | Baseline median time-to-antibiotics (real) | **1.883 h** | `prep.baseline_benchmark.tta_median_h` |
| §4.5.2 | Baseline emergency-phase duration (real, median) | **3.219 h** | `prep.baseline_benchmark.ephase_median_h` |
| §2.10 / §4.5.1 | Fidelity - control-flow distance (CFLD) | **0.0521** | `fidelity_cfld` |
| §2.10 / §4.5.1 | Fidelity - cycle-time relative EMD | **0.1356** | `fidelity_cycletime_distance` |
| §2.10 / §4.5.1 | Fidelity gate passed | **True** | `fidelity_gate_passed` |
| §4.5.2 / abstract | Predicted 1-hour compliance (top candidate) | **80.371%** | `phase3...C1_resequence_abx_before_lactate.predicted_compliance_1h` |
| §4.5.2 | Predicted median TTA (top candidate) | **0.39 h** | `phase3...C1_resequence_abx_before_lactate.predicted_median_tta_h` |
| §4.5.2 | Predicted scoped cycle time (top candidate) [NOT the §7.1 emergency phase] | **1.556 h** | `phase3...C1_resequence_abx_before_lactate.predicted_scoped_cycle_time_h` |
| §4.5.2 / abstract | Composite risk reduction (top candidate) | **53.59%** | `phase3...C1_resequence_abx_before_lactate.predicted_risk_reduction` |
