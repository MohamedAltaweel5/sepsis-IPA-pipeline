# Sepsis IPA Pipeline

Process-mining and simulation pipeline for guideline-based redesign of the emergency-department
sepsis pathway: discovery (Phase I), benchmark diagnosis and candidate generation against the
Surviving Sepsis Campaign 2021 guideline (Phase II), and a fidelity-gated simulation of the
candidates (Phase III).

## Event log

The pipeline uses the public **Sepsis Cases – Event Log** (4TU.ResearchData), which is not
redistributed here:

https://doi.org/10.4121/uuid:915d2bfb-7e84-49ad-a286-dc35f063a460

Download it, decompress it if needed, and place it at:

```
data/Sepsis_Cases_-_Event_Log.xes
```

## Running

Python 3.9.23 is required.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python run_pipeline.py      # writes outputs/results.json, outputs/figures/, outputs/models/
pytest tests                # checkpoint tests
```

All reported numbers are in `outputs/results.json`; deviations from the original design are
documented in `DEVIATIONS.md`.
