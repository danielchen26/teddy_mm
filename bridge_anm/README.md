# ANM × TEDDY CITE bridge (v0 → scaled stat proof)

Minimal but real bridge: **TEDDY phase-1 typed evidence (δu)** → **ANM finite-field schema** with declared **P_f / Q_f**, plus response / attribution / flip-distance on real CITE-derived events.

## Framing

| role | object |
|---|---|
| TEDDY | typed evidence engine: predicted ADT, modality mask `rna_only`, z |
| ANM | declared field + nested readouts P_f (workability) / Q_f (exactness) + attribution / flip distances |

Success ≠ beating Pearson 0.61. Success = events exportable; schema declared; at least one response path runs; criterion editable without retraining TEDDY.

Claim boundary: **local response diagnosis only** — not clinical superiority.

## Criteria

| id | what changes | expected_action vs O0 |
|---|---|---|
| O0 | soft thr=0.12, equal panel mean | baseline |
| O1 | thr=0.28, key_marker_boost=2.0, margin=0.12 | usually **unchanged** (stricter abstain) |
| O2 | key-marker priority + lineage weights B/T>myeloid | **rewrites** labels (~12% disagreement) |

## Layout

```
bridge_anm/
  export_cite_events.py   # full site4/test (+ OOD) → JSONL typed events
  adapt_to_anm.py         # events → ANM instance bundle under O0/O1/O2
  run_demo.py             # small demo field + LOO + flip
  run_bakeoff.py          # early 500-cell 3-arm bakeoff
  run_stat_proof.py       # scaled STAT_PROOF (bootstrap + permutation)
  schemas/cite_lineage_finite_field_v0.json
  readouts/cite_lineage_O{0,1,2}.yaml
  lib/lineage_panels.py   # panels + CRITERIA
```

## Re-run (scaled statistical proof)

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
PYTORCH_ENABLE_MPS_FALLBACK=1 .venv/bin/python bridge_anm/export_cite_events.py --n-cells 0 --ood-n 2000
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_stat_proof.py --attr-n 250
```

Outputs: `outputs/anm_cite_bridge/stat_proof/STAT_PROOF.md`, `stat_proof_results.json`.
Does **not** retrain `best.pt`. Jev arm is a log-loss Choice stand-in.


## Missing-modality × ANM demo

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
# phase-2 ckpt (weak ~0.25–0.28; optional if already present)
# PYTORCH_ENABLE_MPS_FALLBACK=1 .venv/bin/python scripts/06_train_bidirectional.py --processed data/processed/cite --out outputs/cite_phase2 --epochs 6
PYTORCH_ENABLE_MPS_FALLBACK=1 .venv/bin/python bridge_anm/export_missing_modality_events.py --n-cells 0
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_missing_modality_demo.py --attr-n 1200
python3 outputs/anm_cite_bridge/missing_modality/make_mm_viz.py
```

Outputs: `outputs/anm_cite_bridge/missing_modality/MISSING_MODALITY_ANM_DEMO.md`.
Hybrid export when phase-2 RNA pearson < 0.35. Do **not** claim win over phase-1 ~0.61.

## Scope-refine proof (P_f / attribution)

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_scope_refine_proof.py
```

Outputs: `outputs/anm_cite_bridge/scope_refine/SCOPE_REFINE_PROOF.md`, `scope_refine_results.json`.
Reuses `hard_proof/attr_compact.npz` + phase-1 panels / missing-modality adt_only. No `best.pt` retrain.
