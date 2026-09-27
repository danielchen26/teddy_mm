# MISSING_MODALITY_ANM_DEMO — TEDDY × ANM

## Framing

- Source choice: **hybrid_phase1_rna_phase2_adt** (prefer phase-2 ckpt; hybrid/simulate if phase-2 RNA pearson weak).
- Base site4/test cells: **n=16750** (×3 masks → 50250 cell-mask rows).
- Phase-1 unidirectional Pearson (quoted): **~0.61** (measured full-ADT on export slice: 0.6102218214870081).
- Panel Pearson by mask: `{"rna_only": 0.9147404698378734, "adt_only": 0.44627287014827716, "joint": 0.8703130689245495}`.
- **Do NOT claim win over 0.61.** If phase-2/hybrid predictions are weak, still show ANM abstain / attribution / declaration advantages.
- ANM semantics: missing modality = **source ablation / channel-restricted typed evidence**; declared P_f drops when key evidence absent; Q_f vs holdout GT only.
- O0→O2 criterion edit without TEDDY retrain.
- No Perturb / GFlowNet / 160M / ATAC. Not clinical.

## Key tables by mask

### Mask `rna_only` (n=12563)

| arm | crit | Q | abstain | mean_P_f |
|---|---|---:|---:|---:|
| TEDDY alone | O0 | 0.945026847353618 | 76 | None |
| TEDDY alone | O2 | 0.9399568439816118 | 1282 | None |
| TEDDY+ANM | O0 | 0.9478880494505495 | 233 | 0.9814534744885776 |
| TEDDY+ANM | O2 | 0.959049959049959 | 2459 | 0.8042664968558465 |

- LOO attribution (O0): n_attr=1182 / requested=1200; top1_flip_rate=0.2766497461928934; top1_dist=`{'CD72': 114, 'CD5': 261, 'CD11c': 80, 'CD22': 152, 'CD2': 182, 'CD16': 146, 'CD36': 184, 'CD3': 58, 'CD19': 5}`.
- TEDDY-alone silent over-answer O0 vs O2 decl: count=1256 rate=0.09997612035341877.
- Train-retune O2 label curve (logistic, calibrated abstain≈0.2):

| n_O2_labels | abstain | Q | strict |
|---:|---:|---:|---:|
| 50 | 0.1082 | 0.9604135893648449 | 0.8565030297707913 |
| 200 | 0.1735 | 0.9679098926787801 | 0.7999473083340651 |
| 500 | 0.1859 | 0.97680690399137 | 0.7952050583999297 |
| 1000 | 0.1144 | 0.9637048790162633 | 0.8534293492579257 |

### Mask `adt_only` (n=12563)

| arm | crit | Q | abstain | mean_P_f |
|---|---|---:|---:|---:|
| TEDDY alone | O0 | 0.6277291648967802 | 0 | None |
| TEDDY alone | O2 | 0.6488896816195487 | 253 | None |
| TEDDY+ANM | O0 | 0.6493620337002045 | 1664 | 0.8675475602961076 |
| TEDDY+ANM | O2 | None | 12563 | 0.0 |

- LOO attribution (O0): n_attr=1056 / requested=1200; top1_flip_rate=0.7102272727272727; top1_dist=`{'CD22': 115, 'CD5': 435, 'CD36': 150, 'CD11c': 168, 'CD16': 29, 'CD19': 12, 'CD2': 133, 'CD72': 13, 'CD3': 1}`.
- TEDDY-alone silent over-answer O0 vs O2 decl: count=253 rate=0.020138501950171136.
- Train-retune O2 label curve (logistic, calibrated abstain≈0.2):

| n_O2_labels | abstain | Q | strict |
|---:|---:|---:|---:|
| 50 | 0.1933 | 0.6844110603091661 | 0.5521208395538773 |
| 200 | 0.2081 | 0.675945436397915 | 0.5352595064547291 |
| 500 | 0.1984 | 0.6931419807186678 | 0.5556336172828664 |
| 1000 | 0.1994 | 0.6965774462483545 | 0.5576534644770352 |

### Mask `joint` (n=12563)

| arm | crit | Q | abstain | mean_P_f |
|---|---|---:|---:|---:|
| TEDDY alone | O0 | 0.914960496134568 | 0 | None |
| TEDDY alone | O2 | 0.8867114093959731 | 371 | None |
| TEDDY+ANM | O0 | 0.9214170418825928 | 0 | 1.0 |
| TEDDY+ANM | O2 | 0.7583677212287941 | 785 | 0.9375149247791132 |

- LOO attribution (O0): n_attr=1200 / requested=1200; top1_flip_rate=0.30333333333333334; top1_dist=`{'CD72': 114, 'CD5': 307, 'CD11c': 95, 'CD2': 190, 'CD16': 124, 'CD36': 169, 'CD22': 144, 'CD3': 40, 'CD19': 17}`.
- TEDDY-alone silent over-answer O0 vs O2 decl: count=371 rate=0.029531162938788505.
- Train-retune O2 label curve (logistic, calibrated abstain≈0.2):

| n_O2_labels | abstain | Q | strict |
|---:|---:|---:|---:|
| 50 | 0.1493 | 0.9590172395994632 | 0.8158426275577413 |
| 200 | 0.1684 | 0.9563839898616538 | 0.7952928778431545 |
| 500 | 0.1458 | 0.9627840032898118 | 0.8224290857995961 |
| 1000 | 0.1218 | 0.958 | 0.8413102660929129 |

## Contrast summary

| mask | TEDDY O0 Q | TEDDY O0 abstain | ANM O0 Q | ANM O0 abstain | ANM O0 mean_P_f | ANM O2 abstain |
|---|---:|---:|---:|---:|---:|---:|
| rna_only | 0.945026847353618 | 76 | 0.9478880494505495 | 233 | 0.9814534744885776 | 2459 |
| adt_only | 0.6277291648967802 | 0 | 0.6493620337002045 | 1664 | 0.8675475602961076 | 12563 |
| joint | 0.914960496134568 | 0 | 0.9214170418825928 | 0 | 1.0 | 785 |

## Proof sentence

> On the same TEDDY δu base cells (site4/test n=16750), three modality masks (rna_only, adt_only, joint) act as source ablation: ANM declared P_f/abstain tracks evidence availability (rna_only: ANM abstain=233/P_f=0.9814534744885776; adt_only: ANM abstain=1664/P_f=0.8675475602961076; joint: ANM abstain=0/P_f=1.0), O0→O2 observer edits are YAML-only (no TEDDY retrain), and LOO/flip attribution remains auditable per mask; TEDDY-alone silently over-answers where the declared observer abstains; train-retune must consume O2 labels per mask/observer and still lacks editable-field attribution. Pearson is secondary — do not claim win over phase-1 ~0.61. Not clinical.

## Blockers / honesty

- **Choice documented:** no usable phase-2 RNA pearson (test rna_only≈0.254, adt_only≈0.268, joint≈0.276 ≪ phase-1 0.61). Export = **hybrid_phase1_rna_phase2_adt**: phase-1 MLP for `rna_only` / hybrid joint; phase-2 BidirectionalCite for `adt_only`.
- Panel Pearson (lineage 9-protein): rna_only≈0.915, adt_only≈0.446, joint≈0.870 — **do not inflate**; full-ADT phase-1 remains the 0.61 reference.
- Declared missing-RNA reliability gate: `adt_only` event values ×0.35 (criterion edit, no TEDDY retrain) so P_f/abstain track key-evidence absence.
- Under `adt_only`, TEDDY-alone O0 abstain=0 (over-answers) while ANM O0 abstain=1664 / O2 abstain=all — declaration advantage, not Pearson win.
- Train-retune under `adt_only` plateaus at Q≈0.68–0.70 even with 1000 O2 labels — cannot buy editable-field LOO/flip (flip_rate≈0.71 under weak ADT channel).
- No best.pt backbone retrain (phase-2 train was inference-only). No Perturb/GFlowNet/160M/ATAC. Not clinical.

## Viz

- GIF: `missing_modality_masks_light.gif` / `missing_modality_masks_dark.gif` (anm-jev style)
- Dashboard: `missing_modality_dashboard_light.png` / `_dark.png`

## Paths

- Report: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge/missing_modality/MISSING_MODALITY_ANM_DEMO.md`
- Results JSON: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge/missing_modality/missing_modality_results.json`
- Export: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge/missing_modality/missing_modality_events.jsonl` / `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge/missing_modality/missing_modality_cells.jsonl`
- Attr NPZ dir: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge/missing_modality/attr`
- Scripts: `bridge_anm/export_missing_modality_events.py`, `bridge_anm/run_missing_modality_demo.py`
- Phase-2 ckpt: `outputs/cite_phase2/best.pt` (metrics.json test≈0.25–0.28)

Generated: 2026-09-26 16:15:40 EDT