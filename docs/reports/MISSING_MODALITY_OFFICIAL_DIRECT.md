# MISSING_MODALITY_ANM_DEMO — TEDDY × ANM

## Framing

- Source choice: **hybrid_phase1_rna_phase2_adt** (prefer phase-2 ckpt; hybrid/simulate if phase-2 RNA pearson weak).
- Base site4/test cells: **n=16750** (×3 masks → 50250 cell-mask rows).
- Phase-1 unidirectional Pearson (quoted): **~0.61** (measured full-ADT on export slice: 0.6027010897971204).
- Panel Pearson by mask: `{"rna_only": 0.8998163546566806, "adt_only": 0.9448817711440152, "joint": 0.9453124536550299}`.
- **Do NOT claim win over 0.61.** If phase-2/hybrid predictions are weak, still show ANM abstain / attribution / declaration advantages.
- ANM semantics: missing modality = **source ablation / channel-restricted typed evidence**; declared P_f drops when key evidence absent; Q_f vs holdout GT only.
- O0→O2 criterion edit without TEDDY retrain.
- No Perturb / GFlowNet / 160M / ATAC. Not clinical.

## Reading note (TEDDY's role)

- **`adt_only` = RNA missing, so TEDDY does not run.** Its evidence comes from the phase-2 BidirectionalCite model run with its RNA input off: it reads the cell's 134 measured proteins and reconstructs the 9 panel proteins (panel Pearson 0.945). The same measured proteins also set the answer key.
- **Arm names in the tables below.** "TEDDY alone" is the fixed rule with no ANM layer, reading each mask's predicted panel unscaled (O0: average each lineage's 3 markers, call the highest if it reaches 0.12; O2: key marker only × lineage weight B 1.5 / T 1.3 / myeloid 0.5, bar 0.20). "TEDDY+ANM" is ANM on the same evidence (×0.35 for `adt_only`). Under `adt_only` neither arm involves TEDDY, so those rows are labelled "Fixed rule (stand-in)" and "ANM (stand-in)".
- **The `adt_only` abstentions follow the declared trust factor 0.35** (`modality_reliability`), which was chosen, not estimated. Under O2 (key-marker rule) no score can reach the 0.20 bar with that factor: the best case is 0.126 (`scripts/missing_modality_bounds.py`), so all `adt_only` cells are declined by construction.

## Key tables by mask

### Mask `rna_only` (n=11725)

| arm | crit | Q | abstain | mean_P_f |
|---|---|---:|---:|---:|
| TEDDY alone | O0 | 0.9406350667280258 | 196 | None |
| TEDDY alone | O2 | 0.9550879396984925 | 1753 | None |
| TEDDY+ANM | O0 | 0.9406350667280258 | 196 | 0.9832835820895522 |
| TEDDY+ANM | O2 | 0.9550879396984925 | 1753 | 0.8504904051172708 |

- LOO attribution (O0): n_attr=1476 / requested=1500; top1_flip_rate=0.2899728997289973; top1_dist=`{'CD2': 266, 'CD19': 248, 'CD11c': 96, 'CD36': 215, 'CD5': 215, 'CD3': 125, 'CD16': 193, 'CD72': 97, 'CD22': 21}`.
- TEDDY-alone silent over-answer O0 vs O2 decl: count=1559 rate=0.1329637526652452.
- Train-retune O2 label curve (logistic, calibrated abstain≈0.2):

| n_O2_labels | abstain | Q | strict |
|---:|---:|---:|---:|
| 50 | 0.1083 | 0.9205395721361577 | 0.8208814960999906 |
| 200 | 0.1991 | 0.9827505280450598 | 0.7870500892773236 |
| 500 | 0.1200 | 0.9625160187953866 | 0.8470068602574946 |
| 1000 | 0.1181 | 0.965153452685422 | 0.851141809980265 |

### Mask `adt_only` (n=11725)

| arm | crit | Q | abstain | mean_P_f |
|---|---|---:|---:|---:|
| Fixed rule (stand-in) | O0 | 0.9478435305917753 | 21 | None |
| Fixed rule (stand-in) | O2 | 0.9430573374124149 | 985 | None |
| ANM (stand-in) | O0 | 0.9754988111237465 | 1697 | 0.8552665245202559 |
| ANM (stand-in) | O2 | 0.9966514867398875 | 4221 | 0.64 |

- LOO attribution (O0): n_attr=1276 / requested=1500; top1_flip_rate=0.29545454545454547; top1_dist=`{'CD2': 236, 'CD22': 123, 'CD11c': 61, 'CD19': 183, 'CD5': 220, 'CD3': 179, 'CD16': 120, 'CD72': 90, 'CD36': 64}`.
- Fixed rule (stand-in) silent over-answer O0 vs O2 decl: count=965 rate=0.08230277185501066.
- Train-retune O2 label curve (logistic, calibrated abstain≈0.2):

| n_O2_labels | abstain | Q | strict |
|---:|---:|---:|---:|
| 50 | 0.0970 | 0.9349568113227182 | 0.8442815524856686 |
| 200 | 0.1326 | 0.9710725893824486 | 0.8423080537543464 |
| 500 | 0.0954 | 0.9761063785580719 | 0.8829997180716098 |
| 1000 | 0.1190 | 0.9838933333333333 | 0.86683582370078 |

### Mask `joint` (n=11725)

| arm | crit | Q | abstain | mean_P_f |
|---|---|---:|---:|---:|
| TEDDY alone | O0 | 0.9746512264064922 | 15 | None |
| TEDDY alone | O2 | 0.9579317269076305 | 1236 | None |
| TEDDY+ANM | O0 | 0.9746512264064922 | 15 | 0.9987206823027719 |
| TEDDY+ANM | O2 | 0.9579317269076305 | 1236 | 0.8945842217484008 |

- LOO attribution (O0): n_attr=1499 / requested=1500; top1_flip_rate=0.28218812541694466; top1_dist=`{'CD2': 286, 'CD19': 235, 'CD11c': 72, 'CD36': 181, 'CD5': 214, 'CD3': 158, 'CD16': 185, 'CD72': 89, 'CD22': 79}`.
- TEDDY-alone silent over-answer O0 vs O2 decl: count=1221 rate=0.10413646055437101.
- Train-retune O2 label curve (logistic, calibrated abstain≈0.2):

| n_O2_labels | abstain | Q | strict |
|---:|---:|---:|---:|
| 50 | 0.1014 | 0.9346371052081155 | 0.8398646743727093 |
| 200 | 0.1739 | 0.9848708906836537 | 0.8136453340851424 |
| 500 | 0.0998 | 0.9733792671468838 | 0.8762334367070764 |
| 1000 | 0.1282 | 0.982753045165463 | 0.8567803777840428 |

## Contrast summary

| mask | Fixed rule O0 Q | Fixed rule O0 abstain | ANM O0 Q | ANM O0 abstain | ANM O0 mean_P_f | ANM O2 abstain |
|---|---:|---:|---:|---:|---:|---:|
| rna_only | 0.9406350667280258 | 196 | 0.9406350667280258 | 196 | 0.9832835820895522 | 1753 |
| adt_only | 0.9478435305917753 | 21 | 0.9754988111237465 | 1697 | 0.8552665245202559 | 4221 |
| joint | 0.9746512264064922 | 15 | 0.9746512264064922 | 15 | 0.9987206823027719 | 1236 |

## Proof sentence

> On 16,750 site4/test base cells (11,725 scored per mask; 5,025 held for re-tuning), three masks act as source ablation: rna_only = TEDDY phase-1 predictions; adt_only = RNA missing, so TEDDY does not run and the phase-2 protein-only stand-in supplies the evidence; joint = the average of the two. ANM's P_f/abstain follow the declared per-source trust (rna_only/joint 1.0, adt_only 0.35, chosen not estimated) (rna_only: ANM abstain=196/P_f=0.9832835820895522; adt_only: ANM abstain=1697/P_f=0.8552665245202559; joint: ANM abstain=15/P_f=0.9987206823027719), O0→O2 observer edits are YAML-only (no TEDDY retrain), and LOO/flip attribution remains auditable per mask; the fixed rule (TEDDY alone's rule, no trust setting) answers every adt_only cell at O0; train-retune must consume O2 labels per mask/observer and still lacks editable-field attribution. Pearson is secondary — do not claim win over phase-1 ~0.61. Not clinical.

## Blockers / honesty

- export source=hybrid_phase1_rna_phase2_adt: phase-2 RNA path weak or absent; documented hybrid/simulate. Not a claim that phase-2 beats 0.61.
- adt_only / joint may use ADT-channel evidence (phase-2 recon or observed ADT in simulate); Q_f can look strong by construction on ADT-present masks — interpret as modality-availability, not SOTA Pearson.

## Paths

- Report: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_official/missing_modality_direct_e20/MISSING_MODALITY_ANM_DEMO.md`
- Results JSON: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_official/missing_modality_direct_e20/missing_modality_results.json`
- Export: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_official/missing_modality_direct_e20/missing_modality_events.jsonl` / `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_official/missing_modality_direct_e20/missing_modality_cells.jsonl`
- Attr NPZ dir: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_official/missing_modality_direct_e20/attr`

Generated: 2026-09-30 18:11:25 EDT