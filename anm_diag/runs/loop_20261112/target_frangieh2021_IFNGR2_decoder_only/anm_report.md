# ANM diagnosis: loop_target_frangieh2021_IFNGR2_decoder_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR2']` vs `['control']` · 278 treated cells in 28 clusters, 2278 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `bb3dade22cbc5693…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.861 |
| HEAD_LAT | 0.792 |
| HEAD_OUT | 0.751 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.124, -0.002] |
| HEAD_LAT->HEAD_OUT | [-0.059, +0.010] |

## Readout PDL1

Response Delta = -0.844 [-1.054, -0.565]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.02 (ridge) | [-0.22, +0.07] | **support** | ridge -0.02 [-0.22, +0.07]; hgb +0.17 [-0.12, +0.25] |
| HEAD_LAT | chain | +0.23 (ridge) | [+0.01, +0.33] | **unresolved** | ridge +0.23 [+0.01, +0.33]; hgb +0.48 [+0.33, +0.58] |
| HEAD_OUT | chain | +0.51 (hgb) | [+0.38, +0.62] | **reject** | ridge +0.52 [+0.36, +0.60]; hgb +0.51 [+0.38, +0.62] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.13, +0.34] | ridge [+0.13, +0.34]; hgb [+0.25, +0.62] |
| HEAD_LAT->HEAD_OUT | [+0.18, +0.52] | ridge [+0.21, +0.44]; hgb [-0.02, +0.13] |

**Root cause: loss_at.** Representation loss: a stage is rejected after a supported stage. distinction merged between 'Z' and 'HEAD_OUT' (via unresolved ['HEAD_LAT'])

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 95.9 s
