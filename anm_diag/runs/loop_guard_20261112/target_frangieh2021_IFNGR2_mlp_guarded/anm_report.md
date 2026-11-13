# ANM diagnosis: loopguard_target_frangieh2021_IFNGR2_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR2']` vs `['control']` · 278 treated cells in 28 clusters, 2278 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `a0cceba461bb4d52…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.861 |
| HEAD_LAT | 0.854 |
| HEAD_OUT | 0.825 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.028, +0.015] |
| HEAD_LAT->HEAD_OUT | [-0.044, +0.005] |

## Readout PDL1

Response Delta = -0.844 [-1.054, -0.565]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.02 (ridge) | [-0.22, +0.07] | **support** | ridge -0.02 [-0.22, +0.07]; hgb +0.17 [-0.12, +0.25] |
| HEAD_LAT | chain | -0.01 (hgb) | [-0.19, +0.06] | **support** | ridge +0.05 [-0.19, +0.13]; hgb -0.01 [-0.19, +0.06] |
| HEAD_OUT | chain | +0.13 (ridge) | [-0.07, +0.21] | **unresolved** | ridge +0.13 [-0.07, +0.21]; hgb +0.14 [-0.01, +0.21] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.08, +0.10] | ridge [-0.04, +0.12]; hgb [-0.24, +0.02] |
| HEAD_LAT->HEAD_OUT | [+0.07, +0.18] | ridge [+0.04, +0.16]; hgb [+0.09, +0.24] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_OUT'] unresolved after 'HEAD_LAT'

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 94.0 s
