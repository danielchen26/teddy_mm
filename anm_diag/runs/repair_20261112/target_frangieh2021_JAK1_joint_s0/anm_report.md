# ANM diagnosis: repair_target_frangieh2021_JAK1_joint_s0

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK1']` vs `['control']` · 248 treated cells in 23 clusters, 2248 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `a68346f597c04646…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.843 |
| HEAD_LAT | 0.819 |
| HEAD_OUT | 0.796 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.054, -0.002] |
| HEAD_LAT->HEAD_OUT | [-0.034, -0.005] |

## Readout PDL1

Response Delta = -0.714 [-0.833, -0.557]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.04 (ridge) | [-0.24, +0.05] | **support** | ridge -0.04 [-0.24, +0.05]; hgb +0.20 [+0.04, +0.28] |
| HEAD_LAT | chain | -0.01 (ridge) | [-0.17, +0.10] | **unresolved** | ridge -0.01 [-0.17, +0.10]; hgb +0.04 [-0.16, +0.13] |
| HEAD_OUT | chain | +0.10 (hgb) | [-0.05, +0.19] | **unresolved** | ridge +0.17 [+0.00, +0.25]; hgb +0.10 [-0.05, +0.19] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.05, +0.15] | ridge [-0.05, +0.15]; hgb [-0.29, -0.05] |
| HEAD_LAT->HEAD_OUT | [-0.01, +0.23] | ridge [+0.06, +0.28]; hgb [-0.01, +0.19] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 66.5 s
