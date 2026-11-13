# ANM diagnosis: repair_pooled_frangieh2021_joint_s1

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 3171 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `534db15ab832ece7…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, JAK1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.889 |
| HEAD_LAT | 0.305 |
| HEAD_OUT | 0.855 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.670, -0.506] |
| HEAD_LAT->HEAD_OUT | [+0.474, +0.639] |

## Readout MHC1

Response Delta = -0.130 [-0.178, -0.080]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.46 (hgb) | [-1.11, +0.05] | **support** | ridge -0.33 [-0.93, +0.20]; hgb -0.46 [-1.11, +0.05] |
| HEAD_LAT | chain | -0.43 (ridge) | [-1.51, +0.26] | **unresolved** | ridge -0.43 [-1.51, +0.26]; hgb -0.08 [-0.62, +0.29] |
| HEAD_OUT | chain | -0.34 (hgb) | [-0.95, +0.09] | **support** | ridge -0.32 [-1.00, +0.11]; hgb -0.34 [-0.95, +0.09] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.73, +0.58] | ridge [-0.88, +0.38]; hgb [-0.12, +0.84] |
| HEAD_LAT->HEAD_OUT | [-0.46, +0.88] | ridge [-0.33, +0.76]; hgb [-0.69, +0.22] |

**Root cause: no_deficit.** No deficit: every stage supported. every chain stage supported for this readout

## Readout PDL1

Response Delta = -0.861 [-0.944, -0.765]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.03 (ridge) | [-0.03, +0.08] | **support** | ridge +0.03 [-0.03, +0.08]; hgb +0.07 [-0.00, +0.14] |
| HEAD_LAT | chain | +0.04 (ridge) | [-0.06, +0.12] | **unresolved** | ridge +0.04 [-0.06, +0.12]; hgb +0.12 [+0.03, +0.21] |
| HEAD_OUT | chain | +0.14 (hgb) | [+0.08, +0.21] | **unresolved** | ridge +0.14 [+0.07, +0.21]; hgb +0.14 [+0.08, +0.21] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.08, +0.10] | ridge [-0.08, +0.10]; hgb [-0.05, +0.16] |
| HEAD_LAT->HEAD_OUT | [-0.01, +0.21] | ridge [+0.01, +0.20]; hgb [-0.08, +0.11] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout HLAE

Response Delta = -0.073 [-0.114, -0.031]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.06 (ridge) | [-1.21, +0.53] | **unresolved** | ridge -0.06 [-1.21, +0.53]; hgb +0.36 [-0.58, +0.87] |
| HEAD_LAT | chain | +0.94 (ridge) | [+0.50, +1.58] | **reject** | ridge +0.94 [+0.50, +1.58]; hgb +1.41 [+0.94, +2.80] |
| HEAD_OUT | chain | +0.78 (hgb) | [-0.05, +1.27] | **unresolved** | ridge +1.09 [+0.62, +1.54]; hgb +0.78 [-0.05, +1.27] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.24, +2.48] | ridge [+0.24, +2.48]; hgb [+0.27, +3.02] |
| HEAD_LAT->HEAD_OUT | [-1.26, +0.33] | ridge [-0.53, +0.62]; hgb [-2.58, -0.06] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 146.9 s
