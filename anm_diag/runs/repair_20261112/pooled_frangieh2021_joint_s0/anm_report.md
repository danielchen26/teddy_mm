# ANM diagnosis: repair_pooled_frangieh2021_joint_s0

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 3171 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `502b23f2bed1a7e1…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, JAK1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.882 |
| HEAD_LAT | 0.261 |
| HEAD_OUT | 0.832 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.681, -0.578] |
| HEAD_LAT->HEAD_OUT | [+0.522, +0.632] |

## Readout MHC1

Response Delta = -0.136 [-0.184, -0.089]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.43 (hgb) | [-1.03, +0.02] | **support** | ridge -0.33 [-0.80, +0.09]; hgb -0.43 [-1.03, +0.02] |
| HEAD_LAT | chain | -0.10 (hgb) | [-0.59, +0.23] | **unresolved** | ridge -0.04 [-0.81, +0.27]; hgb -0.10 [-0.59, +0.23] |
| HEAD_OUT | chain | -0.12 (ridge) | [-0.52, +0.20] | **unresolved** | ridge -0.12 [-0.52, +0.20]; hgb -0.05 [-0.48, +0.26] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.21, +0.84] | ridge [-0.35, +0.79]; hgb [-0.21, +0.84] |
| HEAD_LAT->HEAD_OUT | [-0.45, +0.45] | ridge [-0.49, +0.49]; hgb [-0.40, +0.56] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout PDL1

Response Delta = -0.870 [-0.953, -0.781]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.06 (ridge) | [-0.00, +0.11] | **unresolved** | ridge +0.06 [-0.00, +0.11]; hgb +0.13 [+0.07, +0.21] |
| HEAD_LAT | chain | +0.13 (ridge) | [+0.05, +0.20] | **unresolved** | ridge +0.13 [+0.05, +0.20]; hgb +0.18 [+0.11, +0.23] |
| HEAD_OUT | chain | +0.17 (hgb) | [+0.10, +0.23] | **unresolved** | ridge +0.22 [+0.14, +0.29]; hgb +0.17 [+0.10, +0.23] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.01, +0.13] | ridge [-0.01, +0.13]; hgb [-0.03, +0.12] |
| HEAD_LAT->HEAD_OUT | [-0.03, +0.13] | ridge [+0.04, +0.16]; hgb [-0.07, +0.03] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout HLAE

Response Delta = -0.076 [-0.116, -0.035]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.36 (ridge) | [-1.84, +0.26] | **unresolved** | ridge -0.36 [-1.84, +0.26]; hgb +0.02 [-1.22, +0.52] |
| HEAD_LAT | chain | +0.71 (ridge) | [+0.37, +1.48] | **reject** | ridge +0.71 [+0.37, +1.48]; hgb +0.88 [+0.55, +1.56] |
| HEAD_OUT | chain | +0.96 (hgb) | [+0.36, +1.50] | **reject** | ridge +1.00 [+0.51, +1.39]; hgb +0.96 [+0.36, +1.50] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.43, +3.12] | ridge [+0.43, +3.12]; hgb [+0.20, +2.62] |
| HEAD_LAT->HEAD_OUT | [-0.87, +0.78] | ridge [-0.70, +0.74]; hgb [-0.98, +0.66] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 162.8 s
