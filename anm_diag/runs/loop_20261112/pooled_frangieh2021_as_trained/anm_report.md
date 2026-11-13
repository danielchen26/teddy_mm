# ANM diagnosis: loop_pooled_frangieh2021_as_trained

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 3171 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `5c834a9c095a4865…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, JAK1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.882 |
| HEAD_LAT | 0.842 |
| HEAD_OUT | 0.813 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.048, -0.030] |
| HEAD_LAT->HEAD_OUT | [-0.036, -0.022] |

## Readout MHC1

Response Delta = -0.136 [-0.184, -0.089]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.43 (hgb) | [-1.03, +0.02] | **support** | ridge -0.33 [-0.80, +0.09]; hgb -0.43 [-1.03, +0.02] |
| HEAD_LAT | chain | -0.06 (ridge) | [-0.40, +0.24] | **unresolved** | ridge -0.06 [-0.40, +0.24]; hgb +0.06 [-0.24, +0.35] |
| HEAD_OUT | chain | +0.19 (ridge) | [-0.07, +0.40] | **unresolved** | ridge +0.19 [-0.07, +0.40]; hgb +0.22 [-0.05, +0.49] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.10, +0.74] | ridge [+0.03, +0.51]; hgb [+0.20, +0.96] |
| HEAD_LAT->HEAD_OUT | [+0.08, +0.43] | ridge [+0.08, +0.43]; hgb [-0.03, +0.32] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout PDL1

Response Delta = -0.870 [-0.953, -0.781]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.06 (ridge) | [-0.00, +0.11] | **unresolved** | ridge +0.06 [-0.00, +0.11]; hgb +0.13 [+0.07, +0.21] |
| HEAD_LAT | chain | +0.22 (ridge) | [+0.15, +0.27] | **reject** | ridge +0.22 [+0.15, +0.27]; hgb +0.40 [+0.34, +0.46] |
| HEAD_OUT | chain | +0.45 (hgb) | [+0.40, +0.49] | **reject** | ridge +0.47 [+0.41, +0.52]; hgb +0.45 [+0.40, +0.49] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.11, +0.21] | ridge [+0.11, +0.21]; hgb [+0.19, +0.35] |
| HEAD_LAT->HEAD_OUT | [+0.19, +0.27] | ridge [+0.22, +0.28]; hgb [+0.01, +0.08] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Readout HLAE

Response Delta = -0.076 [-0.116, -0.035]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.36 (ridge) | [-1.84, +0.26] | **unresolved** | ridge -0.36 [-1.84, +0.26]; hgb +0.02 [-1.22, +0.52] |
| HEAD_LAT | chain | +0.74 (ridge) | [+0.10, +1.08] | **reject** | ridge +0.74 [+0.10, +1.08]; hgb +1.01 [+0.44, +1.53] |
| HEAD_OUT | chain | +1.11 (hgb) | [+0.52, +1.66] | **reject** | ridge +1.28 [+0.93, +1.80]; hgb +1.11 [+0.52, +1.66] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.61, +2.12] | ridge [+0.61, +2.12]; hgb [+0.55, +2.08] |
| HEAD_LAT->HEAD_OUT | [+0.10, +0.96] | ridge [+0.30, +1.24]; hgb [-0.20, +0.47] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 195.7 s
