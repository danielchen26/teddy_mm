# ANM diagnosis: loopguard_pooled_frangieh2021_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 3171 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `8ded67bf4b37b7a5…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, JAK1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.882 |
| HEAD_LAT | 0.327 |
| HEAD_OUT | 0.846 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.634, -0.483] |
| HEAD_LAT->HEAD_OUT | [+0.446, +0.597] |

## Readout MHC1

Response Delta = -0.136 [-0.184, -0.089]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.43 (hgb) | [-1.03, +0.02] | **support** | ridge -0.33 [-0.80, +0.09]; hgb -0.43 [-1.03, +0.02] |
| HEAD_LAT | chain | -0.19 (ridge) | [-1.29, +0.28] | **unresolved** | ridge -0.19 [-1.29, +0.28]; hgb -0.04 [-0.39, +0.12] |
| HEAD_OUT | chain | -0.01 (ridge) | [-0.38, +0.29] | **unresolved** | ridge -0.01 [-0.38, +0.29]; hgb +0.05 [-0.32, +0.39] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.60, +0.81] | ridge [-0.86, +0.78]; hgb [+0.01, +0.81] |
| HEAD_LAT->HEAD_OUT | [-0.38, +1.18] | ridge [-0.38, +1.18]; hgb [-0.23, +0.48] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout PDL1

Response Delta = -0.870 [-0.953, -0.781]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.06 (ridge) | [-0.00, +0.11] | **unresolved** | ridge +0.06 [-0.00, +0.11]; hgb +0.13 [+0.07, +0.21] |
| HEAD_LAT | chain | +0.10 (ridge) | [+0.05, +0.16] | **unresolved** | ridge +0.10 [+0.05, +0.16]; hgb +0.10 [+0.05, +0.15] |
| HEAD_OUT | chain | +0.12 (hgb) | [+0.06, +0.17] | **unresolved** | ridge +0.13 [+0.06, +0.19]; hgb +0.12 [+0.06, +0.17] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.01, +0.09] | ridge [-0.01, +0.09]; hgb [-0.09, +0.01] |
| HEAD_LAT->HEAD_OUT | [-0.05, +0.09] | ridge [-0.03, +0.10]; hgb [-0.03, +0.07] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout HLAE

Response Delta = -0.076 [-0.116, -0.035]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.36 (ridge) | [-1.84, +0.26] | **unresolved** | ridge -0.36 [-1.84, +0.26]; hgb +0.02 [-1.22, +0.52] |
| HEAD_LAT | chain | +0.74 (ridge) | [+0.46, +1.37] | **reject** | ridge +0.74 [+0.46, +1.37]; hgb +1.21 [+0.84, +2.26] |
| HEAD_OUT | chain | +1.04 (hgb) | [+0.45, +1.67] | **reject** | ridge +1.19 [+0.75, +1.79]; hgb +1.04 [+0.45, +1.67] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.52, +2.85] | ridge [+0.52, +2.85]; hgb [+0.60, +3.10] |
| HEAD_LAT->HEAD_OUT | [-0.53, +0.81] | ridge [-0.21, +0.85]; hgb [-1.12, +0.27] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 196.1 s
