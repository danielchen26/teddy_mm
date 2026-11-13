# ANM diagnosis: loop_pooled_frangieh2021_mlp_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 3171 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `201e91470e9d155d…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, JAK1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.882 |
| HEAD_LAT | 0.269 |
| HEAD_OUT | 0.848 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.704, -0.496] |
| HEAD_LAT->HEAD_OUT | [+0.468, +0.668] |

## Readout MHC1

Response Delta = -0.136 [-0.184, -0.089]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.43 (hgb) | [-1.03, +0.02] | **support** | ridge -0.33 [-0.80, +0.09]; hgb -0.43 [-1.03, +0.02] |
| HEAD_LAT | chain | -0.05 (ridge) | [-1.00, +0.34] | **unresolved** | ridge -0.05 [-1.00, +0.34]; hgb +0.05 [-0.43, +0.26] |
| HEAD_OUT | chain | -0.00 (ridge) | [-0.36, +0.30] | **unresolved** | ridge -0.00 [-0.36, +0.30]; hgb +0.01 [-0.40, +0.34] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.29, +0.90] | ridge [-0.54, +0.87]; hgb [-0.00, +0.94] |
| HEAD_LAT->HEAD_OUT | [-0.47, +0.88] | ridge [-0.47, +0.88]; hgb [-0.44, +0.47] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout PDL1

Response Delta = -0.870 [-0.953, -0.781]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.06 (ridge) | [-0.00, +0.11] | **unresolved** | ridge +0.06 [-0.00, +0.11]; hgb +0.13 [+0.07, +0.21] |
| HEAD_LAT | chain | +0.10 (hgb) | [+0.01, +0.18] | **unresolved** | ridge +0.11 [+0.01, +0.20]; hgb +0.10 [+0.01, +0.18] |
| HEAD_OUT | chain | +0.12 (hgb) | [+0.05, +0.17] | **unresolved** | ridge +0.14 [+0.07, +0.20]; hgb +0.12 [+0.05, +0.17] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.04, +0.11] | ridge [-0.04, +0.13]; hgb [-0.12, +0.03] |
| HEAD_LAT->HEAD_OUT | [-0.04, +0.08] | ridge [-0.04, +0.12]; hgb [-0.04, +0.08] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout HLAE

Response Delta = -0.076 [-0.116, -0.035]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.36 (ridge) | [-1.84, +0.26] | **unresolved** | ridge -0.36 [-1.84, +0.26]; hgb +0.02 [-1.22, +0.52] |
| HEAD_LAT | chain | +0.91 (ridge) | [+0.49, +2.04] | **reject** | ridge +0.91 [+0.49, +2.04]; hgb +1.32 [+0.73, +2.48] |
| HEAD_OUT | chain | +0.91 (hgb) | [+0.16, +1.55] | **reject** | ridge +1.12 [+0.62, +1.72]; hgb +0.91 [+0.16, +1.55] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.53, +3.64] | ridge [+0.53, +3.64]; hgb [+0.57, +3.62] |
| HEAD_LAT->HEAD_OUT | [-1.62, +0.70] | ridge [-1.13, +0.81]; hgb [-1.86, +0.23] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 198.5 s
