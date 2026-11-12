# ANM diagnosis: loop_transfer_frangieh2021_decoder_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 5171 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `70092131bfd57712…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, JAK1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.884 |
| HEAD_LAT | 0.843 |
| HEAD_OUT | 0.822 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.056, -0.027] |
| HEAD_LAT->HEAD_OUT | [-0.029, -0.012] |

## Readout MHC1

Response Delta = -0.138 [-0.182, -0.094]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.37 (hgb) | [-0.92, +0.02] | **support** | ridge -0.20 [-0.62, +0.20]; hgb -0.37 [-0.92, +0.02] |
| HEAD_LAT | chain | +0.01 (ridge) | [-0.34, +0.27] | **unresolved** | ridge +0.01 [-0.34, +0.27]; hgb +0.06 [-0.24, +0.29] |
| HEAD_OUT | chain | +0.10 (hgb) | [-0.21, +0.35] | **unresolved** | ridge +0.18 [-0.08, +0.41]; hgb +0.10 [-0.21, +0.35] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.19, +0.65] | ridge [-0.00, +0.40]; hgb [+0.20, +0.80] |
| HEAD_LAT->HEAD_OUT | [-0.02, +0.21] | ridge [+0.05, +0.33]; hgb [-0.07, +0.14] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout PDL1

Response Delta = -0.849 [-0.928, -0.762]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.04 (ridge) | [-0.01, +0.09] | **support** | ridge +0.04 [-0.01, +0.09]; hgb +0.12 [+0.03, +0.20] |
| HEAD_LAT | chain | +0.20 (ridge) | [+0.13, +0.26] | **reject** | ridge +0.20 [+0.13, +0.26]; hgb +0.43 [+0.38, +0.48] |
| HEAD_OUT | chain | +0.44 (hgb) | [+0.38, +0.48] | **reject** | ridge +0.48 [+0.41, +0.52]; hgb +0.44 [+0.38, +0.48] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.12, +0.20] | ridge [+0.12, +0.20]; hgb [+0.24, +0.39] |
| HEAD_LAT->HEAD_OUT | [+0.20, +0.27] | ridge [+0.25, +0.30]; hgb [-0.02, +0.03] |

**Root cause: loss_at.** Representation loss: a stage is rejected after a supported stage. distinction merged between 'Z' and 'HEAD_LAT'

## Readout HLAE

Response Delta = -0.079 [-0.115, -0.046]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.18 (hgb) | [-0.49, +0.57] | **unresolved** | ridge +0.20 [-0.39, +0.61]; hgb +0.18 [-0.49, +0.57] |
| HEAD_LAT | chain | +0.94 (ridge) | [+0.50, +1.27] | **reject** | ridge +0.94 [+0.50, +1.27]; hgb +1.38 [+1.09, +1.83] |
| HEAD_OUT | chain | +1.37 (ridge) | [+1.08, +1.78] | **reject** | ridge +1.37 [+1.08, +1.78]; hgb +1.38 [+1.09, +1.94] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.44, +1.21] | ridge [+0.36, +1.19]; hgb [+0.76, +2.04] |
| HEAD_LAT->HEAD_OUT | [+0.24, +0.88] | ridge [+0.24, +0.88]; hgb [-0.19, +0.33] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 234.5 s
