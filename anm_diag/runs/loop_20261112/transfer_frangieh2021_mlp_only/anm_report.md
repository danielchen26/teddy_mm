# ANM diagnosis: loop_transfer_frangieh2021_mlp_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 5171 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `18219bdd28505b20…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, JAK1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.884 |
| HEAD_LAT | 0.839 |
| HEAD_OUT | 0.824 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.063, -0.028] |
| HEAD_LAT->HEAD_OUT | [-0.019, -0.012] |

## Readout MHC1

Response Delta = -0.138 [-0.182, -0.094]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.37 (hgb) | [-0.92, +0.02] | **support** | ridge -0.20 [-0.62, +0.20]; hgb -0.37 [-0.92, +0.02] |
| HEAD_LAT | chain | -0.07 (ridge) | [-0.40, +0.17] | **unresolved** | ridge -0.07 [-0.40, +0.17]; hgb +0.14 [-0.11, +0.34] |
| HEAD_OUT | chain | +0.15 (ridge) | [-0.13, +0.38] | **unresolved** | ridge +0.15 [-0.13, +0.38]; hgb +0.16 [-0.15, +0.39] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.07, +0.59] | ridge [-0.11, +0.33]; hgb [+0.26, +0.90] |
| HEAD_LAT->HEAD_OUT | [+0.14, +0.34] | ridge [+0.14, +0.34]; hgb [-0.13, +0.16] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout PDL1

Response Delta = -0.849 [-0.928, -0.762]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.04 (ridge) | [-0.01, +0.09] | **support** | ridge +0.04 [-0.01, +0.09]; hgb +0.12 [+0.03, +0.20] |
| HEAD_LAT | chain | +0.17 (ridge) | [+0.08, +0.23] | **unresolved** | ridge +0.17 [+0.08, +0.23]; hgb +0.37 [+0.29, +0.42] |
| HEAD_OUT | chain | +0.43 (ridge) | [+0.36, +0.48] | **reject** | ridge +0.43 [+0.36, +0.48]; hgb +0.44 [+0.38, +0.49] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.07, +0.17] | ridge [+0.07, +0.17]; hgb [+0.18, +0.31] |
| HEAD_LAT->HEAD_OUT | [+0.24, +0.29] | ridge [+0.24, +0.29]; hgb [+0.04, +0.11] |

**Root cause: loss_at.** Representation loss: a stage is rejected after a supported stage. distinction merged between 'Z' and 'HEAD_OUT' (via unresolved ['HEAD_LAT'])

## Readout HLAE

Response Delta = -0.079 [-0.115, -0.046]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.18 (hgb) | [-0.49, +0.57] | **unresolved** | ridge +0.20 [-0.39, +0.61]; hgb +0.18 [-0.49, +0.57] |
| HEAD_LAT | chain | +0.68 (ridge) | [+0.26, +0.95] | **reject** | ridge +0.68 [+0.26, +0.95]; hgb +0.91 [+0.52, +1.24] |
| HEAD_OUT | chain | +1.04 (hgb) | [+0.76, +1.41] | **reject** | ridge +1.19 [+0.92, +1.53]; hgb +1.04 [+0.76, +1.41] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.22, +0.91] | ridge [+0.13, +0.86]; hgb [+0.44, +1.26] |
| HEAD_LAT->HEAD_OUT | [+0.16, +0.83] | ridge [+0.35, +0.91]; hgb [-0.04, +0.42] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 237.1 s
