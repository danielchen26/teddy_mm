# ANM diagnosis: loop_transfer_papalexi2021_decoder_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 7478 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `2d0bbea7ab5b18ab…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, IRF1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.757 |
| HEAD_LAT | 0.636 |
| HEAD_OUT | 0.582 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.158, -0.062] |
| HEAD_LAT->HEAD_OUT | [-0.070, -0.026] |

## Readout PDL1

Response Delta = -0.693 [-0.856, -0.469]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.15 (ridge) | [+0.04, +0.32] | **unresolved** | ridge +0.15 [+0.04, +0.32]; hgb +0.24 [+0.16, +0.38] |
| HEAD_LAT | chain | +0.35 (ridge) | [+0.30, +0.45] | **reject** | ridge +0.35 [+0.30, +0.45]; hgb +0.43 [+0.38, +0.56] |
| HEAD_OUT | chain | +0.47 (ridge) | [+0.43, +0.58] | **reject** | ridge +0.47 [+0.43, +0.58]; hgb +0.50 [+0.47, +0.61] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.12, +0.28] | ridge [+0.12, +0.28]; hgb [+0.16, +0.24] |
| HEAD_LAT->HEAD_OUT | [+0.11, +0.15] | ridge [+0.11, +0.15]; hgb [+0.05, +0.10] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Readout PDL2

Response Delta = -0.287 [-0.338, -0.229]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.13 (ridge) | [-0.11, +0.36] | **unresolved** | ridge +0.13 [-0.11, +0.36]; hgb +0.23 [+0.04, +0.42] |
| HEAD_LAT | chain | +0.40 (ridge) | [+0.22, +0.62] | **reject** | ridge +0.40 [+0.22, +0.62]; hgb +0.51 [+0.35, +0.72] |
| HEAD_OUT | chain | +0.57 (hgb) | [+0.40, +0.80] | **reject** | ridge +0.57 [+0.40, +0.82]; hgb +0.57 [+0.40, +0.80] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.22, +0.35] | ridge [+0.22, +0.35]; hgb [+0.25, +0.34] |
| HEAD_LAT->HEAD_OUT | [+0.13, +0.21] | ridge [+0.13, +0.21]; hgb [+0.03, +0.09] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 218.4 s
