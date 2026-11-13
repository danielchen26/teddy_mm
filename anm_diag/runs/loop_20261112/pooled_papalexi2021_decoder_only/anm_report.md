# ANM diagnosis: loop_pooled_papalexi2021_decoder_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 6285 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `09868eab59c20db7…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, IRF1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.835 |
| HEAD_LAT | 0.755 |
| HEAD_OUT | 0.724 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.103, -0.052] |
| HEAD_LAT->HEAD_OUT | [-0.045, -0.015] |

## Readout PDL1

Response Delta = -0.700 [-0.865, -0.476]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.12 (ridge) | [+0.01, +0.29] | **unresolved** | ridge +0.12 [+0.01, +0.29]; hgb +0.25 [+0.15, +0.41] |
| HEAD_LAT | chain | +0.39 (ridge) | [+0.31, +0.52] | **reject** | ridge +0.39 [+0.31, +0.52]; hgb +0.46 [+0.40, +0.57] |
| HEAD_OUT | chain | +0.47 (hgb) | [+0.41, +0.58] | **reject** | ridge +0.51 [+0.43, +0.64]; hgb +0.47 [+0.41, +0.58] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.22, +0.32] | ridge [+0.22, +0.32]; hgb [+0.16, +0.27] |
| HEAD_LAT->HEAD_OUT | [+0.03, +0.13] | ridge [+0.07, +0.16]; hgb [-0.02, +0.03] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Readout PDL2

Response Delta = -0.313 [-0.367, -0.250]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.16 (ridge) | [-0.08, +0.47] | **unresolved** | ridge +0.16 [-0.08, +0.47]; hgb +0.29 [+0.05, +0.57] |
| HEAD_LAT | chain | +0.47 (ridge) | [+0.25, +0.77] | **reject** | ridge +0.47 [+0.25, +0.77]; hgb +0.55 [+0.34, +0.80] |
| HEAD_OUT | chain | +0.53 (hgb) | [+0.31, +0.80] | **reject** | ridge +0.57 [+0.34, +0.90]; hgb +0.53 [+0.31, +0.80] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.24, +0.38] | ridge [+0.24, +0.38]; hgb [+0.19, +0.35] |
| HEAD_LAT->HEAD_OUT | [-0.00, +0.11] | ridge [+0.04, +0.16]; hgb [-0.09, +0.03] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 198.7 s
