# ANM diagnosis: loop_target_papalexi2021_IRF1_decoder_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IRF1']` vs `['control']` · 1072 treated cells in 4 clusters, 2265 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `cd0308290d3805fe…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.786 |
| HEAD_LAT | 0.585 |
| HEAD_OUT | 0.538 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.245, -0.152] |
| HEAD_LAT->HEAD_OUT | [-0.083, -0.022] |

## Readout PDL1

Response Delta = -0.290 [-0.531, -0.029]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.33 (ridge) | [-0.16, +0.52] | **unresolved** | ridge +0.33 [-0.16, +0.52]; hgb +0.51 [-0.33, +0.72] |
| HEAD_LAT | chain | +0.88 (ridge) | [+0.75, +1.16] | **reject** | ridge +0.88 [+0.75, +1.16]; hgb +1.07 [+0.95, +1.57] |
| HEAD_OUT | chain | +1.01 (hgb) | [+0.87, +1.34] | **reject** | ridge +1.01 [+0.73, +1.14]; hgb +1.01 [+0.87, +1.34] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.38, +1.22] | ridge [+0.38, +1.22]; hgb [+0.33, +1.87] |
| HEAD_LAT->HEAD_OUT | [+0.00, +0.25] | ridge [-0.30, +0.19]; hgb [-0.27, +0.01] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 88.1 s
