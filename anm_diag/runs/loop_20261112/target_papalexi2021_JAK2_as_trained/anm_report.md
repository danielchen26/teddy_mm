# ANM diagnosis: loop_target_papalexi2021_JAK2_as_trained

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK2']` vs `['control']` · 1024 treated cells in 4 clusters, 2217 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `49363e398105008f…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.893 |
| HEAD_LAT | 0.819 |
| HEAD_OUT | 0.763 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.107, -0.041] |
| HEAD_LAT->HEAD_OUT | [-0.078, -0.034] |

## Readout PDL1

Response Delta = -0.908 [-0.999, -0.820]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.15 (ridge) | [+0.08, +0.22] | **unresolved** | ridge +0.15 [+0.08, +0.22]; hgb +0.31 [+0.25, +0.36] |
| HEAD_LAT | chain | +0.45 (ridge) | [+0.39, +0.50] | **reject** | ridge +0.45 [+0.39, +0.50]; hgb +0.54 [+0.49, +0.58] |
| HEAD_OUT | chain | +0.64 (ridge) | [+0.58, +0.68] | **reject** | ridge +0.64 [+0.58, +0.68]; hgb +0.66 [+0.62, +0.72] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.24, +0.36] | ridge [+0.24, +0.36]; hgb [+0.18, +0.29] |
| HEAD_LAT->HEAD_OUT | [+0.13, +0.24] | ridge [+0.13, +0.24]; hgb [+0.09, +0.17] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 90.7 s
