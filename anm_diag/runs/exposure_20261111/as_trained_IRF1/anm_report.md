# ANM diagnosis: exposure_as_trained_IRF1

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IRF1']` vs `['control']` · 539 treated cells in 4 clusters, 1732 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `10cf7446decf6fee…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.768 |
| HEAD_LAT | 0.573 |
| HEAD_OUT | 0.497 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.254, -0.140] |
| HEAD_LAT->HEAD_OUT | [-0.114, -0.038] |

## Readout PDL1

Response Delta = -0.285 [-0.533, -0.014]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.49 (ridge) | [+0.27, +0.81] | **unresolved** | ridge +0.49 [+0.27, +0.81]; hgb +0.59 [-0.08, +0.82] |
| HEAD_LAT | chain | +0.95 (ridge) | [+0.81, +1.68] | **reject** | ridge +0.95 [+0.81, +1.68]; hgb +1.07 [+0.93, +1.85] |
| HEAD_OUT | chain | +1.13 (ridge) | [+0.99, +1.82] | **reject** | ridge +1.13 [+0.99, +1.82]; hgb +1.17 [+1.02, +2.29] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.27, +1.16] | ridge [+0.27, +1.16]; hgb [+0.20, +1.77] |
| HEAD_LAT->HEAD_OUT | [-0.04, +0.32] | ridge [-0.04, +0.32]; hgb [-0.00, +0.50] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 128.7 s
