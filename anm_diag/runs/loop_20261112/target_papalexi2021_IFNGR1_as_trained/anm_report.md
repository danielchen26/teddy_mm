# ANM diagnosis: loop_target_papalexi2021_IFNGR1_as_trained

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1']` vs `['control']` · 1339 treated cells in 4 clusters, 2532 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `b3468ed892fd32bd…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.865 |
| HEAD_LAT | 0.766 |
| HEAD_OUT | 0.691 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.126, -0.067] |
| HEAD_LAT->HEAD_OUT | [-0.090, -0.059] |

## Readout PDL1

Response Delta = -0.745 [-0.968, -0.486]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.01 (ridge) | [-0.14, +0.13] | **unresolved** | ridge +0.01 [-0.14, +0.13]; hgb +0.18 [-0.04, +0.33] |
| HEAD_LAT | chain | +0.38 (ridge) | [+0.19, +0.50] | **reject** | ridge +0.38 [+0.19, +0.50]; hgb +0.43 [+0.20, +0.57] |
| HEAD_OUT | chain | +0.55 (ridge) | [+0.40, +0.65] | **reject** | ridge +0.55 [+0.40, +0.65]; hgb +0.56 [+0.38, +0.68] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.29, +0.44] | ridge [+0.29, +0.44]; hgb [+0.19, +0.32] |
| HEAD_LAT->HEAD_OUT | [+0.12, +0.23] | ridge [+0.12, +0.23]; hgb [+0.07, +0.21] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 96.3 s
