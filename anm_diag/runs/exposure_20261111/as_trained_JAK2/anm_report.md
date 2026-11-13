# ANM diagnosis: exposure_as_trained_JAK2

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK2']` vs `['control']` · 513 treated cells in 4 clusters, 1706 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `67652dd5a2981526…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.897 |
| HEAD_LAT | 0.817 |
| HEAD_OUT | 0.767 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.099, -0.059] |
| HEAD_LAT->HEAD_OUT | [-0.067, -0.026] |

## Readout PDL1

Response Delta = -0.955 [-1.038, -0.872]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.22 (ridge) | [+0.16, +0.29] | **reject** | ridge +0.22 [+0.16, +0.29]; hgb +0.39 [+0.33, +0.44] |
| HEAD_LAT | chain | +0.53 (ridge) | [+0.46, +0.58] | **reject** | ridge +0.53 [+0.46, +0.58]; hgb +0.57 [+0.50, +0.62] |
| HEAD_OUT | chain | +0.66 (ridge) | [+0.61, +0.70] | **reject** | ridge +0.66 [+0.61, +0.70]; hgb +0.66 [+0.62, +0.71] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.24, +0.36] | ridge [+0.24, +0.36]; hgb [+0.12, +0.23] |
| HEAD_LAT->HEAD_OUT | [+0.06, +0.19] | ridge [+0.06, +0.19]; hgb [+0.05, +0.17] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is reject (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 124.7 s
