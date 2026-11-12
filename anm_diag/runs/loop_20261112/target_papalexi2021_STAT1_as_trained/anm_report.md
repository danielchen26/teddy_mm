# ANM diagnosis: loop_target_papalexi2021_STAT1_as_trained

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['STAT1']` vs `['control']` · 444 treated cells in 4 clusters, 1637 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `136073ce2ffc9658…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.815 |
| HEAD_LAT | 0.809 |
| HEAD_OUT | 0.688 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.026, +0.046] |
| HEAD_LAT->HEAD_OUT | [-0.151, -0.039] |

## Readout PDL1

Response Delta = -0.681 [-0.853, -0.589]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.16 (ridge) | [-0.29, +0.19] | **unresolved** | ridge -0.16 [-0.29, +0.19]; hgb +0.01 [-0.09, +0.38] |
| HEAD_LAT | chain | +0.18 (hgb) | [+0.05, +0.59] | **unresolved** | ridge +0.18 [+0.06, +0.45]; hgb +0.18 [+0.05, +0.59] |
| HEAD_OUT | chain | +0.36 (hgb) | [+0.30, +0.58] | **reject** | ridge +0.46 [+0.42, +0.61]; hgb +0.36 [+0.30, +0.58] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.26, +0.50] | ridge [+0.24, +0.48]; hgb [+0.12, +0.28] |
| HEAD_LAT->HEAD_OUT | [-0.04, +0.28] | ridge [+0.10, +0.39]; hgb [-0.04, +0.28] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 81.7 s
