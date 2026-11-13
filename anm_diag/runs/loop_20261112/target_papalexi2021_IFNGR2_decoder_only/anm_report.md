# ANM diagnosis: loop_target_papalexi2021_IFNGR2_decoder_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR2']` vs `['control']` · 1213 treated cells in 4 clusters, 2406 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `5737f02b405ea40e…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.854 |
| HEAD_LAT | 0.789 |
| HEAD_OUT | 0.760 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.103, -0.024] |
| HEAD_LAT->HEAD_OUT | [-0.040, -0.010] |

## Readout PDL1

Response Delta = -0.843 [-0.948, -0.720]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.14 (ridge) | [+0.09, +0.22] | **unresolved** | ridge +0.14 [+0.09, +0.22]; hgb +0.24 [+0.18, +0.35] |
| HEAD_LAT | chain | +0.43 (ridge) | [+0.36, +0.51] | **reject** | ridge +0.43 [+0.36, +0.51]; hgb +0.53 [+0.46, +0.62] |
| HEAD_OUT | chain | +0.48 (hgb) | [+0.41, +0.59] | **reject** | ridge +0.53 [+0.44, +0.63]; hgb +0.48 [+0.41, +0.59] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.21, +0.37] | ridge [+0.21, +0.37]; hgb [+0.23, +0.34] |
| HEAD_LAT->HEAD_OUT | [+0.02, +0.10] | ridge [+0.07, +0.13]; hgb [-0.07, -0.01] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 92.9 s
