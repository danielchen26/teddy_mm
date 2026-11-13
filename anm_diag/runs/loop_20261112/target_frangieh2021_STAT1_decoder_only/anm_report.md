# ANM diagnosis: loop_target_frangieh2021_STAT1_decoder_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['STAT1']` vs `['control']` · 203 treated cells in 17 clusters, 2203 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `42a02111b7e5aa83…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.901 |
| HEAD_LAT | 0.851 |
| HEAD_OUT | 0.830 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.071, +0.002] |
| HEAD_LAT->HEAD_OUT | [-0.042, +0.005] |

## Readout PDL1

Response Delta = -0.977 [-1.102, -0.799]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.04 (ridge) | [-0.12, +0.16] | **unresolved** | ridge +0.04 [-0.12, +0.16]; hgb +0.23 [+0.04, +0.31] |
| HEAD_LAT | chain | +0.28 (ridge) | [+0.10, +0.37] | **reject** | ridge +0.28 [+0.10, +0.37]; hgb +0.43 [+0.21, +0.51] |
| HEAD_OUT | chain | +0.46 (hgb) | [+0.32, +0.56] | **reject** | ridge +0.54 [+0.41, +0.62]; hgb +0.46 [+0.32, +0.56] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.10, +0.33] | ridge [+0.10, +0.33]; hgb [+0.04, +0.31] |
| HEAD_LAT->HEAD_OUT | [+0.09, +0.29] | ridge [+0.20, +0.36]; hgb [-0.03, +0.21] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 93.0 s
