# ANM diagnosis: loop_target_frangieh2021_IFNGR1_as_trained

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1']` vs `['control']` · 198 treated cells in 22 clusters, 2198 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `00151ddad89fbaf4…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.891 |
| HEAD_LAT | 0.818 |
| HEAD_OUT | 0.778 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.119, -0.014] |
| HEAD_LAT->HEAD_OUT | [-0.056, -0.026] |

## Readout PDL1

Response Delta = -0.921 [-1.016, -0.681]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.14 (ridge) | [-0.10, +0.23] | **unresolved** | ridge +0.14 [-0.10, +0.23]; hgb +0.20 [-0.13, +0.27] |
| HEAD_LAT | chain | +0.28 (ridge) | [+0.07, +0.34] | **unresolved** | ridge +0.28 [+0.07, +0.34]; hgb +0.45 [+0.29, +0.52] |
| HEAD_OUT | chain | +0.46 (hgb) | [+0.35, +0.52] | **reject** | ridge +0.49 [+0.36, +0.56]; hgb +0.46 [+0.35, +0.52] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.04, +0.27] | ridge [+0.04, +0.27]; hgb [+0.15, +0.54] |
| HEAD_LAT->HEAD_OUT | [+0.10, +0.36] | ridge [+0.17, +0.34]; hgb [-0.05, +0.12] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'HEAD_LAT' -> 'HEAD_OUT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 93.0 s
