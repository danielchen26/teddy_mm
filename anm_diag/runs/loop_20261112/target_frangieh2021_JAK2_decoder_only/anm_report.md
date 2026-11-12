# ANM diagnosis: loop_target_frangieh2021_JAK2_decoder_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK2']` vs `['control']` · 244 treated cells in 24 clusters, 2244 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `d851d60a560453d1…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.883 |
| HEAD_LAT | 0.844 |
| HEAD_OUT | 0.831 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.048, +0.015] |
| HEAD_LAT->HEAD_OUT | [-0.033, +0.024] |

## Readout PDL1

Response Delta = -0.928 [-1.058, -0.669]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.18 (ridge) | [+0.05, +0.28] | **unresolved** | ridge +0.18 [+0.05, +0.28]; hgb +0.30 [+0.12, +0.38] |
| HEAD_LAT | chain | +0.27 (ridge) | [+0.11, +0.34] | **reject** | ridge +0.27 [+0.11, +0.34]; hgb +0.42 [+0.32, +0.50] |
| HEAD_OUT | chain | +0.50 (hgb) | [+0.40, +0.63] | **reject** | ridge +0.53 [+0.41, +0.61]; hgb +0.50 [+0.40, +0.63] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.05, +0.18] | ridge [-0.05, +0.18]; hgb [+0.03, +0.30] |
| HEAD_LAT->HEAD_OUT | [+0.17, +0.39] | ridge [+0.20, +0.36]; hgb [+0.02, +0.18] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'Z' is unresolved (the readout is not fully explained even by the input view), but 'HEAD_LAT' -> 'HEAD_OUT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 93.2 s
