# ANM diagnosis: loop_target_frangieh2021_JAK1_as_trained

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK1']` vs `['control']` · 248 treated cells in 23 clusters, 2248 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `03a4cbcf7d0ebe2c…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.843 |
| HEAD_LAT | 0.795 |
| HEAD_OUT | 0.751 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.103, +0.001] |
| HEAD_LAT->HEAD_OUT | [-0.065, -0.014] |

## Readout PDL1

Response Delta = -0.714 [-0.833, -0.557]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.04 (ridge) | [-0.24, +0.05] | **support** | ridge -0.04 [-0.24, +0.05]; hgb +0.20 [+0.04, +0.28] |
| HEAD_LAT | chain | +0.20 (ridge) | [+0.01, +0.28] | **unresolved** | ridge +0.20 [+0.01, +0.28]; hgb +0.37 [+0.26, +0.48] |
| HEAD_OUT | chain | +0.41 (hgb) | [+0.29, +0.51] | **reject** | ridge +0.48 [+0.40, +0.59]; hgb +0.41 [+0.29, +0.51] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.14, +0.35] | ridge [+0.14, +0.35]; hgb [+0.10, +0.31] |
| HEAD_LAT->HEAD_OUT | [+0.11, +0.41] | ridge [+0.20, +0.52]; hgb [-0.06, +0.13] |

**Root cause: loss_at.** Representation loss: a stage is rejected after a supported stage. distinction merged between 'Z' and 'HEAD_OUT' (via unresolved ['HEAD_LAT'])

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 90.1 s
