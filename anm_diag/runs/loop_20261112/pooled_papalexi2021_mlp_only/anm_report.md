# ANM diagnosis: loop_pooled_papalexi2021_mlp_only

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 6285 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `3c518609ddb103f3…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, IRF1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.835 |
| HEAD_LAT | 0.149 |
| HEAD_OUT | 0.630 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.732, -0.624] |
| HEAD_LAT->HEAD_OUT | [+0.293, +0.636] |

## Readout PDL1

Response Delta = -0.700 [-0.865, -0.476]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.12 (ridge) | [+0.01, +0.29] | **unresolved** | ridge +0.12 [+0.01, +0.29]; hgb +0.25 [+0.15, +0.41] |
| HEAD_LAT | chain | +0.18 (hgb) | [+0.11, +0.32] | **reject** | ridge +0.25 [+0.14, +0.44]; hgb +0.18 [+0.11, +0.32] |
| HEAD_OUT | chain | +0.16 (hgb) | [+0.05, +0.30] | **unresolved** | ridge +0.20 [+0.07, +0.39]; hgb +0.16 [+0.05, +0.30] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.00, +0.14] | ridge [+0.06, +0.25]; hgb [-0.11, +0.01] |
| HEAD_LAT->HEAD_OUT | [-0.10, +0.01] | ridge [-0.15, +0.01]; hgb [-0.10, +0.01] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout PDL2

Response Delta = -0.313 [-0.367, -0.250]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.16 (ridge) | [-0.08, +0.47] | **unresolved** | ridge +0.16 [-0.08, +0.47]; hgb +0.29 [+0.05, +0.57] |
| HEAD_LAT | chain | -0.12 (hgb) | [-0.45, +0.34] | **unresolved** | ridge +0.10 [-0.21, +0.49]; hgb -0.12 [-0.45, +0.34] |
| HEAD_OUT | chain | +0.15 (ridge) | [-0.19, +0.51] | **unresolved** | ridge +0.15 [-0.19, +0.51]; hgb +0.16 [-0.17, +0.52] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.43, -0.10] | ridge [-0.23, +0.17]; hgb [-0.59, -0.20] |
| HEAD_LAT->HEAD_OUT | [+0.08, +0.39] | ridge [-0.25, +0.26]; hgb [+0.07, +0.44] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 200.1 s
