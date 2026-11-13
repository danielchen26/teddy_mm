# ANM diagnosis: exposure_without_IRF1_JAK2

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK2']` vs `['control']` · 513 treated cells in 4 clusters, 1706 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `f6a7a83355a74a4f…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.897 |
| HEAD_LAT | 0.857 |
| HEAD_OUT | 0.821 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.050, -0.028] |
| HEAD_LAT->HEAD_OUT | [-0.051, -0.015] |

## Readout PDL1

Response Delta = -0.955 [-1.038, -0.872]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.22 (ridge) | [+0.16, +0.29] | **reject** | ridge +0.22 [+0.16, +0.29]; hgb +0.39 [+0.33, +0.44] |
| HEAD_LAT | chain | +0.32 (ridge) | [+0.24, +0.39] | **reject** | ridge +0.32 [+0.24, +0.39]; hgb +0.35 [+0.29, +0.41] |
| HEAD_OUT | chain | +0.36 (hgb) | [+0.30, +0.41] | **reject** | ridge +0.39 [+0.32, +0.44]; hgb +0.36 [+0.30, +0.41] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.04, +0.16] | ridge [+0.04, +0.16]; hgb [-0.10, +0.02] |
| HEAD_LAT->HEAD_OUT | [-0.04, +0.09] | ridge [+0.00, +0.11]; hgb [-0.05, +0.06] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is rejected: the input view itself does not reproduce this response for any declared readout class, so by these classes the response is not in the input (a limit of the data or the readout classes, not a loss inside the system); no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 124.8 s
