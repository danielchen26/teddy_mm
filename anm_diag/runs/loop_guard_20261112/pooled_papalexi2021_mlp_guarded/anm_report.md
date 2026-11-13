# ANM diagnosis: loopguard_pooled_papalexi2021_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 6285 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `282d45fb8472a58a…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, IRF1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.835 |
| HEAD_LAT | 0.186 |
| HEAD_OUT | 0.739 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.708, -0.601] |
| HEAD_LAT->HEAD_OUT | [+0.452, +0.634] |

## Readout PDL1

Response Delta = -0.700 [-0.865, -0.476]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.12 (ridge) | [+0.01, +0.29] | **unresolved** | ridge +0.12 [+0.01, +0.29]; hgb +0.25 [+0.15, +0.41] |
| HEAD_LAT | chain | +0.16 (hgb) | [+0.07, +0.32] | **unresolved** | ridge +0.29 [+0.18, +0.46]; hgb +0.16 [+0.07, +0.32] |
| HEAD_OUT | chain | +0.18 (hgb) | [+0.07, +0.37] | **unresolved** | ridge +0.22 [+0.07, +0.49]; hgb +0.18 [+0.07, +0.37] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.00, +0.11] | ridge [+0.09, +0.28]; hgb [-0.13, -0.03] |
| HEAD_LAT->HEAD_OUT | [-0.05, +0.07] | ridge [-0.21, +0.08]; hgb [-0.05, +0.07] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout PDL2

Response Delta = -0.313 [-0.367, -0.250]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.16 (ridge) | [-0.08, +0.47] | **unresolved** | ridge +0.16 [-0.08, +0.47]; hgb +0.29 [+0.05, +0.57] |
| HEAD_LAT | chain | +0.02 (hgb) | [-0.35, +0.46] | **unresolved** | ridge +0.14 [-0.20, +0.56]; hgb +0.02 [-0.35, +0.46] |
| HEAD_OUT | chain | +0.15 (ridge) | [-0.25, +0.57] | **unresolved** | ridge +0.15 [-0.25, +0.57]; hgb +0.17 [-0.21, +0.61] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.32, +0.08] | ridge [-0.19, +0.20]; hgb [-0.47, -0.03] |
| HEAD_LAT->HEAD_OUT | [-0.02, +0.21] | ridge [-0.24, +0.18]; hgb [-0.08, +0.33] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 198.7 s
