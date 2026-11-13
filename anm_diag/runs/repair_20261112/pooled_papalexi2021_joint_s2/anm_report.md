# ANM diagnosis: repair_pooled_papalexi2021_joint_s2

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 6285 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `e6e7275945871ede…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, IRF1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.826 |
| HEAD_LAT | 0.186 |
| HEAD_OUT | 0.606 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.706, -0.548] |
| HEAD_LAT->HEAD_OUT | [+0.255, +0.612] |

## Readout PDL1

Response Delta = -0.683 [-0.852, -0.455]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.14 (ridge) | [+0.02, +0.34] | **unresolved** | ridge +0.14 [+0.02, +0.34]; hgb +0.23 [+0.14, +0.38] |
| HEAD_LAT | chain | +0.17 (hgb) | [+0.09, +0.22] | **unresolved** | ridge +0.36 [+0.30, +0.42]; hgb +0.17 [+0.09, +0.22] |
| HEAD_OUT | chain | +0.15 (hgb) | [+0.03, +0.36] | **unresolved** | ridge +0.16 [+0.04, +0.30]; hgb +0.15 [+0.03, +0.36] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.15, +0.14] | ridge [+0.05, +0.33]; hgb [-0.20, +0.03] |
| HEAD_LAT->HEAD_OUT | [-0.13, +0.18] | ridge [-0.31, -0.09]; hgb [-0.13, +0.18] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout PDL2

Response Delta = -0.284 [-0.338, -0.222]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.16 (ridge) | [-0.14, +0.54] | **unresolved** | ridge +0.16 [-0.14, +0.54]; hgb +0.27 [+0.01, +0.57] |
| HEAD_LAT | chain | -0.06 (hgb) | [-0.32, +0.21] | **unresolved** | ridge +0.23 [+0.07, +0.37]; hgb -0.06 [-0.32, +0.21] |
| HEAD_OUT | chain | +0.09 (ridge) | [-0.15, +0.41] | **unresolved** | ridge +0.09 [-0.15, +0.41]; hgb +0.17 [-0.22, +0.73] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.41, -0.00] | ridge [-0.27, +0.36]; hgb [-0.47, -0.16] |
| HEAD_LAT->HEAD_OUT | [-0.07, +0.29] | ridge [-0.43, +0.12]; hgb [-0.10, +0.55] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 155.4 s
