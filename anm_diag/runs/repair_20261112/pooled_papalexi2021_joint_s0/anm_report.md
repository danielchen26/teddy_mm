# ANM diagnosis: repair_pooled_papalexi2021_joint_s0

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 6285 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `fed918d4a505acef…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, IRF1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.835 |
| HEAD_LAT | 0.194 |
| HEAD_OUT | 0.672 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.698, -0.585] |
| HEAD_LAT->HEAD_OUT | [+0.290, +0.600] |

## Readout PDL1

Response Delta = -0.700 [-0.865, -0.476]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.12 (ridge) | [+0.01, +0.29] | **unresolved** | ridge +0.12 [+0.01, +0.29]; hgb +0.25 [+0.15, +0.41] |
| HEAD_LAT | chain | +0.18 (hgb) | [+0.07, +0.37] | **unresolved** | ridge +0.29 [+0.19, +0.45]; hgb +0.18 [+0.07, +0.37] |
| HEAD_OUT | chain | +0.17 (hgb) | [+0.03, +0.38] | **unresolved** | ridge +0.19 [+0.03, +0.44]; hgb +0.17 [+0.03, +0.38] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.01, +0.13] | ridge [+0.10, +0.27]; hgb [-0.12, -0.00] |
| HEAD_LAT->HEAD_OUT | [-0.08, +0.04] | ridge [-0.24, +0.01]; hgb [-0.08, +0.04] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout PDL2

Response Delta = -0.313 [-0.367, -0.250]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.16 (ridge) | [-0.08, +0.47] | **unresolved** | ridge +0.16 [-0.08, +0.47]; hgb +0.29 [+0.05, +0.57] |
| HEAD_LAT | chain | +0.04 (hgb) | [-0.33, +0.48] | **unresolved** | ridge +0.15 [-0.15, +0.53]; hgb +0.04 [-0.33, +0.48] |
| HEAD_OUT | chain | +0.24 (hgb) | [-0.14, +0.67] | **unresolved** | ridge +0.30 [-0.11, +0.76]; hgb +0.24 [-0.14, +0.67] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.32, +0.10] | ridge [-0.15, +0.18]; hgb [-0.45, -0.02] |
| HEAD_LAT->HEAD_OUT | [-0.08, +0.39] | ridge [-0.19, +0.33]; hgb [-0.08, +0.39] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 156.8 s
