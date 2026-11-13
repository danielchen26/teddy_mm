# ANM diagnosis: repair_pooled_frangieh2021_joint_s2

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 3171 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `02efbcbd0ec75d88…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, JAK1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.888 |
| HEAD_LAT | 0.277 |
| HEAD_OUT | 0.840 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.698, -0.525] |
| HEAD_LAT->HEAD_OUT | [+0.483, +0.645] |

## Readout MHC1

Response Delta = -0.162 [-0.207, -0.111]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.36 (hgb) | [-0.82, +0.04] | **support** | ridge -0.26 [-0.65, +0.12]; hgb -0.36 [-0.82, +0.04] |
| HEAD_LAT | chain | -0.06 (ridge) | [-0.52, +0.30] | **unresolved** | ridge -0.06 [-0.52, +0.30]; hgb +0.07 [-0.07, +0.17] |
| HEAD_OUT | chain | -0.12 (hgb) | [-0.47, +0.22] | **unresolved** | ridge -0.11 [-0.48, +0.29]; hgb -0.12 [-0.47, +0.22] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.02, +0.56] | ridge [-0.09, +0.46]; hgb [+0.08, +0.81] |
| HEAD_LAT->HEAD_OUT | [-0.30, +0.24] | ridge [-0.28, +0.27]; hgb [-0.44, +0.07] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout PDL1

Response Delta = -0.864 [-0.948, -0.771]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.01 (ridge) | [-0.04, +0.07] | **support** | ridge +0.01 [-0.04, +0.07]; hgb +0.11 [+0.05, +0.18] |
| HEAD_LAT | chain | +0.03 (ridge) | [-0.07, +0.15] | **unresolved** | ridge +0.03 [-0.07, +0.15]; hgb +0.16 [+0.05, +0.24] |
| HEAD_OUT | chain | +0.16 (hgb) | [+0.09, +0.22] | **unresolved** | ridge +0.17 [+0.09, +0.25]; hgb +0.16 [+0.09, +0.22] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.08, +0.14] | ridge [-0.08, +0.14]; hgb [-0.05, +0.15] |
| HEAD_LAT->HEAD_OUT | [+0.01, +0.20] | ridge [+0.03, +0.22]; hgb [-0.08, +0.08] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['HEAD_LAT', 'HEAD_OUT'] unresolved after 'Z'

## Readout HLAE

Response Delta = -0.098 [-0.138, -0.055]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.03 (ridge) | [-0.69, +0.50] | **unresolved** | ridge +0.03 [-0.69, +0.50]; hgb +0.17 [-0.59, +0.63] |
| HEAD_LAT | chain | +0.57 (ridge) | [-0.07, +1.05] | **unresolved** | ridge +0.57 [-0.07, +1.05]; hgb +0.87 [+0.20, +1.68] |
| HEAD_OUT | chain | +0.87 (hgb) | [+0.40, +1.24] | **reject** | ridge +0.96 [+0.60, +1.31]; hgb +0.87 [+0.40, +1.24] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.03, +1.21] | ridge [+0.03, +1.21]; hgb [+0.24, +1.57] |
| HEAD_LAT->HEAD_OUT | [-0.01, +0.73] | ridge [+0.08, +0.92]; hgb [-0.71, +0.44] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 158.5 s
