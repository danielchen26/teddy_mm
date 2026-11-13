# ANM diagnosis: repair_target_papalexi2021_JAK2_joint_s0

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK2']` vs `['control']` · 1024 treated cells in 4 clusters, 2217 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `40cc989fa3f20dcb…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.893 |
| HEAD_LAT | 0.867 |
| HEAD_OUT | 0.839 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.041, -0.008] |
| HEAD_LAT->HEAD_OUT | [-0.036, -0.019] |

## Readout PDL1

Response Delta = -0.908 [-0.999, -0.820]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.15 (ridge) | [+0.08, +0.22] | **unresolved** | ridge +0.15 [+0.08, +0.22]; hgb +0.31 [+0.25, +0.36] |
| HEAD_LAT | chain | +0.21 (hgb) | [+0.17, +0.27] | **reject** | ridge +0.22 [+0.15, +0.29]; hgb +0.21 [+0.17, +0.27] |
| HEAD_OUT | chain | +0.22 (hgb) | [+0.17, +0.28] | **reject** | ridge +0.24 [+0.18, +0.29]; hgb +0.22 [+0.17, +0.28] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.02, +0.11] | ridge [+0.02, +0.13]; hgb [-0.14, -0.04] |
| HEAD_LAT->HEAD_OUT | [-0.02, +0.04] | ridge [-0.01, +0.04]; hgb [-0.02, +0.04] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 59.0 s
