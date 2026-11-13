# ANM diagnosis: repair_pooled_papalexi2021_joint_s1

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 6285 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `d55631e53228a961…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'site' (5 values: IFNGR1, IFNGR2, IRF1, JAK2, STAT1) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each site
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.836 |
| HEAD_LAT | 0.182 |
| HEAD_OUT | 0.711 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.711, -0.564] |
| HEAD_LAT->HEAD_OUT | [+0.365, +0.630] |

## Readout PDL1

Response Delta = -0.714 [-0.878, -0.479]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.13 (ridge) | [+0.03, +0.31] | **unresolved** | ridge +0.13 [+0.03, +0.31]; hgb +0.22 [+0.14, +0.36] |
| HEAD_LAT | chain | +0.24 (hgb) | [+0.16, +0.36] | **reject** | ridge +0.36 [+0.18, +0.57]; hgb +0.24 [+0.16, +0.36] |
| HEAD_OUT | chain | +0.15 (hgb) | [+0.04, +0.35] | **unresolved** | ridge +0.16 [+0.03, +0.37]; hgb +0.15 [+0.04, +0.35] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.01, +0.17] | ridge [+0.06, +0.38]; hgb [-0.06, +0.08] |
| HEAD_LAT->HEAD_OUT | [-0.17, +0.02] | ridge [-0.34, -0.02]; hgb [-0.17, +0.02] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout PDL2

Response Delta = -0.307 [-0.359, -0.244]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.13 (ridge) | [-0.11, +0.44] | **unresolved** | ridge +0.13 [-0.11, +0.44]; hgb +0.24 [+0.02, +0.50] |
| HEAD_LAT | chain | +0.06 (hgb) | [-0.14, +0.41] | **unresolved** | ridge +0.38 [-0.02, +0.73]; hgb +0.06 [-0.14, +0.41] |
| HEAD_OUT | chain | +0.10 (ridge) | [-0.30, +0.57] | **unresolved** | ridge +0.10 [-0.30, +0.57]; hgb +0.17 [-0.21, +0.68] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.16, +0.02] | ridge [-0.10, +0.59]; hgb [-0.27, -0.06] |
| HEAD_LAT->HEAD_OUT | [-0.20, +0.22] | ridge [-0.77, +0.22]; hgb [-0.10, +0.28] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 158.5 s
