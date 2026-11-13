# ANM diagnosis: repair_target_papalexi2021_IFNGR2_joint_s0

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR2']` vs `['control']` · 1213 treated cells in 4 clusters, 2406 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `902585ea68bad0a4…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.854 |
| HEAD_LAT | 0.831 |
| HEAD_OUT | 0.813 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.056, -0.003] |
| HEAD_LAT->HEAD_OUT | [-0.032, -0.008] |

## Readout PDL1

Response Delta = -0.843 [-0.948, -0.720]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.14 (ridge) | [+0.09, +0.22] | **unresolved** | ridge +0.14 [+0.09, +0.22]; hgb +0.24 [+0.18, +0.35] |
| HEAD_LAT | chain | +0.18 (ridge) | [+0.12, +0.27] | **reject** | ridge +0.18 [+0.12, +0.27]; hgb +0.19 [+0.11, +0.30] |
| HEAD_OUT | chain | +0.15 (hgb) | [+0.10, +0.25] | **unresolved** | ridge +0.17 [+0.10, +0.28]; hgb +0.15 [+0.10, +0.25] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.00, +0.09] | ridge [-0.00, +0.09]; hgb [-0.13, +0.02] |
| HEAD_LAT->HEAD_OUT | [-0.06, -0.00] | ridge [-0.03, +0.02]; hgb [-0.08, +0.01] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 63.1 s
