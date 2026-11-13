# ANM diagnosis: loopguard_target_papalexi2021_IFNGR1_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1']` vs `['control']` · 1339 treated cells in 4 clusters, 2532 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `c535bc155472253b…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.865 |
| HEAD_LAT | 0.820 |
| HEAD_OUT | 0.798 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.061, -0.032] |
| HEAD_LAT->HEAD_OUT | [-0.031, -0.011] |

## Readout PDL1

Response Delta = -0.745 [-0.968, -0.486]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.01 (ridge) | [-0.14, +0.13] | **unresolved** | ridge +0.01 [-0.14, +0.13]; hgb +0.18 [-0.04, +0.33] |
| HEAD_LAT | chain | +0.13 (ridge) | [-0.03, +0.24] | **unresolved** | ridge +0.13 [-0.03, +0.24]; hgb +0.14 [+0.02, +0.24] |
| HEAD_OUT | chain | +0.14 (ridge) | [+0.01, +0.24] | **unresolved** | ridge +0.14 [+0.01, +0.24]; hgb +0.17 [+0.06, +0.26] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.08, +0.17] | ridge [+0.08, +0.17]; hgb [-0.12, +0.06] |
| HEAD_LAT->HEAD_OUT | [-0.02, +0.04] | ridge [-0.02, +0.04]; hgb [-0.01, +0.06] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 94.3 s
