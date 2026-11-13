# ANM diagnosis: loopguard_target_papalexi2021_STAT1_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['STAT1']` vs `['control']` · 444 treated cells in 4 clusters, 1637 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `2ec3f4817b1f05a7…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.815 |
| HEAD_LAT | 0.820 |
| HEAD_OUT | 0.780 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.010, +0.055] |
| HEAD_LAT->HEAD_OUT | [-0.056, -0.019] |

## Readout PDL1

Response Delta = -0.681 [-0.853, -0.589]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.16 (ridge) | [-0.29, +0.19] | **unresolved** | ridge -0.16 [-0.29, +0.19]; hgb +0.01 [-0.09, +0.38] |
| HEAD_LAT | chain | -0.02 (hgb) | [-0.13, +0.31] | **unresolved** | ridge +0.04 [-0.06, +0.31]; hgb -0.02 [-0.13, +0.31] |
| HEAD_OUT | chain | +0.01 (hgb) | [-0.07, +0.25] | **unresolved** | ridge +0.02 [-0.08, +0.32]; hgb +0.01 [-0.07, +0.25] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.03, +0.24] | ridge [+0.07, +0.26]; hgb [-0.11, +0.02] |
| HEAD_LAT->HEAD_OUT | [-0.11, +0.09] | ridge [-0.06, +0.04]; hgb [-0.11, +0.09] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 90.3 s
