# ANM diagnosis: loopguard_target_papalexi2021_IRF1_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IRF1']` vs `['control']` · 1072 treated cells in 4 clusters, 2265 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `19d62c5dd2474534…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.786 |
| HEAD_LAT | 0.628 |
| HEAD_OUT | 0.579 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.190, -0.130] |
| HEAD_LAT->HEAD_OUT | [-0.066, -0.028] |

## Readout PDL1

Response Delta = -0.290 [-0.531, -0.029]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.33 (ridge) | [-0.16, +0.52] | **unresolved** | ridge +0.33 [-0.16, +0.52]; hgb +0.51 [-0.33, +0.72] |
| HEAD_LAT | chain | +0.52 (ridge) | [-0.14, +0.67] | **unresolved** | ridge +0.52 [-0.14, +0.67]; hgb +0.60 [+0.24, +0.76] |
| HEAD_OUT | chain | +0.60 (ridge) | [-0.03, +0.74] | **unresolved** | ridge +0.60 [-0.03, +0.74]; hgb +0.65 [+0.27, +0.79] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.18, +0.34] | ridge [-0.18, +0.34]; hgb [-0.06, +0.53] |
| HEAD_LAT->HEAD_OUT | [+0.02, +0.19] | ridge [+0.02, +0.19]; hgb [-0.14, +0.15] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 96.6 s
