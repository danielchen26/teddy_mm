# ANM diagnosis: exposure_without_IRF1_IRF1

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IRF1']` vs `['control']` · 539 treated cells in 4 clusters, 1732 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `0a575a5ef79b2fb2…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.768 |
| HEAD_LAT | 0.611 |
| HEAD_OUT | 0.552 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.190, -0.129] |
| HEAD_LAT->HEAD_OUT | [-0.089, -0.028] |

## Readout PDL1

Response Delta = -0.285 [-0.533, -0.014]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.49 (ridge) | [+0.27, +0.81] | **unresolved** | ridge +0.49 [+0.27, +0.81]; hgb +0.59 [-0.08, +0.82] |
| HEAD_LAT | chain | +0.65 (ridge) | [+0.47, +0.98] | **reject** | ridge +0.65 [+0.47, +0.98]; hgb +0.71 [+0.54, +1.02] |
| HEAD_OUT | chain | +0.72 (ridge) | [+0.54, +0.97] | **reject** | ridge +0.72 [+0.54, +0.97]; hgb +0.74 [+0.53, +1.33] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [+0.01, +0.39] | ridge [+0.01, +0.39]; hgb [-0.11, +1.04] |
| HEAD_LAT->HEAD_OUT | [-0.21, +0.19] | ridge [-0.21, +0.19]; hgb [-0.17, +0.35] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 125.3 s
