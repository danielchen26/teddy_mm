# ANM diagnosis: loopguard_target_papalexi2021_IFNGR2_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR2']` vs `['control']` · 1213 treated cells in 4 clusters, 2406 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `a060ac5eb2008fe5…`; anm_diag 0.1.0

**Checks**
- [S2] warning: the source pools events that differ in 'stratum' (4 values: rep1-tx, rep2-tx, rep3-tx, rep4-tx) without declaring it under source.pooling: the verdicts hold for the pooled source, not for each stratum
- [C5] warning: no positive control node: support verdicts will be marked uninformative
- [K1] warning: only 4 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.854 |
| HEAD_LAT | 0.836 |
| HEAD_OUT | 0.817 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.046, +0.006] |
| HEAD_LAT->HEAD_OUT | [-0.028, -0.011] |

## Readout PDL1

Response Delta = -0.843 [-0.948, -0.720]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.14 (ridge) | [+0.09, +0.22] | **unresolved** | ridge +0.14 [+0.09, +0.22]; hgb +0.24 [+0.18, +0.35] |
| HEAD_LAT | chain | +0.18 (ridge) | [+0.12, +0.27] | **reject** | ridge +0.18 [+0.12, +0.27]; hgb +0.19 [+0.13, +0.31] |
| HEAD_OUT | chain | +0.14 (hgb) | [+0.08, +0.23] | **unresolved** | ridge +0.16 [+0.10, +0.25]; hgb +0.14 [+0.08, +0.23] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.01, +0.09] | ridge [-0.01, +0.09]; hgb [-0.11, +0.01] |
| HEAD_LAT->HEAD_OUT | [-0.08, +0.00] | ridge [-0.04, +0.01]; hgb [-0.10, -0.01] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 98.6 s
