# ANM diagnosis: loopguard_target_frangieh2021_IFNGR1_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1']` vs `['control']` · 198 treated cells in 22 clusters, 2198 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `deacdf5bc720d39b…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.891 |
| HEAD_LAT | 0.877 |
| HEAD_OUT | 0.859 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.053, +0.020] |
| HEAD_LAT->HEAD_OUT | [-0.024, -0.001] |

## Readout PDL1

Response Delta = -0.921 [-1.016, -0.681]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.14 (ridge) | [-0.10, +0.23] | **unresolved** | ridge +0.14 [-0.10, +0.23]; hgb +0.20 [-0.13, +0.27] |
| HEAD_LAT | chain | +0.04 (hgb) | [-0.15, +0.10] | **support** | ridge +0.06 [-0.15, +0.12]; hgb +0.04 [-0.15, +0.10] |
| HEAD_OUT | chain | +0.02 (hgb) | [-0.15, +0.08] | **support** | ridge +0.11 [-0.10, +0.18]; hgb +0.02 [-0.15, +0.08] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.19, +0.02] | ridge [-0.17, +0.05]; hgb [-0.23, +0.05] |
| HEAD_LAT->HEAD_OUT | [-0.08, +0.08] | ridge [-0.04, +0.11]; hgb [-0.08, +0.08] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 90.6 s
