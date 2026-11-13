# ANM diagnosis: loopguard_target_frangieh2021_JAK2_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK2']` vs `['control']` · 244 treated cells in 24 clusters, 2244 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `0c327c917f804a2b…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.883 |
| HEAD_LAT | 0.880 |
| HEAD_OUT | 0.865 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.013, +0.039] |
| HEAD_LAT->HEAD_OUT | [-0.037, -0.000] |

## Readout PDL1

Response Delta = -0.928 [-1.058, -0.669]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.18 (ridge) | [+0.05, +0.28] | **unresolved** | ridge +0.18 [+0.05, +0.28]; hgb +0.30 [+0.12, +0.38] |
| HEAD_LAT | chain | +0.10 (hgb) | [-0.05, +0.18] | **unresolved** | ridge +0.16 [+0.02, +0.24]; hgb +0.10 [-0.05, +0.18] |
| HEAD_OUT | chain | +0.16 (hgb) | [+0.01, +0.24] | **unresolved** | ridge +0.17 [+0.01, +0.24]; hgb +0.16 [+0.01, +0.24] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.20, +0.00] | ridge [-0.10, +0.05]; hgb [-0.27, -0.10] |
| HEAD_LAT->HEAD_OUT | [-0.02, +0.13] | ridge [-0.07, +0.07]; hgb [-0.02, +0.13] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 92.5 s
