# ANM diagnosis: repair_target_frangieh2021_JAK2_joint_s0

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK2']` vs `['control']` · 244 treated cells in 24 clusters, 2244 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `7e8822f0600a364c…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.883 |
| HEAD_LAT | 0.872 |
| HEAD_OUT | 0.862 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.022, +0.040] |
| HEAD_LAT->HEAD_OUT | [-0.024, +0.008] |

## Readout PDL1

Response Delta = -0.928 [-1.058, -0.669]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.18 (ridge) | [+0.05, +0.28] | **unresolved** | ridge +0.18 [+0.05, +0.28]; hgb +0.30 [+0.12, +0.38] |
| HEAD_LAT | chain | +0.15 (ridge) | [+0.00, +0.23] | **unresolved** | ridge +0.15 [+0.00, +0.23]; hgb +0.27 [+0.13, +0.37] |
| HEAD_OUT | chain | +0.24 (hgb) | [+0.12, +0.35] | **reject** | ridge +0.29 [+0.19, +0.39]; hgb +0.24 [+0.12, +0.35] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.15, +0.06] | ridge [-0.15, +0.06]; hgb [-0.13, +0.12] |
| HEAD_LAT->HEAD_OUT | [+0.03, +0.20] | ridge [+0.10, +0.25]; hgb [-0.09, +0.06] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 67.1 s
