# ANM diagnosis: repair_target_frangieh2021_STAT1_joint_s0

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['STAT1']` vs `['control']` · 203 treated cells in 17 clusters, 2203 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `03a8724090499ebe…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.901 |
| HEAD_LAT | 0.882 |
| HEAD_OUT | 0.884 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.038, +0.004] |
| HEAD_LAT->HEAD_OUT | [-0.007, +0.026] |

## Readout PDL1

Response Delta = -0.977 [-1.102, -0.799]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.04 (ridge) | [-0.12, +0.16] | **unresolved** | ridge +0.04 [-0.12, +0.16]; hgb +0.23 [+0.04, +0.31] |
| HEAD_LAT | chain | +0.12 (ridge) | [-0.04, +0.20] | **unresolved** | ridge +0.12 [-0.04, +0.20]; hgb +0.18 [+0.03, +0.25] |
| HEAD_OUT | chain | +0.27 (hgb) | [+0.15, +0.36] | **reject** | ridge +0.31 [+0.17, +0.40]; hgb +0.27 [+0.15, +0.36] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.05, +0.16] | ridge [-0.05, +0.16]; hgb [-0.14, +0.08] |
| HEAD_LAT->HEAD_OUT | [+0.09, +0.27] | ridge [+0.14, +0.29]; hgb [+0.03, +0.18] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 66.5 s
