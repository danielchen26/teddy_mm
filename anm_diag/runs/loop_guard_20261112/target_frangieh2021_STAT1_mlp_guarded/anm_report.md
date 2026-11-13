# ANM diagnosis: loopguard_target_frangieh2021_STAT1_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['STAT1']` vs `['control']` · 203 treated cells in 17 clusters, 2203 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `cf520a68ac619856…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.901 |
| HEAD_LAT | 0.897 |
| HEAD_OUT | 0.881 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.022, +0.035] |
| HEAD_LAT->HEAD_OUT | [-0.025, -0.006] |

## Readout PDL1

Response Delta = -0.977 [-1.102, -0.799]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | +0.04 (ridge) | [-0.12, +0.16] | **unresolved** | ridge +0.04 [-0.12, +0.16]; hgb +0.23 [+0.04, +0.31] |
| HEAD_LAT | chain | +0.09 (hgb) | [-0.04, +0.17] | **unresolved** | ridge +0.10 [-0.09, +0.22]; hgb +0.09 [-0.04, +0.17] |
| HEAD_OUT | chain | +0.17 (ridge) | [+0.06, +0.25] | **unresolved** | ridge +0.17 [+0.06, +0.25]; hgb +0.17 [+0.07, +0.25] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.08, +0.16] | ridge [-0.05, +0.12]; hgb [-0.21, -0.02] |
| HEAD_LAT->HEAD_OUT | [+0.01, +0.17] | ridge [+0.01, +0.17]; hgb [+0.02, +0.18] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'Z' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 89.6 s
