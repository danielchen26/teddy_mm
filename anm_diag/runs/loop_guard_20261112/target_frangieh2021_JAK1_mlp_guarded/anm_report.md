# ANM diagnosis: loopguard_target_frangieh2021_JAK1_mlp_guarded

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['JAK1']` vs `['control']` · 248 treated cells in 23 clusters, 2248 cells in all
- chain: Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `957cea049626af50…`; anm_diag 0.1.0

**Checks**
- [C5] warning: no positive control node: support verdicts will be marked uninformative

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| Z | 0.843 |
| HEAD_LAT | 0.828 |
| HEAD_OUT | 0.800 |

| step | change, 95% interval |
|---|---|
| Z->HEAD_LAT | [-0.036, -0.001] |
| HEAD_LAT->HEAD_OUT | [-0.045, -0.005] |

## Readout PDL1

Response Delta = -0.714 [-0.833, -0.557]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| Z | chain | -0.04 (ridge) | [-0.24, +0.05] | **support** | ridge -0.04 [-0.24, +0.05]; hgb +0.20 [+0.04, +0.28] |
| HEAD_LAT | chain | -0.01 (hgb) | [-0.29, +0.08] | **support** | ridge -0.02 [-0.24, +0.11]; hgb -0.01 [-0.29, +0.08] |
| HEAD_OUT | chain | -0.04 (hgb) | [-0.25, +0.03] | **support** | ridge +0.02 [-0.16, +0.11]; hgb -0.04 [-0.25, +0.03] |

| step | change in u (best class) | same class |
|---|---|---|
| Z->HEAD_LAT | [-0.11, +0.12] | ridge [-0.09, +0.12]; hgb [-0.38, -0.12] |
| HEAD_LAT->HEAD_OUT | [-0.10, +0.08] | ridge [-0.04, +0.16]; hgb [-0.10, +0.08] |

**Root cause: no_deficit.** No deficit: every stage supported. every chain stage supported for this readout

## Provenance
- adapter: `precomputed(perturb_cite features exported on the source machine)`
- python: `3.13.18`
- runtime: 93.9 s
