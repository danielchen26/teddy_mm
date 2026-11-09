# ANM diagnosis: teddy_mm_papalexi2021_ifng_pathway_pdl1_headsplit

- status: **sealed**, exploratory evidence (paper p. 8): Designed after the sealed headfix run, which located a resolved loss at Z -> HEAD_LAT. This run splits that step into teddy_mm's own operations (unit-length input, linear 1, SiLU, linear 2, SiLU, linear 3) to locate it further. A method chosen after observing responses gives exploratory evidence and leaves the headfix classifications unchanged (ANM paper p. 8).

- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 7478 cells in all
- chain: RAW → TOKENS → Z → Z_UNIT → L1 → H1 → L2 → H2 → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `840a4b88a1353a63…`; anm_diag 0.3.0
- verdict rule (two-sided equivalence, paper §2.3 and Table 2): **support** when some declared readout class has its whole 95% interval of u inside (−0.1, 0.1); **reject** when every class lies wholly outside it (this rejects the declared readout classes on that state, Eq. 5, not the representation alone); otherwise **unresolved**. A step adds a resolved loss when its change in |u| (distance from 0) has a lower bound above the step margin.
- u is the overlap-weighted average same-state difference divided by Delta: support tests a necessary condition for mean sufficiency. Positive control in the paper's sense: Delta's interval excludes 0. Method control: the readout itself as a state (u should be 0). Finite contrasts only; no first-order claim.

**Checks**
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

State response, reported beside the readout response Delta; a joint change is association under intervention, not mediation (paper SI S1.3).

| state | AUC |
|---|---|
| RAW | 0.767 |
| TOKENS | 0.752 |
| Z | 0.757 |
| Z_UNIT | 0.755 |
| L1 | 0.728 |
| H1 | 0.727 |
| L2 | 0.651 |
| H2 | 0.652 |
| HEAD_LAT | 0.636 |
| HEAD_OUT | 0.577 |

| step | change, 95% interval |
|---|---|
| RAW->TOKENS | [-0.035, -0.003] |
| TOKENS->Z | [-0.003, +0.016] |
| Z->Z_UNIT | [-0.004, +0.000] |
| Z_UNIT->L1 | [-0.043, -0.006] |
| L1->H1 | [-0.003, -0.001] |
| H1->L2 | [-0.101, -0.040] |
| L2->H2 | [-0.006, +0.007] |
| H2->HEAD_LAT | [-0.025, -0.004] |
| HEAD_LAT->HEAD_OUT | [-0.077, -0.030] |

## Readout PDL1

Response Delta = -0.693 [-0.856, -0.469]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | +0.10 (hgb) | [+0.00, +0.23] | **unresolved** | ridge +0.12 [-0.05, +0.37]; hgb +0.10 [+0.00, +0.23] |
| TOKENS | chain | +0.12 (ridge) | [-0.04, +0.36] | **unresolved** | ridge +0.12 [-0.04, +0.36]; hgb +0.13 [-0.01, +0.32] |
| Z | chain | +0.15 (ridge) | [+0.04, +0.32] | **unresolved** | ridge +0.15 [+0.04, +0.32]; hgb +0.24 [+0.16, +0.38] |
| Z_UNIT | chain | +0.15 (ridge) | [+0.04, +0.33] | **unresolved** | ridge +0.15 [+0.04, +0.33]; hgb +0.24 [+0.15, +0.39] |
| L1 | chain | +0.16 (ridge) | [+0.06, +0.32] | **unresolved** | ridge +0.16 [+0.06, +0.32]; hgb +0.34 [+0.27, +0.50] |
| H1 | chain | +0.17 (ridge) | [+0.09, +0.31] | **unresolved** | ridge +0.17 [+0.09, +0.31]; hgb +0.34 [+0.27, +0.50] |
| L2 | chain | +0.19 (ridge) | [+0.12, +0.34] | **reject** | ridge +0.19 [+0.12, +0.34]; hgb +0.40 [+0.35, +0.53] |
| H2 | chain | +0.22 (ridge) | [+0.14, +0.35] | **reject** | ridge +0.22 [+0.14, +0.35]; hgb +0.40 [+0.35, +0.53] |
| HEAD_LAT | chain | +0.35 (ridge) | [+0.30, +0.45] | **reject** | ridge +0.35 [+0.30, +0.45]; hgb +0.43 [+0.38, +0.56] |
| HEAD_OUT | chain | +0.47 (ridge) | [+0.43, +0.58] | **reject** | ridge +0.47 [+0.43, +0.58]; hgb +0.51 [+0.47, +0.62] |
| ADT_PDL1 | method control | -0.00 (ridge) | [-0.00, +0.00] | **support** | ridge -0.00 [-0.00, +0.00]; hgb +0.00 [-0.00, +0.00] |

| step | change in u (best class) | change in abs(u), used for the root cause | same class |
|---|---|---|---|
| RAW->TOKENS | [-0.06, +0.13] | [-0.05, +0.13] | ridge [-0.02, +0.03]; hgb [-0.02, +0.10] |
| TOKENS->Z | [-0.05, +0.10] | [-0.05, +0.08] | ridge [-0.05, +0.10]; hgb [+0.05, +0.17] |
| Z->Z_UNIT | [+0.00, +0.01] | [+0.00, +0.01] | ridge [+0.00, +0.01]; hgb [-0.01, +0.01] |
| Z_UNIT->L1 | [-0.01, +0.03] | [-0.01, +0.03] | ridge [-0.01, +0.03]; hgb [+0.09, +0.14] |
| L1->H1 | [-0.01, +0.03] | [-0.01, +0.03] | ridge [-0.01, +0.03]; hgb [-0.00, +0.00] |
| H1->L2 | [+0.02, +0.04] | [+0.02, +0.04] | ridge [+0.02, +0.04]; hgb [+0.02, +0.09] |
| L2->H2 | [+0.01, +0.04] | [+0.01, +0.04] | ridge [+0.01, +0.04]; hgb [-0.01, +0.01] |
| H2->HEAD_LAT | [+0.09, +0.17] | [+0.09, +0.17] | ridge [+0.09, +0.17]; hgb [+0.02, +0.04] |
| HEAD_LAT->HEAD_OUT | [+0.11, +0.15] | [+0.11, +0.15] | ridge [+0.11, +0.15]; hgb [+0.06, +0.10] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'RAW' is unresolved (the readout is not fully explained even by the input view), but 'HEAD_LAT' -> 'HEAD_OUT' adds a resolved loss (change in |u| > 0.1)

## Readout PDL2

Response Delta = -0.287 [-0.338, -0.229]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | +0.08 (ridge) | [-0.23, +0.34] | **unresolved** | ridge +0.08 [-0.23, +0.34]; hgb +0.14 [-0.14, +0.38] |
| TOKENS | chain | +0.08 (ridge) | [-0.19, +0.31] | **unresolved** | ridge +0.08 [-0.19, +0.31]; hgb +0.14 [-0.13, +0.37] |
| Z | chain | +0.13 (ridge) | [-0.11, +0.36] | **unresolved** | ridge +0.13 [-0.11, +0.36]; hgb +0.23 [+0.04, +0.42] |
| Z_UNIT | chain | +0.14 (ridge) | [-0.10, +0.36] | **unresolved** | ridge +0.14 [-0.10, +0.36]; hgb +0.22 [+0.03, +0.41] |
| L1 | chain | +0.16 (ridge) | [-0.05, +0.36] | **unresolved** | ridge +0.16 [-0.05, +0.36]; hgb +0.44 [+0.24, +0.68] |
| H1 | chain | +0.17 (ridge) | [-0.02, +0.37] | **unresolved** | ridge +0.17 [-0.02, +0.37]; hgb +0.44 [+0.24, +0.68] |
| L2 | chain | +0.20 (ridge) | [+0.02, +0.40] | **unresolved** | ridge +0.20 [+0.02, +0.40]; hgb +0.52 [+0.34, +0.74] |
| H2 | chain | +0.28 (ridge) | [+0.10, +0.47] | **unresolved** | ridge +0.28 [+0.10, +0.47]; hgb +0.53 [+0.36, +0.74] |
| HEAD_LAT | chain | +0.40 (ridge) | [+0.23, +0.62] | **reject** | ridge +0.40 [+0.23, +0.62]; hgb +0.51 [+0.35, +0.72] |
| HEAD_OUT | chain | +0.55 (ridge) | [+0.39, +0.79] | **reject** | ridge +0.55 [+0.39, +0.79]; hgb +0.56 [+0.39, +0.80] |
| ADT_PDL2 | method control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb +0.00 [-0.00, +0.01] |

| step | change in u (best class) | change in abs(u), used for the root cause | same class |
|---|---|---|---|
| RAW->TOKENS | [-0.05, +0.08] | [-0.08, +0.03] | ridge [-0.05, +0.08]; hgb [-0.02, +0.02] |
| TOKENS->Z | [+0.01, +0.14] | [-0.12, +0.10] | ridge [+0.01, +0.14]; hgb [+0.03, +0.19] |
| Z->Z_UNIT | [-0.00, +0.01] | [-0.01, +0.01] | ridge [-0.00, +0.01]; hgb [-0.03, +0.02] |
| Z_UNIT->L1 | [-0.00, +0.06] | [-0.05, +0.05] | ridge [-0.00, +0.06]; hgb [+0.18, +0.27] |
| L1->H1 | [+0.01, +0.03] | [-0.02, +0.03] | ridge [+0.01, +0.03]; hgb [-0.00, +0.00] |
| H1->L2 | [+0.01, +0.07] | [-0.00, +0.07] | ridge [+0.01, +0.07]; hgb [+0.05, +0.12] |
| L2->H2 | [+0.04, +0.10] | [+0.04, +0.10] | ridge [+0.04, +0.10]; hgb [-0.01, +0.02] |
| H2->HEAD_LAT | [+0.10, +0.16] | [+0.10, +0.16] | ridge [+0.10, +0.16]; hgb [-0.03, +0.01] |
| HEAD_LAT->HEAD_OUT | [+0.13, +0.18] | [+0.13, +0.18] | ridge [+0.13, +0.18]; hgb [+0.02, +0.09] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'RAW' is unresolved (the readout is not fully explained even by the input view), but 'HEAD_LAT' -> 'HEAD_OUT' adds a resolved loss (change in |u| > 0.1)

## Provenance
- adapter: `perturb_cite`
- sha256_rna_h5ad: `03e9602d124c2612`
- sha256_adt_h5ad: `e1e293c75dfebe09`
- sha256_vocab: `e184b550fa2c3971`
- sha256_medians: `3d3637fa24ef81bf`
- sha256_head_ckpt: `b0c543e285123205`
- sha256_teddy_states: `d2b2974da8c700fc`
- python: `3.11.13`
- runtime: 1437.1 s
