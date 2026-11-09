# ANM diagnosis: teddy_mm_frangieh2021_ifng_pathway_headsplit

- status: **sealed**, exploratory evidence (paper p. 8): Designed after the sealed headfix run, which located a resolved loss at Z -> HEAD_LAT. This run splits that step into teddy_mm's own operations (unit-length input, linear 1, SiLU, linear 2, SiLU, linear 3) to locate it further. A method chosen after observing responses gives exploratory evidence and leaves the headfix classifications unchanged (ANM paper p. 8).

- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 5171 cells in all
- chain: RAW → TOKENS → Z → Z_UNIT → L1 → H1 → L2 → H2 → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `795ef90f5a1b0d93…`; anm_diag 0.3.0
- verdict rule (two-sided equivalence, paper §2.3 and Table 2): **support** when some declared readout class has its whole 95% interval of u inside (−0.1, 0.1); **reject** when every class lies wholly outside it (this rejects the declared readout classes on that state, Eq. 5, not the representation alone); otherwise **unresolved**. A step adds a resolved loss when its change in |u| (distance from 0) has a lower bound above the step margin.
- u is the overlap-weighted average same-state difference divided by Delta: support tests a necessary condition for mean sufficiency. Positive control in the paper's sense: Delta's interval excludes 0. Method control: the readout itself as a state (u should be 0). Finite contrasts only; no first-order claim.

**Checks**
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

State response, reported beside the readout response Delta; a joint change is association under intervention, not mediation (paper SI S1.3).

| state | AUC |
|---|---|
| RAW | 0.895 |
| TOKENS | 0.889 |
| Z | 0.884 |
| Z_UNIT | 0.884 |
| L1 | 0.889 |
| H1 | 0.889 |
| L2 | 0.867 |
| H2 | 0.860 |
| HEAD_LAT | 0.843 |
| HEAD_OUT | 0.811 |

| step | change, 95% interval |
|---|---|
| RAW->TOKENS | [-0.010, -0.001] |
| TOKENS->Z | [-0.011, +0.002] |
| Z->Z_UNIT | [-0.001, -0.000] |
| Z_UNIT->L1 | [+0.001, +0.010] |
| L1->H1 | [-0.000, +0.001] |
| H1->L2 | [-0.031, -0.013] |
| L2->H2 | [-0.009, -0.004] |
| H2->HEAD_LAT | [-0.022, -0.013] |
| HEAD_LAT->HEAD_OUT | [-0.040, -0.024] |

## Readout MHC1

Response Delta = -0.138 [-0.182, -0.094]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.23 (hgb) | [-0.71, +0.15] | **unresolved** | ridge -0.27 [-0.71, +0.11]; hgb -0.23 [-0.71, +0.15] |
| TOKENS | chain | -0.25 (hgb) | [-0.70, +0.12] | **unresolved** | ridge -0.29 [-0.70, +0.07]; hgb -0.25 [-0.70, +0.12] |
| Z | chain | -0.20 (ridge) | [-0.62, +0.20] | **unresolved** | ridge -0.20 [-0.62, +0.20]; hgb -0.37 [-0.92, +0.02] |
| Z_UNIT | chain | -0.18 (ridge) | [-0.60, +0.21] | **unresolved** | ridge -0.18 [-0.60, +0.21]; hgb -0.38 [-0.94, +0.07] |
| L1 | chain | -0.09 (hgb) | [-0.53, +0.24] | **unresolved** | ridge -0.29 [-0.72, +0.03]; hgb -0.09 [-0.53, +0.24] |
| H1 | chain | -0.09 (hgb) | [-0.53, +0.24] | **unresolved** | ridge -0.26 [-0.70, +0.06]; hgb -0.09 [-0.53, +0.24] |
| L2 | chain | +0.08 (hgb) | [-0.24, +0.37] | **unresolved** | ridge -0.19 [-0.54, +0.04]; hgb +0.08 [-0.24, +0.37] |
| H2 | chain | +0.05 (hgb) | [-0.28, +0.36] | **unresolved** | ridge -0.21 [-0.62, +0.08]; hgb +0.05 [-0.28, +0.36] |
| HEAD_LAT | chain | +0.01 (ridge) | [-0.34, +0.27] | **unresolved** | ridge +0.01 [-0.34, +0.27]; hgb +0.06 [-0.24, +0.30] |
| HEAD_OUT | chain | +0.19 (ridge) | [-0.10, +0.38] | **unresolved** | ridge +0.19 [-0.10, +0.38]; hgb +0.20 [-0.11, +0.44] |
| ADT_MHC1 | method control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb +0.00 [-0.02, +0.02] |

| step | change in u (best class) | change in abs(u), used for the root cause | same class |
|---|---|---|---|
| RAW->TOKENS | [-0.13, +0.07] | [-0.08, +0.12] | ridge [-0.21, +0.22]; hgb [-0.13, +0.07] |
| TOKENS->Z | [-0.09, +0.23] | [-0.23, +0.11] | ridge [-0.14, +0.31]; hgb [-0.32, +0.05] |
| Z->Z_UNIT | [-0.00, +0.03] | [-0.03, +0.01] | ridge [-0.00, +0.03]; hgb [-0.12, +0.11] |
| Z_UNIT->L1 | [-0.11, +0.29] | [-0.26, +0.15] | ridge [-0.24, -0.01]; hgb [+0.09, +0.54] |
| L1->H1 | [+0.00, +0.00] | [+0.00, +0.00] | ridge [+0.01, +0.06]; hgb [+0.00, +0.00] |
| H1->L2 | [-0.01, +0.44] | [-0.37, +0.24] | ridge [-0.07, +0.20]; hgb [-0.01, +0.44] |
| L2->H2 | [-0.12, +0.05] | [-0.10, +0.08] | ridge [-0.16, +0.11]; hgb [-0.12, +0.05] |
| H2->HEAD_LAT | [-0.24, +0.11] | [-0.19, +0.18] | ridge [+0.12, +0.35]; hgb [-0.15, +0.15] |
| HEAD_LAT->HEAD_OUT | [+0.04, +0.32] | [-0.27, +0.23] | ridge [+0.04, +0.32]; hgb [-0.01, +0.24] |

**Root cause: reference_unresolved.** The reference stage is not supported. Rejected: the input view itself does not reproduce the response for any declared readout class (a limit of the data or the classes). Unresolved: the data cannot tell. Either way nothing is charged to the system. reference stage 'RAW' is unresolved: the data cannot tell whether the input view reproduces this response; no claim about later stages

## Readout PDL1

Response Delta = -0.849 [-0.928, -0.762]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.01 (hgb) | [-0.06, +0.04] | **support** | ridge -0.01 [-0.09, +0.05]; hgb -0.01 [-0.06, +0.04] |
| TOKENS | chain | +0.01 (hgb) | [-0.03, +0.05] | **support** | ridge +0.00 [-0.07, +0.08]; hgb +0.01 [-0.03, +0.05] |
| Z | chain | +0.04 (ridge) | [-0.01, +0.09] | **support** | ridge +0.04 [-0.01, +0.09]; hgb +0.12 [+0.03, +0.20] |
| Z_UNIT | chain | +0.04 (ridge) | [-0.01, +0.09] | **support** | ridge +0.04 [-0.01, +0.09]; hgb +0.12 [+0.03, +0.20] |
| L1 | chain | +0.03 (ridge) | [-0.03, +0.09] | **support** | ridge +0.03 [-0.03, +0.09]; hgb +0.14 [+0.07, +0.21] |
| H1 | chain | +0.04 (ridge) | [-0.02, +0.10] | **support** | ridge +0.04 [-0.02, +0.10]; hgb +0.14 [+0.07, +0.21] |
| L2 | chain | +0.04 (ridge) | [-0.04, +0.10] | **unresolved** | ridge +0.04 [-0.04, +0.10]; hgb +0.35 [+0.30, +0.40] |
| H2 | chain | +0.05 (ridge) | [-0.00, +0.10] | **unresolved** | ridge +0.05 [-0.00, +0.10]; hgb +0.36 [+0.31, +0.41] |
| HEAD_LAT | chain | +0.20 (ridge) | [+0.13, +0.26] | **reject** | ridge +0.20 [+0.13, +0.26]; hgb +0.43 [+0.38, +0.48] |
| HEAD_OUT | chain | +0.46 (hgb) | [+0.41, +0.50] | **reject** | ridge +0.47 [+0.41, +0.51]; hgb +0.46 [+0.41, +0.50] |
| ADT_PDL1 | method control | +0.00 (ridge) | [+0.00, +0.00] | **support** | ridge +0.00 [+0.00, +0.00]; hgb +0.00 [-0.00, +0.00] |

| step | change in u (best class) | change in abs(u), used for the root cause | same class |
|---|---|---|---|
| RAW->TOKENS | [-0.01, +0.05] | [-0.04, +0.03] | ridge [-0.01, +0.05]; hgb [-0.01, +0.05] |
| TOKENS->Z | [-0.00, +0.07] | [-0.03, +0.06] | ridge [-0.01, +0.09]; hgb [+0.05, +0.16] |
| Z->Z_UNIT | [-0.01, +0.00] | [-0.01, +0.00] | ridge [-0.01, +0.00]; hgb [-0.02, +0.01] |
| Z_UNIT->L1 | [-0.04, +0.01] | [-0.03, +0.02] | ridge [-0.04, +0.01]; hgb [-0.02, +0.07] |
| L1->H1 | [+0.00, +0.02] | [-0.01, +0.02] | ridge [+0.00, +0.02]; hgb [+0.00, +0.00] |
| H1->L2 | [-0.03, +0.02] | [-0.02, +0.03] | ridge [-0.03, +0.02]; hgb [+0.18, +0.25] |
| L2->H2 | [-0.02, +0.06] | [-0.03, +0.05] | ridge [-0.02, +0.06]; hgb [-0.01, +0.02] |
| H2->HEAD_LAT | [+0.11, +0.19] | [+0.11, +0.19] | ridge [+0.11, +0.19]; hgb [+0.04, +0.10] |
| HEAD_LAT->HEAD_OUT | [+0.22, +0.29] | [+0.22, +0.29] | ridge [+0.24, +0.29]; hgb [+0.00, +0.06] |

**Root cause: loss_at.** Representation loss: a stage is rejected after a supported stage. distinction merged between 'H1' and 'HEAD_LAT' (via unresolved ['L2', 'H2'])

## Readout HLAE

Response Delta = -0.079 [-0.115, -0.046]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | +0.00 (hgb) | [-0.68, +0.41] | **unresolved** | ridge -0.27 [-1.28, +0.25]; hgb +0.00 [-0.68, +0.41] |
| TOKENS | chain | +0.11 (hgb) | [-0.51, +0.47] | **unresolved** | ridge -0.16 [-1.08, +0.35]; hgb +0.11 [-0.51, +0.47] |
| Z | chain | +0.18 (hgb) | [-0.49, +0.57] | **unresolved** | ridge +0.20 [-0.39, +0.61]; hgb +0.18 [-0.49, +0.57] |
| Z_UNIT | chain | +0.19 (hgb) | [-0.65, +0.65] | **unresolved** | ridge +0.22 [-0.37, +0.62]; hgb +0.19 [-0.65, +0.65] |
| L1 | chain | +0.31 (ridge) | [-0.39, +0.71] | **unresolved** | ridge +0.31 [-0.39, +0.71]; hgb +0.56 [-0.23, +0.96] |
| H1 | chain | +0.31 (ridge) | [-0.41, +0.71] | **unresolved** | ridge +0.31 [-0.41, +0.71]; hgb +0.56 [-0.23, +0.96] |
| L2 | chain | +0.38 (ridge) | [-0.17, +0.73] | **unresolved** | ridge +0.38 [-0.17, +0.73]; hgb +0.87 [+0.42, +1.21] |
| H2 | chain | +0.82 (ridge) | [+0.30, +1.17] | **reject** | ridge +0.82 [+0.30, +1.17]; hgb +0.91 [+0.37, +1.36] |
| HEAD_LAT | chain | +0.94 (ridge) | [+0.51, +1.27] | **reject** | ridge +0.94 [+0.51, +1.27]; hgb +1.37 [+1.08, +1.80] |
| HEAD_OUT | chain | +1.36 (hgb) | [+1.09, +1.88] | **reject** | ridge +1.41 [+1.13, +1.87]; hgb +1.36 [+1.09, +1.88] |
| ADT_HLAE | method control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb +0.01 [-0.02, +0.04] |

| step | change in u (best class) | change in abs(u), used for the root cause | same class |
|---|---|---|---|
| RAW->TOKENS | [-0.09, +0.32] | [-0.27, +0.24] | ridge [-0.23, +0.40]; hgb [-0.09, +0.32] |
| TOKENS->Z | [-0.17, +0.30] | [-0.18, +0.28] | ridge [-0.04, +1.03]; hgb [-0.17, +0.30] |
| Z->Z_UNIT | [-0.25, +0.19] | [-0.14, +0.23] | ridge [-0.01, +0.06]; hgb [-0.25, +0.19] |
| Z_UNIT->L1 | [-0.13, +0.47] | [-0.37, +0.35] | ridge [-0.10, +0.26]; hgb [+0.07, +0.77] |
| L1->H1 | [-0.04, +0.02] | [-0.03, +0.03] | ridge [-0.04, +0.02]; hgb [+0.00, +0.00] |
| H1->L2 | [-0.20, +0.32] | [-0.27, +0.28] | ridge [-0.20, +0.32]; hgb [+0.10, +0.70] |
| L2->H2 | [+0.31, +0.68] | [+0.16, +0.65] | ridge [+0.31, +0.68]; hgb [-0.11, +0.22] |
| H2->HEAD_LAT | [-0.06, +0.28] | [-0.06, +0.28] | ridge [-0.06, +0.28]; hgb [+0.08, +1.08] |
| HEAD_LAT->HEAD_OUT | [+0.16, +0.99] | [+0.16, +0.99] | ridge [+0.25, +1.01]; hgb [-0.19, +0.20] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'RAW' is unresolved (the readout is not fully explained even by the input view), but 'L2' -> 'H2' adds a resolved loss (change in |u| > 0.1)

## Provenance
- adapter: `perturb_cite`
- sha256_rna_h5ad: `cc42ef38bcf703a0`
- sha256_adt_h5ad: `1f85827b5afad11a`
- sha256_vocab: `e184b550fa2c3971`
- sha256_medians: `3d3637fa24ef81bf`
- sha256_head_ckpt: `b0c543e285123205`
- sha256_teddy_states: `216b77a003c62ad4`
- python: `3.11.13`
- runtime: 1533.3 s
