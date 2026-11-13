# ANM diagnosis: teddy_mm_frangieh2021_ifng_pathway

- status: **sealed**
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK1', 'JAK2', 'STAT1']` vs `['control']` · 1171 treated cells in 5 clusters, 5171 cells in all
- chain: RAW → TOKENS → Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `90c55801b45dd5b3…`; anm_diag 0.1.0

**Checks**
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| RAW | 0.895 |
| TOKENS | 0.889 |
| Z | 0.884 |
| HEAD_LAT | 0.832 |
| HEAD_OUT | 0.768 |
| SET | 0.886 |
| MRNA_HLA | 0.810 |
| MRNA_PDL1 | 0.358 |
| MRNA_HLAE | 0.832 |

| step | change, 95% interval |
|---|---|
| RAW->TOKENS | [-0.010, -0.001] |
| TOKENS->Z | [-0.011, +0.002] |
| Z->HEAD_LAT | [-0.070, -0.036] |
| HEAD_LAT->HEAD_OUT | [-0.082, -0.044] |
| TOKENS->SET | [-0.009, +0.002] |
| RAW->MRNA_HLA | [-0.115, -0.053] |
| RAW->MRNA_PDL1 | [-0.652, -0.410] |
| RAW->MRNA_HLAE | [-0.081, -0.046] |

## Readout MHC1

Response Delta = -0.138 [-0.182, -0.094]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.27 (ridge) | [-0.71, +0.11] | **unresolved** | ridge -0.27 [-0.71, +0.11]; hgb -0.23 [-0.71, +0.15] |
| TOKENS | chain | -0.29 (ridge) | [-0.70, +0.07] | **support** | ridge -0.29 [-0.70, +0.07]; hgb -0.25 [-0.70, +0.12] |
| Z | chain | -0.37 (hgb) | [-0.92, +0.02] | **support** | ridge -0.20 [-0.62, +0.20]; hgb -0.37 [-0.92, +0.02] |
| HEAD_LAT | chain | +0.12 (ridge) | [-0.16, +0.34] | **unresolved** | ridge +0.12 [-0.16, +0.34]; hgb +0.33 [+0.04, +0.57] |
| HEAD_OUT | chain | +0.56 (hgb) | [+0.34, +0.70] | **reject** | ridge +0.57 [+0.30, +0.75]; hgb +0.56 [+0.34, +0.70] |
| SET | probe | -0.30 (hgb) | [-0.79, +0.08] | **support** | ridge -0.23 [-0.60, +0.10]; hgb -0.30 [-0.79, +0.08] |
| MRNA_HLA | probe | -0.59 (ridge) | [-1.43, -0.17] | **support** | ridge -0.59 [-1.43, -0.17]; hgb -0.53 [-1.25, -0.15] |
| MRNA_PDL1 | probe | +1.00 (hgb) | [+0.93, +1.04] | **reject** | ridge +1.03 [+0.96, +1.06]; hgb +1.00 [+0.93, +1.04] |
| MRNA_HLAE | probe | +0.08 (hgb) | [-0.35, +0.32] | **unresolved** | ridge +0.17 [-0.37, +0.45]; hgb +0.08 [-0.35, +0.32] |
| ADT_MHC1 | positive_control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb +0.00 [-0.02, +0.02] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [-0.21, +0.22] | ridge [-0.21, +0.22]; hgb [-0.13, +0.07] |
| TOKENS->Z | [-0.42, +0.22] | ridge [-0.14, +0.31]; hgb [-0.32, +0.05] |
| Z->HEAD_LAT | [+0.26, +0.85] | ridge [+0.09, +0.53]; hgb [+0.46, +1.09] |
| HEAD_LAT->HEAD_OUT | [+0.26, +0.60] | ridge [+0.31, +0.58]; hgb [+0.03, +0.39] |
| TOKENS->SET | [-0.23, +0.18] | ridge [-0.09, +0.22]; hgb [-0.16, +0.05] |
| RAW->MRNA_HLA | [-0.90, +0.03] | ridge [-0.90, +0.03]; hgb [-0.80, -0.02] |
| RAW->MRNA_PDL1 | [+0.89, +1.69] | ridge [+0.90, +1.71]; hgb [+0.85, +1.67] |
| RAW->MRNA_HLAE | [-0.01, +0.60] | ridge [+0.07, +0.70]; hgb [+0.01, +0.57] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'RAW' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Readout PDL1

Response Delta = -0.849 [-0.928, -0.762]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.01 (hgb) | [-0.06, +0.04] | **support** | ridge -0.01 [-0.09, +0.05]; hgb -0.01 [-0.06, +0.04] |
| TOKENS | chain | +0.01 (hgb) | [-0.03, +0.05] | **support** | ridge +0.00 [-0.07, +0.08]; hgb +0.01 [-0.03, +0.05] |
| Z | chain | +0.04 (ridge) | [-0.01, +0.09] | **support** | ridge +0.04 [-0.01, +0.09]; hgb +0.12 [+0.03, +0.20] |
| HEAD_LAT | chain | +0.31 (ridge) | [+0.25, +0.35] | **reject** | ridge +0.31 [+0.25, +0.35]; hgb +0.47 [+0.40, +0.52] |
| HEAD_OUT | chain | +0.58 (hgb) | [+0.51, +0.63] | **reject** | ridge +0.60 [+0.53, +0.65]; hgb +0.58 [+0.51, +0.63] |
| SET | probe | +0.01 (ridge) | [-0.08, +0.08] | **support** | ridge +0.01 [-0.08, +0.08]; hgb +0.03 [-0.06, +0.11] |
| MRNA_HLA | probe | +0.47 (hgb) | [+0.34, +0.54] | **reject** | ridge +0.55 [+0.43, +0.63]; hgb +0.47 [+0.34, +0.54] |
| MRNA_PDL1 | probe | +0.93 (hgb) | [+0.92, +0.95] | **reject** | ridge +0.94 [+0.92, +0.95]; hgb +0.93 [+0.92, +0.95] |
| MRNA_HLAE | probe | +0.49 (hgb) | [+0.40, +0.56] | **reject** | ridge +0.54 [+0.45, +0.60]; hgb +0.49 [+0.40, +0.56] |
| ADT_PDL1 | positive_control | +0.00 (ridge) | [+0.00, +0.00] | **support** | ridge +0.00 [+0.00, +0.00]; hgb +0.00 [-0.00, +0.00] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [-0.01, +0.05] | ridge [-0.01, +0.05]; hgb [-0.01, +0.05] |
| TOKENS->Z | [-0.00, +0.07] | ridge [-0.01, +0.09]; hgb [+0.05, +0.16] |
| Z->HEAD_LAT | [+0.22, +0.30] | ridge [+0.22, +0.30]; hgb [+0.29, +0.42] |
| HEAD_LAT->HEAD_OUT | [+0.24, +0.30] | ridge [+0.26, +0.32]; hgb [+0.09, +0.13] |
| TOKENS->SET | [-0.06, +0.04] | ridge [-0.03, +0.03]; hgb [-0.04, +0.07] |
| RAW->MRNA_HLA | [+0.37, +0.56] | ridge [+0.50, +0.62]; hgb [+0.37, +0.56] |
| RAW->MRNA_PDL1 | [+0.90, +0.99] | ridge [+0.88, +1.02]; hgb [+0.90, +0.99] |
| RAW->MRNA_HLAE | [+0.43, +0.57] | ridge [+0.51, +0.60]; hgb [+0.43, +0.57] |

**Root cause: loss_at.** Representation loss: a stage is rejected after a supported stage. distinction merged between 'Z' and 'HEAD_LAT'

## Readout HLAE

Response Delta = -0.079 [-0.115, -0.046]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.27 (ridge) | [-1.28, +0.25] | **unresolved** | ridge -0.27 [-1.28, +0.25]; hgb +0.00 [-0.68, +0.41] |
| TOKENS | chain | -0.16 (ridge) | [-1.08, +0.35] | **unresolved** | ridge -0.16 [-1.08, +0.35]; hgb +0.11 [-0.51, +0.47] |
| Z | chain | +0.18 (hgb) | [-0.49, +0.57] | **unresolved** | ridge +0.20 [-0.39, +0.61]; hgb +0.18 [-0.49, +0.57] |
| HEAD_LAT | chain | +1.19 (hgb) | [+0.80, +1.56] | **reject** | ridge +1.27 [+0.98, +1.63]; hgb +1.19 [+0.80, +1.56] |
| HEAD_OUT | chain | +1.15 (hgb) | [+0.74, +1.57] | **reject** | ridge +1.27 [+0.93, +1.69]; hgb +1.15 [+0.74, +1.57] |
| SET | probe | -0.04 (ridge) | [-0.82, +0.41] | **unresolved** | ridge -0.04 [-0.82, +0.41]; hgb +0.06 [-0.70, +0.45] |
| MRNA_HLA | probe | +0.53 (hgb) | [-0.13, +0.88] | **unresolved** | ridge +0.75 [+0.28, +1.01]; hgb +0.53 [-0.13, +0.88] |
| MRNA_PDL1 | probe | +1.14 (hgb) | [+1.04, +1.36] | **reject** | ridge +1.17 [+1.09, +1.36]; hgb +1.14 [+1.04, +1.36] |
| MRNA_HLAE | probe | +0.23 (hgb) | [-0.42, +0.56] | **unresolved** | ridge +0.28 [-0.37, +0.63]; hgb +0.23 [-0.42, +0.56] |
| ADT_HLAE | positive_control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb +0.01 [-0.02, +0.04] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [-0.23, +0.40] | ridge [-0.23, +0.40]; hgb [-0.09, +0.32] |
| TOKENS->Z | [-0.02, +0.88] | ridge [-0.04, +1.03]; hgb [-0.17, +0.30] |
| Z->HEAD_LAT | [+0.62, +1.66] | ridge [+0.64, +1.72]; hgb [+0.62, +1.66] |
| HEAD_LAT->HEAD_OUT | [-0.27, +0.22] | ridge [-0.17, +0.19]; hgb [-0.27, +0.22] |
| TOKENS->SET | [-0.07, +0.45] | ridge [-0.07, +0.45]; hgb [-0.32, +0.11] |
| RAW->MRNA_HLA | [+0.38, +1.32] | ridge [+0.58, +1.72]; hgb [+0.16, +0.89] |
| RAW->MRNA_PDL1 | [+0.86, +2.63] | ridge [+0.89, +2.63]; hgb [+0.69, +1.97] |
| RAW->MRNA_HLAE | [-0.09, +1.16] | ridge [-0.04, +1.21]; hgb [-0.25, +0.68] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'RAW' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `perturb_cite`
- sha256_rna_h5ad: `cc42ef38bcf703a0`
- sha256_adt_h5ad: `1f85827b5afad11a`
- sha256_vocab: `e184b550fa2c3971`
- sha256_medians: `3d3637fa24ef81bf`
- sha256_head_ckpt: `b0c543e285123205`
- sha256_teddy_states: `216b77a003c62ad4`
- python: `3.11.15`
- runtime: 2059.4 s
