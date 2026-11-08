# ANM diagnosis: teddy_mm_papalexi2021_ifng_pathway_pdl1

- status: **sealed**
- system: TEDDY-G-70M + teddy_mm phase-1 head (cite_phase1_official/best.pt) (frozen) · inputs ['rna_counts']
- source: crispr_ko · treated `['IFNGR1', 'IFNGR2', 'JAK2', 'STAT1', 'IRF1']` vs `['control']` · 5092 treated cells in 5 clusters, 7478 cells in all
- chain: RAW → TOKENS → Z → HEAD_LAT → HEAD_OUT
- margin u = 0.1, alpha = 0.05; declaration sha256 `afd4ed69e1c5db94…`; anm_diag 0.1.0

**Checks**
- [K1] warning: only 5 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| RAW | 0.767 |
| TOKENS | 0.753 |
| Z | 0.757 |
| HEAD_LAT | 0.631 |
| HEAD_OUT | 0.568 |
| SET | 0.749 |
| MRNA_PDL1 | 0.292 |
| MRNA_PDL2 | 0.249 |

| step | change, 95% interval |
|---|---|
| RAW->TOKENS | [-0.034, -0.002] |
| TOKENS->Z | [-0.004, +0.015] |
| Z->HEAD_LAT | [-0.165, -0.061] |
| HEAD_LAT->HEAD_OUT | [-0.082, -0.034] |
| TOKENS->SET | [-0.010, +0.002] |
| RAW->MRNA_PDL1 | [-0.589, -0.317] |
| RAW->MRNA_PDL2 | [-0.629, -0.357] |

## Readout PDL1

Response Delta = -0.693 [-0.856, -0.469]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | +0.10 (hgb) | [+0.00, +0.23] | **unresolved** | ridge +0.12 [-0.05, +0.37]; hgb +0.10 [+0.00, +0.23] |
| TOKENS | chain | +0.12 (ridge) | [-0.04, +0.37] | **unresolved** | ridge +0.12 [-0.04, +0.37]; hgb +0.13 [-0.01, +0.32] |
| Z | chain | +0.15 (ridge) | [+0.04, +0.32] | **unresolved** | ridge +0.15 [+0.04, +0.32]; hgb +0.24 [+0.16, +0.38] |
| HEAD_LAT | chain | +0.36 (ridge) | [+0.31, +0.45] | **reject** | ridge +0.36 [+0.31, +0.45]; hgb +0.43 [+0.40, +0.51] |
| HEAD_OUT | chain | +0.47 (hgb) | [+0.44, +0.56] | **reject** | ridge +0.49 [+0.46, +0.57]; hgb +0.47 [+0.44, +0.56] |
| SET | probe | +0.14 (ridge) | [-0.05, +0.40] | **unresolved** | ridge +0.14 [-0.05, +0.40]; hgb +0.17 [-0.01, +0.41] |
| MRNA_PDL1 | probe | +0.69 (hgb) | [+0.63, +0.78] | **reject** | ridge +0.72 [+0.65, +0.82]; hgb +0.69 [+0.63, +0.78] |
| MRNA_PDL2 | probe | +0.81 (hgb) | [+0.69, +0.91] | **reject** | ridge +0.83 [+0.70, +0.94]; hgb +0.81 [+0.69, +0.91] |
| ADT_PDL1 | positive_control | -0.00 (ridge) | [-0.00, +0.00] | **support** | ridge -0.00 [-0.00, +0.00]; hgb +0.00 [-0.00, +0.00] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [-0.06, +0.13] | ridge [-0.02, +0.02]; hgb [-0.02, +0.10] |
| TOKENS->Z | [-0.05, +0.10] | ridge [-0.05, +0.10]; hgb [+0.05, +0.17] |
| Z->HEAD_LAT | [+0.13, +0.29] | ridge [+0.13, +0.29]; hgb [+0.12, +0.26] |
| HEAD_LAT->HEAD_OUT | [+0.09, +0.18] | ridge [+0.11, +0.18]; hgb [+0.03, +0.07] |
| TOKENS->SET | [-0.01, +0.04] | ridge [-0.01, +0.04]; hgb [-0.01, +0.09] |
| RAW->MRNA_PDL1 | [+0.41, +0.74] | ridge [+0.29, +0.84]; hgb [+0.41, +0.74] |
| RAW->MRNA_PDL2 | [+0.47, +0.87] | ridge [+0.32, +0.97]; hgb [+0.47, +0.87] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'RAW' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Readout PDL2

Response Delta = -0.287 [-0.338, -0.229]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | +0.08 (ridge) | [-0.23, +0.34] | **unresolved** | ridge +0.08 [-0.23, +0.34]; hgb +0.14 [-0.14, +0.38] |
| TOKENS | chain | +0.08 (ridge) | [-0.19, +0.31] | **unresolved** | ridge +0.08 [-0.19, +0.31]; hgb +0.14 [-0.13, +0.37] |
| Z | chain | +0.13 (ridge) | [-0.11, +0.36] | **unresolved** | ridge +0.13 [-0.11, +0.36]; hgb +0.23 [+0.04, +0.42] |
| HEAD_LAT | chain | +0.42 (ridge) | [+0.25, +0.64] | **reject** | ridge +0.42 [+0.25, +0.64]; hgb +0.54 [+0.39, +0.75] |
| HEAD_OUT | chain | +0.51 (hgb) | [+0.35, +0.74] | **reject** | ridge +0.53 [+0.38, +0.74]; hgb +0.51 [+0.35, +0.74] |
| SET | probe | +0.09 (ridge) | [-0.22, +0.37] | **unresolved** | ridge +0.09 [-0.22, +0.37]; hgb +0.19 [-0.15, +0.48] |
| MRNA_PDL1 | probe | +0.76 (hgb) | [+0.72, +0.82] | **reject** | ridge +0.80 [+0.77, +0.86]; hgb +0.76 [+0.72, +0.82] |
| MRNA_PDL2 | probe | +0.80 (hgb) | [+0.77, +0.88] | **reject** | ridge +0.83 [+0.78, +0.91]; hgb +0.80 [+0.77, +0.88] |
| ADT_PDL2 | positive_control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb +0.00 [-0.00, +0.01] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [-0.05, +0.08] | ridge [-0.05, +0.08]; hgb [-0.02, +0.02] |
| TOKENS->Z | [+0.01, +0.14] | ridge [+0.01, +0.14]; hgb [+0.03, +0.19] |
| Z->HEAD_LAT | [+0.24, +0.38] | ridge [+0.24, +0.38]; hgb [+0.28, +0.38] |
| HEAD_LAT->HEAD_OUT | [+0.06, +0.13] | ridge [+0.08, +0.14]; hgb [-0.05, +0.01] |
| TOKENS->SET | [-0.07, +0.07] | ridge [-0.07, +0.07]; hgb [-0.03, +0.11] |
| RAW->MRNA_PDL1 | [+0.45, +0.98] | ridge [+0.47, +1.04]; hgb [+0.42, +0.89] |
| RAW->MRNA_PDL2 | [+0.44, +1.05] | ridge [+0.45, +1.09]; hgb [+0.40, +0.96] |

**Root cause: relative_loss_at.** Additional loss: the reference does not fully explain the readout, but one step of the chain adds a resolved loss on top (a weaker, relative claim). reference 'RAW' is unresolved (the readout is not fully explained even by the input view), but 'Z' -> 'HEAD_LAT' adds a resolved loss (change in u > 0.1)

## Provenance
- adapter: `perturb_cite`
- sha256_rna_h5ad: `03e9602d124c2612`
- sha256_adt_h5ad: `e1e293c75dfebe09`
- sha256_vocab: `e184b550fa2c3971`
- sha256_medians: `3d3637fa24ef81bf`
- sha256_head_ckpt: `b0c543e285123205`
- sha256_teddy_states: `d2b2974da8c700fc`
- python: `3.11.15`
- runtime: 3406.2 s
