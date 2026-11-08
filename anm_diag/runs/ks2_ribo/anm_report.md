# ANM diagnosis: ks2_teddy_ribo_dev

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M (frozen) · inputs ['rna_counts']
- source: crispri · treated `^RP[LS]\d` vs `control` · 480 treated cells in 12 clusters, 2196 cells in all
- chain: RAW → TOKENS → Z
- margin u = 0.1, alpha = 0.05; declaration sha256 `7edd53acc63abe86…`; anm_diag 0.1.0

**Checks**
- [R2] warning: readout 'RP' shares measurement 'rna_counts' with the candidates (declared_shared): counting noise is shared; count splitting would remove it
- [R2] warning: readout 'MRP' shares measurement 'rna_counts' with the candidates (declared_shared): counting noise is shared; count splitting would remove it
- [R2] warning: readout 'MT' shares measurement 'rna_counts' with the candidates (declared_shared): counting noise is shared; count splitting would remove it

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| RAW | 0.855 |
| TOKENS | 0.839 |
| Z | 0.817 |
| SET | 0.844 |

| step | change, 95% interval |
|---|---|
| RAW->TOKENS | [-0.027, -0.006] |
| TOKENS->Z | [-0.034, -0.010] |
| TOKENS->SET | [-0.008, +0.022] |

## Readout RP

Response Delta = -0.128 [-0.181, -0.075]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.03 (ridge) | [-0.13, +0.02] | **support** | ridge -0.03 [-0.13, +0.02]; hgb -0.01 [-0.06, +0.03] |
| TOKENS | chain | -0.11 (ridge) | [-0.26, -0.03] | **support** | ridge -0.11 [-0.26, -0.03]; hgb +0.02 [-0.05, +0.06] |
| Z | chain | -0.02 (ridge) | [-0.09, +0.04] | **support** | ridge -0.02 [-0.09, +0.04]; hgb +0.02 [-0.21, +0.14] |
| SET | probe | -0.05 (ridge) | [-0.43, +0.20] | **unresolved** | ridge -0.05 [-0.43, +0.20]; hgb +0.09 [-0.37, +0.32] |
| READ_RP | positive_control | -0.00 (ridge) | [-0.00, +0.00] | **support** | ridge -0.00 [-0.00, +0.00]; hgb -0.03 [-0.10, +0.02] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [-0.14, -0.04] | ridge [-0.14, -0.04]; hgb [-0.02, +0.07] |
| TOKENS->Z | [+0.02, +0.22] | ridge [+0.02, +0.22]; hgb [-0.17, +0.09] |
| TOKENS->SET | [-0.19, +0.25] | ridge [-0.19, +0.25]; hgb [-0.34, +0.29] |

**Root cause: no_deficit.** No deficit: every stage supported. every chain stage supported for this readout

**Exact same-state test**: cells = 480; tokens_identical = 0.252; median_move_over_delta = 0.086; frac_full_delta_hidden = 0.0; frac_resolved_beyond_noise = 0.0; note = minimal family: readout genes move between rank neighbours, vocabulary total fixed

## Readout MRP

Response Delta = +0.074 [+0.048, +0.096]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.08 (ridge) | [-0.21, +0.00] | **support** | ridge -0.08 [-0.21, +0.00]; hgb -0.02 [-0.21, +0.08] |
| TOKENS | chain | -0.03 (ridge) | [-0.16, +0.06] | **support** | ridge -0.03 [-0.16, +0.06]; hgb +0.03 [-0.11, +0.12] |
| Z | chain | -0.03 (ridge) | [-0.17, +0.05] | **support** | ridge -0.03 [-0.17, +0.05]; hgb +0.02 [-0.16, +0.13] |
| SET | probe | +0.00 (ridge) | [-0.25, +0.14] | **unresolved** | ridge +0.00 [-0.25, +0.14]; hgb +0.30 [+0.07, +0.42] |
| READ_MRP | positive_control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb +0.04 [-0.09, +0.13] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [+0.00, +0.10] | ridge [+0.00, +0.10]; hgb [-0.00, +0.14] |
| TOKENS->Z | [-0.09, +0.09] | ridge [-0.09, +0.09]; hgb [-0.11, +0.07] |
| TOKENS->SET | [-0.10, +0.11] | ridge [-0.10, +0.11]; hgb [+0.14, +0.35] |

**Root cause: no_deficit.** No deficit: every stage supported. every chain stage supported for this readout

**Exact same-state test**: cells = 480; tokens_identical = 0.752; median_move_over_delta = 0.0; frac_full_delta_hidden = 0.0; frac_resolved_beyond_noise = 0.0; note = minimal family: readout genes move between rank neighbours, vocabulary total fixed

## Readout MT

Response Delta = -0.343 [-0.391, -0.296]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.00 (hgb) | [-0.02, +0.01] | **support** | ridge +0.02 [-0.02, +0.06]; hgb -0.00 [-0.02, +0.01] |
| TOKENS | chain | +0.03 (hgb) | [-0.00, +0.05] | **support** | ridge +0.08 [+0.00, +0.14]; hgb +0.03 [-0.00, +0.05] |
| Z | chain | +0.07 (ridge) | [+0.03, +0.11] | **unresolved** | ridge +0.07 [+0.03, +0.11]; hgb +0.27 [+0.20, +0.32] |
| SET | probe | +0.16 (ridge) | [+0.04, +0.25] | **unresolved** | ridge +0.16 [+0.04, +0.25]; hgb +0.35 [+0.27, +0.43] |
| READ_MT | positive_control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb -0.01 [-0.03, +0.00] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [+0.01, +0.06] | ridge [+0.01, +0.10]; hgb [+0.01, +0.06] |
| TOKENS->Z | [-0.00, +0.08] | ridge [-0.06, +0.06]; hgb [+0.18, +0.30] |
| TOKENS->SET | [+0.03, +0.21] | ridge [+0.00, +0.15]; hgb [+0.25, +0.40] |

**Root cause: unresolved.** Unresolved: later stages neither supported nor rejected. stages ['Z'] unresolved after 'TOKENS'

**Exact same-state test**: cells = 480; tokens_identical = 0.24; median_move_over_delta = 0.039; frac_full_delta_hidden = 0.0; frac_resolved_beyond_noise = 0.0; note = minimal family: readout genes move between rank neighbours, vocabulary total fixed

## Provenance
- adapter: `ks2_teddy`
- measurement_sharing_split: `False`
- sha256_counts: `70903aa8f85073c7`
- sha256_prep: `86d7619b2d1cb436`
- sha256_cells: `26334a983f6b91c4`
- sha256_vocab: `e184b550fa2c3971`
- sha256_medians: `3d3637fa24ef81bf`
- sha256_teddy_encoder: `6768e30688f89cdc`
- torch: `2.15.0+cu132`
- python: `3.13.18`
- runtime: 4212.4 s
