# ANM diagnosis: ks2_teddy_mitoribo_dev

- status: **exploratory** (exploratory: margins and design were not sealed before these data)
- system: TEDDY-G-70M (frozen) · inputs ['rna_counts']
- source: crispri · treated `^MRP[LS]\d` vs `control` · 240 treated cells in 6 clusters, 1956 cells in all
- chain: RAW → TOKENS → Z
- margin u = 0.1, alpha = 0.05; declaration sha256 `bfe0a7edc0f1b24b…`; anm_diag 0.1.0

**Checks**
- [R2] warning: readout 'MT' shares measurement 'rna_counts' with the candidates (declared_shared): counting noise is shared; count splitting would remove it
- [R2] warning: readout 'PSM' shares measurement 'rna_counts' with the candidates (declared_shared): counting noise is shared; count splitting would remove it
- [K1] warning: only 6 treated clusters: cluster-bootstrap intervals are unreliable; read directions only

## Forward: does the source reach each state? (held-out-cluster AUC)

| state | AUC |
|---|---|
| RAW | 0.789 |
| TOKENS | 0.717 |
| Z | 0.701 |
| SET | 0.670 |
| DEPTH | 0.370 |

| step | change, 95% interval |
|---|---|
| RAW->TOKENS | [-0.105, -0.037] |
| TOKENS->Z | [-0.105, +0.061] |
| TOKENS->SET | [-0.066, -0.029] |

## Readout MT

Response Delta = +0.124 [+0.052, +0.190]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | -0.01 (hgb) | [-0.08, +0.08] | **support** | ridge +0.39 [+0.11, +0.49]; hgb -0.01 [-0.08, +0.08] |
| TOKENS | chain | +0.01 (hgb) | [-0.14, +0.09] | **support** | ridge +0.63 [+0.31, +0.75]; hgb +0.01 [-0.14, +0.09] |
| Z | chain | +0.47 (ridge) | [+0.29, +0.64] | **reject** | ridge +0.47 [+0.29, +0.64]; hgb +0.78 [+0.64, +0.86] |
| SET | probe | +0.78 (hgb) | [+0.47, +0.90] | **reject** | ridge +0.80 [+0.47, +0.92]; hgb +0.78 [+0.47, +0.90] |
| DEPTH | technical | +0.87 (hgb) | [+0.79, +0.96] | **reject** | ridge +0.92 [+0.84, +0.97]; hgb +0.87 [+0.79, +0.96] |
| READ_MT | positive_control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb -0.00 [-0.07, +0.10] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [-0.21, +0.10] | ridge [+0.11, +0.33]; hgb [-0.21, +0.10] |
| TOKENS->Z | [+nan, +nan] | ridge [-0.33, +0.16]; hgb [+0.66, +0.89] |
| TOKENS->SET | [+0.60, +0.87] | ridge [+0.13, +0.21]; hgb [+0.60, +0.87] |

**Root cause: loss_at.** Representation loss: a stage is rejected after a supported stage. distinction merged between 'TOKENS' and 'Z'

**Exact same-state test**: cells = 240; tokens_identical = 0.171; median_move_over_delta = 0.147; frac_full_delta_hidden = 0.0; frac_resolved_beyond_noise = 0.0; note = minimal family: readout genes move between rank neighbours, vocabulary total fixed

## Readout PSM

Response Delta = -0.072 [-0.088, -0.057]

| state | role | u (best class) | 95% interval | verdict | all classes |
|---|---|---|---|---|---|
| RAW | chain | +0.19 (ridge) | [+0.08, +0.29] | **unresolved** | ridge +0.19 [+0.08, +0.29]; hgb +0.26 [+0.19, +0.34] |
| TOKENS | chain | +0.31 (hgb) | [+0.23, +0.38] | **reject** | ridge +0.38 [+0.31, +0.45]; hgb +0.31 [+0.23, +0.38] |
| Z | chain | +0.43 (hgb) | [+0.31, +0.54] | **reject** | ridge +0.44 [+0.32, +0.56]; hgb +0.43 [+0.31, +0.54] |
| SET | probe | +0.46 (ridge) | [+0.37, +0.55] | **reject** | ridge +0.46 [+0.37, +0.55]; hgb +0.55 [+0.45, +0.65] |
| DEPTH | technical | +0.84 (hgb) | [+0.71, +0.95] | **reject** | ridge +0.89 [+0.80, +0.97]; hgb +0.84 [+0.71, +0.95] |
| READ_PSM | positive_control | +0.00 (ridge) | [-0.00, +0.00] | **support** | ridge +0.00 [-0.00, +0.00]; hgb +0.14 [+0.06, +0.20] |

| step | change in u (best class) | same class |
|---|---|---|
| RAW->TOKENS | [+nan, +nan] | ridge [+0.09, +0.27]; hgb [-0.02, +0.12] |
| TOKENS->Z | [+0.02, +0.22] | ridge [-0.07, +0.20]; hgb [+0.02, +0.22] |
| TOKENS->SET | [+nan, +nan] | ridge [-0.02, +0.17]; hgb [+0.14, +0.35] |

**Root cause: reference_unresolved.** The reference stage is not supported: the readout class or the data cannot resolve this response; nothing is charged to the system. reference stage 'RAW' is unresolved: the readout class or the data cannot resolve this response; no claim about later stages

**Exact same-state test**: cells = 240; tokens_identical = 0.887; median_move_over_delta = 0.0; frac_full_delta_hidden = 0.0; frac_resolved_beyond_noise = 0.0; note = minimal family: readout genes move between rank neighbours, vocabulary total fixed

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
- runtime: 3335.4 s
