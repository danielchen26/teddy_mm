# MUST-SEPARATE autopsy (Mode B)
**Not Mode A.** Decomposes reverse_step1 near-z pairs by protein/lineage.
- pairs analyzed: **5000**
- bio_frac: **1.0000** | noise_like_frac: **0.0000** | bio_minus_noise: **1.0000**
- identical fine cell_type frac: **0.5810**

## Buckets
- `bio_key_marker`: 2467 (0.493)
- `bio_adt_lineage`: 2059 (0.412)
- `bio_coarse_cell_type`: 474 (0.095)

## Protein rank (mean |Δ| among must pairs)
| protein | lineage | key? | mean|Δ| | median|Δ| |
|---|---|---|---:|---:|
| CD5 | t_lineage | False | 0.7024 | 0.4455 |
| CD3 | t_lineage | True | 0.7006 | 0.3810 |
| CD16 | myeloid | True | 0.6943 | 0.3455 |
| CD2 | t_lineage | False | 0.6654 | 0.5225 |
| CD36 | myeloid | False | 0.4979 | 0.3238 |
| CD72 | b_lineage | False | 0.4406 | 0.2963 |
| CD19 | b_lineage | True | 0.3081 | 0.2070 |
| CD11c | myeloid | False | 0.2315 | 0.1038 |
| CD22 | b_lineage | False | 0.2103 | 0.1355 |

## Lineage confusion (ADT labels)
- myeloid|t_lineage: 1490
- b_lineage|t_lineage: 375
- abstain|t_lineage: 213
- b_lineage|myeloid: 194
- abstain|myeloid: 168
- abstain|b_lineage: 24

Do not claim Pearson > ~0.61. Not clinical.
