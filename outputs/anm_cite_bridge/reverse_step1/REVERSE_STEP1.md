# Reverse step 1 — must-separate pairs (Mode B)

**Not Mode A.** This probe asks whether nearby points in the TEDDY evidence
export can carry different true lineage / protein labels
(**sufficiency-for-readout**). It does **not** measure how `z` responds to
gene or protein perturbs (**perturb-response** / reverse closed loop).

## z definition

- `z_512` = mean-pool last-layer tokens @ context **1024**
  (not pretrain 2048, not disease token); L2-normalized full **512-D**.
- This run used `outputs/anm_cite_bridge/z_rna_512.npy`: site4 n=16750, dim=512
  (not compact `z_rna_export.npy` / `z_keep=32`).

## Params (primary pair list)

- k=30, min_cos=0.98, adt_margin=0.05,
  protein_l1_rel=0.35

## Preliminary stats (site4/test n=16750, primary min_cos=0.98)

| Metric | Value |
|---|---:|
| Undirected near pairs (i&lt;j, cos≥min) | 177704 |
| Must-separate pairs (undirected) | 9723 |
| Pairs written | 5000 |
| Cells with ≥1 must neighbor | 3666 (0.2189) |
| Diff ADT lineage (among unique near) | 4051 |
| Diff key-marker lineage | 13359 |
| Diff coarse cell_type | 2094 |
| Diff protein profile (rel L1) | 71180 |

Must-pair cos quantiles: `{'q25': 0.9812130928039551, 'q50': 0.982680469751358, 'q75': 0.984502375125885, 'min': 0.9800000190734863, 'max': 0.993342936038971}`

## Cosine threshold sweep (same kNN)

| min_cos | near pairs | must-separate | frac must | cells w/ must | diff ADT lin | diff key | diff coarse | diff protein |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.995 | 0 | 0 | nan | 0 | 0 | 0 | 0 | 0 |
| 0.99 | 558 | 94 | 0.1685 | 112 | 12 | 207 | 0 | 169 |
| 0.98 | 177704 | 9723 | 0.0547 | 3666 | 4051 | 13359 | 2094 | 71180 |
| 0.97 | 308800 | 25178 | 0.0815 | 7014 | 11925 | 33734 | 5924 | 133475 |
| 0.95 | 365944 | 38690 | 0.1057 | 8878 | 18818 | 52354 | 7407 | 161398 |

Outputs: `must_separate_pairs.jsonl`, `reverse_step1_stats.json`.

Do not claim Pearson > ~0.61. Not clinical.
