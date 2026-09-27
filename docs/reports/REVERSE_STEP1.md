# Reverse step 1 — must-separate pairs (Mode B)

**Not Mode A.** This probe asks whether nearby points in the TEDDY evidence
export can carry different true lineage / protein labels
(**sufficiency-for-readout**). It does **not** measure how `z` responds to
gene or protein perturbs (**perturb-response** / reverse closed loop).

## z definition

- Full loading factor: `z_512` = mean-pool last-layer tokens @ context **1024**
  (not pretrain 2048, not disease token).
- This run used export sidecar `outputs/anm_cite_bridge/z_rna_export.npy`:
  site4 n=16750, dim=32 (compact `z_keep`, not full 512).

## Params (primary pair list)

- k=30, min_cos=0.98, adt_margin=0.05,
  protein_l1_rel=0.35

## Preliminary stats (site4/test n=16750, primary min_cos=0.98)

| Metric | Value |
|---|---:|
| Undirected near pairs (i&lt;j, cos≥min) | 262089 |
| Must-separate pairs (undirected) | 21267 |
| Pairs written | 5000 |
| Cells with ≥1 must neighbor | 6659 (0.3976) |
| Diff ADT lineage (among unique near) | 9023 |
| Diff key-marker lineage | 30046 |
| Diff coarse cell_type | 5460 |
| Diff protein profile (rel L1) | 117726 |

Must-pair cos quantiles: `{'q25': 0.9822211265563965, 'q50': 0.9844339191913605, 'q75': 0.9867914766073227, 'min': 0.9800009727478027, 'max': 0.9950082898139954}`

## Cosine threshold sweep (same kNN)

| min_cos | near pairs | must-separate | frac must | cells w/ must | diff ADT lin | diff key | diff coarse | diff protein |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.995 | 85 | 4 | 0.0471 | 8 | 0 | 12 | 0 | 32 |
| 0.99 | 26122 | 954 | 0.0365 | 871 | 309 | 1960 | 74 | 10977 |
| 0.98 | 262089 | 21267 | 0.0811 | 6659 | 9023 | 30046 | 5460 | 117726 |
| 0.97 | 337102 | 35177 | 0.1044 | 8794 | 16262 | 47382 | 7842 | 155961 |
| 0.95 | 375430 | 47281 | 0.1259 | 9830 | 22775 | 60060 | 9939 | 178812 |

Outputs: `must_separate_pairs.jsonl`, `reverse_step1_stats.json`.

Do not claim Pearson > ~0.61. Not clinical.
