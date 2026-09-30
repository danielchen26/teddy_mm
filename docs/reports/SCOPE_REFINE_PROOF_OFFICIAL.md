# SCOPE_REFINE_PROOF — ANM × TEDDY CITE

## 中文摘要（给 Daniel）

在同一批 TEDDY δu（site4/test）+ holdout GT 上证明两件事：**(1) 用 ANM 连续 workability `soft_P=max(action_scores)`（二元 verifier `P_f` 的连续底层）做弃权门控，收紧 τ 后条件正确率 Q 上升**；**(2) ANM LOO / flip-sensitive 归因划定的 hard 细胞与优先蛋白子集，会改变测到的正确率，并在小标签预算下对 hard holdout 更省标签。**不声称临床；不声称 phase-1 Pearson>0.61 全局，除非数字本身显示。

### Claim 1 — soft_P / abstain 门控（coverage–Q）

- 主切片：`site4_phase1_export` n=16750；二元 `P_f`∈{0,1}（有推荐动作=1）；连续门控用同场 `soft_P`。

#### O0

- TEDDY always Q=0.9381 (n_lab=15711); TEDDY thr Q=0.9428553062929871; ANM default Q=0.9428553062929871 (abstain=272); mean binary P_f=0.9838.
- soft_P 收紧：baseline Q=0.9381325186175291 @cov=1.0; peak Q=1.0 @cov=0.3333731343283582 (ΔQ_peak=0.06186748138247089); best@cov≥0.4 Q=0.9992532855436081 (Δ=0.06112076692607904); tightest Q=1.0 @cov=0.06788059701492537; 非降比例=1.0.

| τ (soft_P) | coverage | Q (TEDDY∨gate) | n_lab_ans |
|---:|---:|---:|---:|
| 0.0000 | 1.0000 | 0.9381 | 15711 |
| 0.3040 | 0.9333 | 0.9564 | 14964 |
| 0.4369 | 0.8000 | 0.9816 | 13035 |
| 0.5106 | 0.6667 | 0.9946 | 11111 |
| 0.6263 | 0.5333 | 0.9983 | 8923 |
| 0.7709 | 0.4000 | 0.9993 | 6696 |
| 0.8790 | 0.2667 | 1.0000 | 4465 |
| 1.0141 | 0.1334 | 1.0000 | 2233 |
| 1.1101 | 0.0679 | 1.0000 | 1136 |

- TEDDY 自身 lineage margin 门控（top1−top2）：baseline Q=0.9381; peak Q=1.0000 @cov=0.3334 (ΔQ_peak=0.0619); best@cov≥0.4 Δ=0.0613.

| coverage | Q gate=ANM soft_P | Q gate=TEDDY margin | Δ (soft_P − margin) |
|---:|---:|---:|---:|
| 1.00 | 0.9381 | 0.9381 | +0.0000 |
| 0.90 | 0.9627 | 0.9696 | -0.0068 |
| 0.80 | 0.9816 | 0.9874 | -0.0058 |
| 0.70 | 0.9926 | 0.9968 | -0.0042 |
| 0.60 | 0.9965 | 0.9978 | -0.0013 |
| 0.50 | 0.9983 | 0.9984 | -0.0001 |
| 0.40 | 0.9993 | 0.9994 | -0.0001 |
| 0.30 | 1.0000 | 1.0000 | +0.0000 |
| 0.20 | 1.0000 | 1.0000 | +0.0000 |
| 0.10 | 1.0000 | 1.0000 | +0.0000 |

#### O2

- TEDDY always Q=0.9108 (n_lab=15180); TEDDY thr Q=0.9544553000808051; ANM default Q=0.9544553000808051 (abstain=2551); mean binary P_f=0.8477.
- soft_P 收紧：baseline Q=0.9108036890645587 @cov=1.0; peak Q=1.0 @cov=0.2 (ΔQ_peak=0.08919631093544134); best@cov≥0.4 Q=0.9988047213506649 (Δ=0.08800103228610623); tightest Q=1.0 @cov=0.1679402985074627; 非降比例=0.9230769230769231.

| τ (soft_P) | coverage | Q (TEDDY∨gate) | n_lab_ans |
|---:|---:|---:|---:|
| 0.0000 | 1.0000 | 0.9108 | 15180 |
| 0.0207 | 1.0000 | 0.9108 | 15180 |
| 0.0295 | 0.9333 | 0.9245 | 14608 |
| 0.0446 | 0.8666 | 0.9470 | 13848 |
| 0.0586 | 0.8000 | 0.9644 | 12994 |
| 0.0986 | 0.7333 | 0.9805 | 12055 |
| 0.1233 | 0.7064 | 0.9872 | 11634 |
| 0.1723 | 0.5333 | 0.9981 | 8913 |
| 0.1937 | 0.4667 | 0.9987 | 7805 |
| 0.2117 | 0.4000 | 0.9988 | 6693 |
| 0.2342 | 0.3334 | 0.9987 | 5578 |
| 0.2764 | 0.2667 | 0.9996 | 4466 |
| 0.3278 | 0.2000 | 1.0000 | 3349 |
| 0.3700 | 0.1679 | 1.0000 | 2812 |

- TEDDY 自身 lineage margin 门控（top1−top2）：baseline Q=0.9108; peak Q=1.0000 @cov=0.2000 (ΔQ_peak=0.0892); best@cov≥0.4 Δ=0.0880.

| coverage | Q gate=ANM soft_P | Q gate=TEDDY margin | Δ (soft_P − margin) |
|---:|---:|---:|---:|
| 1.00 | 0.9108 | 0.9108 | +0.0000 |
| 0.90 | 0.9343 | 0.9364 | -0.0021 |
| 0.80 | 0.9644 | 0.9645 | -0.0001 |
| 0.70 | 0.9880 | 0.9886 | -0.0006 |
| 0.60 | 0.9972 | 0.9973 | -0.0001 |
| 0.50 | 0.9987 | 0.9987 | -0.0000 |
| 0.40 | 0.9988 | 0.9988 | +0.0000 |
| 0.30 | 0.9988 | 0.9988 | +0.0000 |
| 0.20 | 1.0000 | 1.0000 | +0.0000 |
| 0.10 | 1.0000 | 1.0000 | +0.0000 |

### Claim 2 — 归因划定 scope 改变测到的正确率

- LOO attr n=16478; flip-sensitive cells n=4942 (top1_flip_rate=0.2999).

| cell scope | Q | n_decided |
|---|---:|---:|
| full panel cells | 0.9381 | 15711 |
| flip-sensitive (ANM LOO flipped) | 0.7913 | 4068 |
| nonflip | 0.9965 | 11489 |
| random same-n as flip (mean±std, 40 draws) | 0.9388±0.0029 | 4942 |

- ΔQ (flip − random same-n) = **-0.1475** → attributed-hard scope 明显更难。

#### ANM-flagged hard cells — cell_type composition (n_flagged=4942 / 16750)

| coarse cell_type | n_flagged | % of flagged | % of all | flag rate | enrichment |
|---|---:|---:|---:|---:|---:|
| NK / ILC | 1756 | 35.5% | 11.6% | 90.5% | 3.07 |
| erythroid | 1475 | 29.8% | 8.8% | 99.9% | 3.38 |
| T | 729 | 14.8% | 41.9% | 10.4% | 0.35 |
| progenitor | 620 | 12.5% | 9.6% | 38.5% | 1.31 |
| myeloid (mono/DC) | 268 | 5.4% | 10.5% | 15.2% | 0.52 |
| B / plasma | 94 | 1.9% | 17.6% | 3.2% | 0.11 |

| fine cell_type (top 12 by n_flagged) | n_flagged | % of flagged | flag rate | enrichment |
|---|---:|---:|---:|---:|
| NK | 1350 | 27.3% | 97.5% | 3.31 |
| Proerythroblast | 619 | 12.5% | 99.7% | 3.38 |
| gdT CD158b+ | 383 | 7.7% | 88.7% | 3.00 |
| Reticulocyte | 371 | 7.5% | 100.0% | 3.39 |
| Erythroblast | 339 | 6.9% | 100.0% | 3.39 |
| NK CD158e1+ | 303 | 6.1% | 99.0% | 3.36 |
| HSC | 210 | 4.2% | 54.4% | 1.84 |
| G/M prog | 209 | 4.2% | 57.4% | 1.95 |
| cDC2 | 175 | 3.5% | 87.9% | 2.98 |
| Normoblast | 146 | 3.0% | 100.0% | 3.39 |
| MK/E prog | 140 | 2.8% | 71.4% | 2.42 |
| CD8+ T CD57+ CD45RO+ | 128 | 2.6% | 12.3% | 0.42 |

| protein scope | subset | Q |
|---|---|---:|
| full panel | ['CD19', 'CD72', 'CD22', 'CD3', 'CD2', 'CD5', 'CD16', 'CD11c', 'CD36'] | 0.9381 |
| flip-sensitive top-3 | ['CD16', 'CD36', 'CD11c'] | 0.2789 |
| overall LOO top-k | ['CD2', 'CD5', 'CD19'] | 0.6810 |
| balanced 1/lineage by flip | ['CD72', 'CD2', 'CD16'] | 0.8485 |
| random k proteins (mean±std) | k=3 | 0.6510±0.1498 |
| random 1/lineage (mean±std) | — | 0.8505±0.0225 |

- Flip-top-k 蛋白子集 Q=0.2789 vs random-k mean 0.6510（归因优先蛋白改变 panel 正确率；注：flip-top 常偏 myeloid，balanced 对照 Δ=-0.0020）。

### Claim 2b — 标签预算（frozen features + 小 head）

- Train labeled=3925 (hard=1026); eval labeled=11786 (hard=3042). TEDDY Q eval=0.9387, hard-eval=0.7936.
- Frozen TEDDY panel features; small logreg/MLP head only. Eval target = ANM flip-sensitive hard holdout. Positive delta ⇒ ANM-flagged hard labels buy more Q-on-hard per label.

#### logreg → Q on hard holdout

| n_labels | hard→Q_hard | random→Q_hard | Δ |
|---:|---:|---:|---:|
| 50 | 0.7875±0.0168 | 0.7290±0.0721 | +0.0585 |
| 100 | 0.8086±0.0088 | 0.7905±0.0152 | +0.0181 |
| 200 | 0.8137±0.0054 | 0.7982±0.0157 | +0.0155 |
| 400 | 0.8189±0.0035 | 0.8072±0.0070 | +0.0116 |

#### mlp → Q on hard holdout

| n_labels | hard→Q_hard | random→Q_hard | Δ |
|---:|---:|---:|---:|
| 50 | 0.7804±0.0243 | 0.7236±0.0771 | +0.0568 |
| 100 | 0.8089±0.0116 | 0.7906±0.0117 | +0.0183 |
| 200 | 0.8119±0.0081 | 0.7990±0.0140 | +0.0129 |
| 400 | 0.8237±0.0034 | 0.8046±0.0066 | +0.0191 |

### 一句话证明

> On the same TEDDY δu (site4/test n=16750), ANM soft_P gating raises conditional Q O0 0.9381→1.0000 (ΔQ_peak=0.0619; @cov≥0.4 best=0.9993) and O2 0.9108→1.0000 (ΔQ_peak=0.0892; mid operating region — extreme-low coverage can dip). ANM LOO flip-sensitive cells Q=0.7913 vs random same-n 0.9388 (Δ=-0.1475). Label-budget: frozen-feature head on ANM flip-hard cells beats random labels on hard holdout at n=50 (ΔQ=+0.0585). Conditional correctness / label efficiency only — not clinical; no best.pt retrain.

### 声明边界

- Holdout `adt_true` = verifier only；未把 ANM 场律拟合到标签。
- 未重训 `best.pt`。条件正确率 / 标签效率，非临床。
- Verifier `P_f` 二元；τ 扫描用同场连续 `soft_P=max(action_scores)`。
- 不声称 phase-1 Pearson >0.61 全局（除非导出数字本身如此）。
- 无 Perturb / GFlowNet / 160M / ATAC。

### Re-run

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_scope_refine_proof.py
```

Inputs: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_official`. Outputs: `SCOPE_REFINE_PROOF.md`, `scope_refine_results.json` in the chosen `--out-dir` (default `outputs/anm_cite_bridge/scope_refine`). Flags: `--bridge-dir`, `--workers N`, `--skip-adt-only`, `--max-cells N`.

Generated: 2026-09-30 02:22:33 EDT
