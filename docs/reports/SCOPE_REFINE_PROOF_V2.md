# SCOPE_REFINE_PROOF — ANM × TEDDY CITE

> **Corrected rerun (v2), copied from `outputs/anm_cite_bridge_v2/scope_refine/`.** Run on the v2 export (all markers at once, train-median size factor 0.92663, shared scoring) with `--skip-adt-only`, so there is no RNA-missing view. The auto-written one-line proof below compares the soft_P gate with answering every cell; against TEDDY's own margin gate at matched coverage (the tables under each criterion) soft_P is no better. soft_P reaches 1.11, so it is a score, not a probability. The first-run report is [SCOPE_REFINE_PROOF.md](SCOPE_REFINE_PROOF.md).


## 中文摘要（给 Daniel）

在同一批 TEDDY δu（site4/test）+ holdout GT 上证明两件事：**(1) 用 ANM 连续 workability `soft_P=max(action_scores)`（二元 verifier `P_f` 的连续底层）做弃权门控，收紧 τ 后条件正确率 Q 上升**；**(2) ANM LOO / flip-sensitive 归因划定的 hard 细胞与优先蛋白子集，会改变测到的正确率，并在小标签预算下对 hard holdout 更省标签。**不声称临床；不声称 phase-1 Pearson>0.61 全局，除非数字本身显示。

### Claim 1 — soft_P / abstain 门控（coverage–Q）

- 主切片：`site4_phase1_export` n=16750；二元 `P_f`∈{0,1}（有推荐动作=1）；连续门控用同场 `soft_P`。

#### O0

- TEDDY always Q=0.9452 (n_lab=15711); TEDDY thr Q=0.9530461777260381; ANM default Q=0.9530461777260381 (abstain=448); mean binary P_f=0.9733.
- soft_P 收紧：baseline Q=0.945197632232194 @cov=1.0; peak Q=1.0 @cov=0.0666865671641791 (ΔQ_peak=0.054802367767806004); best@cov≥0.4 Q=0.9986555123991634 (Δ=0.053457880166969396); tightest Q=1.0 @cov=0.01635820895522388; 非降比例=0.9375.

| τ (soft_P) | coverage | Q (TEDDY∨gate) | n_lab_ans |
|---:|---:|---:|---:|
| 0.0000 | 1.0000 | 0.9452 | 15711 |
| 0.2684 | 0.9333 | 0.9616 | 14983 |
| 0.4224 | 0.8000 | 0.9881 | 13078 |
| 0.4877 | 0.6667 | 0.9963 | 11131 |
| 0.6563 | 0.5333 | 0.9984 | 8924 |
| 0.7731 | 0.4000 | 0.9987 | 6694 |
| 0.8519 | 0.2667 | 0.9996 | 4466 |
| 0.9413 | 0.1334 | 0.9996 | 2233 |
| 1.1101 | 0.0164 | 1.0000 | 273 |

- TEDDY 自身 lineage margin 门控（top1−top2）：baseline Q=0.9452; peak Q=1.0000 @cov=0.1334 (ΔQ_peak=0.0548); best@cov≥0.4 Δ=0.0536.

| coverage | Q gate=ANM soft_P | Q gate=TEDDY margin | Δ (soft_P − margin) |
|---:|---:|---:|---:|
| 1.00 | 0.9452 | 0.9452 | +0.0000 |
| 0.90 | 0.9662 | 0.9732 | -0.0070 |
| 0.80 | 0.9881 | 0.9855 | +0.0026 |
| 0.70 | 0.9955 | 0.9969 | -0.0014 |
| 0.60 | 0.9976 | 0.9981 | -0.0005 |
| 0.50 | 0.9986 | 0.9986 | +0.0000 |
| 0.40 | 0.9987 | 0.9988 | -0.0001 |
| 0.30 | 0.9996 | 0.9992 | +0.0004 |
| 0.20 | 0.9997 | 0.9997 | +0.0000 |
| 0.10 | 0.9994 | 1.0000 | -0.0006 |

#### O2

- TEDDY always Q=0.9148 (n_lab=15180); TEDDY thr Q=0.9634897695291702; ANM default Q=0.9634897695291702 (abstain=3169); mean binary P_f=0.8108.
- soft_P 收紧：baseline Q=0.9148221343873518 @cov=1.0; peak Q=0.999409681227863 @cov=0.10113432835820896 (ΔQ_peak=0.08458754684051117); best@cov≥0.4 Q=0.9980576721948304 (Δ=0.08323553780747861); tightest Q=0.999409681227863 @cov=0.10113432835820896; 非降比例=1.0.

| τ (soft_P) | coverage | Q (TEDDY∨gate) | n_lab_ans |
|---:|---:|---:|---:|
| 0.0000 | 1.0000 | 0.9148 | 15180 |
| 0.0314 | 0.9333 | 0.9302 | 14618 |
| 0.0511 | 0.8000 | 0.9658 | 13001 |
| 0.1233 | 0.7219 | 0.9858 | 11890 |
| 0.1904 | 0.5333 | 0.9976 | 8915 |
| 0.2214 | 0.4000 | 0.9981 | 6693 |
| 0.2718 | 0.2667 | 0.9991 | 4465 |
| 0.3484 | 0.1334 | 0.9991 | 2234 |
| 0.3700 | 0.1011 | 0.9994 | 1694 |

- TEDDY 自身 lineage margin 门控（top1−top2）：baseline Q=0.9148; peak Q=1.0000 @cov=0.0667 (ΔQ_peak=0.0852); best@cov≥0.4 Δ=0.0837.

| coverage | Q gate=ANM soft_P | Q gate=TEDDY margin | Δ (soft_P − margin) |
|---:|---:|---:|---:|
| 1.00 | 0.9148 | 0.9148 | +0.0000 |
| 0.90 | 0.9399 | 0.9419 | -0.0020 |
| 0.80 | 0.9658 | 0.9665 | -0.0007 |
| 0.70 | 0.9879 | 0.9891 | -0.0012 |
| 0.60 | 0.9957 | 0.9956 | +0.0001 |
| 0.50 | 0.9977 | 0.9977 | +0.0000 |
| 0.40 | 0.9981 | 0.9985 | -0.0004 |
| 0.30 | 0.9988 | 0.9990 | -0.0002 |
| 0.20 | 0.9991 | 0.9991 | +0.0000 |
| 0.10 | 0.9994 | 1.0000 | -0.0006 |

### Claim 2 — 归因划定 scope 改变测到的正确率

- LOO attr n=16302; flip-sensitive cells n=4476 (top1_flip_rate=0.2746).

| cell scope | Q | n_decided |
|---|---:|---:|
| full panel cells | 0.9452 | 15711 |
| flip-sensitive (ANM LOO flipped) | 0.8177 | 3686 |
| nonflip | 0.9954 | 11776 |
| random same-n as flip (mean±std, 40 draws) | 0.9454±0.0034 | 4476 |

- ΔQ (flip − random same-n) = **-0.1277** → attributed-hard scope 明显更难。

#### ANM-flagged hard cells — cell_type composition (n_flagged=4476 / 16750)

| coarse cell_type | n_flagged | % of flagged | % of all | flag rate | enrichment |
|---|---:|---:|---:|---:|---:|
| NK / ILC | 1780 | 39.8% | 11.6% | 91.8% | 3.43 |
| erythroid | 1475 | 33.0% | 8.8% | 99.9% | 3.74 |
| progenitor | 432 | 9.7% | 9.6% | 26.8% | 1.00 |
| T | 421 | 9.4% | 41.9% | 6.0% | 0.22 |
| myeloid (mono/DC) | 275 | 6.1% | 10.5% | 15.6% | 0.58 |
| B / plasma | 93 | 2.1% | 17.6% | 3.2% | 0.12 |

| fine cell_type (top 12 by n_flagged) | n_flagged | % of flagged | flag rate | enrichment |
|---|---:|---:|---:|---:|
| NK | 1373 | 30.7% | 99.2% | 3.71 |
| Proerythroblast | 619 | 13.8% | 99.7% | 3.73 |
| Reticulocyte | 371 | 8.3% | 100.0% | 3.74 |
| Erythroblast | 339 | 7.6% | 100.0% | 3.74 |
| gdT CD158b+ | 331 | 7.4% | 76.6% | 2.87 |
| NK CD158e1+ | 306 | 6.8% | 100.0% | 3.74 |
| cDC2 | 177 | 4.0% | 88.9% | 3.33 |
| G/M prog | 170 | 3.8% | 46.7% | 1.75 |
| Normoblast | 146 | 3.3% | 100.0% | 3.74 |
| HSC | 134 | 3.0% | 34.7% | 1.30 |
| MK/E prog | 108 | 2.4% | 55.1% | 2.06 |
| ILC | 98 | 2.2% | 95.1% | 3.56 |

| protein scope | subset | Q |
|---|---|---:|
| full panel | ['CD19', 'CD72', 'CD22', 'CD3', 'CD2', 'CD5', 'CD16', 'CD11c', 'CD36'] | 0.9452 |
| flip-sensitive top-3 | ['CD16', 'CD36', 'CD11c'] | 0.2789 |
| overall LOO top-k | ['CD5', 'CD36', 'CD16'] | 0.6969 |
| balanced 1/lineage by flip | ['CD72', 'CD2', 'CD16'] | 0.8961 |
| random k proteins (mean±std) | k=3 | 0.6955±0.1393 |
| random 1/lineage (mean±std) | — | 0.8569±0.0226 |

- Flip-top-k 蛋白子集 Q=0.2789 vs random-k mean 0.6955（归因优先蛋白改变 panel 正确率；注：flip-top 常偏 myeloid，balanced 对照 Δ=0.0392）。

### Claim 2b — 标签预算（frozen features + 小 head）

- Train labeled=3925 (hard=940); eval labeled=11786 (hard=2746). TEDDY Q eval=0.9460, hard-eval=0.8201.
- Frozen TEDDY panel features; small logreg/MLP head only. Eval target = ANM flip-sensitive hard holdout. Positive delta ⇒ ANM-flagged hard labels buy more Q-on-hard per label.

#### logreg → Q on hard holdout

| n_labels | hard→Q_hard | random→Q_hard | Δ |
|---:|---:|---:|---:|
| 50 | 0.7968±0.0142 | 0.7419±0.0813 | +0.0549 |
| 100 | 0.8063±0.0079 | 0.8038±0.0116 | +0.0025 |
| 200 | 0.8173±0.0087 | 0.8100±0.0095 | +0.0073 |
| 400 | 0.8228±0.0040 | 0.8186±0.0057 | +0.0042 |

#### mlp → Q on hard holdout

| n_labels | hard→Q_hard | random→Q_hard | Δ |
|---:|---:|---:|---:|
| 50 | 0.7954±0.0172 | 0.7413±0.0770 | +0.0541 |
| 100 | 0.8028±0.0120 | 0.8078±0.0088 | -0.0050 |
| 200 | 0.8191±0.0068 | 0.8084±0.0080 | +0.0107 |
| 400 | 0.8252±0.0045 | 0.8123±0.0075 | +0.0129 |

### 一句话证明

> On the same TEDDY δu (site4/test n=16750), ANM soft_P gating raises conditional Q O0 0.9452→1.0000 (ΔQ_peak=0.0548; @cov≥0.4 best=0.9987) and O2 0.9148→0.9994 (ΔQ_peak=0.0846; mid operating region — extreme-low coverage can dip). ANM LOO flip-sensitive cells Q=0.8177 vs random same-n 0.9454 (Δ=-0.1277). Label-budget: frozen-feature head on ANM flip-hard cells beats random labels on hard holdout at n=50 (ΔQ=+0.0549). Conditional correctness / label efficiency only — not clinical; no best.pt retrain.

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

Inputs: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_v2`. Outputs: `SCOPE_REFINE_PROOF.md`, `scope_refine_results.json` in the chosen `--out-dir` (default `outputs/anm_cite_bridge/scope_refine`). Flags: `--bridge-dir`, `--workers N`, `--skip-adt-only`, `--max-cells N`.

Generated: 2026-09-29 13:30:37 EDT
