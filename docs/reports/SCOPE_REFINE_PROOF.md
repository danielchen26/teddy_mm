# SCOPE_REFINE_PROOF — ANM × TEDDY CITE

## 中文摘要（给 Daniel）

在同一批 TEDDY δu（site4/test）+ holdout GT 上证明两件事：**(1) 用 ANM 连续 workability `soft_P=max(action_scores)`（二元 verifier `P_f` 的连续底层）做弃权门控，收紧 τ 后条件正确率 Q 上升**；**(2) ANM LOO / flip-sensitive 归因划定的 hard 细胞与优先蛋白子集，会改变测到的正确率，并在小标签预算下对 hard holdout 更省标签。**不声称临床；不声称 phase-1 Pearson>0.61 全局，除非数字本身显示。

### Claim 1 — soft_P / abstain 门控（coverage–Q）

- 主切片：`site4_phase1_export` n=16750；二元 `P_f`∈{0,1}（有推荐动作=1）；连续门控用同场 `soft_P`。

#### O0

- TEDDY always Q=0.9453 (n_lab=15711); TEDDY thr Q=0.9471936657940107; ANM default Q=0.9496203834770299 (abstain=317); mean binary P_f=0.9811.
- soft_P 收紧：baseline Q=0.9452612819043982 @cov=1.0; peak Q=1.0 @cov=0.0666865671641791 (ΔQ_peak=0.054738718095601824); best@cov≥0.4 Q=0.9983560005978179 (Δ=0.05309471869341975); tightest Q=1.0 @cov=0.012059701492537314; 非降比例=0.875.

| τ (soft_P) | coverage | Q (TEDDY∨gate) | n_lab_ans |
|---:|---:|---:|---:|
| 0.0000 | 1.0000 | 0.9453 | 15711 |
| 0.2978 | 0.9333 | 0.9614 | 14980 |
| 0.4352 | 0.8000 | 0.9835 | 13114 |
| 0.5838 | 0.6667 | 0.9961 | 11136 |
| 0.6801 | 0.5333 | 0.9980 | 8915 |
| 0.7832 | 0.4000 | 0.9984 | 6691 |
| 0.8109 | 0.2667 | 0.9982 | 4460 |
| 0.8952 | 0.1334 | 0.9996 | 2232 |
| 0.9778 | 0.0121 | 1.0000 | 202 |

#### O2

- TEDDY always Q=0.9149 (n_lab=15180); TEDDY thr Q=0.9428168019767031; ANM default Q=0.960305225836288 (abstain=3373); mean binary P_f=0.7986.
- soft_P 收紧：baseline Q=0.9148880105401844 @cov=1.0; peak Q=0.9703325624052955 @cov=0.7653731343283582 (ΔQ_peak=0.0554445518651111); best@cov≥0.4 Q=0.9703325624052955 (Δ=0.0554445518651111); tightest Q=0.8808167141500475 @cov=0.13886567164179103; 非降比例=0.7142857142857143.

| τ (soft_P) | coverage | Q (TEDDY∨gate) | n_lab_ans |
|---:|---:|---:|---:|
| 0.0000 | 1.0000 | 0.9149 | 15180 |
| 0.0730 | 1.0000 | 0.9149 | 15180 |
| 0.1462 | 0.9333 | 0.9314 | 14601 |
| 0.1797 | 0.8666 | 0.9466 | 13831 |
| 0.1995 | 0.8000 | 0.9624 | 12993 |
| 0.2435 | 0.7654 | 0.9703 | 12539 |
| 0.3062 | 0.5454 | 0.9628 | 8881 |
| 0.3602 | 0.1389 | 0.8808 | 2106 |

#### Missing-modality `adt_only` (n=16750)

> Reading note: with `adt_only` (RNA missing) TEDDY does not run. "TEDDY always" here is TEDDY alone's fixed rule applied to the phase-2 protein-only stand-in (see MISSING_MODALITY_ANM_DEMO.md), so this curve gates the stand-in's calls, not TEDDY's.

- **O0**: TEDDY always Q=0.6288; ΔQ_peak=0.32061029112476147 (Q 0.6288279111224295→0.949438202247191 @cov 0.010626865671641792); best@cov≥0.4 Δ=0.14510126778587418; nondec frac=1.0.

| τ (soft_P) | coverage | Q | n_lab_ans |
|---:|---:|---:|---:|
| 0.0000 | 1.0000 | 0.6288 | 15707 |
| 0.3182 | 0.9333 | 0.6423 | 14716 |
| 0.3653 | 0.8000 | 0.6703 | 12669 |
| 0.4064 | 0.6667 | 0.7046 | 10618 |
| 0.4540 | 0.5333 | 0.7371 | 8556 |
| 0.5101 | 0.4000 | 0.7739 | 6467 |
| 0.5868 | 0.2667 | 0.8117 | 4345 |
| 0.7165 | 0.1334 | 0.8578 | 2187 |
| 0.9778 | 0.0106 | 0.9494 | 178 |
- **O2**: TEDDY always Q=0.6470; ΔQ_peak=0.013795984430029051 (Q 0.6469542311491604→0.6607502155791894 @cov 0.9000597014925373); best@cov≥0.4 Δ=0.013795984430029051; nondec frac=0.7142857142857143.

| τ (soft_P) | coverage | Q | n_lab_ans |
|---:|---:|---:|---:|
| 0.0000 | 1.0000 | 0.6470 | 15185 |
| 0.0901 | 1.0000 | 0.6470 | 15185 |
| 0.2169 | 0.9333 | 0.6565 | 14357 |
| 0.2435 | 0.9001 | 0.6608 | 13916 |
| 0.2525 | 0.8000 | 0.6492 | 12496 |
| 0.2909 | 0.7333 | 0.6541 | 11525 |
| 0.3062 | 0.7054 | 0.6588 | 11120 |
| 0.3602 | 0.0167 | 0.6181 | 254 |

### Claim 2 — 归因划定 scope 改变测到的正确率

- LOO attr n=16433; flip-sensitive cells n=4557 (top1_flip_rate=0.2773).

| cell scope | Q | n_decided |
|---|---:|---:|
| full panel cells | 0.9453 | 15711 |
| flip-sensitive (ANM LOO flipped) | 0.8099 | 3724 |
| nonflip | 0.9948 | 11818 |
| random same-n as flip (mean±std, 40 draws) | 0.9458±0.0027 | 4557 |

- ΔQ (flip − random same-n) = **-0.1359** → attributed-hard scope 明显更难。

| protein scope | subset | Q |
|---|---|---:|
| full panel | ['CD19', 'CD72', 'CD22', 'CD3', 'CD2', 'CD5', 'CD16', 'CD11c', 'CD36'] | 0.9453 |
| flip-sensitive top-3 | ['CD16', 'CD36', 'CD11c'] | 0.2789 |
| overall LOO top-k | ['CD5', 'CD2', 'CD36'] | 0.6280 |
| balanced 1/lineage by flip | ['CD19', 'CD2', 'CD16'] | 0.8470 |
| random k proteins (mean±std) | k=3 | 0.6692±0.1334 |
| random 1/lineage (mean±std) | — | 0.8547±0.0195 |

- Flip-top-k 蛋白子集 Q=0.2789 vs random-k mean 0.6692（归因优先蛋白改变 panel 正确率；注：flip-top 常偏 myeloid，balanced 对照 Δ=-0.0077）。

### Claim 2b — 标签预算（frozen features + 小 head）

- Train labeled=3925 (hard=945); eval labeled=11786 (hard=2779). TEDDY Q eval=0.9461, hard-eval=0.8129.
- Frozen TEDDY panel features; small logreg/MLP head only. Eval target = ANM flip-sensitive hard holdout. Positive delta ⇒ ANM-flagged hard labels buy more Q-on-hard per label.

#### logreg → Q on hard holdout

| n_labels | hard→Q_hard | random→Q_hard | Δ |
|---:|---:|---:|---:|
| 50 | 0.7769±0.0156 | 0.7426±0.0731 | +0.0342 |
| 100 | 0.7893±0.0169 | 0.7927±0.0106 | -0.0034 |
| 200 | 0.8048±0.0061 | 0.8002±0.0136 | +0.0046 |
| 400 | 0.8103±0.0027 | 0.8057±0.0050 | +0.0046 |

#### mlp → Q on hard holdout

| n_labels | hard→Q_hard | random→Q_hard | Δ |
|---:|---:|---:|---:|
| 50 | 0.7719±0.0128 | 0.7287±0.0839 | +0.0432 |
| 100 | 0.7878±0.0143 | 0.7956±0.0080 | -0.0078 |
| 200 | 0.8071±0.0079 | 0.7970±0.0093 | +0.0101 |
| 400 | 0.8143±0.0055 | 0.8025±0.0061 | +0.0118 |

### 一句话证明

> On the same TEDDY δu (site4/test n=16750), ANM soft_P gating raises conditional Q O0 0.9453→1.0000 (ΔQ_peak=0.0547; @cov≥0.4 best=0.9984) and O2 0.9149→0.9703 (ΔQ_peak=0.0554; mid operating region — extreme-low coverage can dip). ANM LOO flip-sensitive cells Q=0.8099 vs random same-n 0.9458 (Δ=-0.1359). Separately, on missing-modality adt_only O0 (RNA missing, so TEDDY does not run; these are TEDDY alone's fixed rule's calls on the phase-2 protein-only stand-in), soft_P gate lifts Q 0.6288→0.9494 (ΔQ_peak=0.3206); not a gain on TEDDY's calls. Label-budget: frozen-feature head on ANM flip-hard cells beats random labels on hard holdout at n=50 (ΔQ=+0.0342). Conditional correctness / label efficiency only — not clinical; no best.pt retrain.

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

Outputs: `outputs/anm_cite_bridge/scope_refine/SCOPE_REFINE_PROOF.md`, `scope_refine_results.json`.

Generated: 2026-09-27 01:06:40 EDT
