# HARD_PROOF — ANM × TEDDY CITE（最强统计档）

## 中文摘要（给 Daniel）

在 **同一批 δu**（site4/test **n=16750**；OOD **n=2000**，val_non_site4）上做三故事对照：**（1）TEDDY multimodal alone**、**（2）TEDDY+ANM 声明场**、**（3）Train+ANM（用标签硬训逼近）**。准确率是次要的；根本变化 = 可编辑观察者 + 弃权/P_f + 可审计归因（GT 仅 verifier）。归因全量 **n_attr=16433**；bootstrap B=200；permutation n_perm=100。

### 准则分歧 / O2 标签改写

| 对比 | both_defined | disagree | disagree_rate |
|---|---:|---:|---:|
| O0_vs_O1 | 14488 | 0 | 0.0000 |
| O0_vs_O2 | 14379 | 1794 | 0.1248 |
| O1_vs_O2 | 13400 | 1414 | 0.1055 |

- O0 vs O1：加严弃权（标签改写率 **0.0000**）。
- O0 vs O2：key-marker 优先 **改写** expected_action（分歧率 **0.1248**）。

### 三臂 Q vs GT（site4，full n）

| 臂 | O0 Q | O0 abstain | O1 Q | O1 abstain | O2 Q | O2 abstain | 准则变更成本 |
|---|---:|---:|---:|---:|---:|---:|---|
| （1）TEDDY alone | 0.9472 | 101 | 0.9735 | 881 | 0.9428 | 1796 | 需标签重调 |
| （2）TEDDY+ANM | 0.9394 | 317 | 0.9578 | 834 | 0.8208 | 3373 | 仅改声明 |

### 全量归因 / Bootstrap / Permutation（ANM，O0 场）

- 归因细胞数：**16433** / requested=16750
- top-1 蛋白分布：`{'CD16': 2087, 'CD3': 668, 'CD72': 1455, 'CD36': 2394, 'CD22': 2120, 'CD5': 3882, 'CD2': 2516, 'CD11c': 1195, 'CD19': 116}`
- top-1 flip rate：**0.27730785614312664**
- flip distance quantiles：`{'q25': 0.5270483050519688, 'q50': 0.75, 'q75': 0.8, 'n_with_flip': 4532, 'n_no_flip': 11901}`
- bootstrap top-1：mode=CD5 frac=0.2362 boot_mean=0.23629951926002554 ci95=[0.23014209213168624, 0.24329550295137833] (B=200)
- permutation：null_max_mean=0.17131827709000183 null_ci95=[0.16605098172593907, 0.17704780643315893] p=0.009900990099009901 (n_perm=100)

### 静默过度回答 / 弃权是否跟随声明

- TEDDY-alone O0 规则 vs O1 声明弃权：`{'count': 780, 'rate': 0.046567164179104475}`
- TEDDY-alone O0 规则 vs O2 声明弃权：`{'count': 1752, 'rate': 0.10459701492537313}`
- ANM abstain O0/O1/O2：[317, 834, 3373]
- Train 臂 O0 模型零样本弃权不跟随：`{'O0': 185, 'O1': 185, 'O2': 185}`

### （3）Train+ANM 标签成本曲线（逼近声明 O2）

目标（ANM 声明 O2）：abstain_rate=0.2014，Q=0.9603。

| 方法 | 匹配 headline |
|---|---|
| logistic_calibrated | `{'n_o2_labels': 500, 'accuracy_strict_labeled': 0.8129117259552042, 'Q_analogue': 0.9633099141295862, 'abstain_rate': 0.20250746268656716, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |
| mlp_calibrated | `{'n_o2_labels': 200, 'accuracy_strict_labeled': 0.8019762845849803, 'Q_analogue': 0.9515397842738784, 'abstain_rate': 0.19761194029850745, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |
| grid_calibrated | `{'n_o2_labels': 50, 'accuracy_strict_labeled': 0.8227272727272728, 'Q_analogue': 0.9665660552588808, 'abstain_rate': 0.20626865671641792, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |

Logistic（Jev 站位）校准弃权后曲线（节选）：

| n_O2_labels | abstain | Q | Q_gap_vs_ANM |
|---:|---:|---:|---:|
| 50 | 0.1959 | 0.9383 | 0.022046523368243864 |
| 100 | 0.2047 | 0.9456 | 0.014700830231892392 |
| 200 | 0.1975 | 0.9480 | 0.012325625993212297 |
| 500 | 0.2025 | 0.9633 | 0.003004688293298252 |
| 1000 | 0.1996 | 0.9656 | 0.00527631405960205 |
| 2000 | 0.1997 | 0.9707 | 0.010426104768391897 |
| 5000 | 0.2030 | 0.9740 | 0.013688582213247646 |
| 10619 | 0.1992 | 0.9741 | 0.013745748214686104 |


### 标签成本数字（中文）

- 要在数值上逼近 ANM 声明 O2 的 **弃权率≈0.201 + strict≈0.821**：
  - thr/boost grid：**约 50** 个 O2 标签；
  - 小 MLP head：**约 200** 个；
  - Logistic（Jev 站位）：**约 500** 个。
- 但这只是数值校准：Train 臂仍 **没有** 可编辑场声明、**没有** LOO/翻转距离归因；O0→O2 观察者切换每次都要重新耗标签。ANM **零标签改 YAML/准则** 即切换弃权与 expected_action（O0→O2 标签改写率 12.48%）。
- O0 训好的模型弃权在 O0/O1/O2 上固定为 185——**不跟随声明**；ANM 弃权 317→834→3373。

**标签成本 headline：** logistic_calibrated matched at n_O2_labels=500 (strict=0.8129117259552042, abstain=0.20250746268656716)；mlp_calibrated matched at n_O2_labels=200 (strict=0.8019762845849803, abstain=0.19761194029850745)；grid_calibrated matched at n_O2_labels=50 (strict=0.8227272727272728, abstain=0.20626865671641792)。关键：即便用满池 O2 标签把弃权率/strict 数值校准到接近 ANM，Train 臂仍无声明场、无 LOO 归因/翻转距离，且 O0→O2 观察者切换必须重新消耗标签；ANM 零标签改声明即切换弃权与 expected_action。


### 一句话证明

> **On the same TEDDY δu (site4/test full-n), ANM declaration edits the observer (O0→O1 tighter abstain; O0→O2 rewrites expected_action) with zero endpoint-label fit, keeping verifiable P_f/Q_f and auditable full-n LOO attribution (n_attr=16433, bootstrap B=200, permutation p=0.009900990099009901). TEDDY-alone silently over-answers where the declared observer abstains; Train+ANM (logistic/MLP/grid) must consume O2 labels to approach ANM's O2 abstain/strict operating point and still cannot buy editable-field attribution or label-free observer edits. Accuracy is secondary; not clinical.**

### 声明边界

- Holdout `adt_true` = verifier only；未拟合 ANM 场律到标签。
- 未重训 `best.pt`。未声称临床优越。Jev/MLP 为训练站位，非 live Jev API。
- 无 Perturb / GFlowNet / 160M / ATAC。

### OOD 摘要

| crit | TEDDY Q | TEDDY abstain | ANM Q | ANM abstain | ANM mean_P_f |
|---|---:|---:|---:|---:|---:|
| O0 | 0.9816980172852059 | 1 | 0.9816980172852059 | 1 | 0.9995 |
| O1 | 0.9904661016949152 | 19 | 0.9899257688229056 | 21 | 0.9895 |
| O2 | 0.9170428893905191 | 111 | 0.9783037475345168 | 431 | 0.7845 |

### Disk

Avail 312.8 GiB / total 1995.2 GiB (84% used). No cache purge needed.

### Re-run

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100
```

Outputs: `outputs/anm_cite_bridge/hard_proof/HARD_PROOF.md`, `hard_proof_results.json`, `attr_compact.npz`.

Generated: 2026-09-26 15:55:04 EDT
