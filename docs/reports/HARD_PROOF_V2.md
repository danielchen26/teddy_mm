# HARD_PROOF — ANM × TEDDY CITE（最强统计档）

> **Corrected rerun (v2), copied from `outputs/anm_cite_bridge_v2/hard_proof/`.** Branch `fix/anm-bridge-corrections`, commit `bf76d00`: all markers enter at once, train-median size factor 0.92663, one scoring function shared by the answer key, the rule and ANM, disjoint label-cost split, within-lineage permutation null. In the three-arm table below, the TEDDY-alone row lists accuracy of calls made (Q) and the TEDDY+ANM row lists exactness (Q_f, a decline counts as a miss). ANM makes the same call as the TEDDY-alone rule on every cell, so its accuracy of calls made is identical (0.9530 / 0.9776 / 0.9635). The first-run report is [HARD_PROOF.md](HARD_PROOF.md).


## 中文摘要（给 Daniel）

在 **同一批 δu**（site4/test **n=16750**；OOD **n=2000**，val_non_site4）上做三故事对照：**（1）TEDDY multimodal alone**、**（2）TEDDY+ANM 声明场**、**（3）Train+ANM（用标签硬训逼近）**。准确率是次要的；根本变化 = 可编辑观察者 + 弃权/P_f + 可审计归因（GT 仅 verifier）。归因全量 **n_attr=16302**；bootstrap B=200；permutation n_perm=100。

### 准则分歧 / O2 标签改写

| 对比 | both_defined | disagree | disagree_rate |
|---|---:|---:|---:|
| O0_vs_O1 | 14433 | 16 | 0.0011 |
| O0_vs_O2 | 14379 | 1794 | 0.1248 |
| O1_vs_O2 | 13732 | 1082 | 0.0788 |

- O0 vs O1：加严弃权（标签改写率 **0.0011**）。
- O0 vs O2：key-marker 优先 **改写** expected_action（分歧率 **0.1248**）。

### 三臂 Q vs GT（site4，full n）

| 臂 | O0 Q | O0 abstain | O1 Q | O1 abstain | O2 Q | O2 abstain | 准则变更成本 |
|---|---:|---:|---:|---:|---:|---:|---|
| （1）TEDDY alone | 0.9530 | 448 | 0.9776 | 2014 | 0.9635 | 3169 | 改声明（无标签），但无场/归因 |
| （2）TEDDY+ANM | 0.9379 | 448 | 0.9273 | 2014 | 0.8345 | 3169 | 仅改声明 |

### 全量归因 / Bootstrap / Permutation（ANM，O0 场）

- 归因细胞数：**16302** / requested=16750
- top-1 蛋白分布：`{'CD16': 2195, 'CD3': 2194, 'CD72': 1232, 'CD36': 2336, 'CD19': 1763, 'CD5': 3178, 'CD2': 1562, 'CD11c': 952, 'CD22': 890}`
- top-1 flip rate：**0.27456753772543246**
- flip distance quantiles：`{'q25': 0.4148986853337806, 'q50': 0.6858404569028007, 'q75': 0.8, 'n_with_flip': 4474, 'n_no_flip': 11828}`
- bootstrap top-1：mode=CD5 frac=0.1949 boot_mean=0.19492853637590482 ci95=[0.1886210280947123, 0.2007775119617225] (B=200)
- permutation：null_max_mean=0.14412648754754018 null_ci95=[0.1421589375536744, 0.14703717335296282] p=0.009900990099009901 (n_perm=100)

### 静默过度回答 / 弃权是否跟随声明

- TEDDY-alone O0 规则 vs O1 声明弃权：`{'count': 1566, 'rate': 0.09349253731343284}`
- TEDDY-alone O0 规则 vs O2 声明弃权：`{'count': 2726, 'rate': 0.16274626865671643}`
- ANM abstain O0/O1/O2：[448, 2014, 3169]
- Train 臂 O0 模型零样本弃权不跟随：`{'O0': 40, 'O1': 40, 'O2': 40}`

### 注释分级（cell_type；范围内 B/T/myeloid 准确率，范围外调用占比）

| 臂 | crit | 范围内 n | 范围内准确率(strict) | 范围外 n | 范围外调用 | 范围外调用率 |
|---|---|---:|---:|---:|---:|---:|
| TEDDY alone | O0 | 11632 | 0.9696 | 5118 | 4678 | 0.9140 |
| TEDDY alone | O1 | 11632 | 0.9555 | 5118 | 3257 | 0.6364 |
| TEDDY alone | O2 | 11632 | 0.8482 | 5118 | 3384 | 0.6612 |
| TEDDY+ANM | O0 | 11632 | 0.9696 | 5118 | 4678 | 0.9140 |
| TEDDY+ANM | O1 | 11632 | 0.9555 | 5118 | 3257 | 0.6364 |
| TEDDY+ANM | O2 | 11632 | 0.8482 | 5118 | 3384 | 0.6612 |


### （3）Train+ANM 标签成本曲线（逼近声明 O2）

划分：mode=disjoint，train_pool=11725，eval=5025，overlap=0，阈值选择集=train_pool。

目标（ANM 声明 O2）：abstain_rate=0.1897，Q=0.9628。

| 方法 | 匹配 headline |
|---|---|
| logistic_calibrated | `{'n_o2_labels': 500, 'accuracy_strict_labeled': 0.8188993641745231, 'Q_analogue': 0.9557318321392017, 'abstain_rate': 0.18865671641791046, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |
| mlp_calibrated | `{'n_o2_labels': 2000, 'accuracy_strict_labeled': 0.8241613681210261, 'Q_analogue': 0.9613810741687979, 'abstain_rate': 0.18766169154228857, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |
| grid_calibrated | `{'n_o2_labels': 50, 'accuracy_strict_labeled': 0.8403858802894102, 'Q_analogue': 0.9606516290726818, 'abstain_rate': 0.1799004975124378, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |

Logistic（Jev 站位）校准弃权后曲线（节选）：

| n_O2_labels | abstain | Q | Q_gap_vs_ANM |
|---:|---:|---:|---:|
| 50 | 0.1906 | 0.9380 | 0.024807028634564854 |
| 100 | 0.1936 | 0.9399 | 0.0229572431327294 |
| 200 | 0.1813 | 0.9438 | 0.019062828069798088 |
| 500 | 0.1887 | 0.9557 | 0.007090626130904609 |
| 1000 | 0.1986 | 0.9626 | 0.00024181310881599138 |
| 2000 | 0.1930 | 0.9674 | 0.004549450983016401 |
| 5000 | 0.1914 | 0.9705 | 0.007645973501787862 |
| 10619 | 0.1932 | 0.9722 | 0.00941391462851826 |

**标签成本 headline：** logistic_calibrated matched at n_O2_labels=500；mlp_calibrated matched at n_O2_labels=2000；grid_calibrated matched at n_O2_labels=50。关键：即便用满池 O2 标签把 Q/弃权率数值校准到接近 ANM，Train 臂仍无声明场、无 LOO 归因/翻转距离，且 O0→O2 观察者切换必须重新消耗标签；ANM 零标签改声明即切换弃权与 expected_action。


### 一句话证明

> **On the same TEDDY δu (site4/test full-n), ANM declaration edits the observer (O0→O1 tighter abstain; O0→O2 rewrites expected_action) with zero endpoint-label fit, keeping verifiable P_f/Q_f and auditable full-n LOO attribution (n_attr=16302, bootstrap B=200, permutation p=0.009900990099009901). TEDDY-alone silently over-answers where the declared observer abstains; Train+ANM (logistic/MLP/grid) must consume O2 labels to approach ANM's O2 Q/abstain numbers and still cannot buy editable-field attribution or label-free observer edits. Accuracy is secondary; not clinical.**

### 声明边界

- Holdout `adt_true` = verifier only；未拟合 ANM 场律到标签。
- 未重训 `best.pt`。未声称临床优越。Jev/MLP 为训练站位，非 live Jev API。
- 无 Perturb / GFlowNet / 160M / ATAC。

### OOD 摘要

| crit | TEDDY Q | TEDDY abstain | ANM Q | ANM abstain | ANM mean_P_f |
|---|---:|---:|---:|---:|---:|
| O0 | 0.9816887080366226 | 2 | 0.9816887080366226 | 2 | 0.999 |
| O1 | 0.9907030796048809 | 79 | 0.9907030796048809 | 79 | 0.9605 |
| O2 | 0.9217040875071963 | 155 | 0.9217040875071963 | 155 | 0.9225 |

### Disk

Avail 41.0 GiB / total 994.6 GiB (96% used). No cache purge needed.

### Re-run

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100 --workers 9  # add --label-cost-split legacy --perm-null within_cell for the old numbers
```

Outputs: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_v2/hard_proof/HARD_PROOF.md`, `hard_proof_results.json`, `attr_compact.npz`.

Generated: 2026-09-29 13:27:36 EDT
