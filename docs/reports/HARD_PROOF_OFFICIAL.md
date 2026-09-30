# HARD_PROOF — ANM × TEDDY CITE（最强统计档）

## 中文摘要（给 Daniel）

在 **同一批 δu**（site4/test **n=16750**；OOD **n=2000**，val_non_site4）上做三故事对照：**（1）TEDDY multimodal alone**、**（2）TEDDY+ANM 声明场**、**（3）Train+ANM（用标签硬训逼近）**。准确率是次要的；根本变化 = 可编辑观察者 + 弃权/P_f + 可审计归因（GT 仅 verifier）。归因全量 **n_attr=16478**；bootstrap B=200；permutation n_perm=100。

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
| （1）TEDDY alone | 0.9429 | 272 | 0.9667 | 1493 | 0.9545 | 2551 | 改声明（无标签），但无场/归因 |
| （2）TEDDY+ANM | 0.9336 | 272 | 0.9312 | 1493 | 0.8559 | 2551 | 仅改声明 |

### 全量归因 / Bootstrap / Permutation（ANM，O0 场）

- 归因细胞数：**16478** / requested=16750
- top-1 蛋白分布：`{'CD16': 2221, 'CD2': 2945, 'CD36': 2500, 'CD19': 2749, 'CD5': 2849, 'CD3': 1011, 'CD11c': 1049, 'CD22': 285, 'CD72': 869}`
- top-1 flip rate：**0.29991503823279525**
- flip distance quantiles：`{'q25': 0.46572335471778076, 'q50': 0.7394980447674924, 'q75': 0.8, 'n_with_flip': 4931, 'n_no_flip': 11547}`
- bootstrap top-1：mode=CD2 frac=0.1787 boot_mean=0.17925476392766115 ci95=[0.17343579317878383, 0.1835781041388518] (B=200)
- permutation：null_max_mean=0.15725755552858356 null_ci95=[0.15481247724238378, 0.16012410486709552] p=0.009900990099009901 (n_perm=100)

### 静默过度回答 / 弃权是否跟随声明

- TEDDY-alone O0 规则 vs O1 声明弃权：`{'count': 1221, 'rate': 0.0728955223880597}`
- TEDDY-alone O0 规则 vs O2 声明弃权：`{'count': 2279, 'rate': 0.13605970149253732}`
- ANM abstain O0/O1/O2：[272, 1493, 2551]
- Train 臂 O0 模型零样本弃权不跟随：`{'O0': 43, 'O1': 43, 'O2': 43}`

### 注释分级（cell_type；范围内 B/T/myeloid 准确率，范围外调用占比）

| 臂 | crit | 范围内 n | 范围内准确率(strict) | 范围外 n | 范围外调用 | 范围外调用率 |
|---|---|---:|---:|---:|---:|---:|
| TEDDY alone | O0 | 11632 | 0.9564 | 5118 | 4853 | 0.9482 |
| TEDDY alone | O1 | 11632 | 0.9442 | 5118 | 3715 | 0.7259 |
| TEDDY alone | O2 | 11632 | 0.8319 | 5118 | 4022 | 0.7859 |
| TEDDY+ANM | O0 | 11632 | 0.9564 | 5118 | 4853 | 0.9482 |
| TEDDY+ANM | O1 | 11632 | 0.9442 | 5118 | 3715 | 0.7259 |
| TEDDY+ANM | O2 | 11632 | 0.8319 | 5118 | 4022 | 0.7859 |


### （3）Train+ANM 标签成本曲线（逼近声明 O2）

划分：mode=disjoint，train_pool=11725，eval=5025，overlap=0，阈值选择集=train_pool。

目标（ANM 声明 O2）：abstain_rate=0.1534，Q=0.9553。

| 方法 | 匹配 headline |
|---|---|
| logistic_calibrated | `{'n_o2_labels': 500, 'accuracy_strict_labeled': 0.8414821311115983, 'Q_analogue': 0.9476543209876543, 'abstain_rate': 0.14925373134328357, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |
| mlp_calibrated | `{'n_o2_labels': 1000, 'accuracy_strict_labeled': 0.838193378645034, 'Q_analogue': 0.9483999007690399, 'abstain_rate': 0.1572139303482587, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |
| grid_calibrated | `{'n_o2_labels': 50, 'accuracy_strict_labeled': 0.85792589344442, 'Q_analogue': 0.954390243902439, 'abstain_rate': 0.15044776119402986, 'matched': True, 'match_metric': 'abstain+accuracy_strict vs ANM O2 declared'}` |

Logistic（Jev 站位）校准弃权后曲线（节选）：

| n_O2_labels | abstain | Q | Q_gap_vs_ANM |
|---:|---:|---:|---:|
| 50 | 0.1475 | 0.9025 | 0.052796322871820456 |
| 100 | 0.1590 | 0.9326 | 0.02265179449497301 |
| 200 | 0.1433 | 0.9362 | 0.019034484213620195 |
| 500 | 0.1493 | 0.9477 | 0.007624271387712267 |
| 1000 | 0.1536 | 0.9532 | 0.0020491937116694503 |
| 2000 | 0.1538 | 0.9601 | 0.004839576014589064 |
| 5000 | 0.1542 | 0.9618 | 0.006487564506725163 |
| 10619 | 0.1508 | 0.9629 | 0.007584172012237933 |

**标签成本 headline：** logistic_calibrated matched at n_O2_labels=500；mlp_calibrated matched at n_O2_labels=1000；grid_calibrated matched at n_O2_labels=50。关键：即便用满池 O2 标签把 Q/弃权率数值校准到接近 ANM，Train 臂仍无声明场、无 LOO 归因/翻转距离，且 O0→O2 观察者切换必须重新消耗标签；ANM 零标签改声明即切换弃权与 expected_action。


### 一句话证明

> **On the same TEDDY δu (site4/test full-n), ANM declaration edits the observer (O0→O1 tighter abstain; O0→O2 rewrites expected_action) with zero endpoint-label fit, keeping verifiable P_f/Q_f and auditable full-n LOO attribution (n_attr=16478, bootstrap B=200, permutation p=0.009900990099009901). TEDDY-alone silently over-answers where the declared observer abstains; Train+ANM (logistic/MLP/grid) must consume O2 labels to approach ANM's O2 Q/abstain numbers and still cannot buy editable-field attribution or label-free observer edits. Accuracy is secondary; not clinical.**

### 声明边界

- Holdout `adt_true` = verifier only；未拟合 ANM 场律到标签。
- 未重训 `best.pt`。未声称临床优越。Jev/MLP 为训练站位，非 live Jev API。
- 无 Perturb / GFlowNet / 160M / ATAC。

### OOD 摘要

| crit | TEDDY Q | TEDDY abstain | ANM Q | ANM abstain | ANM mean_P_f |
|---|---:|---:|---:|---:|---:|
| O0 | 0.9796644636502287 | 1 | 0.9796644636502287 | 1 | 0.9995 |
| O1 | 0.9872979214780601 | 57 | 0.9872979214780601 | 57 | 0.9715 |
| O2 | 0.9158134243458476 | 129 | 0.9158134243458476 | 129 | 0.9355 |

### Disk

Avail 28.5 GiB / total 994.6 GiB (97% used). No cache purge needed.

### Re-run

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100 --workers 9  # add --label-cost-split legacy --perm-null within_cell for the old numbers
```

Outputs: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge_official/hard_proof/HARD_PROOF.md`, `hard_proof_results.json`, `attr_compact.npz`.

Generated: 2026-09-30 02:21:56 EDT
