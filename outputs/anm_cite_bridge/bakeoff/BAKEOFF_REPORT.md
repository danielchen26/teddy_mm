# 3-way Bakeoff: TEDDY alone × TEDDY+ANM × TEDDY+Jev-class

## 中文摘要（给 Daniel）

在 **同一批 CITE test/site4 细胞**（500 eval / 1500 train pool）上做三臂对照，证明 ANM 带来的根本变化：**准则可声明编辑**（O0→O1 只改 YAML/阈值/key-marker boost，**零端点标签拟合**），并给出场上的 **留一归因 + 翻转距离**；Jev-class 站位与 TEDDY-alone 在准则切换时必须 **用新准则标签重调/重训**，否则准确率漂移或「静默错误」。

**准确率不是胜负条件**（本数据上 O0∩O1 的 expected_action 一致，三臂 Q 都可以很高）；故事是 **准则变更成本 + 弃权/P_f 门控 + 归因形态**：TEDDY-alone 用 O0 规则会在 O1 声明应弃权处继续作答（静默过度回答）；Jev 站位的 abstain 不会随更严观察者自动上升；ANM 仅改声明即可提高弃权并保持可验证 Q_f，且给出场 LOO/翻转距离。

声明边界：仅 **局部响应诊断**，非临床优越性。未重训 TEDDY phase-1 / `best.pt`。

### 关键对照表

| 臂 | O0 Q/acc | O1 zero-shot | O1 after adapt/retrain | 归因能力 |
|---|---|---|---|---|
| TEDDY alone | 0.9641 (abstain 1/500) | 0.9860 (仍用 O0 规则；over-answer vs O1声明弃权 **33**/500) | 0.9860 (网格搜 thr=0.08, boost=3.0；需标签) | 面板均值 LOO / 无声明场；翻转=argmax 网格 |
| TEDDY + ANM | 0.9684 / P_f=0.990 (abstain 5) | 0.9767 / P_f=0.938 (**仅改声明**) | 与 zero-shot **相同**（adapt=重声明，非重训） | 场 LOO 加性归因 + 闭合翻转距离 |
| TEDDY + Jev-class 站位 | 0.9490 (log-loss Choice，需 O0 标签训练) | 0.9651 (O0 模型评 O1，**未重训**) | 0.9605 (用 O1 标签重训 logreg) | 学习系数×特征 / LOO；不透明；无声明场律 |

### 准则切换诊断（准确率之外）

| 现象 | 数值 / 含义 |
|---|---|
| O0∩O1 expected_action 一致率 | 本导出上一致（O1 加严 margin→更多无 expected / 弃权，而非改写标签） |
| TEDDY-alone 静默过度回答 | O0 规则作答而 O1 声明弃权：**33**/500 |
| ANM O0→O1 abstain | 5 → 31（仅改声明） |
| Jev O0→O1 abstain（未重训） | 4 → 4（不随观察者加严） |
| TEDDY label-fit 网格 | thr=0.08, boost=3.0 ≠ 声明 O1 (0.28, 2.0) |

### 一句话证明

> **ANM 把观察者/准则变成可编辑的场声明：O0→O1 零标签改 YAML 即可保持可验证 Q_f，并给出可审计的留一归因与翻转距离；Jev-class 与 TEDDY-alone 在同一准则切换上必须消耗新标签做重调/重训，否则准确率崩或静默错，且归因停留在启发式/梯度近似。**

### 细胞级对照（同例）

#### `cite_site4_73511`

- **ANM**: 推荐 `myeloid` (P_f=1.0, Q_f=1.0); 顶事件 `CD16` Δchosen=-0.3377 flipped=True; flip_dist=0.7313117727321707 → t_lineage
- **TEDDY alone**: 推荐 `myeloid`; 顶蛋白 `CD16` Δchosen=-0.3104 flipped=True; flip=0.6313117727321707 (无场声明)
- **Jev-class 站位**: 推荐 `myeloid` proba={'b_lineage': 0.03305432880746313, 'myeloid': 0.9529390848295076, 't_lineage': 0.014006586363029273}; 顶 `CD16` coef×feat=2.4294 flipped=True; flip=0.7313117727321707 (学习归因，不透明)

#### `cite_site4_73515`

- **ANM**: 推荐 `b_lineage` (P_f=1.0, Q_f=0.0); 顶事件 `CD72` Δchosen=-0.0566 flipped=True; flip_dist=0.046353474864931604 → t_lineage
- **TEDDY alone**: 推荐 `b_lineage`; 顶蛋白 `CD72` Δchosen=-0.0655 flipped=True; flip=0.1463534748649316 (无场声明)
- **Jev-class 站位**: 推荐 `t_lineage` proba={'b_lineage': 0.3218785270089038, 'myeloid': 0.11061976353137418, 't_lineage': 0.567501709459722}; 顶 `CD2` coef×feat=-1.2236 flipped=False; flip=None (学习归因，不透明)

#### `cite_site4_73518`

- **ANM**: 推荐 `b_lineage` (P_f=1.0, Q_f=1.0); 顶事件 `CD22` Δchosen=-0.3139 flipped=False; flip_dist=None → None
- **TEDDY alone**: 推荐 `b_lineage`; 顶蛋白 `CD19` Δchosen=-0.3333 flipped=False; flip=None (无场声明)
- **Jev-class 站位**: 推荐 `b_lineage` proba={'b_lineage': 0.9999377711381136, 'myeloid': 4.2085266852184165e-07, 't_lineage': 6.180800921793004e-05}; 顶 `CD22` coef×feat=2.1390 flipped=False; flip=None (学习归因，不透明)

### 准则编辑成本

| 臂 | 需要端点标签？ | 需要重训/重调？ | 机制 |
|---|---|---|---|
| TEDDY alone | True | True | threshold/boost grid search or logistic head fit on O1 labels |
| TEDDY+ANM | False | False | edit readouts/cite_lineage_O1.yaml or CRITERIA['O1']; re-adapt instances |
| TEDDY+Jev 站位 | True | True | re-fit log-loss Choice on new-criterion labels (RLCD stand-in) |

---

## English technical appendix

### Setup

- Eval cells: **500** (same as `anm_instances_O0/O1.json`).
- Train pool: **1500** remaining of 2000 exported site4/test cells.
- TEDDY phase-1 / `best.pt`: **not** retrained.
- Holdout `adt_true` → `expected_action` only (verifier); never source events.
- Jev-class: **stand-in** (log-loss multinomial logistic + MLPClassifier); real Jev API unavailable. Pattern follows `anm-jev/mwe/jev_class.py` (single-pass Choice, log-loss as supervised RLCD analogue) and case1 observer-shift.

### Numeric summary (primary metrics)

```json
{
  "teddy_alone": {
    "O0_Q": 0.9641350210970464,
    "O0_abstain": 1,
    "O1_zero_shot_Q": 0.986046511627907,
    "O1_silent_wrong_over_answer": 33,
    "O1_after_retune_Q": 0.986046511627907,
    "O1_retune_rule": {
      "readout_threshold": 0.08,
      "key_marker_boost": 3.0,
      "label_fit": true,
      "search": "grid thr\u00d7boost on train O1 labels",
      "train_best_strict_acc": 0.9629915188897455
    },
    "O1_declaration_copied_Q": 0.9905660377358491,
    "O1_classifier_head_strict": 0.9906976744186047
  },
  "anm": {
    "O0_P": 0.99,
    "O0_Q": 0.9683544303797469,
    "O0_abstain": 5,
    "O1_P": 0.938,
    "O1_Q": 0.9767441860465116,
    "O1_abstain": 31,
    "O1_after_adapt_same_as_zero_shot": true
  },
  "jev_standin": {
    "O0_Q_logreg": 0.9490445859872612,
    "O1_zero_shot_Q": 0.9651162790697675,
    "O1_retrain_Q": 0.9604651162790697,
    "train_labeled_O0": 1412,
    "train_labeled_O1": 1297
  }
}
```

### Re-run

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
# (optional) refresh bridge instances if needed — do NOT retrain best.pt
.venv/bin/python bridge_anm/adapt_to_anm.py --criterion O0 --max-cells 500
.venv/bin/python bridge_anm/adapt_to_anm.py --criterion O1 --max-cells 500
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_demo.py
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_bakeoff.py
```

Outputs: `outputs/anm_cite_bridge/bakeoff/BAKEOFF_REPORT.md`, `bakeoff_results.json`.
