# STAT_PROOF — ANM × TEDDY CITE（规模化统计证明）

## 中文摘要（给 Daniel）

在 **同一批 δu**（site4/test **n=16750**；OOD **n=2000**，val_non_site4）上做三臂对照 + O0/O1/O2 观察者。**准确率是次要的**；根本变化 = **可编辑观察者** + **弃权/P_f** + **可审计归因**（GT 仅作 verifier）。

### 准则分歧（expected_action）

| 对比 | both_defined | disagree | disagree_rate |
|---|---:|---:|---:|
| O0_vs_O1 | 14488 | 0 | 0.0000 |
| O0_vs_O2 | 14379 | 1794 | 0.1248 |
| O1_vs_O2 | 13400 | 1414 | 0.1055 |

- O0 vs O1：多为加严弃权（标签改写率 **0.0000**）。
- O0 vs O2：key-marker 优先 + B/T 权重 **改写** 标签（分歧率 **0.1248**）。

### 三臂 Q vs GT（site4）

| 臂 | O0 Q | O0 abstain | O1 Q | O1 abstain | O2 Q | O2 abstain | 准则变更成本 |
|---|---:|---:|---:|---:|---:|---:|---|
| TEDDY alone | 0.9472 | 101 | 0.9735 | 881 | 0.9428 | 1796 | 需标签重调 |
| TEDDY+ANM | 0.9394 | 317 | 0.9578 | 834 | 0.8208 | 3373 | 仅改声明 |
| TEDDY+Jev站位 (O1/O2=O0模型零样本) | 0.9552 | 185 | 0.9733 | 185 | 0.8277 | 185 | 需标签重训 |
| TEDDY+Jev站位 (各准则重训后) | 0.9552 | 185 | 0.9740 | 157 | 0.9232 | 50 | 消耗新标签 |

### 归因 / Bootstrap / Permutation（ANM，O0 场）

- 归因细胞数：**246**
- top-1 蛋白分布：`{'CD5': 59, 'CD16': 33, 'CD2': 49, 'CD22': 29, 'CD36': 36, 'CD72': 22, 'CD11c': 12, 'CD3': 6}`
- top-1 flip rate：**0.24390243902439024**
- flip distance quantiles：`{'q25': 0.7, 'q50': 0.75, 'q75': 0.8233382545240459, 'n_with_flip': 59, 'n_no_flip': 187}`
- bootstrap top-1：mode=CD5 frac=0.2398 boot_mean=0.24109756097560975 ci95=[0.2032520325203252, 0.2886178861788618]
- permutation：null_max_mean=0.17875799721756264 p=0.024390243902439025 (n_perm=40)

### 静默过度回答 / 弃权是否跟随声明

- TEDDY-alone O0 规则 vs O1 声明弃权：`{'count': 780, 'rate': 0.046567164179104475}`
- TEDDY-alone O0 规则 vs O2 声明弃权：`{'count': 1752, 'rate': 0.10459701492537313}`
- ANM abstain O0/O1/O2：[317, 834, 3373]
- Jev stand-in abstain track：`{'O0_model_abstain_on_O0': 185, 'O0_model_abstain_on_O1': 185, 'O0_model_abstain_on_O2': 185, 'note': 'Jev stand-in abstain is confidence gate, not observer declaration'}`

### 一句话证明

> **ANM 把观察者变成可编辑场声明：O0→O1 加严弃权、O0→O2 改写 expected_action，均零端点标签；同 δu 上保持可验证 Q_f / P_f，并给出可审计 LOO 归因与翻转距离（bootstrap 稳定、permutation 可检）。TEDDY-alone 与 Jev-class 站位在准则切换上必须消耗新标签重调/重训，否则静默过度回答或弃权不跟随声明。准确率非胜负条件；非临床声明。**


### Blockers / honesty notes

- **Disk pressure** on the Mac forced deleting regenerable `data/raw/*.h5ad` and caches mid-run; processed `cite_arrays.npz` + `z_rna.npy` + `best.pt` retained. Full-scale field_artifact JSON (O0/O1/O2 × 16k) intentionally **not** written — metrics streamed in `run_stat_proof.py`.
- **Real Jev API** not called; arm is sklearn log-loss Choice stand-in (`anm-jev/mwe/jev_class.py` intent).
- Attribution bootstrap n=100 / permutation n=40 on **250** attr cells (not all 16k) for LOO+flip cost; Q/P_f evaluated on **full** 16750.
- No retrain of `best.pt`; no clinical claim; no Perturb/GFlowNet/160M/ATAC.

### 声明边界

- Holdout `adt_true` = verifier only；未拟合 ANM 场律到标签。
- 未重训 `best.pt`。未声称临床优越。Jev 为 log-loss Choice **站位**。
- 无 Perturb / GFlowNet / 160M / ATAC。

### OOD 摘要

```json
{
  "n": 2000,
  "slice": "val_non_site4",
  "Q": {
    "O0": {
      "teddy": 0.9816980172852059,
      "anm": 0.9816980172852059,
      "teddy_abstain": 1,
      "anm_abstain": 1
    },
    "O1": {
      "teddy": 0.9904661016949152,
      "anm": 0.9899257688229056,
      "teddy_abstain": 19,
      "anm_abstain": 21
    },
    "O2": {
      "teddy": 0.9170428893905191,
      "anm": 0.9783037475345168,
      "teddy_abstain": 111,
      "anm_abstain": 431
    }
  }
}
```

### Re-run

```bash
cd /Users/tianchichen/Documents/GitHub/teddy_mm
PYTORCH_ENABLE_MPS_FALLBACK=1 .venv/bin/python bridge_anm/export_cite_events.py --n-cells 0 --ood-n 2000
PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_stat_proof.py
```

Outputs: `/Users/tianchichen/Documents/GitHub/teddy_mm/outputs/anm_cite_bridge/stat_proof/STAT_PROOF.md`, `stat_proof_results.json`.

Generated: 2026-09-26 14:18:34 EDT
