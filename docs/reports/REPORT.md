# ANM × TEDDY CITE Bridge — v0 REPORT

## 中文摘要（给 Daniel）

本桥接把 TEDDY phase-1 CITE（RNA→冻结 TEDDY-G→z→MLP→ADT）的 **类型化证据 δu** 导出为 ANM 有限场实例，并在声明的 P_f / Q_f 读出下跑通场演化、留一归因与翻转距离。 **不是** 临床优越性声明；成功标准是：事件可导出、schema 已声明、真实 CITE 事件上响应/归因路径可跑、准则 O0→O1 可改声明而不重训 TEDDY。

- 导出细胞数：**2000**；事件数：**18000**（test / site4）。
- phase-1 测试集整体 MLP Pearson（metrics.json）：**0.610222**。
- 面板蛋白 MLP Pearson（test_per_protein.json，原样引用）：
  - `CD72`: **0.932699**
  - `CD36`: **0.932459**
  - `CD19`: **0.932114**
  - `CD16`: **0.918683**
  - `CD3`: **0.912537**
  - `CD22`: **0.910817**
  - `CD2`: **0.898578**
  - `CD5`: **0.898245**
  - `CD11c`: **0.896533**
- O0 场跑通：`mean_P_f=0.9900`，`mean_Q_f(defined)=0.9683544303797469`，abstain=5 / 500。
- O1（仅改声明：threshold 0.12→0.28，key_marker_boost 1→2）：`mean_P_f=0.9380`，`mean_Q_f(defined)=0.9767441860465116`，abstain=31 / 500。
- ANM `finite_field_runner` 实际执行：**是**。
- 翻转距离示例（事件 `cite_site4_73511:CD16`）： min |Δvalue| = **0.7313** （myeloid → t_lineage）。
- 留一归因（示例细胞，按 |Δ chosen score|）：
  - `CD16` action=myeloid Δscore=-0.3377 flipped=True
  - `CD11c` action=myeloid Δscore=-0.1589 flipped=False
  - `CD36` action=myeloid Δscore=-0.0348 flipped=False
  - `CD19` action=b_lineage Δscore=0.0000 flipped=False
  - `CD72` action=b_lineage Δscore=0.0000 flipped=False

### 准则编辑（O0→O1）而不重训

只改 `bridge_anm/readouts/cite_lineage_O*.yaml` / `lib/lineage_panels.CRITERIA` 与 实例上的 `readout_threshold` / `key_marker_boost`，重新 `adapt_to_anm.py` + `run_demo.py`。 **不** 触碰 `best.pt`、不重跑 phase-1 训练。

## English technical notes

- TEDDY role: typed evidence engine (δu).
- ANM role: declared field + nested P_f / Q_f + attribution / flip distance.
- Holdout `adt_true` is verifier-only; never admitted as source events.
- Claim boundary: local response diagnosis only.

## Numeric blobs

```json
{
  "O0": {
    "n_instances": 500,
    "mean_P_f": 0.99,
    "mean_Q_f_defined": 0.9683544303797469,
    "n_Q_f_defined": 474,
    "abstain_count": 5,
    "action_counts": {
      "myeloid": 214,
      "b_lineage": 77,
      "t_lineage": 204
    }
  },
  "O1": {
    "n_instances": 500,
    "mean_P_f": 0.938,
    "mean_Q_f_defined": 0.9767441860465116,
    "n_Q_f_defined": 430,
    "abstain_count": 31,
    "action_counts": {
      "myeloid": 210,
      "b_lineage": 70,
      "t_lineage": 189
    }
  },
  "attribution_example": {
    "instance_id": "cite_site4_73511",
    "recommended_action": "myeloid",
    "P_f": 1.0,
    "Q_f": 1.0,
    "action_scores": {
      "b_lineage": 0.10386857495660447,
      "t_lineage": 0.2720092563799766,
      "myeloid": 0.5263200659813642
    },
    "top_events": [
      {
        "event_id": "cite_site4_73511:CD16",
        "protein": "CD16",
        "action": "myeloid",
        "value": 0.9313117727321707,
        "base_recommended": "myeloid",
        "ablated_recommended": "t_lineage",
        "flipped": true,
        "delta_scores": {
          "myeloid": -0.33766991631607557,
          "b_lineage": 0.0,
          "t_lineage": 0.0
        },
        "delta_chosen_score": -0.33766991631607557
      },
      {
        "event_id": "cite_site4_73511:CD11c",
        "protein": "CD11c",
        "action": "myeloid",
        "value": 0.4125626782674178,
        "base_recommended": "myeloid",
        "ablated_recommended": "myeloid",
        "flipped": false,
        "delta_scores": {
          "myeloid": -0.15894008480164024,
          "b_lineage": 0.0,
          "t_lineage": 0.0
        },
        "delta_chosen_score": -0.15894008480164024
      },
      {
        "event_id": "cite_site4_73511:CD36",
        "protein": "CD36",
        "action": "myeloid",
        "value": 0.10359825398825429,
        "base_recommended": "myeloid",
        "ablated_recommended": "myeloid",
        "flipped": false,
        "delta_scores": {
          "myeloid": -0.03483497225919174,
          "b_lineage": 0.008194185228326545,
          "t_lineage": 0.015147337640516567
        },
        "delta_chosen_score": -0.03483497225919174
      },
      {
        "event_id": "cite_site4_73511:CD19",
        "protein": "CD19",
        "action": "b_lineage",
        "value": 0.11109760767958753,
        "base_recommended": "myeloid",
        "ablated_recommended": "myeloid",
        "flipped": false,
        "delta_scores": {
          "myeloid": 0.0,
          "b_lineage": -0.03142614416481457,
          "t_lineage": 0.0
        },
        "delta_chosen_score": 0.0
      },
      {
        "event_id": "cite_site4_73511:CD72",
        "protein": "CD72",
        "action": "b_lineage",
        "value": 0.15325553212810547,
        "base_recommended": "myeloid",
        "ablated_recommended": "myeloid",
        "flipped": false,
        "delta_scores": {
          "myeloid": 0.0,
          "b_lineage": -0.04429885744277525,
          "t_lineage": 0.0
        },
        "delta_chosen_score": 0.0
      },
      {
        "event_id": "cite_site4_73511:CD22",
        "protein": "CD22",
        "action": "b_lineage",
        "value": 0.12751951521193522,
        "base_recommended": "myeloid",
        "ablated_recommended": "myeloid",
        "flipped": false,
        "delta_scores": {
          "myeloid": 0.0,
          "b_lineage": -0.04081361460112008,
          "t_lineage": 0.0
        },
        "delta_chosen_score": 0.0
      },
      {
        "event_id": "cite_site4_73511:CD3",
        "protein": "CD3",
        "action": "t_lineage",
        "value": 0.17236604157159863,
        "base_recommended": "myeloid",
        "ablated_recommended": "myeloid",
        "flipped": false,
        "delta_scores": {
          "myeloid": 0.0,
          "b_lineage": 0.0,
          "t_lineage": -0.059558231634370934
        },
        "delta_chosen_score": 0.0
      },
      {
        "event_id": "cite_site4_73511:CD2",
        "protein": "CD2",
        "action": "t_lineage",
        "value": 0.4769148923214528,
        "base_recommended": "myeloid",
        "ablated_recommended": "myeloid",
        "flipped": false,
        "delta_scores": {
          "myeloid": 0.0,
          "b_lineage": 0.0,
          "t_lineage": -0.15956363149463576
        },
        "delta_chosen_score": 0.0
      },
      {
        "event_id": "cite_site4_73511:CD5",
        "protein": "CD5",
        "action": "t_lineage",
        "value": 0.1839910346538489,
        "base_recommended": "myeloid",
        "ablated_recommended": "myeloid",
        "flipped": false,
        "delta_scores": {
          "myeloid": 0.0,
          "b_lineage": 0.0,
          "t_lineage": -0.07089122152712696
        },
        "delta_chosen_score": 0.0
      }
    ]
  },
  "flip_example": {
    "event_id": "cite_site4_73511:CD16",
    "base_value": 0.9313117727321707,
    "base_recommended": "myeloid",
    "flip_exists": true,
    "min_flip_distance": 0.7313117727321707,
    "flip_to": "t_lineage",
    "flip_value": 0.2
  }
}
```
