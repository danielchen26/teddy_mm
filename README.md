# TEDDY × ANM — CITE bridge

**Keep TEDDY frozen. Make its answers editable, honest and auditable with ANM.**

**Site:** [danielchen26.github.io/teddy_mm](https://danielchen26.github.io/teddy_mm/) — interactive charts for every block · [decision-loop report](https://danielchen26.github.io/teddy_mm/anm-loop.html) · source [`docs/`](docs/)

![Global frame: CITE cell → frozen TEDDY map → ANM editable decisions](docs/assets/infographics/00_global_frame_light.png)

## 1 · The problem

TEDDY (frozen TEDDY-G 70M; `z_512` = mean-pool last layer @ context 1024) maps a cell's RNA to its surface proteins (phase-1 full-ADT Pearson ≈ 0.61). As a readout on its own it leaves six gaps. Each became a block with its own test.

| Gap | TEDDY alone today | Block |
|---|---|---|
| Change the question | Hard-wired rule; keeps answering the old question (1,752 cells that O2 would hold back) | 1 |
| Weak or missing input | Answers every cell: 0 abstentions at Q 0.63 under `adt_only` | 2 |
| Explain a call | No declared decision field → no per-marker leave-one-out or flip distance | 3 |
| Know where to trust it | No graded signal to gate on | 4 |
| Price of a new criterion | A trained readout needs 50–500 new labels per switch | 5 |
| Is the embedding enough? | Myeloid↔T pairs at cos 0.982 in `z_512`; readout falsely agrees 85% | 6 |

**Not a goal:** beating Pearson ≈ 0.61 or retraining TEDDY. TEDDY's predictions are the input to every block, unchanged.

## 2 · How we compare — one protocol, every block

- **Same frozen TEDDY** — `best.pt` is never retrained; every arm reads the same predicted 9-protein panel (typed evidence δu).
- **Same cells** — holdout site4/test, n = 16,750 (out-of-site check: val_non_site4, n = 2,000).
- **Truth only verifies** — holdout true ADT sets the expected answer and never enters ANM's evidence (one labelled exception: Block 6's complement probe).
- **Decision metrics first** — abstention, label cost, attribution, coverage vs Q, false agreement. Accuracy is reported, never the win condition.

| Arm | What it is | Adopt a new criterion | Abstention follows it? | Attribution |
|---|---|---|---|---|
| TEDDY alone | Predicted panel → lineage score → thresholded argmax (O0 rule) | Re-tune on labels | No — keeps its old rule | Panel heuristic |
| **TEDDY + ANM** | Same predictions as typed events → ANM `finite_field` under a declared observer YAML (O0 / O1 / O2) | **Edit YAML · 0 labels** | **Yes — 317 → 834 → 3,373** | **Closed-form LOO + flip distance** |
| Train + Jev | Head trained on labels over the same features (log-loss logistic / MLP / grid; Jev-class stand-in) | Retrain on new labels | No — 185 → 185 → 185 | Learned coefficient × feature |

> **help = (TEDDY + ANM) − (TEDDY alone)** on identical δu, with Train + Jev as the "just train harder" control.

Observers live in [`bridge_anm/readouts/cite_lineage_O{0,1,2}.yaml`](bridge_anm/readouts/): O1 raises the threshold 0.12 → 0.28 and doubles key-marker weight (0 expected calls change); O2 switches to key-marker priority with lineage weights (12.48% of expected calls change).

## 3 · Six blocks

| # | Block | Method (workflow) | Key result | Proves ANM helps TEDDY by… | Report |
|---|---|---|---|---|---|
| 1 | [Edit the question](https://danielchen26.github.io/teddy_mm/#edit) | Declare O0 / O1 / O2 in YAML → run field → compare abstention with TEDDY alone (O0 rule) and an O0-trained head | ANM abstain 317 → 834 → 3,373; head 185 flat; TEDDY alone over-answers 1,752 at O2 | serving a new question with 0 labels, no retrain | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 2 | [Missing modality](https://danielchen26.github.io/teddy_mm/#missing) | `rna_only` / `adt_only` / `joint` masks as source ablation + declared reliability gate (×0.35) | `adt_only`: abstain 0 → 1,664, mean P_f 0.868; O2 abstains on all 12,563 | honest silence when TEDDY's input is weak | [MISSING_MODALITY](docs/reports/MISSING_MODALITY_ANM_DEMO.md) |
| 3 | [Attribution](https://danielchen26.github.io/teddy_mm/#attr) | Leave one typed event out → top marker → flip distance; bootstrap B = 200, permutation n = 100 | 16,433 cells · top-1 flip 0.277 · median flip distance 0.75 · p ≈ 0.0099 | an auditable per-marker "why" for every call | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 4 | [Scope gate](https://danielchen26.github.io/teddy_mm/#scope) | Gate on soft_P (τ sweep); ANM hard scope = LOO flip-sensitive cells vs random same-n | O0 Q 0.945 → 0.998 at 40% coverage; hard-scope Q 0.810 vs 0.946 random | showing where TEDDY's calls can be trusted | [SCOPE_REFINE_PROOF](docs/reports/SCOPE_REFINE_PROOF.md) |
| 5 | [Zero-label transfer](https://danielchen26.github.io/teddy_mm/#labels) | Heads trained on 50 → 10,619 O2 labels, gate calibrated to ANM's O2 abstain rate | Match needs ≈ 50 (grid) / 200 (MLP) / 500 (logistic) labels; ANM needs 0 | new criteria at zero label cost | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 6 | [Decision loop](https://danielchen26.github.io/teddy_mm/anm-loop.html) | Veto `z_512` on must-separate pairs → complement (typed true ADT + O1) → verify on the same pairs | Myeloid↔T false agreement 0.85 → 0.40, Q 0.74 → 0.99; z⁺ feature swap −0.27 | a principled way to patch TEDDY's blind spot | [REVERSE_LOOP_SMALL](docs/reports/REVERSE_LOOP_SMALL.md) · [ANM_HELPS_TEDDY_LOOP](docs/reports/ANM_HELPS_TEDDY_LOOP.md) |

**Accuracy check (Block 1, like for like, Q among decided cells):** TEDDY alone 0.947 / 0.974 / 0.943 vs ANM 0.950 / 0.972 / 0.960 (O0 / O1 / O2). Roughly tied, and not the claim. HARD_PROOF's three-arm table lists ANM's mean Q_f, which counts abstentions as misses (0.939 / 0.958 / 0.821).

<details>
<summary>Posters (one per block)</summary>

| | |
|---|---|
| ![Edit the question](docs/assets/infographics/01_edit_question_light.png) | ![Missing modality](docs/assets/infographics/02_missing_modality_light.png) |
| ![Attribution](docs/assets/infographics/03_attribution_light.png) | ![Scope gate](docs/assets/infographics/04_scope_gate_light.png) |
| ![Label cost](docs/assets/infographics/05_label_cost_light.png) | ![Decision loop](docs/assets/infographics/06_anm_loop_poster_light.png) |
| ![Myeloid↔T before/after](docs/assets/infographics/07_myeloid_t_before_after_light.png) | ![ANM vs feature swap](docs/assets/infographics/08_anm_vs_feature_swap_light.png) |

Dark variants: `docs/assets/infographics/*_dark.png`. GIFs: `docs/assets/hard_proof/`, `docs/assets/missing_modality/`.
</details>

## 4 · Claim boundary

| | |
|---|---|
| **Claimed (Mode B)** | TEDDY typed evidence → ANM `finite_field` decisions · editable observer YAML (O0/O1/O2), no TEDDY retrain · missing-modality abstention · LOO/flip attribution on decisions · scope gate and zero-label transfer vs Train + Jev · decision loop veto → complement → verify on must-separate pairs |
| **Not claimed** | Mode A (residual stream as field, layer Jacobian, in-silico gene perturbs) · reverse loop as representation response to perturbs · gated ±ε G1–G4 · Pearson > phase-1 full-ADT ≈ 0.61 · fusion audit · ATAC/chromatin · clinical superiority · live Jev API |
| **`z_512`** | Mean-pool of last-layer tokens at context length **1024** (not pretrain 2048, not a disease token), L2-normalized 512-D (`z_rna_512.npy`). Compact `z_rna_export.npy` (`z_keep=32`) is not the probe space. |
| **Honest weak spots** | ADT-only is a weak channel (panel Pearson 0.446). Block 4 trades coverage for Q and O2 dips at very low coverage. Block 6 uses true ADT as an explicit complement probe; production keeps holdout ADT verifier-only. |

Full scope: [`docs/MODE_B_SCOPE.md`](docs/MODE_B_SCOPE.md).

## 5 · Reproduce (no TEDDY retrain)

| Block | Script | Report |
|---|---|---|
| Export (once) | `bridge_anm/export_cite_events.py --n-cells 0 --ood-n 2000` | [STAT_PROOF](docs/reports/STAT_PROOF.md) |
| 1 · 3 · 5 | `bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100` | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 2 | `bridge_anm/export_missing_modality_events.py` → `bridge_anm/run_missing_modality_demo.py` | [MISSING_MODALITY](docs/reports/MISSING_MODALITY_ANM_DEMO.md) |
| 4 | `bridge_anm/run_scope_refine_proof.py` | [SCOPE_REFINE_PROOF](docs/reports/SCOPE_REFINE_PROOF.md) |
| 6 | `bridge_anm/reverse_step1_must_separate.py` → `bridge_anm/reverse_loop_small.py` | [REVERSE_LOOP_SMALL](docs/reports/REVERSE_LOOP_SMALL.md) |

```bash
cd teddy_mm
PYTHONPATH=/path/to/ANM:. .venv/bin/python bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100
PYTHONPATH=/path/to/ANM:. .venv/bin/python bridge_anm/run_scope_refine_proof.py
PYTHONPATH=/path/to/ANM:. .venv/bin/python bridge_anm/reverse_loop_small.py
```

Figures: `python scripts/make_story_infographics.py` · loop posters: `/usr/bin/python3 scripts/make_anm_helps_teddy_loop_figures.py`.

### 中文摘要

TEDDY 冻结不动（RNA→蛋白，Pearson ≈ 0.61）。ANM 在同一份冻结预测上加一层可编辑、会弃权、可审计的决策层。六个模块（改问题 · 缺模态 · 归因 · 范围门控 · 零标签迁移 · 决策闭环）都用同一套对照：TEDDY alone vs TEDDY + ANM vs 训练头，只证明“帮助 TEDDY”，不证明 Pearson 胜利，不重训，非临床。

---

# TEDDY multimodal phase-1

冻结 TEDDY 编 RNA，预测 NeurIPS 2021 BMMC CITE 的表面蛋白。
对照：MLP vs latent flow matching。

为 Mac M4 Max 128GB（MPS）写的。不要在这台机器上预训练第二模态。

## 一次性安装

```bash
cd teddy_mm
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export PYTORCH_ENABLE_MPS_FALLBACK=1
```

TEDDY-G 70M 权重应在 `../teddy_mwe/ckpt/teddy_g_70M/`（先跑 `teddy_mwe/setup.sh`）。

## 正式数据（约 587 MB）

```bash
bash scripts/01_download_cite.sh
python scripts/02_prepare_cite.py
python scripts/03_embed_rna.py --device mps --seq-len 1024 --batch-size 16
python scripts/04_train.py --device mps
python scripts/05_eval.py --device mps
```

编码 9 万细胞、seq=1024，在 M4 Max 上大概数小时，只做一次。  
想先冒烟：

```bash
python scripts/02_prepare_cite.py --max-cells 4000
python scripts/03_embed_rna.py --seq-len 256 --batch-size 8 --device mps
python scripts/04_train.py --epochs 8 --device mps
```

没有 GEO 文件时可用合成数据测训练环：

```bash
python scripts/00_smoke_synthetic.py
python scripts/04_train.py --processed data/processed/cite_smoke --out outputs/smoke --epochs 5 --device cpu
```

## 切分

- test = `site4`（数据集里的held-out site）
- val = 从剩余 donor 里抽 10%
- train = 其余

## 现在做什么 / 还不做什么

做：RNA → ADT，冻结 TEDDY，MLP 基线 + FM。  
不做：ATAC foundation、解冻 400M、三模态联合生成。

双向 / 模态缺失是下一期：给 ADT 加 encoder，训练时 drop RNA。

## Phase-2 bidirectional CITE (scaffold)

Independent of full phase-1 90k embed. ADT encoder (CLR/log → shallow MLP) +
modality dropout (~18% drop RNA) + latent FM with `cond = z_obs`. No ATAC /
perturb-seq / 160M.

**Do not** write into `data/processed/cite/` or `outputs/cite_phase1/` while
phase-1 is running. Use a separate pack:

```bash
# Option A — tiny synthetic (immediate)
python scripts/00_smoke_synthetic.py
python scripts/06_train_bidirectional.py \
  --processed data/processed/cite_smoke --out outputs/cite_phase2 \
  --epochs 5 --device cpu

# Option B — 4k real CITE in a NEW dir (does not touch data/processed/cite)
python scripts/02_prepare_cite.py --max-cells 4000 --out data/processed/cite_phase2_4k
# Prefer CPU while phase-1 owns MPS:
python scripts/03_embed_rna.py --processed data/processed/cite_phase2_4k \
  --device cpu --seq-len 256 --batch-size 4
python scripts/06_train_bidirectional.py \
  --processed data/processed/cite_phase2_4k --out outputs/cite_phase2 \
  --epochs 8 --device mps   # or cpu if MPS still busy

# After full phase-1 z_rna exists, same script on the full pack:
# python scripts/06_train_bidirectional.py \
#   --processed data/processed/cite --out outputs/cite_phase2 --epochs 40 --device mps
```

Config: `configs/phase2_cite.yaml`. Code: `teddy_mm/bidirectional.py`,
`teddy_mm/modality_dropout.py`, `scripts/06_train_bidirectional.py`.
