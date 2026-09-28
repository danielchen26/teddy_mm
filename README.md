<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/readme/hero-dark.png">
  <img alt="TEDDY × ANM — keep TEDDY frozen; make its answers editable, honest and auditable with ANM" src="docs/assets/readme/hero-light.png" width="100%">
</picture>

<p>
  <a href="https://danielchen26.github.io/teddy_mm/"><b>Live site</b></a> ·
  <a href="https://danielchen26.github.io/teddy_mm/anm-loop.html"><b>Decision-loop report</b></a> ·
  <a href="#quickstart">Quickstart</a> ·
  <a href="docs/reports/">Proof reports</a> ·
  <a href="#claim-boundary">Claim boundary</a>
</p>

<img alt="TEDDY: frozen, never retrained" src="https://img.shields.io/badge/TEDDY-frozen%2C%20never%20retrained-c8641f?style=flat-square">
<img alt="scope: decision layer (Mode B)" src="https://img.shields.io/badge/scope-decision%20layer%20(Mode%20B)-008c7e?style=flat-square">
<img alt="data: CITE-seq, site4, 16,750 cells" src="https://img.shields.io/badge/data-CITE--seq%20%C2%B7%20site4%20%C2%B7%2016%2C750%20cells-6554c9?style=flat-square">
<img alt="numbers: from committed reports" src="https://img.shields.io/badge/numbers-from%20committed%20reports-4c566a?style=flat-square">

</div>

## What this is

**TEDDY** ([Merck TEDDY-G 70M](https://huggingface.co/Merck/TEDDY) · [paper](https://arxiv.org/abs/2503.03485)) is a single-cell foundation model. Kept frozen, it maps a cell's RNA to 9 surface proteins on CITE-seq (NeurIPS 2021 BMMC, holdout site4, phase-1 full-ADT Pearson ≈ 0.61).

**ANM** (Active Neural Matter) is a decision layer. It reads TEDDY's predictions as *typed evidence* and decides per cell under a *declared observer*: call a lineage or abstain, with a workability score P<sub>f</sub>, an exactness check Q<sub>f</sub>, and a per-marker explanation.

**This repository** is the bridge between the two, plus six head-to-head tests of where ANM *helps* TEDDY. TEDDY's weights and predictions never change; only the decision layer does.

## Why ANM

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/readme/why-dark.png">
  <img alt="Six results: 0 labels to change the question; 0 → 1,664 honest abstentions with RNA missing; 16,433 explained cells; Q 0.945 → 0.998 when gating TEDDY's calls; 0 vs 50–500 labels to adopt a new criterion; false agreement 0.85 → 0.40" src="docs/assets/readme/why-light.png" width="100%">
</picture>

<details>
<summary>The same results as a table</summary>

| A biologist asks… | TEDDY alone | TEDDY + ANM |
|---|---|---|
| Can I change the question (O0 → O2)? | Keeps its old rule; answers 1,752 cells the new criterion holds back | A declaration edit with 0 labels; abstentions follow: 317 → 834 → 3,373 |
| What if RNA is missing? | 0 abstentions at Q 0.63 | 1,664 abstentions; mean P<sub>f</sub> 0.868 |
| Why this call? | No declared field to explain | Leave-one-out + flip distance on 16,433 cells (permutation p ≈ 0.0099) |
| Which calls can I trust? | Every call looks equally sure | Q of TEDDY's calls 0.945 → 0.998 keeping the top 40% by soft_P |
| What does a new criterion cost? | A trained head needs 50–500 new labels | 0 labels |
| Can it tell apart cells that look identical to TEDDY? | False agreement 0.85 on myeloid↔T pairs | 0.40 after veto → complement → verify |

Accuracy is not the claim: compared like for like (Q among decided cells), ANM is within 0.02 of TEDDY alone.

</details>

## How it works

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/assets/readme/how-dark.png">
  <img alt="Step 1: TEDDY stays frozen and exports typed evidence. Step 2: declare the question as an observer. Step 3: ANM decides per cell and explains. Compared against TEDDY alone and a trained head." src="docs/assets/readme/how-light.png" width="100%">
</picture>

An observer is a short declaration, not a trained model. O1 is O0 with three numbers changed:

```python
# bridge_anm/lib/lineage_panels.py → CRITERIA  (readable mirror: bridge_anm/readouts/cite_lineage_O1.yaml)
"O1": {
    "score_mode": "equal_panel_mean",
    "readout_threshold": 0.28,   # O0: 0.12  → stricter calls
    "key_marker_boost": 2.0,     # O0: 1.0   → CD19 / CD3 / CD16 count double
    "expected_margin": 0.12,     # O0: 0.05
    ...
},
```

## How we test "ANM helps TEDDY"

Every block follows the same protocol:

- **Same frozen TEDDY.** `best.pt` is never retrained; every arm reads the same predicted 9-protein panel.
- **Same cells.** Holdout site4/test, n = 16,750, with an out-of-site check on val_non_site4 (n = 2,000).
- **Truth only verifies.** Holdout true ADT scores the calls and never enters ANM's evidence. The one labelled exception is Block 6's complement probe.
- **Decision metrics first.** Abstention, label cost, attribution, coverage vs Q and false agreement. Accuracy is reported, never the win condition.

| Arm | What it is | Adopting a new criterion |
|---|---|---|
| TEDDY alone | Predicted panel → lineage score → thresholded argmax | Re-tune the rule on new labels |
| **TEDDY + ANM** | Same predictions as typed events → ANM `finite_field` under a declared observer | **Edit the declaration: 0 labels, no retrain** |
| Train + Jev | A head trained on labels over the same features (logistic / MLP / grid; Jev-class stand-in) | Collect labels and retrain |

## Results by block

| # | Block | Method | Key result | Report |
|---|---|---|---|---|
| 1 | [Edit the question](https://danielchen26.github.io/teddy_mm/#edit) | Declare O0 / O1 / O2 → run field → compare abstention with TEDDY alone and an O0-trained head | ANM 317 → 834 → 3,373; head 185 flat; O0 → O2 rewrites 12.48% of calls | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 2 | [Missing modality](https://danielchen26.github.io/teddy_mm/#missing) | `rna_only` / `adt_only` / `joint` masks as source ablation + a declared reliability gate | `adt_only`: abstain 0 → 1,664; O2 abstains on all 12,563 cells | [MISSING_MODALITY](docs/reports/MISSING_MODALITY_ANM_DEMO.md) |
| 3 | [Attribution](https://danielchen26.github.io/teddy_mm/#attr) | Leave one typed event out → top marker → flip distance; bootstrap B = 200, permutation n = 100 | 16,433 cells · top-1 flip 0.277 · p ≈ 0.0099 | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 4 | [Scope gate](https://danielchen26.github.io/teddy_mm/#scope) | Keep TEDDY's call only where ANM soft_P ≥ τ; ANM hard scope vs random same-n | Q 0.945 → 0.998 at 40% coverage; hard cells 0.810 vs 0.946 | [SCOPE_REFINE_PROOF](docs/reports/SCOPE_REFINE_PROOF.md) |
| 5 | [Zero-label transfer](https://danielchen26.github.io/teddy_mm/#labels) | Heads on 50 → 10,619 O2 labels, gate calibrated to ANM's O2 abstain rate | Match needs ≈ 50 / 200 / 500 labels; ANM needs 0 | [HARD_PROOF](docs/reports/HARD_PROOF.md) |
| 6 | [Decision loop](https://danielchen26.github.io/teddy_mm/anm-loop.html) | Veto `z_512` on must-separate pairs → complement with typed true ADT + O1 → verify | False agreement 0.85 → 0.40; Q 0.74 → 0.99 | [REVERSE_LOOP_SMALL](docs/reports/REVERSE_LOOP_SMALL.md) |

Charts, data tables and every number's source: **[danielchen26.github.io/teddy_mm](https://danielchen26.github.io/teddy_mm/)**.

## Quickstart

**You need**

- Python 3.11 and the packages in `requirements.txt`
- TEDDY-G 70M weights at `../teddy_mwe/ckpt/teddy_g_70M/` (from [Merck/TEDDY](https://huggingface.co/Merck/TEDDY))
- ANM on `PYTHONPATH` (Active Neural Matter; the research code is currently private, so ask the authors for access)
- About 600 MB for the CITE-seq data (GEO GSE194122)

**1 · Environment**

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
export PYTORCH_ENABLE_MPS_FALLBACK=1          # Apple silicon
```

**2 · Build TEDDY's evidence (once)** — TEDDY stays frozen; only a small RNA → protein head is trained

```bash
bash scripts/01_download_cite.sh               # GEO GSE194122, ~587 MB
python scripts/02_prepare_cite.py              # holdout = site4
python scripts/03_embed_rna.py --device mps --seq-len 1024 --batch-size 16
python scripts/04_train.py --device mps
python scripts/05_eval.py --device mps
```

**3 · Export the predictions as typed evidence (δu)**

```bash
python bridge_anm/export_cite_events.py --n-cells 0 --ood-n 2000
```

**4 · Run the tests**

```bash
export PYTHONPATH=/path/to/ANM:.
python bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100   # blocks 1 · 3 · 5
python bridge_anm/run_scope_refine_proof.py                                # block 4
python bridge_anm/export_missing_modality_events.py --n-cells 0            # block 2 (needs phase-2 ckpt)
python bridge_anm/run_missing_modality_demo.py
python bridge_anm/reverse_step1_must_separate.py                           # block 6
python bridge_anm/reverse_loop_small.py
```

Each script writes a Markdown report and a JSON file under `outputs/anm_cite_bridge/<block>/`; the published copies live in [`docs/reports/`](docs/reports/).

**5 · Ask a different question**

Edit an observer's entry in `CRITERIA` ([`bridge_anm/lib/lineage_panels.py`](bridge_anm/lib/lineage_panels.py)): its threshold, key-marker boost or lineage weights. Keep its YAML mirror in [`bridge_anm/readouts/`](bridge_anm/readouts/) in sync, then re-run the proof script. A new observer id also needs adding to `CRIT_IDS` at the top of the proof script. TEDDY is not retrained and no labels are used; holdout truth only scores the result.

## Repository layout

```text
teddy_mm/
├── scripts/                   # phase 1: download → prepare → embed (frozen TEDDY) → train head → eval
├── teddy_mm/                  # TEDDY encoder wrapper, RNA → protein models, phase-2 scaffold
├── bridge_anm/                # TEDDY → ANM bridge
│   ├── export_*.py            #   predictions → typed evidence events
│   ├── run_*.py, reverse_*.py #   the six tests
│   ├── lib/lineage_panels.py  #   marker panels + CRITERIA (the observers the code runs)
│   ├── readouts/              #   readable YAML mirror of each observer
│   └── schemas/               #   ANM finite-field schema
├── configs/                   # phase-1 / phase-2 training configs
├── outputs/anm_cite_bridge/   # generated reports + JSON
└── docs/                      # GitHub Pages site, proof reports, figures
```

## Claim boundary

**Claimed (Mode B: TEDDY as a typed evidence source for ANM decisions)**

- TEDDY typed evidence → ANM `finite_field` decisions, with editable observer declarations (O0 / O1 / O2) and no TEDDY retrain
- Missing-modality abstention; leave-one-out and flip attribution on decisions
- Scope gate and zero-label observer transfer vs Train + Jev
- A decision loop (veto → complement → verify) on must-separate pairs

**Not claimed**

- Mode A: residual stream as field, layer Jacobian, in-silico gene perturbations
- The reverse loop as a representation response to perturbations
- Pearson above phase-1 full-ADT ≈ 0.61 · fusion audit · ATAC / chromatin · clinical superiority · the live Jev API

**Honest weak spots.** ADT-only is a weak channel (panel Pearson 0.446). Block 4 trades coverage for Q, and O2 dips at very low coverage. Block 6 uses true ADT as an explicit complement probe; production keeps holdout ADT verifier-only.

**`z_512`** is the mean-pool of last-layer tokens at context length 1024 (not pretrain 2048, not a disease token), L2-normalized 512-D (`z_rna_512.npy`).

Full scope: [`docs/MODE_B_SCOPE.md`](docs/MODE_B_SCOPE.md).

## Reports

[HARD_PROOF](docs/reports/HARD_PROOF.md) · [SCOPE_REFINE_PROOF](docs/reports/SCOPE_REFINE_PROOF.md) · [MISSING_MODALITY_ANM_DEMO](docs/reports/MISSING_MODALITY_ANM_DEMO.md) · [REVERSE_LOOP_SMALL](docs/reports/REVERSE_LOOP_SMALL.md) · [REVERSE_STEP1](docs/reports/REVERSE_STEP1.md) · [ANM_HELPS_TEDDY_LOOP](docs/reports/ANM_HELPS_TEDDY_LOOP.md) · [CHEAP_PROBES_RANKING](docs/reports/CHEAP_PROBES_RANKING.md) · [STAT_PROOF](docs/reports/STAT_PROOF.md) · [BAKEOFF_REPORT](docs/reports/BAKEOFF_REPORT.md)

Regenerate figures: `python scripts/make_story_infographics.py` (site posters) · `node scripts/render_readme_figures.mjs` (README figures, needs Playwright).

## References

- TEDDY: Merck TEDDY-G 70M on [Hugging Face](https://huggingface.co/Merck/TEDDY); paper [arXiv:2503.03485](https://arxiv.org/abs/2503.03485)
- Data: NeurIPS 2021 BMMC CITE-seq (Open Problems multimodal), GEO [GSE194122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122)
- ANM: Active Neural Matter (research code, currently private)

> **中文简介**：TEDDY 冻结不动（RNA→蛋白，Pearson ≈ 0.61）。ANM 在同一份冻结预测上加一层可编辑、会弃权、可审计的决策层。六个模块（改问题 · 缺模态 · 归因 · 范围门控 · 零标签迁移 · 决策闭环）都用同一套对照：TEDDY alone vs TEDDY + ANM vs 训练头，只证明“帮助 TEDDY”，不证明 Pearson 胜利，不重训，非临床。

<details>
<summary><b>Phase-1 / phase-2 training notes (中文)</b></summary>

### TEDDY multimodal phase-1

冻结 TEDDY 编 RNA，预测 NeurIPS 2021 BMMC CITE 的表面蛋白。
对照：MLP vs latent flow matching。

为 Mac M4 Max 128GB（MPS）写的。不要在这台机器上预训练第二模态。

#### 一次性安装

```bash
cd teddy_mm
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
export PYTORCH_ENABLE_MPS_FALLBACK=1
```

TEDDY-G 70M 权重应在 `../teddy_mwe/ckpt/teddy_g_70M/`（先跑 `teddy_mwe/setup.sh`）。

#### 正式数据（约 587 MB）

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

#### 切分

- test = `site4`（数据集里的held-out site）
- val = 从剩余 donor 里抽 10%
- train = 其余

#### 现在做什么 / 还不做什么

做：RNA → ADT，冻结 TEDDY，MLP 基线 + FM。  
不做：ATAC foundation、解冻 400M、三模态联合生成。

双向 / 模态缺失是下一期：给 ADT 加 encoder，训练时 drop RNA。

#### Phase-2 bidirectional CITE (scaffold)

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

</details>
