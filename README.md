# TEDDY × ANM CITE bridge

Public landing (GitHub Pages): **[docs/index.html](docs/index.html)** — three-arm comparison, observer shift, label cost, missing-modality reverse diagnosis, scope refine. Real numbers only; **do not claim Pearson win over ~0.61**.

| What | Where |
|---|---|
| Landing page | [`docs/`](docs/) |
| Bridge code | [`bridge_anm/`](bridge_anm/) |
| Proof reports | [`docs/reports/`](docs/reports/) · also under `outputs/anm_cite_bridge/` |
| Key GIFs | `docs/assets/hard_proof/`, `docs/assets/missing_modality/` |

### Claim highlights (not accuracy)

- **Criterion edit** O0→O1 (stricter abstain, disagree=0) / O0→O2 (expected_action rewrite **12.48%**) — YAML only, no `best.pt` retrain
- **Abstain** follows declaration: ANM 317 → 834 → 3,373; TEDDY-alone silently over-answers; Train abstain stuck at 185 without new labels
- **Attribution** full-n LOO n_attr=16,433 · top-1 flip rate 0.277 · permutation p≈0.0099
- **Scope refine** soft_P gating ΔQ_peak≈0.055; flip-sensitive cells Q=0.8099 vs random 0.9458
- **Missing modality** `adt_only`: TEDDY abstain=0 vs ANM abstain=1,664 (O2 abstain=all) — reverse diagnosis, not Pearson

Claim boundary: local response diagnosis only · holdout ADT = verifier · Train+Jev = log-loss stand-in · not clinical.

### Re-run (no TEDDY retrain)

```bash
cd teddy_mm
PYTHONPATH=/path/to/ANM:. .venv/bin/python bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100
PYTHONPATH=/path/to/ANM:. .venv/bin/python bridge_anm/run_scope_refine_proof.py
# missing modality: see bridge_anm/README.md
```

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
