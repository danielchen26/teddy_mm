#!/usr/bin/env python3
"""Re-score a saved phase-2 checkpoint on the site4 test cells, two ways.

The phase-2 training script (and the Block 2 export) score predict_adt(..., use_fm=True): one random
flow-matching sample per cell. The model also has a direct decode from the fused latent
(use_fm=False), which is the path its main loss trains. This prints the per-protein mean Pearson for
both, per input mode, next to the recorded metrics.

Size factor (--size-factor): train-median (default) hands the decoder one constant, the median ADT
size factor over the training cells (the value recorded in the checkpoint when present, else the
median over split=='train'); measured (legacy) hands it each test cell's own measured ADT depth.
ADT input transform: taken from the checkpoint; checkpoints without one were trained with the legacy
'clr' encoder and are rebuilt with it (--adt-input-transform overrides).

In adt_only / joint the test cell's measured ADT is the encoder input and also the answer key, so those
rows are reconstructions, not predictions.

Usage (from the teddy_mm repo root; needs torch):
    python scripts/phase2_decode_check.py
    python scripts/phase2_decode_check.py --size-factor measured   # legacy numbers
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from teddy_mm.bidirectional import ADT_INPUT_TRANSFORMS, load_bidirectional  # noqa: E402
from teddy_mm.data import SIZE_FACTOR_MODES, eval_size_factor, load_prepared  # noqa: E402
from teddy_mm.device import resolve_device  # noqa: E402


def pearson_mean(a: np.ndarray, b: np.ndarray) -> float:
    s = [np.corrcoef(a[:, j], b[:, j])[0, 1] for j in range(a.shape[1]) if a[:, j].std() > 1e-8 and b[:, j].std() > 1e-8]
    return float(np.mean(s)) if s else float("nan")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/cite_phase2/best.pt")
    p.add_argument("--out", type=Path, default=None, help="default: <ckpt dir>/decode_check.json")
    p.add_argument("--size-factor", choices=SIZE_FACTOR_MODES, default="train-median")
    p.add_argument("--adt-input-transform", choices=ADT_INPUT_TRANSFORMS, default=None,
                   help="override the transform recorded in the checkpoint")
    p.add_argument("--device", default="auto")
    p.add_argument("--fm-steps", type=int, default=20)
    args = p.parse_args()

    dev = resolve_device(args.device)
    pack = load_prepared(args.processed, load_rna=False)
    z = np.load(args.processed / "z_rna.npy").astype(np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)  # same L2 step as scripts/06_train_bidirectional.py
    split = pack["split"]
    test_idx = np.where(split == "test")[0]
    train_idx = np.where(split == "train")[0]
    ck = torch.load(args.ckpt, map_location=dev, weights_only=False)
    model, info = load_bidirectional(ck, dev, adt_input_transform=args.adt_input_transform)
    sf_all = np.asarray(pack["adt_size_factor"], dtype=np.float32)
    sf_np, sf_info = eval_size_factor(args.size_factor, sf_all, train_idx, test_idx)
    if args.size_factor == "train-median" and info["train_median_size_factor"] is not None:
        sf_np = np.full(test_idx.size, info["train_median_size_factor"], dtype=np.float32)
        sf_info = {"mode": "train-median", "value": float(info["train_median_size_factor"]), "source": "checkpoint"}
    zt = torch.from_numpy(z[test_idx]).to(dev)
    y = torch.from_numpy(pack["adt"][test_idx]).to(dev)
    sf = torch.from_numpy(sf_np).to(dev)
    Y = y.cpu().numpy()
    metrics_path = args.ckpt.parent / "metrics.json"
    recorded = json.loads(metrics_path.read_text())["test"] if metrics_path.exists() else None
    out = {"test_cells": int(test_idx.size), "recorded": recorded, "size_factor": sf_info, "model": info,
           "answer_key_is_encoder_input": {"rna_only": False, "adt_only": True, "joint": True}, "by_mode": {}}
    with torch.no_grad():
        for mode in ("rna_only", "adt_only", "joint"):
            torch.manual_seed(0)
            one = model.predict_adt(zt, y, sf, mode=mode, fm_steps=args.fm_steps, use_fm=True).cpu().numpy()
            five = np.mean([model.predict_adt(zt, y, sf, mode=mode, fm_steps=args.fm_steps, use_fm=True).cpu().numpy()
                            for _ in range(5)], axis=0)
            direct = model.predict_adt(zt, y, sf, mode=mode, use_fm=False).cpu().numpy()
            out["by_mode"][mode] = {"fm_one_sample": round(pearson_mean(Y, one), 4), "fm_mean_of_5": round(pearson_mean(Y, five), 4),
                                    "direct_decode": round(pearson_mean(Y, direct), 4)}
    out_path = args.out or (args.ckpt.parent / "decode_check.json")
    out_path.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
