#!/usr/bin/env python3
"""Phase-2 bidirectional CITE train scaffold (modality dropout + latent FM).

Writes only to --out (default outputs/cite_phase2). Does not touch
data/processed/cite or outputs/cite_phase1.

Defaults (fixed):
  --adt-input-transform auto   the ADT encoder transforms the stored ADT exactly
                               once: 'none' when the pack already holds normalized
                               (non-integer) ADT, 'log1p' for raw integer counts.
  --eval-size-factor train-median
                               val/test decoding uses one constant = median ADT
                               size factor over the training cells used; never the
                               evaluation cells' own measured ADT depth.
  --legacy                     old behaviour: clr on top of the stored ADT and the
                               evaluation cells' measured size factor.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from teddy_mm.bidirectional import (
    ADT_INPUT_TRANSFORMS,
    BidirectionalCite,
    adt_is_count_like,
    resolve_adt_input_transform,
)
from teddy_mm.data import SIZE_FACTOR_MODES, eval_size_factor, load_prepared
from teddy_mm.device import resolve_device

# In adt_only / joint the measured ADT of the scored cell is the encoder input,
# so those scores are reconstructions of the answer key, not predictions.
ANSWER_KEY_IS_ENCODER_INPUT = {"rna_only": False, "adt_only": True, "joint": True}


def _pearson_mean(a: np.ndarray, b: np.ndarray) -> float:
    scores = []
    for j in range(a.shape[1]):
        if a[:, j].std() < 1e-8 or b[:, j].std() < 1e-8:
            continue
        scores.append(np.corrcoef(a[:, j], b[:, j])[0, 1])
    return float(np.mean(scores)) if scores else float("nan")


def split_indices(split, name: str, max_cells: int, rng: np.random.Generator) -> np.ndarray:
    """Global indices of one split; a sorted random subset when max_cells > 0."""
    idx = np.where(np.asarray(split) == name)[0]
    if max_cells and max_cells > 0 and idx.size > max_cells:
        idx = np.sort(rng.choice(idx, size=max_cells, replace=False))
    return idx


def subset(pack, z, idx, sf):
    """Tensors for the cells idx; sf is the size factor handed to the decoder."""
    return {
        "z": torch.from_numpy(z[idx]),
        "adt": torch.from_numpy(pack["adt"][idx]),
        "sf": torch.from_numpy(np.asarray(sf, dtype=np.float32)),
    }


@torch.no_grad()
def eval_split(name, data, model, device, fm_steps):
    """Per-protein mean Pearson per input mode.

    {name}_{mode}_pearson        : one flow-matching sample per cell (use_fm=True);
                                   the key used for checkpoint selection.
    {name}_{mode}_direct_pearson : direct decode of the fused latent (use_fm=False).
    """
    model.eval()
    z = data["z"].to(device)
    y = data["adt"].to(device)
    sf = data["sf"].to(device)
    y_np = y.cpu().numpy()
    out = {}
    for mode in ("rna_only", "adt_only", "joint"):
        pred = model.predict_adt(z, y, sf, mode=mode, fm_steps=fm_steps, use_fm=True)
        out[f"{name}_{mode}_pearson"] = _pearson_mean(y_np, pred.cpu().numpy())
        direct = model.predict_adt(z, y, sf, mode=mode, use_fm=False)
        out[f"{name}_{mode}_direct_pearson"] = _pearson_mean(y_np, direct.cpu().numpy())
    print(out)
    return out


def resolve_processed(explicit: Path | None) -> Path:
    """Prefer phase2_4k, then smoke; never require full cite z_rna."""
    if explicit is not None:
        return explicit
    for cand in (
        ROOT / "data/processed/cite_phase2_4k",
        ROOT / "data/processed/cite_smoke",
    ):
        if (cand / "z_rna.npy").exists() and (cand / "cite_arrays.npz").exists():
            return cand
    raise SystemExit(
        "no ready processed pack with z_rna.npy. "
        "Run: python scripts/02_prepare_cite.py --max-cells 4000 "
        "--out data/processed/cite_phase2_4k && "
        "python scripts/03_embed_rna.py --processed data/processed/cite_phase2_4k "
        "--device cpu --seq-len 256 --batch-size 4\n"
        "Or use smoke: python scripts/00_smoke_synthetic.py"
    )


def main():
    p = argparse.ArgumentParser(
        description="Phase-2 bidirectional CITE (ADT encoder + modality dropout + latent FM)"
    )
    p.add_argument("--processed", type=Path, default=None, help="pack with cite_arrays.npz + z_rna.npy")
    p.add_argument("--out", type=Path, default=ROOT / "outputs/cite_phase2")
    p.add_argument("--epochs", type=int, default=8)
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--hidden", type=int, default=512)
    p.add_argument("--adt-hidden", type=int, default=256)
    p.add_argument("--adt-layers", type=int, default=2)
    p.add_argument("--p-drop-rna", type=float, default=0.18)
    p.add_argument("--p-drop-adt", type=float, default=0.15)
    p.add_argument("--lambda-nll", type=float, default=1.0)
    p.add_argument("--lambda-fm", type=float, default=1.0)
    p.add_argument("--lambda-rna-side", type=float, default=0.5)
    p.add_argument("--fm-steps", type=int, default=20)
    p.add_argument("--device", default="auto")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--select-on", choices=("joint_fm", "joint_direct"), default="joint_fm",
                   help="val metric for checkpoint selection: joint_fm (default, legacy) = one flow-matching sample; "
                   "joint_direct = direct decode (less noisy)")
    p.add_argument(
        "--adt-input-transform",
        choices=("auto",) + ADT_INPUT_TRANSFORMS,
        default=None,
        help="transform applied once inside the ADT encoder (default auto: 'none' for an "
        "already-normalized pack, 'log1p' for integer counts; legacy: clr)",
    )
    p.add_argument(
        "--eval-size-factor",
        choices=SIZE_FACTOR_MODES,
        default=None,
        help="decoder size factor for val/test cells (default train-median over the "
        "training cells used; legacy: measured = each evaluated cell's own ADT depth)",
    )
    p.add_argument(
        "--legacy",
        action="store_true",
        help="reproduce the pre-fix behaviour: --adt-input-transform clr --eval-size-factor measured",
    )
    p.add_argument("--max-train-cells", type=int, default=0, help="0 = all split=='train' cells")
    p.add_argument("--max-val-cells", type=int, default=0, help="0 = all split=='val' cells")
    p.add_argument("--max-test-cells", type=int, default=0, help="0 = all split=='test' cells")
    p.add_argument("--subset-seed", type=int, default=None, help="seed for the cell subsets (default --seed)")
    p.add_argument("--threads", type=int, default=0, help="torch CPU threads (0 = torch default)")
    args = p.parse_args()

    if args.legacy:
        for flag, val, legacy_val in (
            ("--adt-input-transform", args.adt_input_transform, "clr"),
            ("--eval-size-factor", args.eval_size_factor, "measured"),
        ):
            if val is not None and val != legacy_val:
                p.error(f"--legacy conflicts with {flag} {val}")
        args.adt_input_transform, args.eval_size_factor = "clr", "measured"
    args.adt_input_transform = args.adt_input_transform or "auto"
    args.eval_size_factor = args.eval_size_factor or "train-median"
    if args.threads > 0:
        torch.set_num_threads(args.threads)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = resolve_device(args.device)
    processed = resolve_processed(args.processed)
    print("device", device)
    print("processed", processed)

    pack = load_prepared(processed, load_rna=False)
    z_path = processed / "z_rna.npy"
    if not z_path.exists():
        raise SystemExit(f"missing {z_path}")
    z = np.load(z_path).astype(np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)

    split = pack["split"]
    sf_all = np.asarray(pack["adt_size_factor"], dtype=np.float32)
    sub_rng = np.random.default_rng(args.seed if args.subset_seed is None else args.subset_seed)
    train_idx = split_indices(split, "train", args.max_train_cells, sub_rng)
    val_idx = split_indices(split, "val", args.max_val_cells, sub_rng)
    test_idx = split_indices(split, "test", args.max_test_cells, sub_rng)

    adt_count_like = adt_is_count_like(pack["adt"][train_idx])
    adt_transform = resolve_adt_input_transform(args.adt_input_transform, pack["adt"][train_idx])
    print(f"adt_input_transform requested={args.adt_input_transform} resolved={adt_transform} "
          f"(stored ADT count-like={adt_count_like})")
    if adt_transform != "none" and not adt_count_like:
        print("WARNING: stored ADT is already normalized (non-integer); "
              f"'{adt_transform}' transforms it a second time")

    # Training cells: their own measured size factor (their ADT is the training target).
    # Val/test cells: resolved by --eval-size-factor (default: train-median constant).
    val_sf, val_sf_info = eval_size_factor(args.eval_size_factor, sf_all, train_idx, val_idx)
    test_sf, test_sf_info = eval_size_factor(args.eval_size_factor, sf_all, train_idx, test_idx)
    train_median_sf = float(np.median(sf_all[train_idx].astype(np.float64)))
    print(f"eval size factor: {test_sf_info}; train-median over {train_idx.size} training cells = {train_median_sf}")

    train = subset(pack, z, train_idx, sf_all[train_idx])
    val = subset(pack, z, val_idx, val_sf)
    test = subset(pack, z, test_idx, test_sf)
    z_dim = int(z.shape[1])
    n_adt = int(pack["adt"].shape[1])
    print(f"train={len(train['z'])} val={len(val['z'])} test={len(test['z'])} z={z_dim} adt={n_adt}")

    model = BidirectionalCite(
        z_dim,
        n_adt,
        hidden=args.hidden,
        adt_hidden=args.adt_hidden,
        adt_layers=args.adt_layers,
        p_drop_rna=args.p_drop_rna,
        p_drop_adt=args.p_drop_adt,
        adt_input_transform=adt_transform,
    ).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr)

    n_train = len(train["z"])
    if n_train == 0:
        raise SystemExit("train split empty")
    batch_size = min(args.batch_size, n_train)
    drop_last = n_train >= batch_size * 2
    loader = DataLoader(
        TensorDataset(train["z"], train["adt"], train["sf"]),
        batch_size=batch_size,
        shuffle=True,
        drop_last=drop_last,
    )

    args.out.mkdir(parents=True, exist_ok=True)
    history = []
    best = -1e9
    for epoch in range(1, args.epochs + 1):
        model.train()
        losses = []
        for zb, yb, sfb in loader:
            zb, yb, sfb = zb.to(device), yb.to(device), sfb.to(device)
            out = model.forward_train(
                zb,
                yb,
                sfb,
                fm_steps=max(args.fm_steps // 2, 4),
                lambda_nll=args.lambda_nll,
                lambda_fm=args.lambda_fm,
                lambda_rna_side=args.lambda_rna_side,
            )
            opt.zero_grad()
            out["loss"].backward()
            opt.step()
            losses.append(float(out["loss"].detach()))
        metrics = {
            "epoch": epoch,
            "train_loss": float(np.mean(losses)),
            "frac_rna_kept": float(out["frac_rna_kept"]),
            "frac_adt_kept": float(out["frac_adt_kept"]),
        }
        metrics.update(eval_split("val", val, model, device, args.fm_steps))
        history.append(metrics)
        sel_key = "val_joint_pearson" if args.select_on == "joint_fm" else "val_joint_direct_pearson"
        score = metrics.get(sel_key, float("nan"))
        if score == score and score > best:  # not NaN
            best = score
            torch.save(
                {
                    "model": model.state_dict(), "metrics": metrics, "z_dim": z_dim, "n_adt": n_adt,
                    "config": model.config, "train_median_size_factor": train_median_sf,
                    "n_train_cells": int(train_idx.size), "eval_size_factor": test_sf_info,
                },
                args.out / "best.pt",
            )
            print(f"epoch {epoch} saved best {sel_key}={score:.4f}")

    ckpt_path = args.out / "best.pt"
    if ckpt_path.exists():
        ckpt = torch.load(ckpt_path, map_location=device, weights_only=False)
        model.load_state_dict(ckpt["model"])
    test_metrics = eval_split("test", test, model, device, args.fm_steps)
    run_info = {
        "adt_input_transform_requested": args.adt_input_transform,
        "adt_input_transform": adt_transform,
        "stored_adt_count_like": adt_count_like,
        "legacy": bool(args.legacy),
        "train_size_factor": "measured per training cell (training target)",
        "val_size_factor": val_sf_info,
        "test_size_factor": test_sf_info,
        "train_median_size_factor": train_median_sf,
        "n_train_cells": int(train_idx.size),
        "n_val_cells": int(val_idx.size),
        "n_test_cells": int(test_idx.size),
        "epochs": args.epochs,
        "seed": args.seed,
        "answer_key_is_encoder_input": ANSWER_KEY_IS_ENCODER_INPUT,
    }
    (args.out / "metrics.json").write_text(
        json.dumps({"history": history, "test": test_metrics, "processed": str(processed), "run": run_info}, indent=2)
    )
    print("wrote", args.out / "metrics.json")


if __name__ == "__main__":
    main()
