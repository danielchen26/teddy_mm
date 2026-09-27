#!/usr/bin/env python3
"""Phase-2 bidirectional CITE train scaffold (modality dropout + latent FM).

Writes only to --out (default outputs/cite_phase2). Does not touch
data/processed/cite or outputs/cite_phase1.
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

from teddy_mm.bidirectional import BidirectionalCite
from teddy_mm.data import load_prepared
from teddy_mm.device import resolve_device


def _pearson_mean(a: np.ndarray, b: np.ndarray) -> float:
    scores = []
    for j in range(a.shape[1]):
        if a[:, j].std() < 1e-8 or b[:, j].std() < 1e-8:
            continue
        scores.append(np.corrcoef(a[:, j], b[:, j])[0, 1])
    return float(np.mean(scores)) if scores else float("nan")


def subset(pack, z, name):
    m = pack["split"] == name
    return {
        "z": torch.from_numpy(z[m]),
        "adt": torch.from_numpy(pack["adt"][m]),
        "sf": torch.from_numpy(pack["adt_size_factor"][m]),
    }


@torch.no_grad()
def eval_split(name, data, model, device, fm_steps):
    model.eval()
    z = data["z"].to(device)
    y = data["adt"].to(device)
    sf = data["sf"].to(device)
    out = {}
    for mode in ("rna_only", "adt_only", "joint"):
        pred = model.predict_adt(z, y, sf, mode=mode, fm_steps=fm_steps, use_fm=True)
        out[f"{name}_{mode}_pearson"] = _pearson_mean(y.cpu().numpy(), pred.cpu().numpy())
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
    args = p.parse_args()

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = resolve_device(args.device)
    processed = resolve_processed(args.processed)
    print("device", device)
    print("processed", processed)

    pack = load_prepared(processed)
    z_path = processed / "z_rna.npy"
    if not z_path.exists():
        raise SystemExit(f"missing {z_path}")
    z = np.load(z_path).astype(np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)

    train, val, test = subset(pack, z, "train"), subset(pack, z, "val"), subset(pack, z, "test")
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
        score = metrics.get("val_joint_pearson", float("nan"))
        if score == score and score > best:  # not NaN
            best = score
            torch.save({"model": model.state_dict(), "metrics": metrics, "z_dim": z_dim, "n_adt": n_adt}, args.out / "best.pt")
            print(f"epoch {epoch} saved best val_joint_pearson={score:.4f}")

    ckpt_path = args.out / "best.pt"
    if ckpt_path.exists():
        ckpt = torch.load(ckpt_path, map_location=device, weights_only=False)
        model.load_state_dict(ckpt["model"])
    test_metrics = eval_split("test", test, model, device, args.fm_steps)
    (args.out / "metrics.json").write_text(
        json.dumps({"history": history, "test": test_metrics, "processed": str(processed)}, indent=2)
    )
    print("wrote", args.out / "metrics.json")


if __name__ == "__main__":
    main()
