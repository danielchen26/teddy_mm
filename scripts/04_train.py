#!/usr/bin/env python3
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

from teddy_mm.data import load_prepared
from teddy_mm.device import resolve_device
from teddy_mm.models import AdtDecoder, MLP, VelocityNet, integrate_fm, nb_nll


def _r2(y, yhat) -> float:
    ss_res = ((y - yhat) ** 2).sum()
    ss_tot = ((y - y.mean(axis=0)) ** 2).sum()
    return float(1.0 - ss_res / max(ss_tot, 1e-8))


def subset(pack, z, name):
    m = pack["split"] == name
    return {
        "z": torch.from_numpy(z[m]),
        "adt": torch.from_numpy(pack["adt"][m]),
        "sf": torch.from_numpy(pack["adt_size_factor"][m]),
    }


@torch.no_grad()
def eval_split(name, data, mlp, v_net, dec, device, fm_steps):
    mlp.eval(); v_net.eval(); dec.eval()
    z = data["z"].to(device)
    y = data["adt"].to(device)
    sf = data["sf"].to(device)
    pred_mlp = dec(mlp(z), sf)[0]
    z1 = integrate_fm(v_net, torch.randn_like(z), z, fm_steps)
    pred_fm = dec(z1, sf)[0]
    y_np, mlp_np, fm_np = y.cpu().numpy(), pred_mlp.cpu().numpy(), pred_fm.cpu().numpy()
    def pearson_mean(a, b):
        scores = []
        for j in range(a.shape[1]):
            if a[:, j].std() < 1e-8 or b[:, j].std() < 1e-8:
                continue
            scores.append(np.corrcoef(a[:, j], b[:, j])[0, 1])
        return float(np.mean(scores)) if scores else float("nan")
    out = {
        f"{name}_mlp_pearson": pearson_mean(y_np, mlp_np),
        f"{name}_fm_pearson": pearson_mean(y_np, fm_np),
        f"{name}_mlp_r2": _r2(y_np, mlp_np),
        f"{name}_fm_r2": _r2(y_np, fm_np),
    }
    print(out)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--out", type=Path, default=ROOT / "outputs/cite_phase1")
    p.add_argument("--epochs", type=int, default=40)
    p.add_argument("--batch-size", type=int, default=512)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--hidden", type=int, default=512)
    p.add_argument("--lambda-nll", type=float, default=1.0)
    p.add_argument("--lambda-fm", type=float, default=1.0)
    p.add_argument("--fm-steps", type=int, default=20)
    p.add_argument("--device", default="auto")
    args = p.parse_args()

    device = resolve_device(args.device)
    print("device", device)
    pack = load_prepared(args.processed)
    z_path = args.processed / "z_rna.npy"
    if not z_path.exists():
        raise SystemExit(f"missing {z_path}; run scripts/03_embed_rna.py first")
    z = np.load(z_path).astype(np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)

    train, val, test = subset(pack, z, "train"), subset(pack, z, "val"), subset(pack, z, "test")
    z_dim = z.shape[1]
    n_adt = pack["adt"].shape[1]
    print(f"train={len(train['z'])} val={len(val['z'])} test={len(test['z'])} z={z_dim} adt={n_adt}")

    mlp = MLP(z_dim, z_dim, hidden=args.hidden).to(device)
    v_net = VelocityNet(z_dim, cond_dim=z_dim, hidden=args.hidden).to(device)
    dec = AdtDecoder(z_dim, n_adt, hidden=args.hidden).to(device)
    params = list(mlp.parameters()) + list(v_net.parameters()) + list(dec.parameters())
    opt = torch.optim.AdamW(params, lr=args.lr)

    n_train = len(train["z"])
    if n_train == 0:
        raise SystemExit("train split is empty")
    batch_size = min(args.batch_size, n_train)
    drop_last = n_train >= batch_size * 2  # keep drop_last only when >=2 full batches
    if batch_size != args.batch_size:
        print(f"capping batch_size {args.batch_size} -> {batch_size} (n_train={n_train})")
    loader = DataLoader(
        TensorDataset(train["z"], train["adt"], train["sf"]),
        batch_size=batch_size,
        shuffle=True,
        drop_last=drop_last,
    )
    if len(loader) == 0:
        raise SystemExit(f"DataLoader empty: n_train={n_train} batch_size={batch_size} drop_last={drop_last}")

    history = []
    best = -1e9
    args.out.mkdir(parents=True, exist_ok=True)
    for epoch in range(1, args.epochs + 1):
        mlp.train(); v_net.train(); dec.train()
        losses = []
        for zb, yb, sfb in loader:
            zb, yb, sfb = zb.to(device), yb.to(device), sfb.to(device)
            z_mlp = mlp(zb)
            mu_mlp, th = dec(z_mlp, sfb)
            loss_mlp = nb_nll(yb, mu_mlp, th)

            x0 = torch.randn_like(zb)
            t = torch.rand(zb.size(0), device=device)
            # target latent: map RNA embedding through MLP (paired, no OT)
            x1 = z_mlp.detach()
            xt = (1.0 - t[:, None]) * x0 + t[:, None] * x1
            loss_fm = ((v_net(xt, t, zb) - (x1 - x0)) ** 2).mean()
            z1 = integrate_fm(v_net, x0, zb, max(args.fm_steps // 2, 4))
            mu_fm, th_fm = dec(z1, sfb)
            loss_fm_nll = nb_nll(yb, mu_fm, th_fm)

            loss = args.lambda_nll * (loss_mlp + 0.5 * loss_fm_nll) + args.lambda_fm * loss_fm
            opt.zero_grad()
            loss.backward()
            opt.step()
            losses.append(float(loss.detach()))

        metrics = {"epoch": epoch, "train_loss": float(np.mean(losses))}
        metrics.update(eval_split("val", val, mlp, v_net, dec, device, args.fm_steps))
        history.append(metrics)
        score = metrics["val_fm_pearson"]
        if score > best:
            best = score
            torch.save(
                {"mlp": mlp.state_dict(), "v_net": v_net.state_dict(), "dec": dec.state_dict(), "metrics": metrics},
                args.out / "best.pt",
            )
            print(f"epoch {epoch} saved best val_fm_pearson={score:.4f}")

    ckpt = torch.load(args.out / "best.pt", map_location=device, weights_only=False)
    mlp.load_state_dict(ckpt["mlp"])
    v_net.load_state_dict(ckpt["v_net"])
    dec.load_state_dict(ckpt["dec"])
    test_metrics = eval_split("test", test, mlp, v_net, dec, device, args.fm_steps)
    (args.out / "metrics.json").write_text(json.dumps({"history": history, "test": test_metrics}, indent=2))
    print("wrote", args.out / "metrics.json")


if __name__ == "__main__":
    main()
