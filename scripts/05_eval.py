#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from teddy_mm.data import load_prepared
from teddy_mm.device import resolve_device
from teddy_mm.models import AdtDecoder, MLP, VelocityNet, integrate_fm


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite")
    p.add_argument("--ckpt", type=Path, default=ROOT / "outputs/cite_phase1/best.pt")
    p.add_argument("--out", type=Path, default=ROOT / "outputs/cite_phase1")
    p.add_argument("--hidden", type=int, default=512)
    p.add_argument("--fm-steps", type=int, default=20)
    p.add_argument("--device", default="auto")
    args = p.parse_args()

    device = resolve_device(args.device)
    pack = load_prepared(args.processed)
    z = np.load(args.processed / "z_rna.npy").astype(np.float32)
    z = z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)
    mask = pack["split"] == "test"
    zt = torch.from_numpy(z[mask]).to(device)
    y = pack["adt"][mask]
    sf = torch.from_numpy(pack["adt_size_factor"][mask]).to(device)

    z_dim, n_adt = z.shape[1], y.shape[1]
    mlp = MLP(z_dim, z_dim, hidden=args.hidden).to(device)
    v_net = VelocityNet(z_dim, cond_dim=z_dim, hidden=args.hidden).to(device)
    dec = AdtDecoder(z_dim, n_adt, hidden=args.hidden).to(device)
    blob = torch.load(args.ckpt, map_location=device, weights_only=False)
    mlp.load_state_dict(blob["mlp"])
    v_net.load_state_dict(blob["v_net"])
    dec.load_state_dict(blob["dec"])
    mlp.eval(); v_net.eval(); dec.eval()

    with torch.no_grad():
        pred_mlp = dec(mlp(zt), sf)[0].cpu().numpy()
        pred_fm = dec(integrate_fm(v_net, torch.randn_like(zt), zt, args.fm_steps), sf)[0].cpu().numpy()

    names = pack["adt_names"]
    rows = []
    for j, name in enumerate(names):
        if y[:, j].std() < 1e-8:
            continue
        rows.append(
            {
                "protein": str(name),
                "pearson_mlp": float(np.corrcoef(y[:, j], pred_mlp[:, j])[0, 1]),
                "pearson_fm": float(np.corrcoef(y[:, j], pred_fm[:, j])[0, 1]),
            }
        )
    rows.sort(key=lambda r: -r["pearson_fm"])
    out = args.out / "test_per_protein.json"
    out.write_text(json.dumps(rows, indent=2))
    print(f"wrote {out}  n_proteins={len(rows)}")
    print("top 8 FM proteins:")
    for r in rows[:8]:
        print(f"  {r['protein']:20s}  mlp={r['pearson_mlp']:.3f}  fm={r['pearson_fm']:.3f}")


if __name__ == "__main__":
    main()
