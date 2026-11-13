#!/usr/bin/env python3
"""Repair the phase-1 head where the ANM audit locates the loss: its first network (the MLP).

The audit (docs/ANM_DIAGNOSIS.md) finds that the head loses the measured protein response to a knockout between TEDDY's
embedding and the head's latent layer. This script makes the repair that was audited in the loop (docs/anm-loop.html):

  - only the MLP changes; the decoder stays exactly as trained;
  - it is fitted on perturbation cells (knockout and control) whose proteins were measured, to the measured protein,
    rescaled to the mean and spread of the head's own output on the control cells;
  - the head's other outputs on those cells are held where they are;
  - the head's outputs on bone-marrow training cells are held where they are (`--hold-weight`), so its original task is kept.

Two forms of it were audited on knockouts of genes held out of the fit, each against a rule written before it ran:

  one data set, PD-L1 only (anm_diag/runs/loop_guard_20261112)
      Passed its rule on both cell types: the head adds no loss for PD-L1, and its bone-marrow score stays within 0.02.
      The head is for the cell type it was repaired on; a head repaired on one was still found to lose the response on
      the other (anm_diag/runs/loop_20261112).
  both data sets and every measured protein in one fit (anm_diag/runs/repair_20261112)
      Did not pass its rule on melanoma. Not recommended as it stands.

Perturbation cells are given as folders of exported features, one per data set, as written by the anm_diag adapters
(`dense.npz` with A, target, guide, ids, readout__<ID>, and `node__Z.npy` with TEDDY's embedding). The form that passed:

  python scripts/07_repair_head.py \\
      --ckpt outputs/cite_phase1/best.pt --processed data/processed/cite \\
      --data frangieh2021=<features>/frangieh2021:PDL1=CD274 \\
      --out outputs/cite_phase1_repaired_melanoma

The joint form repeats `--data`, with more readouts mapped to head outputs:

      --data frangieh2021=<features>/frangieh2021:PDL1=CD274,MHC1=HLA-A-B-C,HLAE=HLA-E \\
      --data papalexi2021=<features>/papalexi2021:PDL1=CD274

The checkpoint also holds a flow-matching branch (`v_net`), which was trained to reproduce the latent of the MLP as it was.
It is copied unchanged, so only the direct route (MLP, then decoder) is repaired; that route is the one the audit tests.

`--hold-out GENE` leaves one target gene out of the fit (for an audit on it). Half of the control cells are never fitted
(`--seed` draws them, the same way as the audited runs), so they stay available for an audit. The fit is full batch and
has no randomness beyond that split. Writes <out>/best.pt in the format of scripts/04_train.py and <out>/repair_log.json.
A head fitted with no gene held out has not itself been audited: the evidence for it is the audits of the held-out fits.
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

from teddy_mm.data import load_prepared  # noqa: E402
from teddy_mm.models import MLP, AdtDecoder  # noqa: E402


def halves(groups, rng, share=0.5):
    first = np.zeros(len(groups), bool)
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        first[idx[: int(round(share * len(idx)))]] = True
    return first


def unit(z):
    z = np.asarray(z, np.float32)
    return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)


def head_out(mlp, dec, x):
    """The head's outputs as the audit reads them: log of the predicted mean, before the size factor."""
    return torch.log(torch.nn.functional.softplus(dec.mu(mlp(x))).clamp(min=1e-6))


def load_cells(folder: Path, seed: int):
    d = np.load(folder / "dense.npz")
    D = {k: d[k] for k in d.files if not k.startswith("node__")}
    Z = d["node__Z"] if "node__Z" in d.files else np.load(folder / "node__Z.npy")
    A, tg = D["A"], D["target"].astype(str)
    st = D["stratum"].astype(str) if "stratum" in D else np.full(len(A), "all")
    rng = np.random.default_rng(seed)
    ctl = A == 0
    audit_ctl = ctl & halves(np.where(ctl, st, "treated"), rng)                       # never fitted
    grp = np.where(A == 1, np.char.add(np.char.add(D["guide"].astype(str), "|"), st), np.char.add("control|", st))
    val = halves(grp, rng, 0.2) & ~audit_ctl                                           # set aside for the stopping rule
    return D, unit(Z), A, tg, ctl & ~audit_ctl, val, audit_ctl


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ckpt", type=Path, required=True, help="phase-1 checkpoint written by scripts/04_train.py")
    p.add_argument("--processed", type=Path, default=ROOT / "data/processed/cite", help="bone-marrow data of the checkpoint (z_rna.npy and the prepared arrays)")
    p.add_argument("--data", action="append", required=True, metavar="NAME=FOLDER:READOUT=PROTEIN,...", help="perturbation cells and which measured readout trains which head output")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--hold-out", default=None, help="a target gene to leave out of the fit")
    p.add_argument("--hold-weight", type=float, default=1.0, help="weight of the term that holds the bone-marrow outputs")
    p.add_argument("--n-hold", type=int, default=3000, help="bone-marrow training cells in that term")
    p.add_argument("--keep-weight", type=float, default=1.0, help="weight of the term that holds the untrained outputs on the perturbation cells")
    p.add_argument("--lr", type=float, default=3e-4)
    p.add_argument("--steps", type=int, default=400)
    p.add_argument("--every", type=int, default=20, help="steps between checks on the set-aside cells")
    p.add_argument("--patience", type=int, default=3)
    p.add_argument("--seed", type=int, default=0, help="draws the control cells that are never fitted and the set-aside cells")
    args = p.parse_args()
    torch.manual_seed(args.seed)
    torch.set_num_threads(max(torch.get_num_threads(), 1))

    ck = torch.load(args.ckpt, map_location="cpu", weights_only=False)
    pack = load_prepared(args.processed, load_rna=False)
    names = [str(n) for n in pack["adt_names"]]
    z_dim, hidden = ck["mlp"]["net.0.weight"].shape[1], ck["mlp"]["net.0.weight"].shape[0]
    mlp, dec = MLP(z_dim, z_dim, hidden=hidden), AdtDecoder(z_dim, len(names), hidden=hidden)
    mlp.load_state_dict(ck["mlp"]); dec.load_state_dict(ck["dec"])
    mlp.eval(); dec.eval()                                                             # no dropout: the fit is deterministic
    for q in dec.parameters():
        q.requires_grad_(False)
    dec0 = {k: v.clone() for k, v in dec.state_dict().items()}

    # bone-marrow cells whose outputs are held: a seeded sample of the training split, a fifth of it set aside
    zb = np.load(args.processed / "z_rna.npy", mmap_mode="r")
    pick = np.sort(np.random.default_rng(0).choice(np.where(pack["split"] == "train")[0], args.n_hold, replace=False))
    xb = torch.from_numpy(unit(zb[pick]))
    bval = np.zeros(len(pick), bool); bval[np.random.default_rng(args.seed).permutation(len(pick))[: len(pick) // 5]] = True
    with torch.no_grad():
        ob0 = head_out(mlp, dec, xb)

    # perturbation cells: targets, masks and per-entry weights
    X, T, M, O0, DSI, TRAIN, VAL, log_data = [], [], [], [], [], [], [], {}
    for i, spec in enumerate(args.data):
        name, rest = spec.split("=", 1); folder, mapping = rest.rsplit(":", 1)
        D, x, A, tg, train_ctl, val, audit_ctl = load_cells(Path(folder), args.seed)
        xt = torch.from_numpy(x)
        with torch.no_grad():
            out0 = head_out(mlp, dec, xt).numpy()
        t, m = out0.copy(), np.zeros(out0.shape, bool)
        for pair in mapping.split(","):
            rid, prot = pair.split("="); j = names.index(prot); y = D["readout__" + rid].astype(np.float32); c = train_ctl
            t[:, j] = out0[c, j].mean() + out0[c, j].std() * (y - y[c].mean()) / y[c].std(); m[:, j] = True
        train = train_ctl | ((A == 1) & (tg != args.hold_out))
        X.append(x); T.append(t); M.append(m); O0.append(out0); DSI.append(np.full(len(A), i)); TRAIN.append(train); VAL.append(val)
        log_data[name] = {"cells": int(len(A)), "fitted": int((train & ~val).sum()), "set_aside": int((train & val).sum()), "controls_never_fitted": int(audit_ctl.sum()), "trained_outputs": mapping}
    X, T, M, O0, DSI, TRAIN, VAL = (np.concatenate(v) for v in (X, T, M, O0, DSI, TRAIN, VAL))

    def tensors(rows, brows):
        m, ds = M[rows], DSI[rows]
        wt, wk = np.zeros(m.shape, np.float32), np.zeros(m.shape, np.float32)
        for d in np.unique(ds):
            r = ds == d
            wt[r] = m[r] / r.sum()                                                    # one mean squared error per (data set, protein), summed
            wk[r] = (~m[r]) * (args.keep_weight / (r.sum() * max((~m[r][0]).sum(), 1)))
        return tuple(torch.from_numpy(a) for a in (X[rows], T[rows], O0[rows], wt, wk)) + (xb[brows], ob0[brows])

    def loss(batch):
        x, t, o0, wt, wk, xh, oh = batch
        out = head_out(mlp, dec, x)
        return (wt * (out - t) ** 2).sum() + (wk * (out - o0) ** 2).sum() + args.hold_weight * ((head_out(mlp, dec, xh) - oh) ** 2).mean()

    fit_b, val_b = tensors(TRAIN & ~VAL, ~bval), tensors(TRAIN & VAL, bval)
    opt = torch.optim.Adam(mlp.parameters(), lr=args.lr, betas=(0.9, 0.999), eps=1e-8)
    best, best_state, bad, curve = float("inf"), None, 0, []
    for step in range(1, args.steps + 1):
        opt.zero_grad()
        L = loss(fit_b)
        L.backward()
        opt.step()
        if step % args.every == 0:
            with torch.no_grad():
                Lv = float(loss(val_b))
            curve.append([step, float(L), Lv])
            print(f"step {step}: loss {float(L):.4f}, set-aside {Lv:.4f}", flush=True)
            if Lv < best:
                best, best_state, bad = Lv, {k: v.clone() for k, v in mlp.state_dict().items()}, 0
            else:
                bad += 1
                if bad >= args.patience:
                    break
    mlp.load_state_dict(best_state)
    assert all(torch.equal(dec.state_dict()[k], v) for k, v in dec0.items())           # the decoder is exactly the head as trained
    with torch.no_grad():
        moved = float(((head_out(mlp, dec, xb) - ob0) ** 2).mean().sqrt())
    args.out.mkdir(parents=True, exist_ok=True)
    new = dict(ck); new["mlp"] = best_state
    new["repair"] = {"from": args.ckpt.name, "held_out": args.hold_out, "hold_weight": args.hold_weight, "seed": args.seed}
    torch.save(new, args.out / "best.pt")
    json.dump({"data": log_data, "held_out": args.hold_out, "hold_weight": args.hold_weight, "keep_weight": args.keep_weight, "lr": args.lr, "max_steps": args.steps,
               "kept_step": int(min(curve, key=lambda c: c[2])[0]), "curve": curve, "bone_marrow_outputs_moved_rms": moved, "bone_marrow_cells_held": int((~bval).sum())},
              open(args.out / "repair_log.json", "w"), indent=1)
    print(f"kept step {min(curve, key=lambda c: c[2])[0]}; bone-marrow outputs moved by {moved:.4f} (root mean square); wrote {args.out / 'best.pt'}")
    print("score it on the bone-marrow test split with scripts/05_eval.py, and audit it with anm_diag before using it")


if __name__ == "__main__":
    main()
