#!/usr/bin/env python3
"""Fit one head of the round.

The loop is that of scripts/04_train.py, line for line: the same models built in the same order, the same optimizer, the
same batches, the same loss, the same choice of checkpoint. What an arm changes is which cells are fitted on, and for
the `function` arm the width and the number of epochs. With all bone-marrow training cells and nothing else it gives
the numbers of scripts/04_train.py for the same seed and machine (tests/test_round1_train_arm.py).

Two things are added, as DESIGN.md says:

  cells of a second study  Their unmeasured proteins are left out of the loss: the registered loss per protein, averaged
                           over the proteins a cell has, then over cells. Their size factor is the cell's total over its
                           measured proteins divided by the median of that total over the arm's cells of that study.
  stopping and resuming    After every epoch the whole state is written (weights, optimizer, every random generator).
                           A fit that is stopped and resumed gives the same weights as one that is not.

  python train_arm.py --head s0_coverage_M1 --bm-arrays cite_arrays.npz --bm-z z_rna.npy \\
      --pbmc-arrays ext_cite_arrays.npz --pbmc-z ext_z_rna.npy --out OUT [--time-budget-sec 6000] [--threads 2]

Exit code 75: the time budget ended before the fit did; run the same command again.
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from teddy_mm.models import AdtDecoder, MLP, VelocityNet, integrate_fm, nb_nll  # noqa: E402

EXIT_BUDGET = 75
BATCH, LR, FM_STEPS, LAMBDA_NLL, LAMBDA_FM = 512, 1e-3, 20, 1.0, 1.0


# ------------------------------------------------------------------ as in scripts/04_train.py


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def unit_length(z):
    z = np.asarray(z).astype(np.float32)
    return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)


def _r2(y, yhat) -> float:
    ss_res = ((y - yhat) ** 2).sum()
    ss_tot = ((y - y.mean(axis=0)) ** 2).sum()
    return float(1.0 - ss_res / max(ss_tot, 1e-8))


def pearson_mean(a, b):
    scores = []
    for j in range(a.shape[1]):
        if a[:, j].std() < 1e-8 or b[:, j].std() < 1e-8:
            continue
        scores.append(np.corrcoef(a[:, j], b[:, j])[0, 1])
    return float(np.mean(scores)) if scores else float("nan")


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
    return {
        f"{name}_mlp_pearson": pearson_mean(y_np, mlp_np),
        f"{name}_fm_pearson": pearson_mean(y_np, fm_np),
        f"{name}_mlp_r2": _r2(y_np, mlp_np),
        f"{name}_fm_r2": _r2(y_np, fm_np),
    }


# ------------------------------------------------------------------ what the round adds


def nb_nll_cells(x, mu, theta, m):
    """The terms of teddy_mm.models.nb_nll, averaged over the measured proteins of each cell, then over cells."""
    eps = 1e-6
    mu = mu.clamp(min=eps)
    theta = theta.clamp(min=eps)
    t1 = torch.lgamma(theta + x) - torch.lgamma(theta) - torch.lgamma(x + 1.0)
    t2 = theta * torch.log(theta / (theta + mu))
    t3 = x * torch.log(mu / (theta + mu))
    return ((-(t1 + t2 + t3) * m).sum(dim=1) / m.sum(dim=1)).mean()


def load_data(bm_arrays, bm_z, pbmc_arrays=None, pbmc_z=None):
    a = np.load(bm_arrays, allow_pickle=False)
    bm = {"z": unit_length(np.load(bm_z)), "adt": a["adt"], "sf": a["adt_size_factor"], "split": a["split"].astype(str),
          "donors": a["donors"].astype(str), "names": a["adt_names"].astype(str)}
    if len(bm["z"]) != len(bm["adt"]):
        raise SystemExit(f"bone marrow: {len(bm['z'])} states for {len(bm['adt'])} cells")
    pb = None
    if pbmc_arrays is not None:
        e = np.load(pbmc_arrays, allow_pickle=False)
        if [str(x) for x in e["adt_names"]] != [str(x) for x in bm["names"]]:
            raise SystemExit("the two studies do not list the same proteins in the same order")
        pb = {"z": unit_length(np.load(pbmc_z)), "adt": e["adt"], "donors": e["donors"].astype(str),
              "barcodes": e["barcodes"].astype(str) if "barcodes" in e.files else np.arange(len(e["adt"])).astype(str)}
        if len(pb["z"]) != len(pb["adt"]):
            raise SystemExit(f"second study: {len(pb['z'])} states for {len(pb['adt'])} cells")
    return bm, pb


def build_train(cfg, bm, pb):
    """The cells an arm is fitted on. Returns tensors (z, adt, size factor, mask or None) and a record of the cells.

    Bone marrow: the first cells of one random order of the training cells per seed, so arms are nested; all of them
    when cfg["n_bm"] is None, in the order of the pack (the order scripts/04_train.py reads them in).
    Second study: for every fit donor its first n cells of one random order per seed and donor.
    """
    seed = int(cfg["seed"])
    tr = np.where(bm["split"] == "train")[0]
    n_bm = len(tr) if cfg.get("n_bm") is None else int(cfg["n_bm"])
    if n_bm > len(tr):
        raise SystemExit(f"{n_bm} bone-marrow cells asked for, {len(tr)} training cells in the pack")
    sel = tr if n_bm == len(tr) else np.sort(tr[np.random.default_rng([seed, 1]).permutation(len(tr))[:n_bm]])
    z, y, sf = bm["z"][sel], bm["adt"][sel], bm["sf"][sel]
    rec = {"bm_index": sel, "pbmc_index": np.zeros(0, np.int64), "pbmc_row_from": np.zeros(0, np.int64)}
    if not cfg.get("pbmc"):
        return (torch.from_numpy(z), torch.from_numpy(y), torch.from_numpy(sf), None), rec
    if pb is None:
        raise SystemExit("the arm needs cells of the second study; give --pbmc-arrays and --pbmc-z")
    all_donors = sorted(set(pb["donors"].tolist()))
    idx, src = [], []
    for donor, n in sorted(cfg["pbmc"].items()):
        if donor not in all_donors:
            raise SystemExit(f"unknown donor {donor}")
        k = all_donors.index(donor)
        cells = np.where(pb["donors"] == donor)[0]
        if n > len(cells):
            raise SystemExit(f"{n} cells asked of {donor}, which has {len(cells)}")
        pick = cells if n == len(cells) else np.sort(cells[np.random.default_rng([seed, 2, k]).permutation(len(cells))[:n]])
        take = pick[np.random.default_rng([seed, 3, k]).permutation(len(pick))] if cfg.get("shuffle_labels") else pick
        idx.append(pick); src.append(take)
    idx, src = np.concatenate(idx), np.concatenate(src)
    yp = pb["adt"][src].astype(np.float32)                 # the protein row a cell is fitted to (its own unless shuffled)
    m = ~np.isnan(yp)
    tot = np.where(m, yp, 0.0).sum(axis=1, dtype=np.float64)
    sfp = (tot / max(float(np.median(tot)), 1e-6)).astype(np.float32)
    sfp[sfp <= 0] = 1.0
    yp = np.where(m, yp, 0.0).astype(np.float32)
    mask = np.vstack([np.ones(y.shape, np.float32), m.astype(np.float32)])
    rec.update(pbmc_index=idx, pbmc_row_from=src)
    return (torch.from_numpy(np.vstack([z, pb["z"][idx]])), torch.from_numpy(np.vstack([y, yp])),
            torch.from_numpy(np.concatenate([sf, sfp])), torch.from_numpy(mask)), rec


def _save(obj, path: Path):
    tmp = path.with_suffix(path.suffix + ".tmp")
    torch.save(obj, tmp)
    os.replace(tmp, path)


def export_head(ckpt, names, path: Path):
    W = {f"mlp__{k}": v.cpu().numpy() for k, v in ckpt["mlp"].items()}
    W.update({f"dec__{k}": v.cpu().numpy() for k, v in ckpt["dec"].items()})
    np.savez(path, adt_names=np.asarray(names).astype("U64"), **W)


def fit(cfg, bm, pb, out: Path, *, time_budget=None, threads=None, max_epochs_this_call=None, log=print):
    """Fit the head `cfg` describes. Returns True when it is finished, False when the budget ended first."""
    out.mkdir(parents=True, exist_ok=True)
    if (out / "DONE").exists():
        return True
    t_call = time.time()
    if threads:
        torch.set_num_threads(int(threads))
    device = torch.device("cpu")
    seed, hidden, epochs = int(cfg["seed"]), int(cfg.get("hidden", 512)), int(cfg.get("epochs", 40))
    seed_everything(seed)
    (tz, ty, tsf, tm), rec = build_train(cfg, bm, pb)
    val = {k: torch.from_numpy(v[bm["split"] == "val"]) for k, v in (("z", bm["z"]), ("adt", bm["adt"]), ("sf", bm["sf"]))}
    test = {k: torch.from_numpy(v[bm["split"] == "test"]) for k, v in (("z", bm["z"]), ("adt", bm["adt"]), ("sf", bm["sf"]))}
    z_dim, n_adt = tz.shape[1], ty.shape[1]

    mlp = MLP(z_dim, z_dim, hidden=hidden).to(device)
    v_net = VelocityNet(z_dim, cond_dim=z_dim, hidden=hidden).to(device)
    dec = AdtDecoder(z_dim, n_adt, hidden=hidden).to(device)
    params = list(mlp.parameters()) + list(v_net.parameters()) + list(dec.parameters())
    opt = torch.optim.AdamW(params, lr=LR)

    n_train = len(tz)
    batch_size = min(BATCH, n_train)
    drop_last = n_train >= batch_size * 2
    gen = torch.Generator().manual_seed(seed)
    tensors = (tz, ty, tsf) if tm is None else (tz, ty, tsf, tm)
    loader = DataLoader(TensorDataset(*tensors), batch_size=batch_size, shuffle=True, drop_last=drop_last, generator=gen)

    history, best, start, spent = [], -1e9, 1, 0.0
    state_path = out / "state.pt"
    if state_path.exists():
        st = torch.load(state_path, map_location=device, weights_only=False)
        if st["cfg"] != cfg:
            raise SystemExit(f"{out}: the saved state belongs to another head")
        mlp.load_state_dict(st["mlp"]); v_net.load_state_dict(st["v_net"]); dec.load_state_dict(st["dec"]); opt.load_state_dict(st["opt"])
        torch.set_rng_state(st["torch_rng"]); gen.set_state(st["loader_rng"])
        np.random.set_state(st["numpy_rng"]); random.setstate(st["python_rng"])
        history, best, start, spent = st["history"], st["best"], st["epoch"] + 1, st["seconds"]
    else:
        np.savez(out / "fit_cells.npz", **rec)

    done_here, per_epoch = 0, None
    for epoch in range(start, epochs + 1):
        if max_epochs_this_call is not None and done_here >= max_epochs_this_call:
            return False
        if time_budget is not None and per_epoch is not None and (time.time() - t_call) + 1.3 * per_epoch > time_budget:
            return False
        t0 = time.time()
        mlp.train(); v_net.train(); dec.train()
        losses = []
        for batch in loader:
            zb, yb, sfb = batch[0].to(device), batch[1].to(device), batch[2].to(device)
            mb = batch[3].to(device) if tm is not None else None
            z_mlp = mlp(zb)
            mu_mlp, th = dec(z_mlp, sfb)
            loss_mlp = nb_nll(yb, mu_mlp, th) if mb is None else nb_nll_cells(yb, mu_mlp, th, mb)

            x0 = torch.randn_like(zb)
            t = torch.rand(zb.size(0), device=device)
            x1 = z_mlp.detach()
            xt = (1.0 - t[:, None]) * x0 + t[:, None] * x1
            loss_fm = ((v_net(xt, t, zb) - (x1 - x0)) ** 2).mean()
            z1 = integrate_fm(v_net, x0, zb, max(FM_STEPS // 2, 4))
            mu_fm, th_fm = dec(z1, sfb)
            loss_fm_nll = nb_nll(yb, mu_fm, th_fm) if mb is None else nb_nll_cells(yb, mu_fm, th_fm, mb)

            loss = LAMBDA_NLL * (loss_mlp + 0.5 * loss_fm_nll) + LAMBDA_FM * loss_fm
            opt.zero_grad()
            loss.backward()
            opt.step()
            losses.append(float(loss.detach()))

        metrics = {"epoch": epoch, "train_loss": float(np.mean(losses))}
        metrics.update(eval_split("val", val, mlp, v_net, dec, device, FM_STEPS))
        history.append(metrics)
        score = metrics["val_fm_pearson"]
        if score > best:
            best = score
            _save({"mlp": mlp.state_dict(), "v_net": v_net.state_dict(), "dec": dec.state_dict(), "metrics": metrics}, out / "best.pt")
        per_epoch = time.time() - t0
        spent += per_epoch
        _save({"cfg": cfg, "epoch": epoch, "mlp": mlp.state_dict(), "v_net": v_net.state_dict(), "dec": dec.state_dict(),
               "opt": opt.state_dict(), "torch_rng": torch.get_rng_state(), "loader_rng": gen.get_state(),
               "numpy_rng": np.random.get_state(), "python_rng": random.getstate(), "history": history, "best": best,
               "seconds": spent}, state_path)
        done_here += 1
        log(f"  {cfg['name']} epoch {epoch}/{epochs} loss {metrics['train_loss']:.4f} val_mlp {metrics['val_mlp_pearson']:.4f} "
            f"val_fm {score:.4f} ({per_epoch:.0f} s)")

    ckpt = torch.load(out / "best.pt", map_location=device, weights_only=False)
    mlp.load_state_dict(ckpt["mlp"]); v_net.load_state_dict(ckpt["v_net"]); dec.load_state_dict(ckpt["dec"])
    test_metrics = eval_split("test", test, mlp, v_net, dec, device, FM_STEPS)
    export_head(ckpt, bm["names"], out / "head.npz")
    (out / "metrics.json").write_text(json.dumps({
        "head": cfg, "history": history, "test": test_metrics, "checkpoint_epoch": int(ckpt["metrics"]["epoch"]),
        "n_fit": int(n_train), "n_fit_bone_marrow": int(len(rec["bm_index"])), "n_fit_second_study": int(len(rec["pbmc_index"])),
        "masked_loss": bool(tm is not None), "seconds": spent,
        "machine": {"device": "cpu", "threads": int(torch.get_num_threads()), "torch": torch.__version__, "platform": sys.platform,
                    "arch": os.uname().machine},
    }, indent=1))
    state_path.unlink()
    (out / "DONE").write_text("")
    return True


def main():
    from plan import head_by_name
    p = argparse.ArgumentParser()
    p.add_argument("--head", required=True)
    for k in ("bm-arrays", "bm-z", "out"):
        p.add_argument(f"--{k}", type=Path, required=True)
    p.add_argument("--pbmc-arrays", type=Path)
    p.add_argument("--pbmc-z", type=Path)
    p.add_argument("--time-budget-sec", type=float)
    p.add_argument("--threads", type=int, default=2)
    a = p.parse_args()
    bm, pb = load_data(a.bm_arrays, a.bm_z, a.pbmc_arrays, a.pbmc_z)
    cfg = head_by_name(a.head, n_train=int((bm["split"] == "train").sum()))
    ok = fit(cfg, bm, pb, a.out / cfg["name"], time_budget=a.time_budget_sec, threads=a.threads)
    sys.exit(0 if ok else EXIT_BUDGET)


if __name__ == "__main__":
    main()
