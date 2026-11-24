#!/usr/bin/env python3
"""Train the last layer of TEDDY-G 70M by the rows of the TEG, on one half of the cells: one arm, one fold.

  python train.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cache CACHE --cells THP1_ENC.npz --arm ARM --fold 0|1 --out OUT
                  [--steps 600] [--bench N] [--pilot] [--lr RATE]

CACHE: the folder of cache.py (what enters the last layer). THP1_ENC.npz: ids, A (1 = knockout), unit (the gene),
y_PDL1 (the measured protein). The cells are split into the two halves of the earlier rounds (anm_diag.cycle.halves,
seed 0). Fold 0 fits on the first half, fold 1 on the second. **The measured protein of the other half is set to NaN
before anything is computed**, and no cell of the other half enters a loss, a stopping check or the readout.

What is trained: every weight of the last encoder layer. Everything else is frozen, the readout included.

The readout r, fixed before training: a ridge readout of the state (the mean over the gene tokens, at unit length)
fitted on the control cells of the fit half, at the state as released. It reads the protein in units of the spread
among those control cells. It is the readout of the round's arm `controls`: it holds no record of any event.

The loss, DESIGN.md Section 4:
  rows   a row of the TEG written as a term. For every event e of the arm's record: (the shift of r from the control
         cells to the cells of e, minus the measured shift)^2, averaged over the events. It uses the protein of a
         knockout cell only through the mean of its event. In a step the shift of e is taken as: the shift as
         released, read once on all cells of the fit, plus the mean change of r among the step's cells of e, minus
         the mean change of r among the step's control cells.
  hold   the control cells stay where they are: the squared distance of a control cell's state from its state as
         released, in units of the spread of the control cells' states.
  loss = rows + hold  (weights 1 and 1; nothing is tuned)

Arms (what the record holds about the knockout cells of the fit half):
  rows4   the rows of IFNGR1, IFNGR2, JAK2, STAT1. No cell of IRF1 is used: IRF1 is an event without a record
  rows5   the same and the row of IRF1
  mixed   as rows5, with the event labels permuted among the knockout cells: the rows say every event shifts alike
  cells5  the measured protein of every knockout cell, cell by cell, in place of the rows (ordinary fine-tuning)
  prop    rows4, and for the cells of IRF1 a row that was not measured: the mean of the four measured shifts

Steps: Adam without momentum (betas 0 and 0.999), learning rate 1e-5; a step is 64 cells (24 control cells and the
rest in equal parts from the events). The cells of a step are first read as the layer stands, which gives the exact
gradient of the step's loss with respect to every cell's r; then four passes of 16 carry it to the weights. Every 25
steps all cells of the fit half are read: training stops when every row of the record is read to within 0.03
(control SD), or after --steps steps. The weights of the last step are kept.

(The loss is the declared one. How its gradient is taken was changed three times before any state was read, each
time by the training curve of a pilot: DESIGN.md, Section 7. With a running value of the shifts and Adam with
momentum the rows swung about their values. Two forms that ask every knockout cell, one by one, to move by what its
row is missing were tried and dropped: the first did not move the rows, the second held them short of their values,
because it also asks the cells of an event to move alike.)

Writes OUT/<arm>_<fold>.npz: the state of EVERY cell under the trained layer (ZT), the state as released (Z0), the
readout, and the log. --bench N: time N steps and stop (nothing is written). --pilot: train as declared and write the
training curve only (OUT/pilot_<arm>_<fold>.json); no state of a trained layer is computed or kept.
"""
import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np

TYPE_A, ALL5 = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"], ["IFNGR1", "IFNGR2", "IRF1", "JAK2", "STAT1"]
ARMS = ("rows4", "rows5", "mixed", "cells5", "prop")
GRID, FOLDS = [10.0 ** (k / 2) for k in range(-6, 5)], 5        # the grid and folds of anm_diag.cycle.LinearRows
LR, N_CTRL, STEP_CELLS, PASS, EVERY, TOL, MIX_SEED = 1e-5, 24, 64, 16, 25, 0.03, 20261124
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def halves(groups, seed=0, share=0.5):
    """anm_diag.cycle.halves, copied: the halves of the earlier rounds."""
    rng, first = np.random.default_rng(seed), np.zeros(len(groups), bool)
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        first[idx[: int(round(share * len(idx)))]] = True
    return first


def unit(z):
    z = z.astype(np.float64)
    return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)


def ridge(x, t):
    """The readout of anm_diag.cycle.LinearRows, arm `controls`: features centred on the cells given, one ridge weight by
    5-fold cross-validation, the intercept not penalised."""
    xc = x.mean(0)
    F = np.hstack([x - xc, np.ones((len(x), 1))])
    pen = np.r_[np.ones(x.shape[1]), 0.0] * x.var(0).mean()
    fold = np.random.default_rng(0).permutation(len(F)) % FOLDS
    err = np.zeros(len(GRID))
    for k in range(FOLDS):
        tr, te_ = fold != k, fold == k
        M, v = F[tr].T @ F[tr] / tr.sum(), F[tr].T @ t[tr] / tr.sum()
        for i, lam in enumerate(GRID):
            err[i] += ((F[te_] @ np.linalg.solve(M + lam * np.diag(pen), v) - t[te_]) ** 2).sum()
    lam = GRID[int(err.argmin())]
    th = np.linalg.solve(F.T @ F / len(F) + lam * np.diag(pen), F.T @ t / len(F))
    return xc, th[:-1], float(th[-1]), lam, float(err.min() / len(F))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, required=True)
    ap.add_argument("--ckpt", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--arm", choices=ARMS, required=True)
    ap.add_argument("--fold", type=int, choices=(0, 1), required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--steps", type=int, default=600)
    ap.add_argument("--bench", type=int, default=0)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--lr", type=float, default=LR, help="for the pilot only; the trial runs with the declared value")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    out_file = a.out / f"{a.arm}_{a.fold}.npz"
    if out_file.exists() and not (a.bench or a.pilot):
        raise SystemExit(f"{out_file.name} exists: this arm and fold were trained once")
    import torch
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    ix = np.load(a.cache / "index.npz")
    ids, lens = ix["ids"].astype(str), ix["L"].astype(np.int64)
    H = np.load(a.cache / "h.npy", mmap_mode="r")
    K = np.load(a.cells)
    assert np.array_equal(K["ids"].astype(str), ids), "the cells file is of other cells, or in another order"
    A, un = K["A"].astype(int), K["unit"].astype(str)
    first = halves(np.where(A == 1, un, "control"), 0)
    fit = first if a.fold == 0 else ~first
    y = K["y_PDL1"].astype(np.float64).copy()
    y[~fit] = np.nan                                             # the protein of the other half is not available from here on
    C = np.where(fit & (A == 0))[0]
    t = (y - y[C].mean()) / y[C].std()                           # the protein in units of the spread among the control cells of the fit
    label = un.copy()
    if a.arm == "mixed":                                         # the event labels permuted among the knockout cells of the fit
        ko = np.where(fit & (A == 1))[0]
        label[ko] = label[ko][np.random.default_rng(MIX_SEED + a.fold).permutation(len(ko))]
    names = TYPE_A if a.arm in ("rows4", "prop") else ALL5
    ev = {g: np.where(fit & (A == 1) & (label == g))[0] for g in names}
    target = {g: float(t[T].mean()) for g, T in ev.items()}       # the measured shift: the row
    if a.arm == "prop":                                          # a row that was not measured, for the cells of IRF1
        ev["IRF1"] = np.where(fit & (A == 1) & (un == "IRF1"))[0]
        target["IRF1"] = float(np.mean([target[g] for g in TYPE_A]))
    groups = list(ev)
    assert all(np.isfinite(t[T]).all() for T in ev.values()) and np.isfinite(t[C]).all()

    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    for p in m.parameters():
        p.requires_grad_(False)
    last = m.encoder.layers[-1]
    last.eval()                                                  # no dropout; the weights still take gradients
    par = list(last.parameters())
    for p in par:
        p.requires_grad_(True)

    def states(idx, grad=False, batch=32):
        """The mean over the gene tokens of the last layer's output, for the cells idx (in their order)."""
        outs = []
        with torch.set_grad_enabled(grad):
            for s in range(0, len(idx), batch):
                j = np.asarray(idx[s:s + batch])
                o = np.argsort(j)
                Lb = int(lens[j].max())
                h = torch.from_numpy(np.ascontiguousarray(H[j[o], :Lb]).astype(np.float32)).to(dev)
                am = (torch.arange(Lb)[None, :] < torch.from_numpy(lens[j[o]])[:, None]).long().to(dev)
                z = te.pool_hidden(last(src=h, src_key_padding_mask=~am.bool()).float(), am, "gene-mean")
                outs.append(z[torch.from_numpy(np.argsort(o)).to(dev)])
        return torch.cat(outs)

    say(f"arm {a.arm}, fold {a.fold}: {int(fit.sum())} cells of the fit ({len(C)} control), events {json.dumps({g: int(len(T)) for g, T in ev.items()})}")
    Z0 = states(np.arange(len(ids))).cpu().numpy().astype(np.float32)   # the state as released, every cell
    rec = float(np.abs(Z0 - K["L12"]).max())
    assert rec <= 5e-3, f"the state from what is kept differs from the recorded state by {rec}"
    x0 = unit(Z0)
    xc, w, b, lam, cv = ridge(x0[C], t[C])
    v0 = float(((x0[C] - x0[C].mean(0)) ** 2).sum(1).mean())      # the spread of the control cells' states
    W, XC, X0 = torch.tensor(w, dtype=torch.float32, device=dev), torch.tensor(xc, dtype=torch.float32, device=dev), torch.tensor(x0, dtype=torch.float32, device=dev)
    at_ = lambda i: torch.from_numpy(np.asarray(i, np.int64)).to(dev)
    r_of = lambda z: ((z / (z.norm(dim=1, keepdim=True) + 1e-6)) - XC) @ W + b

    def full():
        """Every row of the record read on all cells of the fit, and how far the control cells have moved."""
        with torch.no_grad():
            zc = states(C)
            rc = float(r_of(zc).mean())
            hold = float((((zc / (zc.norm(dim=1, keepdim=True) + 1e-6)) - X0[at_(C)]) ** 2).sum(1).mean()) / v0
            s = {g: float(r_of(states(T)).mean()) - rc for g, T in ev.items()}
        return s, hold

    s_hat, hold = full()
    err = max(abs(s_hat[g] - target[g]) for g in list(ev))       # read on all cells of the fit; only this value stops the training
    r0 = (x0 - xc) @ w + b                                       # the readout at the released state, every cell
    missing = np.zeros(len(ids))                                 # for a knockout cell of the record: what the row says is missing
    for g, T in ev.items():
        missing[T] = target[g] - s_hat[g]
    tt = np.where(np.isfinite(t), t, 0.0)
    log = {"arm": a.arm, "fold": a.fold, "cells_of_the_fit": int(fit.sum()), "control_cells_of_the_fit": int(len(C)), "events": {g: int(len(T)) for g, T in ev.items()}, "rows": target,
           "readout": {"ridge_weight": lam, "cv_error": cv, "spread_of_the_control_states": v0}, "released_state_against_the_record": rec,
           "curve": [{"step": 0, "shift": dict(s_hat), "hold": hold, "largest_error_of_a_row": max(abs(s_hat[g] - target[g]) for g in groups)}]}
    say("step 0:", json.dumps({g: round(s_hat[g], 3) for g in groups}), "rows", json.dumps({g: round(target[g], 3) for g in groups}), f"hold {hold:.4f}")
    opt = torch.optim.Adam(par, lr=a.lr, betas=(0.0, 0.999))
    rng = np.random.default_rng(1000 + a.fold)
    per = (STEP_CELLS - N_CTRL) // len(groups)
    pools = {"control": rng.permutation(C), **{g: rng.permutation(T) for g, T in ev.items()}}
    at = {k: 0 for k in pools}

    def draw(k, n):
        if at[k] + n > len(pools[k]):
            pools[k], at[k] = rng.permutation(pools[k]), 0
        at[k] += n
        return pools[k][at[k] - n:at[k]]

    steps, stopped, tb = 0, "the limit of steps", time.time()
    limit = a.bench or a.steps
    while steps < limit and (a.bench or err > TOL):
        ic = draw("control", N_CTRL)
        ie = {g: draw(g, per) for g in groups}
        idx = np.concatenate([ic] + [ie[g] for g in groups])
        grp = np.concatenate([np.full(len(ic), -1)] + [np.full(per, k) for k, _g in enumerate(groups)])
        with torch.no_grad():                                    # the cells of the step as the layer stands: the exact gradient of the step's loss in r
            rr = r_of(states(idx, batch=PASS)).cpu().numpy().astype(np.float64)
        ctl = grp == -1
        coef = np.zeros(len(idx))
        if a.arm == "cells5":                                    # every cell of the step against its own measured protein
            coef[:] = 2.0 * (rr - tt[idx]) / len(idx)
        else:
            d = rr - r0[idx]                                     # the change of r, cell by cell
            res = np.array([(d - missing[idx])[grp == k].mean() for k in range(len(groups))]) - d[ctl].mean()   # for every event: its shift as it stands, minus the row
            for k in range(len(groups)):
                coef[grp == k] = 2.0 * res[k] / (len(groups) * per)
            coef[ctl] = -2.0 * res.mean() / ctl.sum()
        cf = torch.tensor(coef, dtype=torch.float32, device=dev)
        opt.zero_grad(set_to_none=True)
        o = rng.permutation(len(idx))
        for s in range(0, len(idx), PASS):
            j = o[s:s + PASS]
            z = states(idx[j], grad=True, batch=PASS)
            r = r_of(z)
            g_ = torch.from_numpy(grp[j]).to(dev)
            isc = g_ == -1
            loss = (((z[isc] / (z[isc].norm(dim=1, keepdim=True) + 1e-6)) - X0[at_(idx[j])][isc]) ** 2).sum() / (v0 * len(ic)) if bool(isc.any()) else z.sum() * 0.0
            loss = loss + (cf[at_(j)] * r).sum()
            loss.backward()
        opt.step()
        steps += 1
        if steps % EVERY == 0 and not a.bench:
            s_hat, hold = full()
            err = max(abs(s_hat[g] - target[g]) for g in groups)
            log["curve"].append({"step": steps, "shift": dict(s_hat), "hold": hold, "largest_error_of_a_row": err})
            say(f"step {steps}:", json.dumps({g: round(s_hat[g], 3) for g in groups}), f"largest error of a row {err:.3f}, hold {hold:.4f}, {STEP_CELLS * steps / (time.time() - tb):.1f} cells/s")
            if err <= TOL:
                stopped = "every row read to within the tolerance"
    if a.bench:
        say(f"bench: {steps} steps of {STEP_CELLS} cells in {time.time() - tb:.1f} s: {STEP_CELLS * steps / (time.time() - tb):.1f} cells/s on {dev.type}")
        return
    if steps % EVERY:
        s_hat, hold = full()
        log["curve"].append({"step": steps, "shift": dict(s_hat), "hold": hold, "largest_error_of_a_row": max(abs(s_hat[g] - target[g]) for g in groups)})
    if steps == 0:
        stopped = "every row was read to within the tolerance before the first step"
    if a.pilot:
        log.update({"steps": steps, "stopped_by": stopped, "tolerance": TOL, "learning_rate": a.lr, "seconds": round(time.time() - t0, 1), "cells_per_second": round(STEP_CELLS * steps / max(time.time() - tb, 1e-9), 1)})
        (a.out / f"pilot_{a.arm}_{a.fold}_lr{a.lr:g}.json").write_text(json.dumps(log, indent=1) + "\n")
        say("pilot done:", json.dumps({k: log[k] for k in ("steps", "stopped_by", "seconds")}))
        return
    ZT = states(np.arange(len(ids))).cpu().numpy().astype(np.float32)
    log.update({"steps": steps, "stopped_by": stopped, "tolerance": TOL, "learning_rate": a.lr, "weights_trained": int(sum(p.numel() for p in par)), "device": dev.type,
                "torch": torch.__version__, "numpy": np.__version__, "python": sys.version.split()[0], "seconds": round(time.time() - t0, 1),
                "largest_move_of_a_state": float(np.abs(ZT - Z0).max())})
    np.savez(a.out / f"tmp_{a.arm}_{a.fold}.npz", ids=ids, ZT=ZT, Z0=Z0, fit=fit, w=w.astype(np.float32), xc=xc.astype(np.float32), b=np.float32(b), log=json.dumps(log))
    (a.out / f"tmp_{a.arm}_{a.fold}.npz").rename(out_file)
    say("done:", json.dumps({k: log[k] for k in ("steps", "stopped_by", "seconds")}), json.dumps(log["curve"][-1]))


if __name__ == "__main__":
    main()
