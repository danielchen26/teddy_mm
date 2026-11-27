#!/usr/bin/env python3
"""A state that is not one mean: a pooling over the gene positions, trained by the rows of the TEG. One arm, one fold.

  python pool.py --repo TEDDY_MM --ckpt TEDDY_G_70M --cache CACHE --cells THP1_ENC.npz --arm ARM --fold 0|1 --out OUT
                 [--steps 600] [--pilot] [--lr RATE]
  python pool.py --selftest          made-up cells and a made-up last layer: the code, with no model and no data

NOT RUN. This is the script of DESIGN.md, to be run on the owner's word. CACHE and THP1_ENC.npz are those of
runs/last_layer_trial_20261124 (what enters the last layer of TEDDY-G 70M for the 7,478 cells of the THP-1 stratum;
ids, A, unit, y_PDL1, L12). The halves are those of the earlier rounds. **The measured protein of the other half is
set to NaN before anything is computed.**

What is trained: K vectors of the size of the model's width (K = 4: 2,048 numbers). The encoder is frozen, all of it.
TEDDY gives a cell one vector, the mean of the last layer's output over its gene positions. Here a cell gets K: block k
is the same output averaged with weights softmax_i(q_k . h_i / sqrt(d)) over the cell's gene positions i. With q = 0
every block is the mean, so training starts from the released state exactly. The state of a cell is its K blocks,
each at unit length, side by side.

The loss is one term, a row of the TEG written as a constraint on the state: **a readout taught by control cells alone
has to read the recorded events.** At every step a ridge readout r is fitted on 256 control cells of the fit half at
the state as it stands (in closed form, so the gradient passes through the fit), and for every event e of the arm's
record: (the shift of r from other control cells to the cells of e, minus the measured shift)^2, averaged over the
events. The readout holds no record of any event at any step, so the rows can be met only by a state in which the
events lie where control cells that differ in the protein lie. That is what the loop's own first system reads (S0).
There is no second term: a readout fitted again at every step leaves nothing to hold in place.
Arms: pool4 (the rows of IFNGR1, IFNGR2, JAK2, STAT1; no cell of IRF1 is used), pool5 (and the row of IRF1), mixed (as
pool5 with the event labels permuted among the knockout cells), none (no row: nothing is trained and nothing may move).

Writes OUT/<arm>_<fold>.npz: the state of EVERY cell under the trained pooling (PT, K blocks), the mean as released
(Z0), the vectors q and the log. --pilot: the training curve only; no trained state is kept.
"""
import argparse
import importlib.util
import json
import os
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "last_layer_trial_20261124" / "scripts"))
from train import ALL5, EVERY, MIX_SEED, N_CTRL, PASS, STEP_CELLS, TOL, TYPE_A, halves, ridge  # noqa: E402  (the halves, the readout and the step are those of the trial)

ARMS = ("pool4", "pool5", "mixed", "none")
K, LR, TEACHERS = 4, 1e-2, 256                                 # the number of blocks; the learning rate, to be confirmed by a pilot before any trained state is read; the control cells the readout is fitted on
t0 = time.time()
say = lambda *a: print(f"[{time.time() - t0:7.1f}s]", *a, flush=True)


def run(a, H, ids, lens, cells, last, dev, pooled_mean=None):
    import torch
    A, un = cells["A"].astype(int), cells["unit"].astype(str)
    first = halves(np.where(A == 1, un, "control"), 0)
    fit = first if a.fold == 0 else ~first
    y = cells["y_PDL1"].astype(np.float64).copy()
    y[~fit] = np.nan                                             # the protein of the other half is not available from here on
    C = np.where(fit & (A == 0))[0]
    t = (y - y[C].mean()) / y[C].std()
    label = un.copy()
    if a.arm == "mixed":
        ko = np.where(fit & (A == 1))[0]
        label[ko] = label[ko][np.random.default_rng(MIX_SEED + a.fold).permutation(len(ko))]
    names = [] if a.arm == "none" else [g for g in (TYPE_A if a.arm == "pool4" else ALL5) if ((label == g) & fit & (A == 1)).any()]
    ev = {g: np.where(fit & (A == 1) & (label == g))[0] for g in names}
    target = {g: float(t[T].mean()) for g, T in ev.items()}
    groups = list(ev)
    d = H.shape[2]
    Q = torch.zeros(K, d, device=dev, requires_grad=True)

    def tokens(idx):
        """The last layer's output at every position, for the cells idx; nothing here takes a gradient."""
        j = np.asarray(idx)
        Lb = int(lens[j].max())
        hs, ms = [], []
        with torch.no_grad():
            for s in range(0, len(j), PASS):
                jj = j[s:s + PASS]
                o = np.argsort(jj)
                h = torch.from_numpy(np.ascontiguousarray(H[jj[o], :Lb]).astype(np.float32)).to(dev)
                am = (torch.arange(Lb)[None, :] < torch.from_numpy(lens[jj[o]])[:, None]).to(dev)
                back = torch.from_numpy(np.argsort(o)).to(dev)
                hs.append(last(src=h, src_key_padding_mask=~am).float()[back]); ms.append(am[back])
        return torch.cat(hs), torch.cat(ms)

    def state(h, am, q):
        """K blocks side by side, each at unit length, the whole at unit length."""
        w = torch.einsum("bld,kd->bkl", h, q) / d ** 0.5
        w = torch.softmax(w.masked_fill(~am[:, None, :], float("-inf")), dim=2)
        blocks = torch.einsum("bkl,bld->bkd", w, h)
        return (blocks / (blocks.norm(dim=2, keepdim=True) + 1e-6)).reshape(len(h), -1) / K ** 0.5, blocks

    def states(idx, q, batch=64):
        out, bl = [], []
        with torch.no_grad():
            for s in range(0, len(idx), batch):
                x, b = state(*tokens(idx[s:s + batch]), q)
                out.append(x); bl.append(b)
        return torch.cat(out), torch.cat(bl)

    everyone = np.arange(len(ids))
    X0, B0 = states(everyone, torch.zeros(K, d, device=dev))
    Z0 = B0[:, 0].cpu().numpy().astype(np.float32)               # with q = 0 a block is the mean over the gene positions: the state as released
    if pooled_mean is not None:
        rec = float(np.abs(Z0 - pooled_mean).max())
        assert rec <= 5e-3, f"the mean from what is kept differs from the recorded state by {rec}"
    x0 = X0.cpu().numpy().astype(np.float64)
    order = np.random.default_rng(2000 + a.fold).permutation(C)
    R, base = order[:min(TEACHERS, len(C) // 2)], order[min(TEACHERS, len(C) // 2):]       # the control cells that teach the readout; the control cells an event is read against
    _xc, _w, _b, lam, cv = ridge(x0[R], t[R])                    # the ridge weight, chosen once at the state as released
    pen = float(lam * x0[R].var(0).mean())
    hR, mR = tokens(R)                                           # the last layer's output for the teachers, kept: nothing in it is trained
    tR = torch.tensor(t[R], dtype=torch.float32, device=dev)

    def readout(q):
        """The ridge readout of the teachers at the state q gives, in closed form: r(x) = (x - xc) . w + b."""
        xr = torch.cat([state(hR[s:s + 64], mR[s:s + 64], q)[0] for s in range(0, len(R), 64)])
        xc = xr.mean(0)
        F = xr - xc
        alpha = torch.linalg.solve(F @ F.T / len(R) + pen * torch.eye(len(R), device=dev), (tR - tR.mean()) / len(R))
        return xc, F.T @ alpha, tR.mean()

    def full(q):
        with torch.no_grad():
            xc, w, b = readout(q)
            r_of = lambda idx: (states(idx, q)[0] - xc) @ w + b
            rb = float(r_of(base).mean())
            return {g: float(r_of(T).mean()) - rb for g, T in ev.items()}

    zero = torch.zeros(K, d, device=dev)
    with torch.no_grad():
        xc0, w0, b0 = readout(zero)
        r0 = ((X0 - xc0) @ w0 + b0).cpu().numpy()                # the readout of the teachers at the state as released, every cell
    s0 = full(zero)
    err = max([abs(s0[g] - target[g]) for g in groups] or [0.0])
    log = {"arm": a.arm, "fold": a.fold, "blocks": K, "numbers_trained": int(Q.numel()), "cells_of_the_fit": int(fit.sum()), "control_cells_of_the_fit": int(len(C)), "teachers": int(len(R)),
           "events": {g: int(len(T)) for g, T in ev.items()}, "rows": target, "readout": {"ridge_weight": lam, "cv_error_at_the_start": cv}, "curve": [{"step": 0, "shift": dict(s0), "largest_error_of_a_row": err}]}
    say(f"arm {a.arm}, fold {a.fold}: {int(fit.sum())} cells of the fit ({len(C)} control, {len(R)} teach the readout); step 0:", json.dumps({g: round(s0[g], 3) for g in groups}), "rows", json.dumps({g: round(v, 3) for g, v in target.items()}))
    opt = torch.optim.Adam([Q], lr=a.lr)
    rng = np.random.default_rng(1000 + a.fold)
    per = (STEP_CELLS - N_CTRL) // max(len(groups), 1)
    pools = {"control": rng.permutation(base), **{g: rng.permutation(T) for g, T in ev.items()}}
    at = {k: 0 for k in pools}

    def draw(k, n):
        if at[k] + n > len(pools[k]):
            pools[k], at[k] = rng.permutation(pools[k]), 0
        at[k] += n
        return pools[k][at[k] - n:at[k]]

    steps, stopped, tb = 0, "the limit of steps" if groups else "no row: nothing to train", time.time()
    while groups and steps < a.steps and err > TOL:
        ic = draw("control", N_CTRL)
        ie = {g: draw(g, per) for g in groups}
        idx = np.concatenate([ic] + [ie[g] for g in groups])
        grp = torch.from_numpy(np.concatenate([np.full(len(ic), -1)] + [np.full(per, k) for k, _g in enumerate(groups)])).to(dev)
        xc, w, b = readout(Q)
        x, _b = state(*tokens(idx), Q)
        dr = (x - xc) @ w + b - torch.tensor(r0[idx], dtype=torch.float32, device=dev)     # the change of r, cell by cell, from the state as released
        loss = sum((s0[g] + dr[grp == k].mean() - dr[grp == -1].mean() - target[g]) ** 2 for k, g in enumerate(groups)) / len(groups)
        opt.zero_grad(set_to_none=True)
        loss.backward()
        opt.step()
        steps += 1
        if steps % EVERY == 0 or steps == a.steps:
            s_hat = full(Q.detach())
            err = max(abs(s_hat[g] - target[g]) for g in groups)
            log["curve"].append({"step": steps, "shift": dict(s_hat), "largest_error_of_a_row": err, "loss_of_the_step": float(loss)})
            say(f"step {steps}:", json.dumps({g: round(s_hat[g], 3) for g in groups}), f"largest error of a row {err:.3f}, {len(idx) * steps / (time.time() - tb):.1f} cells/s")
            if err <= TOL:
                stopped = "every row read to within the tolerance"
    log.update({"steps": steps, "stopped_by": stopped, "tolerance": TOL, "learning_rate": a.lr, "device": dev.type, "seconds": round(time.time() - t0, 1)})
    if a.pilot:
        return log, None
    PT = states(everyone, Q.detach())[1].cpu().numpy().astype(np.float32)
    return log, {"ids": ids, "PT": PT, "Z0": Z0, "fit": fit, "teachers": R, "q": Q.detach().cpu().numpy()}


def selftest():
    """Made-up cells: the code from the cache to the written file, with a made-up last layer. Three positions of every cell
    hold a quantity u that the protein follows, and look alike; an event raises u. Averaged over all positions, u is
    one part in eight beside the noise of the others, so a readout of the mean reads a part of an event's shift. A
    pooling that weighs those positions reads more of it: the rows can be met only by finding them."""
    import torch
    torch.manual_seed(0)
    rng = np.random.default_rng(0)
    n, L, d = 1400, 24, 16
    lens = rng.integers(20, L + 1, n)
    H = rng.normal(0, 1, (n, L, d))
    mark, eff = np.linalg.qr(rng.normal(0, 1, (d, 2)))[0].T
    unit_ = np.array(["control"] * 900 + [g for g in ALL5 for _ in range(100)])
    A = (unit_ != "control").astype(int)
    size = dict(zip(ALL5, [1.5, 1.2, 0.9, 1.8, 1.4]))
    u = rng.normal(0, 1, n) + np.array([size.get(g, 0.0) for g in unit_])
    H[:, :3] += 4 * mark + u[:, None, None] * eff                # what the three positions hold
    y = -u + rng.normal(0, 0.3, n)
    for i, l in enumerate(lens):
        H[i, l:] = 0
    H = H.astype(np.float16)
    cells = {"A": A, "unit": unit_, "y_PDL1": y}
    ids = np.array([f"c{i}" for i in range(n)])
    last = lambda src, src_key_padding_mask: src                 # the made-up last layer passes its input on
    out = {}
    for arm in ("none", "pool4", "mixed"):
        out[arm] = run(argparse.Namespace(arm=arm, fold=0, steps=200, lr=0.05, pilot=False), H, ids, lens, cells, last, torch.device("cpu"))
    log, res = out["none"]
    assert log["steps"] == 0 and np.abs(res["q"]).max() == 0 and np.abs(res["PT"] - res["Z0"][:, None]).max() < 1e-5, "with no row nothing may move"
    log, res = out["pool4"]
    c0, c1 = log["curve"][0], log["curve"][-1]
    assert c0["largest_error_of_a_row"] > 0.3 and c1["largest_error_of_a_row"] < 0.8 * c0["largest_error_of_a_row"], (c0, c1)
    assert res["PT"].shape == (n, K, d) and np.isfinite(res["PT"]).all() and np.abs(res["PT"][:, 0] - res["Z0"]).max() > 1e-3
    mean_ok = np.abs(res["Z0"] - np.array([H[i, :l].astype(np.float32).mean(0) for i, l in enumerate(lens)])).max()
    assert mean_ok < 1e-3, mean_ok                               # with q = 0 a block is the mean over the cell's positions
    # the judge is not the readout the pooling was trained under: a readout fitted again on control cells, at each state,
    # reads the event that had no row on the cells of the other half
    other, fitc = ~res["fit"], res["fit"] & (A == 0)
    tt = (y - y[fitc].mean()) / y[fitc].std()
    measured = float(tt[other & (unit_ == "IRF1")].mean() - tt[other & (A == 0)].mean())

    def read(x):
        xc, w, b, _lam, _cv = ridge(x[fitc], tt[fitc])
        r = (x - xc) @ w + b
        return float(r[other & (unit_ == "IRF1")].mean() - r[other & (A == 0)].mean())
    flat = lambda B: (B / (np.linalg.norm(B, axis=2, keepdims=True) + 1e-6)).reshape(n, -1).astype(np.float64) / K ** 0.5
    before, after = read(flat(np.repeat(res["Z0"][:, None], K, 1))), read(flat(res["PT"]))
    assert abs(after - measured) < abs(before - measured) - 0.1, (before, after, measured)
    with tempfile.TemporaryDirectory() as tmp:
        np.savez(Path(tmp) / "x.npz", **res, log=json.dumps(log))
        assert set(np.load(Path(tmp) / "x.npz").files) >= {"PT", "Z0", "q", "fit"}
    print("selftest passed:", {arm: {"largest error of a row": [round(l["curve"][0]["largest_error_of_a_row"], 3), round(l["curve"][-1]["largest_error_of_a_row"], 3)], "steps": l["steps"]} for arm, (l, _r) in out.items()},
          "| IRF1, no row, cells of the other half, a readout fitted again on control cells: at the mean", round(before, 2), "at the trained state", round(after, 2), "measured", round(measured, 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--repo", type=Path)
    ap.add_argument("--ckpt", type=Path)
    ap.add_argument("--cache", type=Path)
    ap.add_argument("--cells", type=Path)
    ap.add_argument("--arm", choices=ARMS)
    ap.add_argument("--fold", type=int, choices=(0, 1))
    ap.add_argument("--out", type=Path)
    ap.add_argument("--steps", type=int, default=600)
    ap.add_argument("--pilot", action="store_true")
    ap.add_argument("--lr", type=float, default=LR, help="for the pilot only; the experiment runs with the declared value")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    assert all(v is not None for v in (a.repo, a.ckpt, a.cache, a.cells, a.arm, a.fold, a.out))
    a.out.mkdir(parents=True, exist_ok=True)
    out_file = a.out / f"{a.arm}_{a.fold}.npz"
    if out_file.exists() and not a.pilot:
        raise SystemExit(f"{out_file.name} exists: this arm and fold were trained once")
    import torch
    spec = importlib.util.spec_from_file_location("teddy_encoder", a.repo / "teddy_mm" / "teddy_encoder.py")
    te = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(te)
    ix = np.load(a.cache / "index.npz")
    ids, lens = ix["ids"].astype(str), ix["L"].astype(np.int64)
    H = np.load(a.cache / "h.npy", mmap_mode="r")
    cells = np.load(a.cells)
    assert np.array_equal(cells["ids"].astype(str), ids), "the cells file is of other cells, or in another order"
    dev = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    m = te.load_teddy(a.ckpt, dev)
    for p in m.parameters():
        p.requires_grad_(False)
    last = m.encoder.layers[-1]
    last.eval()
    log, res = run(a, H, ids, lens, cells, last, dev, pooled_mean=cells["L12"])
    log.update({"torch": torch.__version__, "numpy": np.__version__, "python": sys.version.split()[0]})
    if a.pilot:
        (a.out / f"pilot_{a.arm}_{a.fold}_lr{a.lr:g}.json").write_text(json.dumps(log, indent=1) + "\n")
        return say("pilot done:", json.dumps({k: log[k] for k in ("steps", "stopped_by", "seconds")}))
    np.savez(a.out / f"tmp_{a.arm}_{a.fold}.npz", **res, log=json.dumps(log))
    (a.out / f"tmp_{a.arm}_{a.fold}.npz").rename(out_file)
    say("done:", json.dumps({k: log[k] for k in ("steps", "stopped_by", "seconds")}), json.dumps(log["curve"][-1]))


if __name__ == "__main__":
    main()
