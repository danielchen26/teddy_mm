"""Repair the head's MLP as in runs/loop_20261112, with one more term in the loss that holds the head's outputs on
bone-marrow training cells where they are (DESIGN.md).

  python fit_guarded.py pilot <data set>   one fit per weight of the grid, first target held out; chooses the weight on
                                           bone-marrow validation cells and writes heads/pilot_<data set>.json
  python fit_guarded.py fit <data set>     the fits for the other four held-out targets at the chosen weight
  python fit_guarded.py features           after both data sets: the features of the 12 audits, results/fits.json
  python fit_guarded.py --check            gradient against finite differences on a few cells; one 20-step fit; writes nothing
"""
import json
import shutil
import sys
import time

import numpy as np
from scipy.special import expit

from _paths import ARM, BM_FIT, FEATURES, LOOP, OUT_FEATURES, RES, TARGETS

MLP = ["mlp__net.0", "mlp__net.3", "mlp__net.6"]
DEC = ["dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"]
LAYERS, ACT = MLP + DEC, [True, True, False, True, True, False]
LR, STEPS, EVERY, PATIENCE, KEEP, SEED = 3e-4, 400, 20, 3, 1.0, 0      # as in runs/loop_20261112
GRID, VAL_DROP = (1.0, 10.0, 100.0), 0.01                               # weights of the holding term; allowed fall on validation cells
F = np.float32

silu = lambda a: a * expit(a)


def dsilu(a):
    s = expit(a)
    return s * (1 + a * (1 - s))


def forward(P, x):
    h, ins, pre, lat = x, [], [], None
    for name, act in zip(LAYERS, ACT):
        ins.append(h)
        a = h @ P[name + ".weight"].T + P[name + ".bias"]
        pre.append(a)
        h = silu(a) if act else a
        if name == "mlp__net.6":
            lat = h
    sp = np.logaddexp(F(0), h)
    return lat, np.log(np.maximum(sp, F(1e-6))), (ins, pre, h, sp)


def backward(P, cache, g):
    ins, pre, m, sp = cache
    g = (g * np.where(sp > 1e-6, expit(m) / sp, 0)).astype(g.dtype)
    G = {}
    for i in range(len(LAYERS) - 1, -1, -1):
        if ACT[i]:
            g = g * dsilu(pre[i])
        G[LAYERS[i] + ".weight"], G[LAYERS[i] + ".bias"] = g.T @ ins[i], g.sum(0)
        g = g @ P[LAYERS[i] + ".weight"]
    return G


def loss_grad(P, x, t, out0, j, xb, ob0, lam, grad=True):
    """The loss of runs/loop_20261112 on perturbation cells (squared error of the PD-L1 output against its target, plus KEEP
    times the mean squared change of the other outputs), plus lam times the mean squared change of all outputs on
    bone-marrow cells."""
    _, out, cache = forward(P, x)
    _, ob, cache_b = forward(P, xb)
    n, k = out.shape
    other = np.ones(k, bool); other[j] = False
    L = float(((out[:, j] - t) ** 2).mean() + KEEP * ((out[:, other] - out0[:, other]) ** 2).mean() + lam * ((ob - ob0) ** 2).mean())
    if not grad:
        return L, None
    g = (2 * KEEP / (n * (k - 1))) * (out - out0)
    g[:, j] = 2 * (out[:, j] - t) / n
    G = backward(P, cache, g)
    Gb = backward(P, cache_b, (2 * lam / ob.size) * (ob - ob0))
    return L, {key: G[key] + Gb[key] for key in G}


def fit(P0, x, t, out0, j, use, val, xb, ob0, bval, lam, steps=STEPS):
    """Adam on the three MLP layers only; the weights kept are those with the lowest loss on the set-aside cells."""
    keys = [l + s for l in MLP for s in (".weight", ".bias")]
    P = {k: v.copy() for k, v in P0.items()}
    m = {k: np.zeros_like(P[k]) for k in keys}; v2 = {k: np.zeros_like(P[k]) for k in keys}
    best, best_P, curve, bad = np.inf, None, [], 0
    for step in range(1, steps + 1):
        L, G = loss_grad(P, x[use], t[use], out0[use], j, xb[~bval], ob0[~bval], lam)
        for k in keys:
            m[k] = 0.9 * m[k] + 0.1 * G[k]; v2[k] = 0.999 * v2[k] + 0.001 * G[k] ** 2
            P[k] = (P[k] - LR * (m[k] / (1 - 0.9 ** step)) / (np.sqrt(v2[k] / (1 - 0.999 ** step)) + 1e-8)).astype(F)
        if step % EVERY == 0:
            Lv = loss_grad(P, x[val], t[val], out0[val], j, xb[bval], ob0[bval], lam, grad=False)[0]
            curve.append([step, L, Lv])
            if Lv < best:
                best, best_P, bad = Lv, {k: v.copy() for k, v in P.items()}, 0
            else:
                bad += 1
                if bad >= PATIENCE:
                    break
    assert all(np.array_equal(best_P[k], P0[k]) for k in P0 if k not in keys)     # the decoder is bit-for-bit the head as trained
    return best_P, curve


def halves(groups, rng, share=0.5):
    first = np.zeros(len(groups), bool)
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        first[idx[: int(round(share * len(idx)))]] = True
    return first


def load(ds):
    """The cells and the split of runs/loop_20261112 (same code, same seed); the audit controls are checked against its record."""
    d = np.load(FEATURES / ds / "dense.npz"); D = {k: d[k] for k in d.files}
    Z = (D.pop("node__Z") if "node__Z" in D else np.load(FEATURES / ds / "node__Z.npy")).astype(F)
    D = {k: v for k, v in D.items() if not k.startswith("node__")}
    x = (Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-6)).astype(F)
    A, tg = D["A"], D["target"].astype(str)
    st = D["stratum"].astype(str) if "stratum" in D else np.full(len(A), "all")
    rng = np.random.default_rng(SEED)
    ctl = A == 0
    audit_ctl = ctl & halves(np.where(ctl, st, "treated"), rng)
    grp = np.where(A == 1, np.char.add(np.char.add(D["guide"].astype(str), "|"), st), np.char.add("control|", st))
    val = halves(grp, rng, 0.2) & ~audit_ctl
    ctl_group = np.full(len(A), -1); idx = rng.permutation(np.where(audit_ctl)[0]); ctl_group[idx] = np.arange(len(idx)) % 5
    assert sorted(D["ids"].astype(str)[audit_ctl].tolist()) == json.load(open(LOOP / "results" / f"audit_controls_{ds}.json"))["audit_controls"]
    return dict(D=D, Z=Z, x=x, A=A, tg=tg, y=D["readout__PDL1"].astype(F), audit_ctl=audit_ctl, train_ctl=ctl & ~audit_ctl, val=val, ctl_group=ctl_group)


def unit(z):
    z = z.astype(np.float64)
    return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)


def bm_score(P, x, y):
    """The head's own evaluation: per protein, the correlation of measured counts with the predicted mean; mean over proteins."""
    h = x
    for name, act in zip(LAYERS, ACT):
        h = h @ P[name + ".weight"].T.astype(np.float64) + P[name + ".bias"]
        h = silu(h) if act else h
    pred = np.logaddexp(0.0, h)
    return float(np.mean([np.corrcoef(y[:, c], pred[:, c])[0, 1] for c in range(y.shape[1]) if y[:, c].std() > 1e-8]))


def write(name, D, Z, lat, out, keep):
    p = OUT_FEATURES / name; p.mkdir(parents=True, exist_ok=True)
    np.savez(p / "dense.npz", **{k: v[keep] for k, v in D.items()})
    for k, X in (("Z", Z), ("HEAD_LAT", lat), ("HEAD_OUT", out)):
        np.save(p / f"node__{k}.npy", X[keep].astype(F))


def main():
    check = "--check" in sys.argv
    W = dict(np.load(FEATURES / "head_weights.npz")); names = [str(n) for n in W["adt_names"]]
    j = names.index("CD274")
    P0 = {l + s: W[l + s].astype(F) for l in LAYERS for s in (".weight", ".bias")}
    hdir = OUT_FEATURES / "heads"; hdir.mkdir(parents=True, exist_ok=True)
    stage = "check" if check else sys.argv[1]
    data = {ds: load(ds) for ds in (TARGETS if stage in ("check", "features") else [sys.argv[2]])}
    for ds, d in data.items():
        d["lat0"], d["out0"], _ = forward(P0, d["x"]); c = d["train_ctl"]; y = d["y"]
        d["t"] = (d["out0"][c, j].mean() + d["out0"][c, j].std() * (y - y[c].mean()) / y[c].std()).astype(F)
    if stage != "features":
        b = np.load(BM_FIT); assert [str(n) for n in b["adt_names"]] == names
        xb = unit(b["z_train"]).astype(F); ob0 = forward(P0, xb)[1]
        bval = np.zeros(len(xb), bool); bval[np.random.default_rng(SEED).permutation(len(xb))[: len(xb) // 5]] = True
        xv, yv = unit(b["z_val"]), b["adt_val"].astype(np.float64); s0 = bm_score(P0, xv, yv)
    if check:
        d = data["frangieh2021"]; sel = np.arange(40); P64 = {k: v.astype(np.float64) for k, v in P0.items()}
        P64["mlp__net.3.weight"] = P64["mlp__net.3.weight"] + 0.01 * np.random.default_rng(1).standard_normal(P64["mlp__net.3.weight"].shape)
        a = (d["x"][sel].astype(np.float64), d["t"][sel].astype(np.float64), d["out0"][sel].astype(np.float64), j, xb[:30].astype(np.float64), ob0[:30].astype(np.float64), 10.0)
        _, G = loss_grad(P64, *a); worst = 0.0
        for key, idx in (("mlp__net.0.weight", (3, 5)), ("mlp__net.6.weight", (7, 2)), ("mlp__net.3.bias", (4,))):
            Pp = {k: v.copy() for k, v in P64.items()}; Pm = {k: v.copy() for k, v in P64.items()}
            Pp[key][idx] += 1e-5; Pm[key][idx] -= 1e-5
            num = (loss_grad(Pp, *a, grad=False)[0] - loss_grad(Pm, *a, grad=False)[0]) / 2e-5
            worst = max(worst, abs(num - G[key][idx]) / (abs(num) + 1e-12))
        print("gradient against finite differences, worst relative difference:", f"{worst:.2e}")
        train = d["train_ctl"] | ((d["A"] == 1) & (d["tg"] != "IFNGR1")); t0 = time.time()
        P, curve = fit(P0, d["x"], d["t"], d["out0"], j, train & ~d["val"], train & d["val"], xb, ob0, bval, 10.0, 20)
        print(f"20 steps: {time.time() - t0:.0f}s; set-aside loss {curve[0][2]:.4f}; validation score as trained {s0:.4f}, after {bm_score(P, xv, yv):.4f}")
        return
    if stage in ("pilot", "fit"):
        ds = sys.argv[2]; d = data[ds]; x, A, tg = d["x"], d["A"], d["tg"]

        def one(held, lam):
            train = d["train_ctl"] | ((A == 1) & (tg != held)); t0 = time.time()
            P, curve = fit(P0, x, d["t"], d["out0"], j, train & ~d["val"], train & d["val"], xb, ob0, bval, lam)
            ob = forward(P, xb)[1]
            r = {"held_out": held, "weight": lam, "n_fit": int((train & ~d["val"]).sum()), "n_bone_marrow_fit": int((~bval).sum()), "curve": curve,
                 "kept_step": int(min(curve, key=lambda c: c[2])[0]), "seconds": time.time() - t0,
                 "bone_marrow_outputs_moved": float(np.sqrt(((ob - ob0) ** 2).mean())), "validation_score": bm_score(P, xv, yv), "validation_score_as_trained": s0}
            print(f"{ds}|{held}|weight {lam:g}: kept step {r['kept_step']}, set-aside loss {curve[0][2]:.4f} to {min(c[2] for c in curve):.4f}, "
                  f"validation score {r['validation_score']:.4f} (as trained {s0:.4f}), {r['seconds']:.0f}s", flush=True)
            return P, r

        first = TARGETS[ds][0]
        if stage == "pilot":
            rows = []
            for lam in GRID:
                P, r = one(first, lam); rows.append(r)
                np.savez(hdir / f"{ds}__{first}__weight{lam:g}.npz", **P)
            ok = [r["weight"] for r in rows if s0 - r["validation_score"] <= VAL_DROP]
            chosen = min(ok) if ok else max(GRID)
            json.dump({"grid": list(GRID), "allowed_fall": VAL_DROP, "rows": rows, "chosen_weight": chosen, "met_on_validation": bool(ok)}, open(hdir / f"pilot_{ds}.json", "w"))
            print(ds, "chosen weight", chosen, "| met on validation cells:", bool(ok), flush=True)
        else:
            pilot = json.load(open(hdir / f"pilot_{ds}.json")); lam = pilot["chosen_weight"]; fits = {}
            shutil.copy(hdir / f"{ds}__{first}__weight{lam:g}.npz", hdir / f"{ds}__{first}__{ARM}.npz")
            fits[first] = next(r for r in pilot["rows"] if r["weight"] == lam)
            for held in TARGETS[ds][1:]:
                P, fits[held] = one(held, lam)
                np.savez(hdir / f"{ds}__{held}__{ARM}.npz", **P)
            json.dump(fits, open(hdir / f"fits_{ds}.json", "w"))
        return
    log = {"learning_rate": LR, "max_steps": STEPS, "patience_checks": PATIENCE, "keep_weight": KEEP, "split_seed": SEED, "pdl1_output": j, "pilot": {}, "fits": {}, "side": {}}
    for ds, genes in TARGETS.items():
        d = data[ds]; x, A, tg, y, D, Z = d["x"], d["A"], d["tg"], d["y"], d["D"], d["Z"]
        log["pilot"][ds] = json.load(open(hdir / f"pilot_{ds}.json")); log["fits"][ds] = json.load(open(hdir / f"fits_{ds}.json"))
        lat_p, out_p = np.zeros_like(d["lat0"]), np.zeros_like(d["out0"]); side = {}
        for k, g in enumerate(genes):
            lat, out, _ = forward({kk: v.astype(F) for kk, v in np.load(hdir / f"{ds}__{g}__{ARM}.npz").items()}, x)
            T = (A == 1) & (tg == g); C = d["audit_ctl"]
            write(f"target_{ds}_{g}_{ARM}", D, Z, lat, out, T | C)
            own = T | (d["ctl_group"] == k)
            lat_p[own], out_p[own] = lat[own], out[own]
            side[g] = {"measured_shift_sd": float((y[T].mean() - y[C].mean()) / y[C].std()), "own_output_shift_sd": float((out[T, j].mean() - out[C, j].mean()) / out[C, j].std()),
                       "other_outputs_moved": float(np.sqrt(((np.delete(out, j, 1) - np.delete(d["out0"], j, 1))[T | C] ** 2).mean()))}
        write(f"pooled_{ds}_{ARM}", D, Z, lat_p, out_p, (A == 1) | d["audit_ctl"])
        log["side"][ds] = side
    json.dump(log, open(RES / "fits.json", "w"), indent=1)
    print("features written for", sum(1 for _ in OUT_FEATURES.glob("*/dense.npz")), "audits")


main()
