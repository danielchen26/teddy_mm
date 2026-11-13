"""Round 3 of the loop: one repair of the head's MLP for both cell types and every measured protein the head predicts,
with the head's outputs held on bone-marrow training cells (DESIGN.md).

  python repair.py pilot             seed 0, first gene held out: one fit per weight of the grid, smallest first, until the
                                     bone-marrow validation score is kept; writes heads/pilot.json
  python repair.py fit <seed> <k>    the held-out fits of one split; k = 0 or 1 takes every second gene, so two can run at once
  python repair.py final             seed 0: the head fitted on all genes (the head to use)
  python repair.py features <seed>   the features of that split's audits; with seed 0 also results/fits.json
  python repair.py --check           gradient against finite differences on a few cells; one 20-step fit; writes nothing
"""
import json
import sys
import time

import numpy as np
from scipy.special import expit

from _paths import ARM, BM_FIT, FEATURES, GENES, LOOP, OUT_FEATURES, RES, TARGETS, TRAINED

MLP = ["mlp__net.0", "mlp__net.3", "mlp__net.6"]
DEC = ["dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"]
LAYERS, ACT = MLP + DEC, [True, True, False, True, True, False]
LR, STEPS, EVERY, PATIENCE, KEEP = 3e-4, 400, 20, 3, 1.0               # as in the earlier rounds
GRID, VAL_DROP = (1.0, 10.0, 100.0), 0.01                              # weights of the holding term; allowed fall on validation cells
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


def weights(M, ds_of):
    """Per-entry weights. Trained outputs: each (data set, protein) pair is one mean squared error, all pairs summed.
    The other outputs: KEEP times their mean squared change, per data set."""
    Wt, Wk = np.zeros(M.shape, F), np.zeros(M.shape, F)
    for d in np.unique(ds_of):
        r = ds_of == d
        Wt[r] = M[r] / r.sum()
        Wk[r] = (~M[r]) * (KEEP / (r.sum() * max((~M[r][0]).sum(), 1)))
    return Wt, Wk


def loss_grad(P, x, T, out0, Wt, Wk, xb, ob0, lam, grad=True):
    _, out, cache = forward(P, x)
    _, ob, cache_b = forward(P, xb)
    L = float((Wt * (out - T) ** 2).sum() + (Wk * (out - out0) ** 2).sum() + lam * ((ob - ob0) ** 2).mean())
    if not grad:
        return L, None
    G = backward(P, cache, 2 * (Wt * (out - T) + Wk * (out - out0)))
    Gb = backward(P, cache_b, (2 * lam / ob.size) * (ob - ob0))
    return L, {key: G[key] + Gb[key] for key in G}


def fit(P0, C, use, val, xb, ob0, bval, lam, steps=STEPS):
    """Adam on the three MLP layers only; the weights kept are those with the lowest loss on the set-aside cells."""
    keys = [l + s for l in MLP for s in (".weight", ".bias")]
    P = {k: v.copy() for k, v in P0.items()}
    m = {k: np.zeros_like(P[k]) for k in keys}; v2 = {k: np.zeros_like(P[k]) for k in keys}
    a = (C["x"][use], C["T"][use], C["out0"][use], *weights(C["M"][use], C["ds"][use]), xb[~bval], ob0[~bval], lam)
    b = (C["x"][val], C["T"][val], C["out0"][val], *weights(C["M"][val], C["ds"][val]), xb[bval], ob0[bval], lam)
    best, best_P, curve, bad = np.inf, None, [], 0
    for step in range(1, steps + 1):
        L, G = loss_grad(P, *a)
        for k in keys:
            m[k] = 0.9 * m[k] + 0.1 * G[k]; v2[k] = 0.999 * v2[k] + 0.001 * G[k] ** 2
            P[k] = (P[k] - LR * (m[k] / (1 - 0.9 ** step)) / (np.sqrt(v2[k] / (1 - 0.999 ** step)) + 1e-8)).astype(F)
        if step % EVERY == 0:
            Lv = loss_grad(P, *b, grad=False)[0]
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


def load(ds, seed):
    """The cells and the split code of the earlier rounds. With seed 0 the audit controls are checked against their record."""
    d = np.load(FEATURES / ds / "dense.npz"); D = {k: d[k] for k in d.files}
    Z = (D.pop("node__Z") if "node__Z" in D else np.load(FEATURES / ds / "node__Z.npy")).astype(F)
    D = {k: v for k, v in D.items() if not k.startswith("node__")}
    x = (Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-6)).astype(F)
    A, tg = D["A"], D["target"].astype(str)
    st = D["stratum"].astype(str) if "stratum" in D else np.full(len(A), "all")
    rng = np.random.default_rng(seed)
    ctl = A == 0
    audit_ctl = ctl & halves(np.where(ctl, st, "treated"), rng)
    grp = np.where(A == 1, np.char.add(np.char.add(D["guide"].astype(str), "|"), st), np.char.add("control|", st))
    val = halves(grp, rng, 0.2) & ~audit_ctl
    ctl_group = np.full(len(A), -1); idx = rng.permutation(np.where(audit_ctl)[0]); ctl_group[idx] = np.arange(len(idx)) % 5
    ids = sorted(D["ids"].astype(str)[audit_ctl].tolist())
    if seed == 0:
        assert ids == json.load(open(LOOP / "results" / f"audit_controls_{ds}.json"))["audit_controls"]
    return dict(D=D, Z=Z, x=x, A=A, tg=tg, audit_ctl=audit_ctl, train_ctl=ctl & ~audit_ctl, val=val, ctl_group=ctl_group, audit_ids=ids)


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


def setup(seed, names, P0):
    """Both data sets with their targets, and one table of all their cells for the fit."""
    data = {ds: load(ds, seed) for ds in TARGETS}
    rows = {k: [] for k in ("x", "T", "M", "out0", "ds", "train", "val", "A", "tg")}
    for i, (ds, d) in enumerate(data.items()):
        d["lat0"], d["out0"], _ = forward(P0, d["x"]); c = d["train_ctl"]
        T, M = d["out0"].copy(), np.zeros(d["out0"].shape, bool)
        for rid, prot in TRAINED[ds].items():
            j = names.index(prot); y = d["D"]["readout__" + rid].astype(F)
            T[:, j] = d["out0"][c, j].mean() + d["out0"][c, j].std() * (y - y[c].mean()) / y[c].std(); M[:, j] = True
        for k, v in (("x", d["x"]), ("T", T.astype(F)), ("M", M), ("out0", d["out0"]), ("ds", np.full(len(c), i)), ("train", c), ("val", d["val"]), ("A", d["A"]), ("tg", d["tg"])):
            rows[k].append(v)
    return data, {k: np.concatenate(v) for k, v in rows.items()}


def main():
    check = "--check" in sys.argv
    stage = "check" if check else sys.argv[1]
    seed = int(sys.argv[2]) if stage in ("fit", "features") else 0
    W = dict(np.load(FEATURES / "head_weights.npz")); names = [str(n) for n in W["adt_names"]]
    P0 = {l + s: W[l + s].astype(F) for l in LAYERS for s in (".weight", ".bias")}
    hdir = OUT_FEATURES / "heads"; hdir.mkdir(parents=True, exist_ok=True)
    data, C = setup(seed, names, P0)
    if stage != "features":
        b = np.load(BM_FIT); assert [str(n) for n in b["adt_names"]] == names
        xb = unit(b["z_train"]).astype(F); ob0 = forward(P0, xb)[1]
        bval = np.zeros(len(xb), bool); bval[np.random.default_rng(seed).permutation(len(xb))[: len(xb) // 5]] = True
        xv, yv = unit(b["z_val"]), b["adt_val"].astype(np.float64); s0 = bm_score(P0, xv, yv)

        def one(held, lam, steps=STEPS):
            train = C["train"] | ((C["A"] == 1) & (C["tg"] != held)) if held != "all" else C["train"] | (C["A"] == 1)
            t0 = time.time()
            P, curve = fit(P0, C, train & ~C["val"], train & C["val"], xb, ob0, bval, lam, steps)
            ob = forward(P, xb)[1]
            r = {"held_out": held, "seed": seed, "weight": lam, "n_fit": int((train & ~C["val"]).sum()), "n_bone_marrow_fit": int((~bval).sum()), "curve": curve,
                 "kept_step": int(min(curve, key=lambda c: c[2])[0]), "seconds": time.time() - t0,
                 "bone_marrow_outputs_moved": float(np.sqrt(((ob - ob0) ** 2).mean())), "validation_score": bm_score(P, xv, yv), "validation_score_as_trained": s0}
            print(f"seed {seed}|{held}|weight {lam:g}: kept step {r['kept_step']}, set-aside loss {curve[0][2]:.4f} to {min(c[2] for c in curve):.4f}, "
                  f"validation score {r['validation_score']:.4f} (as trained {s0:.4f}), {r['seconds']:.0f}s", flush=True)
            return P, r

    if check:
        sel = np.r_[0:25, len(data["frangieh2021"]["A"]):len(data["frangieh2021"]["A"]) + 25]
        P64 = {k: v.astype(np.float64) for k, v in P0.items()}
        P64["mlp__net.3.weight"] = P64["mlp__net.3.weight"] + 0.01 * np.random.default_rng(1).standard_normal(P64["mlp__net.3.weight"].shape)
        Wt, Wk = weights(C["M"][sel], C["ds"][sel])
        a = (C["x"][sel].astype(np.float64), C["T"][sel].astype(np.float64), C["out0"][sel].astype(np.float64), Wt.astype(np.float64), Wk.astype(np.float64), xb[:30].astype(np.float64), ob0[:30].astype(np.float64), 10.0)
        _, G = loss_grad(P64, *a); worst = 0.0
        for key, idx in (("mlp__net.0.weight", (3, 5)), ("mlp__net.6.weight", (7, 2)), ("mlp__net.3.bias", (4,))):
            Pp = {k: v.copy() for k, v in P64.items()}; Pm = {k: v.copy() for k, v in P64.items()}
            Pp[key][idx] += 1e-5; Pm[key][idx] -= 1e-5
            num = (loss_grad(Pp, *a, grad=False)[0] - loss_grad(Pm, *a, grad=False)[0]) / 2e-5
            worst = max(worst, abs(num - G[key][idx]) / (abs(num) + 1e-12))
        print("gradient against finite differences, worst relative difference:", f"{worst:.2e}")
        print("trained entries per data set:", {ds: int(C["M"][C["ds"] == i][0].sum()) for i, ds in enumerate(TARGETS)}, "| cells:", {ds: int((C["ds"] == i).sum()) for i, ds in enumerate(TARGETS)})
        one(GENES[0], 1.0, 20)
        return
    if stage == "pilot":
        rows, chosen = [], None
        for lam in GRID:
            P, r = one(GENES[0], lam); rows.append(r)
            np.savez(hdir / f"s0__{GENES[0]}__weight{lam:g}.npz", **P)
            if s0 - r["validation_score"] <= VAL_DROP:
                chosen = lam; break
        met = chosen is not None; chosen = chosen if met else max(GRID)
        np.savez(hdir / f"s0__{GENES[0]}.npz", **dict(np.load(hdir / f"s0__{GENES[0]}__weight{chosen:g}.npz")))
        json.dump({"grid": list(GRID), "allowed_fall": VAL_DROP, "rows": rows, "chosen_weight": chosen, "met_on_validation": met}, open(hdir / "pilot.json", "w"))
        json.dump({GENES[0]: next(r for r in rows if r["weight"] == chosen)}, open(hdir / "fits_s0_pilot.json", "w"))
        print("chosen weight", chosen, "| met on validation cells:", met, flush=True)
    elif stage == "fit":
        k = int(sys.argv[3]); lam = json.load(open(hdir / "pilot.json"))["chosen_weight"]; fits = {}
        for held in [g for i, g in enumerate(GENES) if i % 2 == k and not (seed == 0 and g == GENES[0])]:
            P, fits[held] = one(held, lam)
            np.savez(hdir / f"s{seed}__{held}.npz", **P)
        json.dump(fits, open(hdir / f"fits_s{seed}_{k}.json", "w"))
    elif stage == "final":
        lam = json.load(open(hdir / "pilot.json"))["chosen_weight"]
        P, r = one("all", lam)
        np.savez(hdir / "s0__all.npz", **P); json.dump({"all": r}, open(hdir / "fits_s0_final.json", "w"))
    else:
        side = {}
        for ds, genes in TARGETS.items():
            d = data[ds]; x, A, tg, D, Z = d["x"], d["A"], d["tg"], d["D"], d["Z"]
            lat_p, out_p = np.zeros_like(d["lat0"]), np.zeros_like(d["out0"]); side[ds] = {}
            for k, g in enumerate(genes):
                lat, out, _ = forward({kk: v.astype(F) for kk, v in np.load(hdir / f"s{seed}__{g}.npz").items()}, x)
                T, Cc = (A == 1) & (tg == g), d["audit_ctl"]
                if seed == 0:
                    write(f"target_{ds}_{g}_{ARM}_s0", D, Z, lat, out, T | Cc)
                own = T | (d["ctl_group"] == k)
                lat_p[own], out_p[own] = lat[own], out[own]
                side[ds][g] = {rid: {"measured_shift_sd": float((D["readout__" + rid][T].mean() - D["readout__" + rid][Cc].mean()) / D["readout__" + rid][Cc].std()),
                                     "own_output_shift_sd": float((out[T, names.index(p)].mean() - out[Cc, names.index(p)].mean()) / out[Cc, names.index(p)].std())} for rid, p in TRAINED[ds].items()}
                side[ds][g]["untrained_outputs_moved"] = float(np.sqrt((((out - d["out0"])[T | Cc][:, [j for j in range(out.shape[1]) if names[j] not in TRAINED[ds].values()]]) ** 2).mean()))
            write(f"pooled_{ds}_{ARM}_s{seed}", D, Z, lat_p, out_p, (A == 1) | d["audit_ctl"])
            json.dump({"split_seed": seed, "audit_controls": d["audit_ids"]}, open(RES / f"audit_controls_{ds}_s{seed}.json", "w"))
        fits = {}
        for f in sorted(hdir.glob(f"fits_s{seed}_*.json")):
            fits.update(json.load(open(f)))
        log = {"split_seed": seed, "learning_rate": LR, "max_steps": STEPS, "patience_checks": PATIENCE, "keep_weight": KEEP, "trained": TRAINED, "fits": fits, "side": side}
        if seed == 0:
            log["pilot"] = json.load(open(hdir / "pilot.json"))
        json.dump(log, open(RES / f"fits_s{seed}.json", "w"), indent=1)
        print("features written for seed", seed)


main()
