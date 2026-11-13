"""Retrain one part of the head as DESIGN.md says and write the features every audit needs.

  python fit_heads.py fit <data set>   the 12 fits of one data set ((5 held-out targets + all targets) x 2 parts); heads are saved
  python fit_heads.py features         after both data sets are fitted: the features of every audit, and results/fits.json
  python fit_heads.py --check          one short fit per part on training cells: timing and the frozen-part checks; writes nothing
"""
import json
import sys
import time

import numpy as np
from scipy.special import expit

from _paths import ARMS, FEATURES, OTHER, OUT_FEATURES, RES, TARGETS

MLP = ["mlp__net.0", "mlp__net.3", "mlp__net.6"]
DEC = ["dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"]
LAYERS, ACT = MLP + DEC, [True, True, False, True, True, False]
PART = {"mlp_only": MLP, "decoder_only": DEC}
LR, STEPS, EVERY, PATIENCE, KEEP, SEED = 3e-4, 400, 20, 3, 1.0, 0
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


def loss_grad(P, x, t, out0, j, grad=True):
    """Squared error of the PD-L1 output against its target, plus KEEP times the mean squared change of the other outputs."""
    _, out, (ins, pre, m, sp) = forward(P, x)
    n, k = out.shape
    other = np.ones(k, bool); other[j] = False
    L = float(((out[:, j] - t) ** 2).mean() + KEEP * ((out[:, other] - out0[:, other]) ** 2).mean())
    if not grad:
        return L, None
    g = (2 * KEEP / (n * (k - 1))) * (out - out0)
    g[:, j] = 2 * (out[:, j] - t) / n
    g = (g * np.where(sp > 1e-6, expit(m) / sp, 0)).astype(F)
    G = {}
    for i in range(len(LAYERS) - 1, -1, -1):
        if ACT[i]:
            g = g * dsilu(pre[i])
        G[LAYERS[i] + ".weight"], G[LAYERS[i] + ".bias"] = g.T @ ins[i], g.sum(0)
        g = g @ P[LAYERS[i] + ".weight"]
    return L, G


def fit(P0, x, t, out0, j, use, val, part, steps=STEPS):
    """Adam on the layers of one part only; the weights kept are those with the lowest loss on the set-aside cells."""
    keys = [l + s for l in PART[part] for s in (".weight", ".bias")]
    P = {k: v.copy() for k, v in P0.items()}
    m = {k: np.zeros_like(P[k]) for k in keys}; v2 = {k: np.zeros_like(P[k]) for k in keys}
    best, best_P, curve, bad = np.inf, None, [], 0
    for step in range(1, steps + 1):
        L, G = loss_grad(P, x[use], t[use], out0[use], j)
        for k in keys:
            m[k] = 0.9 * m[k] + 0.1 * G[k]; v2[k] = 0.999 * v2[k] + 0.001 * G[k] ** 2
            P[k] = (P[k] - LR * (m[k] / (1 - 0.9 ** step)) / (np.sqrt(v2[k] / (1 - 0.999 ** step)) + 1e-8)).astype(F)
        if step % EVERY == 0:
            Lv = loss_grad(P, x[val], t[val], out0[val], j, grad=False)[0]
            curve.append([step, L, Lv])
            if Lv < best:
                best, best_P, bad = Lv, {k: v.copy() for k, v in P.items()}, 0
            else:
                bad += 1
                if bad >= PATIENCE:
                    break
    frozen = [k for k in P0 if k not in keys]
    assert all(np.array_equal(best_P[k], P0[k]) for k in frozen)          # the other part is bit-for-bit the head as trained
    return best_P, curve


def halves(groups, rng, share=0.5):
    first = np.zeros(len(groups), bool)
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        first[idx[: int(round(share * len(idx)))]] = True
    return first


def load(ds):
    d = np.load(FEATURES / ds / "dense.npz"); D = {k: d[k] for k in d.files}
    Z = (D.pop("node__Z") if "node__Z" in D else np.load(FEATURES / ds / "node__Z.npy")).astype(F)
    D = {k: v for k, v in D.items() if not k.startswith("node__")}
    x = (Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-6)).astype(F)
    A, tg = D["A"], D["target"].astype(str)
    st = D["stratum"].astype(str) if "stratum" in D else np.full(len(A), "all")
    rng = np.random.default_rng(SEED)
    ctl = A == 0
    audit_ctl = ctl & halves(np.where(ctl, st, "treated"), rng)            # control cells no head is ever fitted on
    grp = np.where(A == 1, np.char.add(np.char.add(D["guide"].astype(str), "|"), st), np.char.add("control|", st))
    val = halves(grp, rng, 0.2) & ~audit_ctl                               # set-aside fifth, used only for stopping
    ctl_group = np.full(len(A), -1); idx = rng.permutation(np.where(audit_ctl)[0]); ctl_group[idx] = np.arange(len(idx)) % 5
    return dict(D=D, Z=Z, x=x, A=A, tg=tg, y=D["readout__PDL1"].astype(F), audit_ctl=audit_ctl, train_ctl=ctl & ~audit_ctl, val=val, ctl_group=ctl_group)


def write(name, D, Z, lat, out, keep):
    p = OUT_FEATURES / name; p.mkdir(parents=True, exist_ok=True)
    np.savez(p / "dense.npz", **{k: v[keep] for k, v in D.items()})
    for k, X in (("Z", Z), ("HEAD_LAT", lat), ("HEAD_OUT", out)):
        np.save(p / f"node__{k}.npy", X[keep].astype(F))


def main():
    check = "--check" in sys.argv
    W = dict(np.load(FEATURES / "head_weights.npz"))
    j = [str(n) for n in W["adt_names"]].index("CD274")
    P0 = {l + s: W[l + s].astype(F) for l in LAYERS for s in (".weight", ".bias")}
    todo = list(TARGETS) if check or sys.argv[1] == "features" else [sys.argv[2]]
    data = {ds: load(ds) for ds in TARGETS}
    for ds in TARGETS:
        d = data[ds]; d["lat0"], d["out0"], _ = forward(P0, d["x"]); c = d["train_ctl"]; y = d["y"]
        d["t"] = (d["out0"][c, j].mean() + d["out0"][c, j].std() * (y - y[c].mean()) / y[c].std()).astype(F)
    hdir = OUT_FEATURES / "heads"; hdir.mkdir(parents=True, exist_ok=True)
    if check or sys.argv[1] == "fit":
        for ds in todo:
            d = data[ds]; x, A, tg = d["x"], d["A"], d["tg"]; fits = {}
            for held in TARGETS[ds] + [None]:
                train = d["train_ctl"] | ((A == 1) & (tg != held))
                for part in PART:
                    t0 = time.time(); tag = f"{ds}|{held or 'all'}|{part}"
                    P, curve = fit(P0, x, d["t"], d["out0"], j, train & ~d["val"], train & d["val"], part, 20 if check else STEPS)
                    fits[tag] = {"n_fit": int((train & ~d["val"]).sum()), "n_stop": int((train & d["val"]).sum()), "curve": curve,
                                 "kept_step": int(min(curve, key=lambda r: r[2])[0]), "seconds": time.time() - t0}
                    print(f"{tag}: {fits[tag]['n_fit']} cells, kept step {fits[tag]['kept_step']}, set-aside loss {curve[0][2]:.4f} to {min(r[2] for r in curve):.4f}, {time.time() - t0:.0f}s", flush=True)
                    if check:
                        print("   latent unchanged:", bool(np.array_equal(forward(P, x[:200])[0], d["lat0"][:200])), "(expected for decoder_only only)")
                    else:
                        np.savez(hdir / (tag.replace("|", "__") + ".npz"), **P)
                if check:
                    break
            if not check:
                json.dump(fits, open(hdir / f"fits_{ds}.json", "w"))
        return
    heads = {(ds, held, part): dict(np.load(hdir / f"{ds}__{held or 'all'}__{part}.npz")) for ds in TARGETS for held in TARGETS[ds] + [None] for part in PART}
    log = {"learning_rate": LR, "max_steps": STEPS, "patience_checks": PATIENCE, "keep_weight": KEEP, "split_seed": SEED, "pdl1_output": j, "fits": {}, "side": {}}
    for ds in TARGETS:
        log["fits"].update(json.load(open(hdir / f"fits_{ds}.json")))
    for ds, genes in TARGETS.items():
        d = data[ds]; x, A, tg, y, D, Z = d["x"], d["A"], d["tg"], d["y"], d["D"], d["Z"]
        for arm in ARMS:
            head = lambda held: P0 if arm == "as_trained" else heads[ds, held, arm]
            lat_p, out_p = np.zeros_like(d["lat0"]), np.zeros_like(d["out0"]); side = {}
            for k, g in enumerate(genes):
                lat, out, _ = forward(head(g), x)
                T = (A == 1) & (tg == g); C = d["audit_ctl"]
                write(f"target_{ds}_{g}_{arm}", D, Z, lat, out, T | C)                       # one head: its held-out target and all audit controls
                own = T | (d["ctl_group"] == k)                                              # pooled audit: every cell through a head not fitted on it
                lat_p[own], out_p[own] = lat[own], out[own]
                side[g] = {"measured_shift_sd": float((y[T].mean() - y[C].mean()) / y[C].std()), "own_output_shift_sd": float((out[T, j].mean() - out[C, j].mean()) / out[C, j].std()),
                           "other_outputs_moved": float(np.sqrt(((np.delete(out, j, 1) - np.delete(d["out0"], j, 1))[T | C] ** 2).mean())),
                           "latent_identical_to_as_trained": bool(np.array_equal(lat, d["lat0"]))}
            write(f"pooled_{ds}_{arm}", D, Z, lat_p, out_p, (A == 1) | d["audit_ctl"])
            log["side"][f"{ds}|{arm}"] = side
            if arm != "as_trained":                                                          # a head fitted on the other data set, on all cells of this one
                lat, out, _ = forward(heads[OTHER[ds], None, arm], x)
                write(f"transfer_{ds}_{arm}", D, Z, lat, out, np.ones(len(A), bool))
        json.dump({"split_seed": SEED, "audit_controls": sorted(D["ids"].astype(str)[d["audit_ctl"]].tolist())}, open(RES / f"audit_controls_{ds}.json", "w"))
    json.dump(log, open(RES / "fits.json", "w"), indent=1)
    print("features written for", sum(1 for _ in OUT_FEATURES.glob("*/dense.npz")), "audits")


main()
