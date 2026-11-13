"""Fine-tune the head as DESIGN.md says and write, for each of the three heads, the features of the audit half.

  python finetune_head.py            split seeds 0 to 2; features are written for seed 0
  python finetune_head.py --check    gradient check and twenty training steps on the training half, nothing written
"""
import json
import sys

import numpy as np
from scipy.special import expit

from _paths import FEATURES, OUT_FEATURES, RES

LAYERS = ["mlp__net.0", "mlp__net.3", "mlp__net.6", "dec__mu.net.0", "dec__mu.net.3", "dec__mu.net.6"]
ACT = [True, True, False, True, True, False]
UPSTREAM = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
TRAIN_ON = {"without_IRF1": UPSTREAM, "with_IRF1": UPSTREAM + ["IRF1"]}
LR, STEPS, EVERY, KEEP = 3e-4, 200, 20, 1.0

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
    sp = np.logaddexp(0.0, h)
    return lat, np.log(np.maximum(sp, 1e-6)), (ins, pre, h, sp)


def loss_grad(P, x, t, out0, j, grad=True):
    _, out, (ins, pre, m, sp) = forward(P, x)
    n, k = out.shape
    other = np.ones(k, bool); other[j] = False
    L = float(((out[:, j] - t) ** 2).mean() + KEEP * ((out[:, other] - out0[:, other]) ** 2).mean())
    if not grad:
        return L, None
    g = 2 * KEEP * (out - out0) / (n * (k - 1))
    g[:, j] = 2 * (out[:, j] - t) / n
    g = g * np.where(sp > 1e-6, expit(m) / sp, 0.0)
    G = {}
    for i in range(len(LAYERS) - 1, -1, -1):
        if ACT[i]:
            g = g * dsilu(pre[i])
        G[LAYERS[i] + ".weight"], G[LAYERS[i] + ".bias"] = g.T @ ins[i], g.sum(0)
        g = g @ P[LAYERS[i] + ".weight"]
    return L, G


def halves(groups, rng, share=0.5):
    """True for the cells that go to the first part: within every group, a random share of its cells."""
    first = np.zeros(len(groups), bool)
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        first[idx[: int(round(share * len(idx)))]] = True
    return first


def finetune(P0, x, t, out0, j, fit, val, check=False):
    P = {k: v.copy() for k, v in P0.items()}
    m = {k: np.zeros_like(v) for k, v in P.items()}; v2 = {k: np.zeros_like(v) for k, v in P.items()}
    best, best_P, curve = np.inf, None, []
    for step in range(1, (20 if check else STEPS) + 1):
        L, G = loss_grad(P, x[fit], t[fit], out0[fit], j)
        for k in P:
            m[k] = 0.9 * m[k] + 0.1 * G[k]; v2[k] = 0.999 * v2[k] + 0.001 * G[k] ** 2
            P[k] -= LR * (m[k] / (1 - 0.9 ** step)) / (np.sqrt(v2[k] / (1 - 0.999 ** step)) + 1e-8)
        if step % EVERY == 0:
            Lv = loss_grad(P, x[val], t[val], out0[val], j, grad=False)[0]
            curve.append([step, L, Lv])
            if Lv < best:
                best, best_P = Lv, {k: v.copy() for k, v in P.items()}
    return best_P, curve


def main():
    check = "--check" in sys.argv
    W = dict(np.load(FEATURES / "head_weights.npz"))
    j = [str(n) for n in W["adt_names"]].index("CD274")
    P0 = {k: W[k].astype(np.float64) for name in LAYERS for k in (name + ".weight", name + ".bias")}
    d = np.load(FEATURES / "papalexi2021" / "dense.npz"); D = {k: d[k] for k in d.files}
    Z = np.load(FEATURES / "papalexi2021" / "node__Z.npy").astype(np.float64)
    x = Z / (np.linalg.norm(Z, axis=1, keepdims=True) + 1e-6)
    A, tg, guide, st, y = D["A"], D["target"].astype(str), D["guide"].astype(str), D["stratum"].astype(str), D["readout__PDL1"]
    groups = np.array([f"{a}|{b}|{c}" for a, b, c in zip(np.where(A == 1, guide, "control"), tg, st)])
    lat0, out0, _ = forward(P0, x)
    if check:                                           # the hand-written gradient against finite differences
        rng = np.random.default_rng(0); s = rng.choice(len(x), 64, replace=False); t = rng.standard_normal(64)
        _, G = loss_grad(P0, x[s], t, out0[s] + 0.1, j)
        for name in ("mlp__net.0.weight", "dec__mu.net.3.bias", "dec__mu.net.6.weight"):
            idx = tuple(rng.integers(0, n) for n in P0[name].shape); P = {k: v.copy() for k, v in P0.items()}; e = 1e-5
            P[name][idx] += e; up = loss_grad(P, x[s], t, out0[s] + 0.1, j, grad=False)[0]
            P[name][idx] -= 2 * e; dn = loss_grad(P, x[s], t, out0[s] + 0.1, j, grad=False)[0]
            print("gradient", name, f"{G[name][idx]:.6e} against {(up - dn) / (2 * e):.6e}")
    log = {"learning_rate": LR, "max_steps": STEPS, "keep_weight": KEEP, "pdl1_output": j, "seeds": {}}
    for seed in range(1 if check else 3):
        rng = np.random.default_rng(seed)
        train = halves(groups, rng)                     # training half; the rest is the audit half
        val = train & halves(np.where(train, groups, "audit"), rng, 0.2)
        ctl = train & (A == 0)
        target = out0[ctl, j].mean() + out0[ctl, j].std() * (y - y[ctl].mean()) / y[ctl].std()
        res = {"n_train": int(train.sum()), "n_audit": int((~train).sum()), "arms": {}}
        heads = {"as_trained": P0}
        for arm, genes in TRAIN_ON.items():
            use = train & ((A == 0) | np.isin(tg, genes))
            heads[arm], curve = finetune(P0, x, target, out0, j, use & ~val, use & val, check)
            res["arms"][arm] = {"n_fit": int((use & ~val).sum()), "n_stop": int((use & val).sum()), "curve": curve,
                                "kept_step": int(min(curve, key=lambda c: c[2])[0])}
            print(f"seed {seed} {arm}: fitted on {res['arms'][arm]['n_fit']} cells, kept step {res['arms'][arm]['kept_step']}, "
                  f"loss on training cells {curve[0][1]:.4f} to {curve[-1][1]:.4f}, on the set-aside fifth {curve[0][2]:.4f} to {min(c[2] for c in curve):.4f}", flush=True)
        if check:
            return
        au = ~train; C = au & (A == 0); res["heads"] = {}
        for arm, P in heads.items():
            lat, out, _ = forward(P, x)
            r = {}
            for g in UPSTREAM + ["IRF1"]:
                T = au & (A == 1) & (tg == g)
                r[g] = {"measured_shift_sd": float((y[T].mean() - y[C].mean()) / y[C].std()),
                        "own_output_shift_sd": float((out[T, j].mean() - out[C, j].mean()) / out[C, j].std()),
                        "gap_over_spread": {k: float(np.linalg.norm(X[T].mean(0) - X[C].mean(0)) / np.sqrt(X[C].var(0).sum())) for k, X in (("HEAD_LAT", lat), ("HEAD_OUT", out))},
                        "other_outputs_moved": float(np.sqrt(((np.delete(out, j, 1) - np.delete(out0, j, 1))[au] ** 2).mean()))}
            res["heads"][arm] = r
            if seed == 0:
                (OUT_FEATURES / arm).mkdir(parents=True, exist_ok=True)
                np.savez(OUT_FEATURES / arm / "dense.npz", **{k: v[au] for k, v in D.items()})
                for k, X in (("Z", Z), ("HEAD_LAT", lat), ("HEAD_OUT", out)):
                    np.save(OUT_FEATURES / arm / f"node__{k}.npy", X[au].astype(np.float32))
                np.savez(OUT_FEATURES / arm / "head_weights.npz", **P)
        if seed == 0:
            json.dump({"split_seed": 0, "audit_ids": sorted(D["ids"].astype(str)[au].tolist())}, open(RES / "audit_half_seed0.json", "w"))
        log["seeds"][str(seed)] = res
        json.dump(log, open(RES / "finetune.json", "w"), indent=1)


main()
