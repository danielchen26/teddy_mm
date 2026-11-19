"""Development, before DESIGN.md was written. Nothing here is a test.

Question: can the ridge weight be chosen from the rows themselves? For each fold of the training controls and each row of the
fit: fit on the other folds, update with the other rows, read the left-out event against the held-out controls; take the weight
with the smallest mean |u| (`rows_x`). Against the weight by cross-validation on the controls (`rows`). Cells already read only
(both data sets, one gene left out at a time, the split of the gate).

  python develop_weight.py     (environment of test.py; writes results/develop_weight.json)
"""
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(RUN.parent / "teg_repair_gate_20261119" / "scripts"))
import gate as G  # noqa: E402
from gate import BM_FIT, forward_np, interval, load, own_output_u, reading, unit, weights  # noqa: E402

GRID = [10.0 ** (k / 2) for k in range(-6, 5)]


def run(ds):
    P0, names, j = weights()
    xb = unit(np.load(BM_FIT)["z_train"]).astype(np.float64)
    d = load(ds); genes = G.TARGETS[ds]
    x, y, A, tg, C = d["x"].astype(np.float64), d["y"], d["A"], d["tg"], d["train_ctl"]
    o = forward_np(P0, d["x"])[1][:, j]
    xc = x[C].mean(0)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(np.vstack([x[C], xb]), np.r_[np.ones(C.sum()), np.zeros(len(xb))])
    gx = clf.predict_proba(x)[:, 1]
    Phi = np.hstack([x - xc, np.ones((len(x), 1))]) * gx[:, None]
    base = (1 - gx) * o
    scale = o[C].std() / y[C].std()
    t = o[C].mean() + scale * (y - y[C].mean()) - base
    pen = np.r_[np.ones(x.shape[1]), 0.0] * x[C].var(0).mean()
    ci = np.where(C)[0]; nc = len(ci)
    rng = np.random.default_rng(0); fold = rng.permutation(nc) % 5

    def fit_controls(idx):
        """Ridge on the control cells idx at every weight, with its own five-fold error (for the uncertainty of the rows)."""
        f2 = np.random.default_rng(1).permutation(len(idx)) % 5
        cv = np.zeros(len(GRID))
        for k in range(5):
            M, v = G.moments(Phi[idx[f2 != k]], t[idx[f2 != k]])
            for i, lam in enumerate(GRID):
                cv[i] += ((Phi[idx[f2 == k]] @ np.linalg.solve(M + lam * np.diag(pen), v) - t[idx[f2 == k]]) ** 2).sum()
        cv /= len(idx)
        M, v = G.moments(Phi[idx], t[idx])
        Minv = [np.linalg.inv(M + lam * np.diag(pen)) for lam in GRID]
        return Minv, [Mi @ v for Mi in Minv], cv

    def rows_of(fit, idx):
        pc, bc, yc, vc = Phi[idx].mean(0), base[idx].mean(), y[idx].mean(), y[idx].var(ddof=1) / len(idx)
        ev = [(A == 1) & (tg == g) for g in fit]
        D = np.array([Phi[T].mean(0) - pc for T in ev]); dy = np.array([y[T].mean() - yc for T in ev])
        se2 = np.array([y[T].var(ddof=1) / T.sum() + vc for T in ev])
        return ev, D, scale * dy - np.array([base[T].mean() - bc for T in ev]), se2

    def upd(Minv, th0, cv, n, D, r, se2):
        return G.update(th0, Minv, D, r, n * scale ** 2 * se2 / cv)[0]

    full = fit_controls(ci)
    icv = int(full[2].argmin())
    parts = [fit_controls(ci[fold != k]) for k in range(5)]
    logs = {k: {} for k in ["controls", "rows", "rows_x", "cells"]}
    picks, Cc = [], d["audit_ctl"]
    curve = np.zeros((len(GRID), len(genes)))
    for kg, held in enumerate(genes):
        fit = [g for g in genes if g != held]
        sc = np.zeros(len(GRID))
        for k in range(5):
            tr_i, te_i = ci[fold != k], ci[fold == k]
            ev, D, r, se2 = rows_of(fit, tr_i)
            Minv, th0, cv = parts[k]
            for i in range(len(GRID)):
                for q, T in enumerate(ev):
                    keep = np.arange(len(ev)) != q
                    th = upd(Minv[i], th0[i], cv[i], len(tr_i), D[keep], r[keep], se2[keep])
                    use = np.r_[np.where(T)[0], te_i]
                    sc[i] += abs(G.u_point(y[use], (base + Phi @ th)[use], A[use])[0])
        sc /= 5 * len(fit)
        ix = int(sc.argmin()); picks.append(GRID[ix])
        ev, D, r, se2 = rows_of(fit, ci)
        Minv, th0, cv = full
        th = {"controls": th0[icv], "rows": upd(Minv[icv], th0[icv], cv[icv], nc, D, r, se2), "rows_x": upd(Minv[ix], th0[ix], cv[ix], nc, D, r, se2)}
        tr = C | ((A == 1) & np.isin(tg, fit))
        lam2, _e, _c = G.cv_weight(Phi[tr], t[tr], pen)
        M2, v2 = G.moments(Phi[tr], t[tr]); th["cells"] = np.linalg.solve(M2 + lam2 * np.diag(pen), v2)
        for k_, v_ in th.items():
            logs[k_][held] = base + Phi @ v_
        use = ((A == 1) & (tg == held)) | Cc
        for i in range(len(GRID)):
            curve[i, kg] = G.u_point(y[use], (base + Phi @ upd(Minv[i], th0[i], cv[i], nc, D, r, se2))[use], A[use])[0]
        print(" ", held, "inner score by weight", [round(float(s), 3) for s in sc], "->", GRID[ix])
    use = (A == 1) | Cc
    out = {"cv_weight": GRID[icv], "picks": picks, "arms": {}}
    for k_, lg in logs.items():
        p = np.zeros(len(A))
        for kg, g in enumerate(genes):
            own = ((A == 1) & (tg == g)) | (d["ctl_group"] == kg)
            p[own] = lg[g][own]
        u, slope = own_output_u(y[use], p[use], A[use], tg[use]); iv = interval(u[1:])
        ev = {g: G.u_point(y[((A == 1) & (tg == g)) | Cc], lg[g][((A == 1) & (tg == g)) | Cc], A[((A == 1) & (tg == g)) | Cc])[0] for g in genes}
        out["arms"][k_] = {"u": float(u[0]), "ci": iv, "reading": reading(iv), "slope": slope, "events": ev}
        print(f" {k_:9s} {u[0]:6.3f} [{iv[0]:6.3f}, {iv[1]:6.3f}] {reading(iv):10s} slope {slope:6.2f} | " + " ".join(f"{g}:{v:.2f}" for g, v in ev.items()))
    print("  held-out event by weight, mean |u|:", {f"{GRID[i]:.3g}": round(float(np.abs(curve[i]).mean()), 3) for i in range(len(GRID))}, "| cv on controls ->", GRID[icv])
    out["held_out_event_by_weight"] = {f"{GRID[i]:.3g}": float(np.abs(curve[i]).mean()) for i in range(len(GRID))}
    return out


if __name__ == "__main__":
    res = {}
    for ds in G.TARGETS:
        print(G.CELLS[ds]); res[ds] = run(ds)
    (RUN / "results" / "develop_weight.json").write_text(json.dumps(res, indent=1))
