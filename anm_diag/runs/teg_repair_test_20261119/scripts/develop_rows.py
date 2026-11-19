"""Development, before DESIGN.md was written. Nothing here is a test.

Question: does the repair by rows get better when the table also holds the rows of the other genes of the library (200 genes
of the IFN-gamma condition outside the pathway, most of them with no change of PD-L1), and with a finer grid of ridge weights?
Read on cells already read (melanoma, IFN-gamma, the five genes, one gene left out at a time; the split of the gate). The cells
of the other 200 genes enter only as rows (their mean displacement and their measured change) or, in `cells_wide`, as cells of a
fit; no output is read on them. CD274 is kept out: it is the probe of the test.

  python develop_rows.py     (environment of test.py; writes results/develop_rows.json)

Kinds: controls; rows (the four pathway rows, weight by cross-validation on the controls); rows_loro (weight by leaving one row
out, read against the training controls); wide and wide_loro (the same with the 200 other rows added); cells; cells_wide.
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(RUN.parent / "teg_repair_gate_20261119" / "scripts"))
import gate as G  # noqa: E402
from gate import BM_FIT, BM_TEST, forward_np, interval, load, own_output_u, reading, unit, weights  # noqa: E402

GRID = [10.0 ** (k / 2) for k in range(-6, 5)]
PROBE = ["CD274"]


def main():
    P0, names, j = weights()
    b = np.load(BM_FIT); t_ = np.load(BM_TEST)
    xb, xt, yt = unit(b["z_train"]).astype(np.float64), unit(t_["z"]).astype(np.float64), t_["adt"].astype(np.float64)
    d = load("frangieh2021"); genes = G.TARGETS["frangieh2021"]
    x, y, A, tg, C = d["x"].astype(np.float64), d["y"], d["A"], d["tg"], d["train_ctl"]
    lib = np.load(Path(os.environ["ANM_TEST_CELLS"]) / "library.npz")
    new = ~np.isin(lib["ids"], d["ids"])
    xl, yl, tl = unit(lib["Z"][new]).astype(np.float64), lib["readout__PDL1"][new], lib["target"][new].astype(str)
    others = sorted(g for g in np.unique(tl) if g not in PROBE)
    o, ol = forward_np(P0, d["x"])[1][:, j], forward_np(P0, xl)[1][:, j]
    xc = x[C].mean(0)
    Fm = lambda z: np.hstack([z - xc, np.ones((len(z), 1))])
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(np.vstack([x[C], xb]), np.r_[np.ones(C.sum()), np.zeros(len(xb))])
    gx, gl = clf.predict_proba(x)[:, 1], clf.predict_proba(xl)[:, 1]
    Phi, Phil = Fm(x) * gx[:, None], Fm(xl) * gl[:, None]
    base, basel = (1 - gx) * o, (1 - gl) * ol
    scale = o[C].std() / y[C].std()
    t = o[C].mean() + scale * (y - y[C].mean()) - base
    tlib = o[C].mean() + scale * (yl - y[C].mean()) - basel
    pen = np.r_[np.ones(x.shape[1]), 0.0] * x[C].var(0).mean()
    pc, bc, yc, nc = Phi[C].mean(0), base[C].mean(), y[C].mean(), C.sum()
    vc = y[C].var(ddof=1) / nc
    # the other rows of the table: one per gene of the library outside the pathway (CD274 kept out, as a probe)
    Do = np.array([Phil[tl == g].mean(0) - pc for g in others])
    dyo = np.array([yl[tl == g].mean() - yc for g in others])
    se2o = np.array([yl[tl == g].var(ddof=1) / (tl == g).sum() + vc for g in others])
    ro = scale * dyo - np.array([basel[tl == g].mean() - bc for g in others])
    zs = dyo / np.sqrt(se2o)
    print("other rows:", len(others), "| measured change: median |.|", round(float(np.median(np.abs(dyo))), 3), "| largest", [(others[i], round(float(dyo[i]), 2), round(float(zs[i]), 1)) for i in np.argsort(-np.abs(dyo))[:8]])
    # ridge on the training controls at every weight of the grid
    rng = np.random.default_rng(0); fold = rng.permutation(nc) % 5
    PC, tC = Phi[C], t[C]
    cv = np.zeros(len(GRID))
    for k in range(5):
        M, v = G.moments(PC[fold != k], tC[fold != k])
        for i, lam in enumerate(GRID):
            cv[i] += ((PC[fold == k] @ np.linalg.solve(M + lam * np.diag(pen), v) - tC[fold == k]) ** 2).sum()
    cv /= nc
    M, v = G.moments(PC, tC)
    Minv = [np.linalg.inv(M + lam * np.diag(pen)) for lam in GRID]
    th0 = [Mi @ v for Mi in Minv]
    icv = int(cv.argmin())
    print("ridge weight by cross-validation on the training controls:", GRID[icv], [round(float(c), 4) for c in cv])

    def rows_of(fit):
        ev = [(g, (A == 1) & (tg == g)) for g in fit]
        D = np.array([Phi[T].mean(0) - pc for _g, T in ev]); dy = np.array([y[T].mean() - yc for _g, T in ev])
        se2 = np.array([y[T].var(ddof=1) / T.sum() + vc for _g, T in ev])
        r = scale * dy - np.array([base[T].mean() - bc for _g, T in ev])
        return ev, D, dy, se2, r

    def upd(i, D, r, se2):
        th, _c, _g = G.update(th0[i], Minv[i], D, r, nc * scale ** 2 * se2 / cv[i])
        return th

    def loro(fit, wide):
        """For every weight: leave each pathway row of the fit out in turn, update with the rest, read the left-out event (its cells and the training controls)."""
        ev, D, dy, se2, r = rows_of(fit)
        sc = np.zeros(len(GRID))
        for i in range(len(GRID)):
            us = []
            for q, (_g, T) in enumerate(ev):
                keep = np.arange(len(ev)) != q
                Dq, rq, sq = D[keep], r[keep], se2[keep]
                if wide:
                    Dq, rq, sq = np.vstack([Dq, Do]), np.r_[rq, ro], np.r_[sq, se2o]
                th = upd(i, Dq, rq, sq)
                use = T | C
                us.append(G.u_point(y[use], (base + Phi @ th)[use], A[use])[0])
            sc[i] = np.mean(np.abs(us))
        return sc

    logs = {k: {} for k in ["controls", "rows", "rows_loro", "wide", "wide_loro", "cells", "cells_wide"]}
    picks = {"rows_loro": [], "wide_loro": []}
    curves = {"rows": np.zeros((len(GRID), 5)), "wide": np.zeros((len(GRID), 5))}
    Cc = d["audit_ctl"]
    for kg, held in enumerate(genes):
        fit = [g for g in genes if g != held]
        ev, D, dy, se2, r = rows_of(fit)
        Dw, rw, sw = np.vstack([D, Do]), np.r_[r, ro], np.r_[se2, se2o]
        s_rows, s_wide = loro(fit, False), loro(fit, True)
        i1, i2 = int(s_rows.argmin()), int(s_wide.argmin())
        picks["rows_loro"].append(GRID[i1]); picks["wide_loro"].append(GRID[i2])
        th = {"controls": th0[icv], "rows": upd(icv, D, r, se2), "rows_loro": upd(i1, D, r, se2), "wide": upd(icv, Dw, rw, sw), "wide_loro": upd(i2, Dw, rw, sw)}
        tr = C | ((A == 1) & np.isin(tg, fit))
        lam2, _e, _c = G.cv_weight(Phi[tr], t[tr], pen)
        M2, v2 = G.moments(Phi[tr], t[tr]); th["cells"] = np.linalg.solve(M2 + lam2 * np.diag(pen), v2)
        Pw, tw = np.vstack([Phi[tr], Phil[np.isin(tl, others)]]), np.r_[t[tr], tlib[np.isin(tl, others)]]
        lam3, _e, _c = G.cv_weight(Pw, tw, pen)
        M3, v3 = G.moments(Pw, tw); th["cells_wide"] = np.linalg.solve(M3 + lam3 * np.diag(pen), v3)
        for k, v_ in th.items():
            logs[k][held] = base + Phi @ v_
        use = ((A == 1) & (tg == held)) | Cc
        for i in range(len(GRID)):                                        # for the record: the held-out event at every weight
            curves["rows"][i, kg] = G.u_point(y[use], (base + Phi @ upd(i, D, r, se2))[use], A[use])[0]
            curves["wide"][i, kg] = G.u_point(y[use], (base + Phi @ upd(i, Dw, rw, sw))[use], A[use])[0]
        print(held, "| inner scores rows", [round(float(s), 2) for s in s_rows], "-> ", GRID[i1], "| wide", [round(float(s), 2) for s in s_wide], "->", GRID[i2], "| cells ridge", lam2, lam3)
    use = (A == 1) | Cc
    res = {"grid": GRID, "cv_weight": GRID[icv], "picks": picks, "arms": {}}
    for k, lg in logs.items():
        p = np.zeros(len(A))
        for kg, g in enumerate(genes):
            own = ((A == 1) & (tg == g)) | (d["ctl_group"] == kg)
            p[own] = lg[g][own]
        u, slope = own_output_u(y[use], p[use], A[use], tg[use]); iv = interval(u[1:])
        ev = {g: G.u_point(y[((A == 1) & (tg == g)) | Cc], lg[g][((A == 1) & (tg == g)) | Cc], A[((A == 1) & (tg == g)) | Cc])[0] for g in genes}
        res["arms"][k] = {"u": float(u[0]), "ci": iv, "reading": reading(iv), "slope": slope, "events": ev}
        print(f"{k:11s} {u[0]:6.3f} [{iv[0]:6.3f}, {iv[1]:6.3f}] {reading(iv):10s} slope {slope:6.2f} | " + " ".join(f"{g}:{v:.2f}" for g, v in ev.items()))
    print("held-out event by weight, mean |u| over the five genes:")
    for nm in ("rows", "wide"):
        print(" ", nm, {f"{GRID[i]:.3g}": round(float(np.abs(curves[nm][i]).mean()), 3) for i in range(len(GRID))})
    res["curves"] = {k: v.tolist() for k, v in curves.items()}
    res["other_rows"] = {"n": len(others), "median_abs_measured_change": float(np.median(np.abs(dyo))), "largest": [[others[i], float(dyo[i])] for i in np.argsort(-np.abs(dyo))[:8]]}
    (RUN / "results" / "develop_rows.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
