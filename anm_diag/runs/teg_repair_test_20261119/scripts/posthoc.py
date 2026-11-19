"""After the reading (post hoc; the co-culture cells have now been read, so nothing here is a test).

  python posthoc.py     (environment of test.py; writes results/posthoc.json)

1. Setting 1 (same events, new cells), fit half against read half: is the remainder there on the cells of the fit too?
2. Setting 1 by ridge weight (the rule chose 0.32 and 0.1 in the two halves).
3. Setting 1 with five folds in place of two halves (each fit sees four fifths of the cells).
4. Why the repair does not carry from the IFN-gamma condition: the displacement of each gene in the two conditions, the readouts
   of the two conditions, and a mixed repair: the control cells of the co-culture condition with the rows of the IFN-gamma condition.
"""
import json
import os
from pathlib import Path

import numpy as np

import test as T
from test import BM_FIT, GRID, KINDS, SEED, TARGETS, forward_np, halves, load, unit, weights

RES = T.RES


class Repair(T.Repair):
    """The repair of test.py with two options: a given ridge weight, and rows brought from other cells (xs, ys, os_, Cs, events)."""

    def __init__(self, x, y, o, C, events, xb, lam=None, rows_from=None):
        if lam is None and rows_from is None:
            super().__init__(x, y, o, C, events, xb)
            return
        from sklearn.linear_model import LogisticRegression
        self.xc = x[C].mean(0)
        self.clf = LogisticRegression(C=1.0, max_iter=2000).fit(np.vstack([x[C], xb]), np.r_[np.ones(len(C)), np.zeros(len(xb))])
        Phi, base = self.features(x, o)
        scale = o[C].std() / y[C].std()
        t = o[C].mean() + scale * (y - y[C].mean()) - base
        pen = np.r_[np.ones(x.shape[1]), 0.0] * x[C].var(0).mean()
        lam0, cv = T.cv_weight(Phi[C], t[C], pen)
        if lam is not None:                                              # the cross-validated error at the given weight
            fold = np.random.default_rng(SEED).permutation(len(C)) % T.FOLDS
            cv = 0.0
            for k in range(T.FOLDS):
                tr, te = C[fold != k], C[fold == k]
                th = np.linalg.solve(Phi[tr].T @ Phi[tr] / len(tr) + lam * np.diag(pen), Phi[tr].T @ t[tr] / len(tr))
                cv += ((Phi[te] @ th - t[te]) ** 2).sum()
            cv /= len(C)
        lam = lam0 if lam is None else lam
        Minv = np.linalg.inv(Phi[C].T @ Phi[C] / len(C) + lam * np.diag(pen))
        self.th = {"controls": Minv @ (Phi[C].T @ t[C] / len(C))}
        if rows_from is None:
            P_, b_, y_, C_, ev = Phi, base, y, C, events
        else:                                                            # the rows of other cells, seen through this repair's gate
            xs, ys, os_, C_, ev = rows_from
            P_, b_ = self.features(xs, os_)
            y_ = ys
        pc, bc, yc, vc = P_[C_].mean(0), b_[C_].mean(), y_[C_].mean(), y_[C_].var(ddof=1) / len(C_)
        D = np.array([P_[q].mean(0) - pc for _n, q in ev])
        dy = np.array([y_[q].mean() - yc for _n, q in ev])
        se2 = np.array([y_[q].var(ddof=1) / len(q) + vc for _n, q in ev])
        r = scale * dy - np.array([b_[q].mean() - bc for _n, q in ev])
        K = Minv @ D.T @ np.linalg.inv(D @ Minv @ D.T + np.diag(len(C) * scale ** 2 * se2 / cv))
        self.th["rows"] = self.th["controls"] + K @ (r - D @ self.th["controls"])
        tr = np.r_[C, np.concatenate([q for _n, q in events])]
        self.th["cells"] = np.linalg.solve(Phi[tr].T @ Phi[tr] / len(tr) + lam * np.diag(pen), Phi[tr].T @ t[tr] / len(tr))


def u_of(y, p, A):
    c, t = A == 0, A == 1
    a = np.cov(p[c], y[c])[0, 1] / p[c].var(ddof=1)
    return float(1 - a * (p[t].mean() - p[c].mean()) / (y[t].mean() - y[c].mean()))


def folds_of(A, tg, k, seed=SEED):
    rng = np.random.default_rng(seed)
    if k == 2:
        return np.where(halves(np.where(A == 1, tg, "control"), rng), 0, 1)
    f = np.zeros(len(A), int)
    for g in np.unique(np.where(A == 1, tg, "control")):
        idx = rng.permutation(np.where(np.where(A == 1, tg, "control") == g)[0])
        f[idx] = np.arange(len(idx)) % k
    return f


def crossfit(x, y, o, A, tg, genes, xb, k=2, lam=None, rows_from=None, inside=False):
    """Setting 1 with k folds: every cell read by the repair fitted on the other folds (or, with inside, on its own fold's fit)."""
    f = folds_of(A, tg, k)
    p = {kind: np.full(len(A), np.nan) for kind in KINDS}
    for q in range(k):
        fit = (f != q) if k > 2 else (f == q)                            # two halves: fit on one, read the other (as test.py)
        rd = fit if inside else ~fit
        rep = Repair(x, y, o, np.where(fit & (A == 0))[0], [(g, np.where(fit & (A == 1) & (tg == g))[0]) for g in genes], xb, lam, rows_from)
        for kind in KINDS:
            if inside and k == 2:
                p[kind][rd] = rep.output(kind, x[rd], o[rd])
            elif not inside:
                p[kind][rd] = rep.output(kind, x[rd], o[rd])
    return p


def main():
    P0, names, j = weights()
    xb = unit(np.load(BM_FIT)["z_train"]).astype(np.float64)
    head = lambda z: forward_np(P0, z)[1][:, j]
    genes = TARGETS["frangieh2021"]
    co = np.load(Path(os.environ["ANM_TEST_CELLS"]) / "coculture.npz")
    x, y, A, tg = unit(co["Z"]).astype(np.float64), co["readout__PDL1"].astype(np.float64), co["A"], co["target"].astype(str)
    o = head(x)
    read = json.loads((RES / "reading.json").read_text())
    res = {"note": "post hoc: the co-culture cells had been read when this was computed; nothing here is a test"}
    # 1. fit half against read half
    out_, in_ = crossfit(x, y, o, A, tg, genes, xb), crossfit(x, y, o, A, tg, genes, xb, inside=True)
    for kind in KINDS:
        assert abs(T.pooled(y, out_[kind], A, tg)[0]["u"] - read["coculture"]["same_events"][kind]["u"]) < 1e-9, kind
    res["fit_half_against_read_half"] = {kind: {"read_half": T.pooled(y, out_[kind], A, tg)[0]["u"], "fit_half": T.pooled(y, in_[kind], A, tg)[0]["u"]} for kind in KINDS}
    # 2. by ridge weight
    res["by_ridge_weight"] = {}
    for lam in GRID:
        p = crossfit(x, y, o, A, tg, genes, xb, lam=lam)
        res["by_ridge_weight"][f"{lam:.3g}"] = {kind: T.pooled(y, p[kind], A, tg)[0]["u"] for kind in KINDS}
    # 3. five folds
    p5 = crossfit(x, y, o, A, tg, genes, xb, k=5)
    res["five_folds"] = {}
    for kind in KINDS:
        r, _d = T.pooled(y, p5[kind], A, tg)
        res["five_folds"][kind] = {q: r[q] for q in ("u", "ci", "reading", "ci_by_cell", "reading_by_cell")} | {"events": {g: T.one_event(y[((A == 1) & (tg == g)) | (A == 0)], p5[kind][((A == 1) & (tg == g)) | (A == 0)], A[((A == 1) & (tg == g)) | (A == 0)]) for g in genes}}
    # 4. the two conditions
    d = load("frangieh2021")
    xi, yi, Ai, ti = d["x"].astype(np.float64), d["y"], d["A"], d["tg"]
    oi, Ci = head(xi), np.where(d["train_ctl"])[0]
    cos = lambda a, b: float(a @ b / np.linalg.norm(a) / np.linalg.norm(b))
    geo = {}
    for g in genes:
        di = xi[(Ai == 1) & (ti == g)].mean(0) - xi[Ai == 0].mean(0)
        dc = x[(A == 1) & (tg == g)].mean(0) - x[A == 0].mean(0)
        geo[g] = {"cosine_of_the_displacements": cos(di, dc), "length_ifn": float(np.linalg.norm(di)), "length_coculture": float(np.linalg.norm(dc)),
                  "measured_change_ifn": float(yi[(Ai == 1) & (ti == g)].mean() - yi[Ai == 0].mean()), "measured_change_coculture": float(y[(A == 1) & (tg == g)].mean() - y[A == 0].mean())}
    shift = x[A == 0].mean(0) - xi[Ai == 0].mean(0)
    rep_i = T.Repair(xi, yi, oi, Ci, [(g, np.where((Ai == 1) & (ti == g))[0]) for g in genes], xb)
    rep_c = T.Repair(x, y, o, np.where(A == 0)[0], [(g, np.where((A == 1) & (tg == g))[0]) for g in genes], xb)
    res["two_conditions"] = {"by_gene": geo, "distance_between_the_control_means": float(np.linalg.norm(shift)),
                             "mean_length_of_a_displacement": float(np.mean([v["length_coculture"] for v in geo.values()])),
                             "control_sd_of_pdl1": {"ifn": float(yi[Ai == 0].std(ddof=1)), "coculture": float(y[A == 0].std(ddof=1))},
                             "cosine_of_the_readouts": {kind: cos(rep_i.th[kind][:-1], rep_c.th[kind][:-1]) for kind in KINDS}}
    # the mixed repair: control cells of the co-culture condition (two halves), rows of the IFN-gamma condition
    pm = crossfit(x, y, o, A, tg, genes, xb, rows_from=(xi, yi, oi, Ci, [(g, np.where((Ai == 1) & (ti == g))[0]) for g in genes]))
    r, _d = T.pooled(y, pm["rows"], A, tg)
    res["mixed_repair"] = {"what": "control cells of the co-culture condition, rows of the IFN-gamma condition", **{q: r[q] for q in ("u", "ci", "reading", "ci_by_cell", "reading_by_cell", "valid")},
                           "events": {g: u_of(y[((A == 1) & (tg == g)) | (A == 0)], pm["rows"][((A == 1) & (tg == g)) | (A == 0)], A[((A == 1) & (tg == g)) | (A == 0)]) for g in genes}}
    (RES / "posthoc.json").write_text(json.dumps(res, indent=1))
    print("1. fit half / read half:", {k: (round(v["fit_half"], 3), round(v["read_half"], 3)) for k, v in res["fit_half_against_read_half"].items()})
    print("2. by ridge weight:")
    for lam, v in res["by_ridge_weight"].items():
        print(f"   {lam:>8s} " + " ".join(f"{k} {v[k]:6.3f}" for k in KINDS))
    print("3. five folds:", {k: (round(v["u"], 3), [round(c, 3) for c in v["ci"]], v["reading"], [round(c, 3) for c in v["ci_by_cell"]], v["reading_by_cell"]) for k, v in res["five_folds"].items()})
    print("   events rows:", {g: (round(e["u"], 2), e["reading"]) for g, e in res["five_folds"]["rows"]["events"].items()})
    print("4. two conditions:", json.dumps({k: v for k, v in res["two_conditions"].items() if k != "by_gene"}))
    for g, v in geo.items():
        print("   ", g, {k: round(q, 3) for k, q in v.items()})
    print("   mixed repair:", round(res["mixed_repair"]["u"], 3), [round(c, 3) for c in res["mixed_repair"]["ci"]], res["mixed_repair"]["reading"], {g: round(v, 2) for g, v in res["mixed_repair"]["events"].items()})


if __name__ == "__main__":
    main()
