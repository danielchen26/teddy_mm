"""After the gate was read (post hoc, development): why place B fails, and what a repair scoped to the cell type of its rows does.
Every corrected head is read on two scales: its log output (the scale of the measured protein) and its predicted mean (as the gate).

  python posthoc.py     (same environment variables as gate.py; writes results/posthoc.json)

1. The cost of place A on bone marrow: the head's PD-L1 score there if its PD-L1 output were the readout of place A everywhere.
2. Place B with the bone-marrow term at weights 0, 0.01, 0.1 and 1 (the gate used 1).
3. A correction scoped to the cell type of its rows. A gate g(x) between 0 and 1 is fitted on embeddings alone (training control
   cells of the data set against the bone-marrow training cells; no protein, no audit cell). The corrected log output is
   o(x) + g(x) * c(x), with c linear as in place B and no bone-marrow term (`gate_add`), or (1 - g(x)) * o(x) + g(x) * c(x)
   with c fitted to the target itself (`gate_replace`).
Kinds as in gate.py: controls, rows (with their uncertainty), cells. Pooled readings with intervals; point values by event.
"""
import json

import numpy as np
from sklearn.linear_model import LogisticRegression

import gate as G
from gate import BM_FIT, BM_TEST, CELLS, FEATURES, MARGIN, TARGETS, forward_np, interval, load, own_output_u, reading, unit, weights

KINDS = ["controls", "rows", "cells"]


def fit_kinds(d, guide, held, genes, Phi, t, C, extra, scale, o=None, gfun=None):
    """Ridge on `controls`, the update by rows, and `cells`, for features Phi and targets t. Returns the three parameter vectors."""
    x, y, A, tg = d["x"].astype(np.float64), d["y"], d["A"], d["tg"]
    fitg = [g for g in genes if g != held]
    pen = np.r_[np.ones(Phi.shape[1] - 1), 0.0] * x[C].var(0).mean()
    lam, cv, _ = G.cv_weight(Phi[C], t[C], pen, extra)
    M, v = G.moments(Phi[C], t[C], extra)
    Minv = np.linalg.inv(M + lam * np.diag(pen))
    th = {"controls": Minv @ v}
    ev = G.events(d, guide, fitg, False)
    D = np.array([Phi[T].mean(0) - Phi[C].mean(0) for _n, T in ev])
    dy = np.array([y[T].mean() - y[C].mean() for _n, T in ev])
    se2 = np.array([y[T].var(ddof=1) / T.sum() + y[C].var(ddof=1) / C.sum() for _n, T in ev])
    r = scale * dy - (0 if o is None else np.array([o[T].mean() - o[C].mean() for _n, T in ev]))
    th["rows"], _closed, _g = G.update(th["controls"], Minv, D, r, C.sum() * scale ** 2 * se2 / cv)
    tr = C | ((A == 1) & np.isin(tg, fitg))
    lam2, _cv2, _ = G.cv_weight(Phi[tr], t[tr], pen, extra)
    M2, v2 = G.moments(Phi[tr], t[tr], extra)
    th["cells"] = np.linalg.solve(M2 + lam2 * np.diag(pen), v2)
    return th


def pooled(d, genes, logs):
    """logs: for each held-out gene, the corrected log output on every cell. Read on the log scale (the scale of the measurement)
    and, as the gate did, as the predicted mean."""
    A, tg, y = d["A"], d["tg"], d["y"]
    use = (A == 1) | d["audit_ctl"]
    T, Cc = A == 1, d["audit_ctl"]
    out = {}
    for name, f in (("log", lambda v: v), ("mean", np.exp)):
        preds = {g: f(v) for g, v in logs.items()}
        p = np.zeros(len(A))
        for k, g in enumerate(genes):
            own = ((A == 1) & (tg == g)) | (d["ctl_group"] == k)
            p[own] = preds[g][own]
        u, slope = own_output_u(y[use], p[use], A[use], tg[use])
        iv = interval(u[1:])
        ev = {g: G.u_point(y[((A == 1) & (tg == g)) | Cc], preds[g][((A == 1) & (tg == g)) | Cc], A[((A == 1) & (tg == g)) | Cc])[0] for g in genes}
        out[name] = {"u": float(u[0]), "ci": iv, "reading": reading(iv), "slope": slope, "valid": bool(slope > 0 and np.sign(p[T].mean() - p[Cc].mean()) == np.sign(y[T].mean() - y[Cc].mean())), "events": ev}
    return out


def bm(out_all, j, yt, new_j):
    """Mean over proteins of the correlation with the measured counts on the bone-marrow test cells, and PD-L1 alone, with a new PD-L1 output (log scale)."""
    def score(col):
        pred = np.exp(out_all); pred[:, j] = np.exp(col)
        rs = [np.corrcoef(yt[:, c], pred[:, c])[0, 1] for c in range(out_all.shape[1]) if yt[:, c].std() > 1e-8]
        return float(np.mean(rs)), float(np.corrcoef(yt[:, j], pred[:, j])[0, 1])
    (m0, p0), (m1, p1) = score(out_all[:, j]), score(new_j)
    return {"mean_fall": m0 - m1, "pdl1_as_trained": p0, "pdl1_now": p1}


def main():
    P0, names, j = weights()
    b = np.load(BM_FIT); t_ = np.load(BM_TEST)
    xb, xt, yt = unit(b["z_train"]).astype(np.float64), unit(t_["z"]).astype(np.float64), t_["adt"].astype(np.float64)
    out_t = forward_np(P0, xt)[1]
    res = {"note": "post hoc, development on cells already read", "data_sets": {}}
    for ds, genes in TARGETS.items():
        d = load(ds); guide = np.load(FEATURES / ds / "dense.npz")["guide"].astype(str)
        x, y, A, tg, C = d["x"].astype(np.float64), d["y"], d["A"], d["tg"], d["train_ctl"]
        o = forward_np(P0, d["x"])[1][:, j]
        xc = x[C].mean(0)
        F, Fb, Ft = (np.hstack([z - xc, np.ones((len(z), 1))]) for z in (x, xb, xt))
        scale = o[C].std() / y[C].std()
        target = o[C].mean() + scale * (y - y[C].mean())                    # the measured value on the head's scale
        entry = {}
        # 1. place A applied to bone marrow: the PD-L1 output replaced everywhere by the readout fitted on this cell type
        thA = {g: fit_kinds(d, guide, g, genes, F, y, C, None, 1.0) for g in genes}
        entry["A_on_bone_marrow"] = {k: {"pdl1_as_trained": bm(out_t, j, yt, out_t[:, j])["pdl1_as_trained"],
                                         "pdl1_if_replaced": float(np.mean([np.corrcoef(yt[:, j], Ft @ thA[g][k])[0, 1] for g in genes]))} for k in KINDS}
        # 2. place B, the bone-marrow term at four weights
        entry["B_by_guard_weight"] = {}
        for w in (0.0, 0.01, 0.1, 1.0):
            extra = w * Fb.T @ Fb / len(Fb)
            th = {g: fit_kinds(d, guide, g, genes, F, target - o, C, extra, scale, o) for g in genes}
            row = {}
            for k in KINDS:
                r = pooled(d, genes, {g: o + F @ th[g][k] for g in genes})
                bmr = [bm(out_t, j, yt, out_t[:, j] + Ft @ th[g][k]) for g in genes]
                r["bone_marrow"] = {"mean_fall": float(np.mean([q["mean_fall"] for q in bmr])), "pdl1_as_trained": bmr[0]["pdl1_as_trained"], "pdl1_now": float(np.mean([q["pdl1_now"] for q in bmr]))}
                row[k] = r
            entry["B_by_guard_weight"][str(w)] = row
        # 3. the correction scoped to the cell type of its rows
        clf = LogisticRegression(C=1.0, max_iter=2000).fit(np.vstack([x[C], xb]), np.r_[np.ones(C.sum()), np.zeros(len(xb))])
        gx, gt = clf.predict_proba(x)[:, 1], clf.predict_proba(xt)[:, 1]
        entry["gate"] = {"training_controls": float(gx[C].mean()), "audit_controls": float(gx[d["audit_ctl"]].mean()), "knockout_cells": float(gx[A == 1].mean()),
                         "knockout_cells_lowest": float(gx[A == 1].min()), "bone_marrow_test_cells": float(gt.mean()), "bone_marrow_test_cells_highest": float(gt.max())}
        Fg, Ftg = F * gx[:, None], Ft * gt[:, None]
        thadd = {g: fit_kinds(d, guide, g, genes, Fg, target - o, C, None, scale, o) for g in genes}
        threp = {g: fit_kinds(d, guide, g, genes, Fg, target - (1 - gx) * o, C, None, scale, (1 - gx) * o) for g in genes}
        entry["scoped"] = {"gate_add": {}, "gate_replace": {}}
        for k in KINDS:
            r = pooled(d, genes, {g: o + Fg @ thadd[g][k] for g in genes})
            bmr = [bm(out_t, j, yt, out_t[:, j] + Ftg @ thadd[g][k]) for g in genes]
            r["bone_marrow"] = {"mean_fall": float(np.mean([q["mean_fall"] for q in bmr])), "pdl1_as_trained": bmr[0]["pdl1_as_trained"], "pdl1_now": float(np.mean([q["pdl1_now"] for q in bmr]))}
            entry["scoped"]["gate_add"][k] = r
            r = pooled(d, genes, {g: (1 - gx) * o + Fg @ threp[g][k] for g in genes})
            bmr = [bm(out_t, j, yt, (1 - gt) * out_t[:, j] + Ftg @ threp[g][k]) for g in genes]
            r["bone_marrow"] = {"mean_fall": float(np.mean([q["mean_fall"] for q in bmr])), "pdl1_as_trained": bmr[0]["pdl1_as_trained"], "pdl1_now": float(np.mean([q["pdl1_now"] for q in bmr]))}
            entry["scoped"]["gate_replace"][k] = r
        res["data_sets"][ds] = entry
        one = lambda q: f"{q['u']:6.3f} [{q['ci'][0]:6.3f}, {q['ci'][1]:6.3f}] {q['reading'][:5]} {'ok' if q['valid'] else 'NO'}"
        f = lambda r: f"log {one(r['log'])} | mean {one(r['mean'])}" + (f" | bm fall {r['bone_marrow']['mean_fall']:.4f}, PD-L1 {r['bone_marrow']['pdl1_as_trained']:.3f} -> {r['bone_marrow']['pdl1_now']:.3f}" if "bone_marrow" in r else "")
        print(f"\n{CELLS[ds]}")
        print(" A applied to bone marrow, PD-L1 score there:", {k: (round(v["pdl1_as_trained"], 3), round(v["pdl1_if_replaced"], 3)) for k, v in entry["A_on_bone_marrow"].items()})
        for w, row in entry["B_by_guard_weight"].items():
            for k in KINDS:
                print(f" B guard {w:5s} {k:9s} {f(row[k])}")
        print(" gate:", {k: round(v, 4) for k, v in entry["gate"].items()})
        for name in ("gate_add", "gate_replace"):
            for k in KINDS:
                print(f" {name:13s} {k:9s} {f(entry['scoped'][name][k])} | events (log) " + " ".join(f"{g}:{v:.2f}" for g, v in entry["scoped"][name][k]["log"]["events"].items()))
    (G.RES / "posthoc.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
