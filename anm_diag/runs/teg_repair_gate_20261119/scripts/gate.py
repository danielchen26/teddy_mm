"""The gate of DESIGN.md: a linear repair computed from the rows of the events of the fit, against one that sees every
knockout cell. Development on cells already read; nothing here is a test.

  ANM_FEATURES=... ANM_BM_FIT=... ANM_BM_TEST=... python gate.py

Writes results/gate.json and prints the tables. Everything is closed form (ridge regressions and one update); it runs
in about a minute on a laptop.

Places:  A  the PD-L1 output is a linear readout of TEDDY's embedding (unit length)
         B  the head as trained, plus a linear correction of its PD-L1 output (log scale), with the bone-marrow guard
Kinds:   controls  the measured PD-L1 of the training control cells
         rows      controls, updated by one row per event of the fit: w = w0 + H^-1 D'(D H^-1 D' + S)^-1 (r - D w0)
         rows_hard the same with S = 0 (every row met exactly); a variant
         rows_guides  rows by guide instead of by gene (guides with at least 20 cells); a variant
         cells     controls and the measured PD-L1 of every knockout cell of the fit (the reference)
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(RUN.parent / "teg_training_20261119" / "scripts"))
from _common import BM_FIT, BM_TEST, CELLS, FEATURES, MARGIN, TARGETS, forward_np, interval, load, own_output_u, reading, unit, weights  # noqa: E402

RES = RUN / "results"
GRID = [10.0 ** k for k in range(-3, 3)]          # ridge weight, in units of the mean variance of a coordinate
GUARD_W, FOLDS, MIN_GUIDE = 1.0, 5, 20
KINDS = ["controls", "rows", "rows_hard", "rows_guides", "cells"]


# --- ridge with a chosen weight ------------------------------------------------------------------------------------------
def moments(F, t, extra=None):
    """Second moment and cross moment of features F with targets t (means over rows); `extra` is added to the second moment."""
    M = F.T @ F / len(F)
    return (M if extra is None else M + extra), F.T @ t / len(F)


def cv_weight(F, t, pen, extra=None, seed=0):
    """Five-fold cross-validation of a ridge fit over GRID. pen: the diagonal that the weight multiplies (0 for an intercept)."""
    rng = np.random.default_rng(seed)
    fold = rng.permutation(len(F)) % FOLDS
    err = np.zeros(len(GRID))
    for k in range(FOLDS):
        tr, te = fold != k, fold == k
        M, v = moments(F[tr], t[tr], extra)
        for i, lam in enumerate(GRID):
            th = np.linalg.solve(M + lam * np.diag(pen), v)
            err[i] += ((F[te] @ th - t[te]) ** 2).sum()
    err /= len(F)
    return GRID[int(err.argmin())], float(err.min()), [float(e) for e in err]


def update(th0, Minv, D, r, S):
    """The rows as observations of linear functionals of the fit: the smallest change, in the metric of M, that brings each row
    inside its own uncertainty. Returns the new fit and, per row, the share of its residual that was closed."""
    G = D @ Minv @ D.T
    K = Minv @ D.T @ np.linalg.inv(G + np.diag(S))
    th = th0 + K @ (r - D @ th0)
    before, after = r - D @ th0, r - D @ th
    return th, [float(1 - abs(a) / abs(b)) if abs(b) > 1e-12 else 1.0 for a, b in zip(after, before)], [float(g) for g in np.diag(G)]


def u_point(y, pred, A):
    """The point value of the statistic: 1 - a * shift(pred) / shift(y), a the slope of y on pred among the control cells."""
    c, t = A == 0, A == 1
    a = np.cov(pred[c], y[c])[0, 1] / pred[c].var(ddof=1)
    return float(1 - a * (pred[t].mean() - pred[c].mean()) / (y[t].mean() - y[c].mean())), float(a)


# --- one data set ----------------------------------------------------------------------------------------------------
def events(d, guide, genes, by_guide):
    """The rows of the fit for a set of genes: name, cells."""
    out = []
    for g in genes:
        T = (d["A"] == 1) & (d["tg"] == g)
        if not by_guide:
            out.append((g, T))
        else:
            for q in np.unique(guide[T]):
                Tq = T & (guide == q)
                if Tq.sum() >= MIN_GUIDE:
                    out.append((f"{g}|{q}", Tq))
    return out


def fit_place(d, guide, held, genes, place, base, bm):
    """All kinds for one held-out gene, in one place. Returns predictions for every cell (the own PD-L1 output, on the scale the
    reading uses) and a log. base: for place B the head's PD-L1 output as trained (log scale) on these cells and on bone marrow."""
    x, y, A, tg = d["x"].astype(np.float64), d["y"], d["A"], d["tg"]
    C = d["train_ctl"]
    fitg = [g for g in genes if g != held]
    xc = x[C].mean(0)
    F = np.hstack([x - xc, np.ones((len(x), 1))])                       # features: the embedding centred on the training controls, and 1
    pen = np.r_[np.ones(x.shape[1]), 0.0] * x[C].var(0).mean()             # the ridge weight is in units of the mean variance of a coordinate
    if place == "A":
        t, extra, scale = y.copy(), None, 1.0
    else:                                                                  # B: targets on the head's scale, as residuals from the head as trained
        o, ob, Fb = base["o"], base["ob"], np.hstack([base["xb"] - xc, np.ones((len(base["xb"]), 1))])
        scale = o[C].std() / y[C].std()
        t = (o[C].mean() + scale * (y - y[C].mean())) - o                  # what has to be added to the head's output
        extra = GUARD_W * Fb.T @ Fb / len(Fb)                              # the correction stays near zero on the bone-marrow training cells
    log, th = {}, {}
    # controls
    lam, cv, curve = cv_weight(F[C], t[C], pen, extra)
    M, v = moments(F[C], t[C], extra)
    Minv = np.linalg.inv(M + lam * np.diag(pen))
    th["controls"] = Minv @ v
    log["controls"] = {"ridge_weight": lam, "cv_error": cv, "cv_curve": curve, "n_controls": int(C.sum())}
    # rows: by gene (with and without their uncertainty), and by guide
    for kind, by_guide, hard in (("rows", False, False), ("rows_hard", False, True), ("rows_guides", True, False)):
        ev = events(d, guide, fitg, by_guide)
        D = np.array([np.r_[x[T].mean(0) - xc, 0.0] for _n, T in ev])
        dy = np.array([y[T].mean() - y[C].mean() for _n, T in ev])
        se2 = np.array([y[T].var(ddof=1) / T.sum() + y[C].var(ddof=1) / C.sum() for _n, T in ev])
        if place == "A":
            r = dy
        else:                                                              # the change the head still owes for this event, on its own scale
            r = scale * dy - np.array([base["o"][T].mean() - base["o"][C].mean() for _n, T in ev])
        S = np.zeros(len(ev)) if hard else C.sum() * scale ** 2 * se2 / cv
        th[kind], closed, G = update(th["controls"], Minv, D, r, S)
        log[kind] = {"events": [n for n, _T in ev], "n_cells": [int(T.sum()) for _n, T in ev], "measured_change": [float(v_) for v_ in dy], "standard_error": [float(np.sqrt(s)) for s in se2],
                     "share_of_the_residual_closed": closed, "distance_of_the_event_from_the_controls": G, "row_uncertainty_on_the_same_scale": [float(s) for s in S]}
    # cells
    tr = C | ((A == 1) & np.isin(tg, fitg))
    lam2, cv2, curve2 = cv_weight(F[tr], t[tr], pen, extra)
    M2, v2 = moments(F[tr], t[tr], extra)
    th["cells"] = np.linalg.solve(M2 + lam2 * np.diag(pen), v2)
    log["cells"] = {"ridge_weight": lam2, "cv_error": cv2, "cv_curve": curve2, "n_cells": int(tr.sum())}
    if place == "A":
        pred = {k: F @ v_ for k, v_ in th.items()}
    else:                                                                  # the own output is the predicted mean: exp of the corrected log output
        pred = {k: np.exp(base["o"] + F @ v_) for k, v_ in th.items()}
        log["bone_marrow_correction_rms"] = {k: float(np.sqrt(((Fb @ v_) ** 2).mean())) for k, v_ in th.items()}
    return pred, th, log, (xc, Minv, pen, lam)


def bm_fall(P0, names, j, xt, yt, xc, th):
    """The head's own score on the bone-marrow test cells with the PD-L1 output corrected: the mean over proteins and PD-L1 alone."""
    out = forward_np(P0, xt)[1]
    out2 = out.copy(); out2[:, j] = out[:, j] + np.hstack([xt.astype(np.float64) - xc, np.ones((len(xt), 1))]) @ th
    def score(o):
        pred = np.exp(o)
        r = [np.corrcoef(yt[:, c], pred[:, c])[0, 1] for c in range(len(names)) if yt[:, c].std() > 1e-8]
        return float(np.mean(r)), float(np.corrcoef(yt[:, j], pred[:, j])[0, 1])
    (m0, p0), (m1, p1) = score(out), score(out2)
    return {"mean_fall": m0 - m1, "pdl1_as_trained": p0, "pdl1_corrected": p1}


def main():
    P0, names, j = weights()
    b = np.load(BM_FIT); t = np.load(BM_TEST)
    xb, xt, yt = unit(b["z_train"]).astype(np.float64), unit(t["z"]), t["adt"].astype(np.float64)
    ob = forward_np(P0, xb)[1][:, j]
    res = {"margin": MARGIN, "note": "development on cells already read; nothing here is a test", "grid": GRID, "data_sets": {}}
    for ds, genes in TARGETS.items():
        d = load(ds)
        guide = np.load(FEATURES / ds / "dense.npz")["guide"].astype(str)
        A, tg, y = d["A"], d["tg"], d["y"]
        base = {"o": forward_np(P0, d["x"])[1][:, j], "ob": ob, "xb": xb}
        use = (A == 1) | d["audit_ctl"]
        entry = {"cells": CELLS[ds], "n_knockout": int((A == 1).sum()), "n_audit_controls": int(d["audit_ctl"].sum()), "n_training_controls": int(d["train_ctl"].sum()), "places": {}, "fits": {}}
        fits = {(pl, g): fit_place(d, guide, g, genes, pl, base, None) for pl in "AB" for g in genes}
        for pl in "AB":
            arms, draws = {}, {}
            for kind in KINDS:
                pooled = np.zeros(len(A))
                for k, g in enumerate(genes):
                    own = ((A == 1) & (tg == g)) | (d["ctl_group"] == k)
                    pooled[own] = fits[pl, g][0][kind][own]
                u, slope = own_output_u(y[use], pooled[use], A[use], tg[use])
                draws[kind] = u
                T, Cc = A == 1, d["audit_ctl"]
                valid = bool(slope > 0 and np.sign(pooled[T].mean() - pooled[Cc].mean()) == np.sign(y[T].mean() - y[Cc].mean()))
                iv = interval(u[1:])
                arms[kind] = {"u": float(u[0]), "ci": iv, "reading": reading(iv), "slope_among_audit_controls": slope, "valid": valid, "events": {}}
                for g in genes:
                    pred = fits[pl, g][0][kind]
                    ug_use = ((A == 1) & (tg == g)) | d["audit_ctl"]
                    ug, sg = own_output_u(y[ug_use], pred[ug_use], A[ug_use], None)
                    up, _a = u_point(y[ug_use], pred[ug_use], A[ug_use]); assert abs(up - ug[0]) < 1e-9
                    Tg = (A == 1) & (tg == g)
                    arms[kind]["events"][g] = {"u": float(ug[0]), "ci": interval(ug[1:]), "reading": reading(interval(ug[1:])), "slope": sg,
                                               "valid": bool(sg > 0 and np.sign(pred[Tg].mean() - pred[Cc].mean()) == np.sign(y[Tg].mean() - y[Cc].mean())),
                                               "correlation_with_measured": float(np.corrcoef(pred[ug_use], y[ug_use])[0, 1]),
                                               "shift_over_measured_shift_in_control_sd": float(((pred[Tg].mean() - pred[Cc].mean()) / pred[Cc].std()) / ((y[Tg].mean() - y[Cc].mean()) / y[Cc].std()))}
                if pl == "B":
                    falls = [bm_fall(P0, names, j, xt, yt, fits[pl, g][3][0], fits[pl, g][1][kind]) for g in genes]
                    arms[kind]["bone_marrow"] = {"mean_fall": float(np.mean([f["mean_fall"] for f in falls])), "pdl1_as_trained": falls[0]["pdl1_as_trained"],
                                                 "pdl1_corrected": float(np.mean([f["pdl1_corrected"] for f in falls])),
                                                 "correction_rms_on_bone_marrow_training_cells": float(np.mean([fits[pl, g][2]["bone_marrow_correction_rms"][kind] for g in genes]))}
            for kind in KINDS:                                            # every kind against `cells`, on the same resamples
                diff = np.abs(draws[kind]) - np.abs(draws["cells"])
                arms[kind]["abs_u_minus_cells"] = {"difference": float(diff[0]), "ci": interval(diff[1:])}
            entry["places"][pl] = arms
            entry["fits"][pl] = {g: fits[pl, g][2] for g in genes}
        # the reading of `rows` (place A) with 1 to 4 events of the fit: every subset, point values on the held-out event
        curve = {k: [] for k in range(0, 5)}
        x = d["x"].astype(np.float64)
        for g in genes:
            xc, Minv, pen, lam = fits["A", g][3]
            F = np.hstack([x - xc, np.ones((len(x), 1))]); C = d["train_ctl"]
            th0 = fits["A", g][1]["controls"]; cv = fits["A", g][2]["controls"]["cv_error"]
            ug_use = ((A == 1) & (tg == g)) | d["audit_ctl"]
            curve[0].append(u_point(y[ug_use], (F @ th0)[ug_use], A[ug_use])[0])
            others = [e for e in genes if e != g]
            for k in range(1, 5):
                for sub in itertools.combinations(others, k):
                    ev = events(d, guide, list(sub), False)
                    D = np.array([np.r_[x[T].mean(0) - xc, 0.0] for _n, T in ev]); r = np.array([y[T].mean() - y[C].mean() for _n, T in ev])
                    S = np.array([C.sum() * (y[T].var(ddof=1) / T.sum() + y[C].var(ddof=1) / C.sum()) / cv for _n, T in ev])
                    th, _c, _g = update(th0, Minv, D, r, S)
                    curve[k].append(u_point(y[ug_use], (F @ th)[ug_use], A[ug_use])[0])
        entry["rows_by_number_of_events"] = {str(k): {"mean_abs_u": float(np.mean(np.abs(v))), "mean_u": float(np.mean(v)), "min": float(np.min(v)), "max": float(np.max(v)), "n": len(v)} for k, v in curve.items()}
        # the displacement of the held-out event against the span of the displacements of the events of the fit
        geo = {}
        for g in genes:
            xc = fits["A", g][3][0]
            dg = x[(A == 1) & (tg == g)].mean(0) - xc
            Dm = np.array([x[(A == 1) & (tg == e)].mean(0) - xc for e in genes if e != g])
            coef, *_ = np.linalg.lstsq(Dm.T, dg, rcond=None)
            proj = Dm.T @ coef
            geo[g] = {"length": float(np.linalg.norm(dg)), "cosine_with_the_span": float(np.linalg.norm(proj) / np.linalg.norm(dg)),
                      "cosine_with_each": {e: float(dg @ v_ / np.linalg.norm(dg) / np.linalg.norm(v_)) for e, v_ in zip([e for e in genes if e != g], Dm)}}
        entry["geometry"] = geo
        res["data_sets"][ds] = entry
    # the decision of DESIGN.md: rows against cells, place A, pooled, point values
    dec = {}
    for ds in TARGETS:
        a = res["data_sets"][ds]["places"]["A"]
        gap = abs(a["rows"]["u"]) - abs(a["cells"]["u"])
        dec[ds] = {"rows_valid": a["rows"]["valid"], "abs_u_rows_minus_cells": gap, "go": bool(a["rows"]["valid"] and gap <= 0.05), "no_go": bool((not a["rows"]["valid"]) or gap >= 0.10),
                   "cells_inside_the_margin": bool(a["cells"]["reading"] == "supported")}
    res["decision"] = {"by_data_set": dec, "outcome": "go" if all(v["go"] for v in dec.values()) else "no-go" if any(v["no_go"] for v in dec.values()) else "undecided"}
    RES.mkdir(exist_ok=True)
    (RES / "gate.json").write_text(json.dumps(res, indent=1))
    f = lambda v: f"{v:6.3f}"
    for ds in TARGETS:
        e = res["data_sets"][ds]
        print(f"\n{CELLS[ds]}: {e['n_knockout']} knockout cells, {e['n_audit_controls']} audit controls, {e['n_training_controls']} training controls")
        for pl in "AB":
            print(f" place {pl}   {'kind':12s} {'u, own PD-L1 output':26s} {'reading':11s} {'slope':>7s} valid {'|u| - cells':>24s}  events" + ("   bone marrow fall, PD-L1 r" if pl == "B" else ""))
            for kind in KINDS:
                r = e["places"][pl][kind]; dd = r["abs_u_minus_cells"]
                print(f"           {kind:12s} {f(r['u'])} [{f(r['ci'][0])}, {f(r['ci'][1])}]   {r['reading']:11s} {r['slope_among_audit_controls']:7.2f} {'ok ' if r['valid'] else 'NO '}  {f(dd['difference'])} [{f(dd['ci'][0])}, {f(dd['ci'][1])}]  "
                      + " ".join(f"{g}:{v['u']:.2f}{'' if v['valid'] else '!'}" for g, v in r["events"].items())
                      + (f"   {r['bone_marrow']['mean_fall']:.4f}, {r['bone_marrow']['pdl1_as_trained']:.3f} -> {r['bone_marrow']['pdl1_corrected']:.3f}" if pl == "B" else ""))
        print("  rows (A) by number of events of the fit, mean |u| on the held-out event:", {k: round(v["mean_abs_u"], 3) for k, v in e["rows_by_number_of_events"].items()})
        print("  held-out displacement, cosine with the span of the events of the fit:", {g: round(v["cosine_with_the_span"], 3) for g, v in e["geometry"].items()})
        g0 = list(TARGETS[ds])[0]
        print("  ridge weight (controls / cells), place A:", e["fits"]["A"][g0]["controls"]["ridge_weight"], e["fits"]["A"][g0]["cells"]["ridge_weight"],
              "| share of each row's residual closed (rows, A, first held-out gene):", [round(v, 2) for v in e["fits"]["A"][g0]["rows"]["share_of_the_residual_closed"]])
    print("\ndecision:", json.dumps(res["decision"]))


if __name__ == "__main__":
    main()
