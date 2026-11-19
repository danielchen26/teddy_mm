"""The test of DESIGN.md: a repair of the PD-L1 output computed from the rows of the table, scoped by a gate, read once on
cells that no run of this repository had read (melanoma, co-culture condition).

  python test.py --development     the same code on cells already read (IFN-gamma melanoma, THP-1); writes results/development.json
  python test.py                   the one reading; writes results/reading.json and refuses to run a second time

Environment: ANM_FEATURES, ANM_BM_FIT, ANM_BM_TEST as for runs/teg_training_20261119, and ANM_TEST_CELLS, the folder that
cells.py wrote (coculture.npz, library.npz).

The repair (frozen in DESIGN.md). For a set of control cells C and a set of events (a gene, its knockout cells):
  gate      g(x), a logistic regression on embeddings: the control cells C against the bone-marrow training cells.
  output    (1 - g(x)) * head(x) + g(x) * c(x), on the log scale; head(x) is the PD-L1 output as trained, c(x) a linear readout
            of the embedding.
  controls  c fitted by ridge on the measured PD-L1 of the control cells C (weight by five-fold cross-validation).
  rows      `controls`, then one update by the rows: w = w0 + H^-1 D'(D H^-1 D' + S)^-1 (r - D w0). A row is one event: the mean
            displacement of its knockout cells from the controls (D), the measured change of PD-L1 (r) and its standard error (S).
            No knockout cell enters one by one.
  cells     ridge on the control cells and on every knockout cell of the events, cell by cell (the reference).
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
from sklearn.linear_model import LogisticRegression

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(RUN.parent / "teg_training_20261119" / "scripts"))
from _common import BM_FIT, BM_TEST, CELLS, FEATURES, MARGIN, TARGETS, forward_np, halves, interval, load, own_output_u, reading, unit, weights  # noqa: E402

RES = RUN / "results"
GRID = [10.0 ** (k / 2) for k in range(-6, 5)]     # ridge weight, in units of the mean variance of a coordinate
FOLDS, SEED, GUARD = 5, 0, 0.02
KINDS = ["controls", "rows", "cells"]
PROBE = "CD274"


# --- the repair -------------------------------------------------------------------------------------------------------
def cv_weight(Phi, t, pen):
    """Five-fold cross-validation of a ridge fit over GRID; pen is the diagonal that the weight multiplies (0 for the intercept)."""
    fold = np.random.default_rng(SEED).permutation(len(Phi)) % FOLDS
    err = np.zeros(len(GRID))
    for k in range(FOLDS):
        tr, te = fold != k, fold == k
        M, v = Phi[tr].T @ Phi[tr] / tr.sum(), Phi[tr].T @ t[tr] / tr.sum()
        for i, lam in enumerate(GRID):
            err[i] += ((Phi[te] @ np.linalg.solve(M + lam * np.diag(pen), v) - t[te]) ** 2).sum()
    err /= len(Phi)
    return GRID[int(err.argmin())], float(err.min())


class Repair:
    def __init__(self, x, y, o, C, events, xb):
        """x: embeddings (unit length) of the cells of the fit, y: their measured PD-L1, o: the head's PD-L1 output as trained (log).
        C: indices of the control cells of the fit. events: [(name, indices of its knockout cells)]. xb: bone-marrow training cells."""
        self.xc = x[C].mean(0)
        self.clf = LogisticRegression(C=1.0, max_iter=2000).fit(np.vstack([x[C], xb]), np.r_[np.ones(len(C)), np.zeros(len(xb))])
        Phi, base = self.features(x, o)
        scale = o[C].std() / y[C].std()
        t = o[C].mean() + scale * (y - y[C].mean()) - base              # the measured value on the head's scale, less what the head keeps
        pen = np.r_[np.ones(x.shape[1]), 0.0] * x[C].var(0).mean()
        lam, cv = cv_weight(Phi[C], t[C], pen)
        Minv = np.linalg.inv(Phi[C].T @ Phi[C] / len(C) + lam * np.diag(pen))
        self.th = {"controls": Minv @ (Phi[C].T @ t[C] / len(C))}
        pc, bc, yc, vc = Phi[C].mean(0), base[C].mean(), y[C].mean(), y[C].var(ddof=1) / len(C)
        D = np.array([Phi[T].mean(0) - pc for _n, T in events])
        dy = np.array([y[T].mean() - yc for _n, T in events])
        se2 = np.array([y[T].var(ddof=1) / len(T) + vc for _n, T in events])
        r = scale * dy - np.array([base[T].mean() - bc for _n, T in events])
        S = len(C) * scale ** 2 * se2 / cv                              # the uncertainty of a row, on the scale of the fit
        G = D @ Minv @ D.T
        K = Minv @ D.T @ np.linalg.inv(G + np.diag(S))
        before = r - D @ self.th["controls"]
        self.th["rows"] = self.th["controls"] + K @ before
        after = r - D @ self.th["rows"]
        tr = np.r_[C, np.concatenate([T for _n, T in events])]
        lam2, _cv2 = cv_weight(Phi[tr], t[tr], pen)
        self.th["cells"] = np.linalg.solve(Phi[tr].T @ Phi[tr] / len(tr) + lam2 * np.diag(pen), Phi[tr].T @ t[tr] / len(tr))
        self.log = {"n_controls": int(len(C)), "ridge_weight": lam, "cv_error": cv, "ridge_weight_cells": lam2,
                    "rows": [{"event": n, "n_cells": int(len(T)), "measured_change": float(d_), "standard_error": float(np.sqrt(s)),
                              "share_of_the_residual_closed": float(1 - abs(a_) / abs(b_)) if abs(b_) > 1e-12 else 1.0}
                             for (n, T), d_, s, a_, b_ in zip(events, dy, se2, after, before)]}

    def gate(self, x):
        return self.clf.predict_proba(x)[:, 1]

    def features(self, x, o):
        g = self.gate(x)
        return np.hstack([x - self.xc, np.ones((len(x), 1))]) * g[:, None], (1 - g) * o

    def output(self, kind, x, o):
        """The PD-L1 output after the repair, log scale."""
        Phi, base = self.features(x, o)
        return base + Phi @ self.th[kind]


# --- the reading ------------------------------------------------------------------------------------------------------
def pooled(y, p, A, tg):
    """u of the own output over the knockout cells of all genes: the knockout cells are resampled by gene (as every pooled reading
    since runs/teg_training_20261119), the control cells by cell. `by_cell` resamples the knockout cells one by one."""
    u, slope = own_output_u(y, p, A, tg)
    u2, _s = own_output_u(y, p, A, None)
    T, C = A == 1, A == 0
    iv = interval(u[1:])
    return {"u": float(u[0]), "ci": iv, "reading": reading(iv), "ci_by_cell": interval(u2[1:]), "reading_by_cell": reading(interval(u2[1:])), "slope": slope,
            "valid": bool(slope > 0 and np.sign(p[T].mean() - p[C].mean()) == np.sign(y[T].mean() - y[C].mean()))}, u


def one_event(y, p, A):
    u, slope = own_output_u(y, p, A, None)
    iv = interval(u[1:])
    T, C = A == 1, A == 0
    return {"u": float(u[0]), "ci": iv, "reading": reading(iv), "slope": slope, "n_knockout": int(T.sum()), "n_control": int(C.sum()),
            "valid": bool(slope > 0 and np.sign(p[T].mean() - p[C].mean()) == np.sign(y[T].mean() - y[C].mean())),
            "measured_change": float(y[T].mean() - y[C].mean())}


def setting(y, A, tg, genes, preds, pool):
    """preds[kind][gene]: the output for every cell under the repair that reads that gene (nan where it does not read).
    pool[kind]: one output per cell for the pooled reading. Returns the table of a setting, with the contrasts on the same resamples."""
    out, draws = {}, {}
    use = np.isfinite(next(iter(pool.values())))
    for kind, p in pool.items():
        out[kind], draws[kind] = pooled(y[use], p[use], A[use], tg[use])
        out[kind]["events"] = {}
        for g in genes:
            q = preds[kind][g]
            m = np.isfinite(q) & (((A == 1) & (tg == g)) | (A == 0))
            out[kind]["events"][g] = one_event(y[m], q[m], A[m])
    for a, b in (("rows", "controls"), ("rows", "cells"), ("rows", "as_trained")):
        if a in draws and b in draws:
            d = np.abs(draws[a]) - np.abs(draws[b])
            out[a][f"abs_u_minus_{b}"] = {"difference": float(d[0]), "ci": interval(d[1:])}
    return out


def bm_score(out_t, j, yt, col):
    pred = np.exp(out_t)
    pred[:, j] = np.exp(col)
    rs = [np.corrcoef(yt[:, c], pred[:, c])[0, 1] for c in range(out_t.shape[1]) if yt[:, c].std() > 1e-8]
    return float(np.mean(rs)), float(np.corrcoef(yt[:, j], pred[:, j])[0, 1])


def guard(reps, xt, out_t, j, yt):
    """The head's own score on the bone-marrow test cells with the repaired PD-L1 output, mean over the repairs of a setting."""
    m0, p0 = bm_score(out_t, j, yt, out_t[:, j])
    res = {}
    for kind in KINDS:
        sc = [bm_score(out_t, j, yt, r.output(kind, xt, out_t[:, j])) for r in reps]
        res[kind] = {"pdl1_as_trained": p0, "pdl1_now": float(np.mean([s[1] for s in sc])), "mean_fall": float(np.mean([m0 - s[0] for s in sc])),
                     "pdl1_fall_largest": float(max(p0 - s[1] for s in sc)), "ok": bool(max(p0 - s[1] for s in sc) <= GUARD and max(m0 - s[0] for s in sc) <= GUARD)}
    res["gate_on_bone_marrow_test_cells"] = {"mean": float(np.mean([r.gate(xt).mean() for r in reps])), "highest": float(max(r.gate(xt).max() for r in reps))}
    return res


def judge(table, grd):
    """The rule of DESIGN.md for one setting, on the `rows` repair."""
    r = table["rows"]
    if r["reading"] == "supported" and r["valid"] and grd["rows"]["ok"]:
        return "holds"
    if r["reading"] == "rejected" or not r["valid"] or not grd["rows"]["ok"]:
        return "does not hold"
    return "undecided"


# --- one condition: same events on new cells, and a new gene --------------------------------------------------------------
def within(x, y, o, A, tg, genes, xb, bm, seed=SEED):
    """Two halves of the cells (knockout cells by gene, control cells at random): fit on one, read the other, and swap.
    same_events: the rows of all genes, read on the other half's cells of the same genes.
    new_gene: for each gene, the rows of the other genes, read on the other half's cells of that gene.
    seed: the split into halves (the resamples of the reading do not change with it)."""
    rng = np.random.default_rng(seed)
    first = halves(np.where(A == 1, tg, "control"), rng)
    grp = np.full(len(A), -1)
    ci_ = rng.permutation(np.where(A == 0)[0])
    grp[ci_] = np.arange(len(ci_)) % len(genes)
    nan = lambda: np.full(len(A), np.nan)
    P = {s: {k: {g: nan() for g in genes} for k in KINDS} for s in ("same_events", "new_gene")}
    pool = {s: {k: nan() for k in KINDS} for s in P}
    reps, logs = {s: [] for s in P}, {s: [] for s in P}
    gates = []
    for fit in (first, ~first):
        rd = ~fit
        C = np.where(fit & (A == 0))[0]
        ev = [(g, np.where(fit & (A == 1) & (tg == g))[0]) for g in genes]
        rep = Repair(x, y, o, C, ev, xb)
        reps["same_events"].append(rep); logs["same_events"].append(rep.log)
        gates.append({"controls_read": float(rep.gate(x[rd & (A == 0)]).mean()), "knockout_read": float(rep.gate(x[rd & (A == 1)]).mean()), "knockout_read_lowest": float(rep.gate(x[rd & (A == 1)]).min())})
        for k in KINDS:
            q = rep.output(k, x[rd], o[rd])
            pool["same_events"][k][rd] = q
            for g in genes:
                P["same_events"][k][g][rd] = q
        for kg, g in enumerate(genes):
            rep = Repair(x, y, o, C, [e for e in ev if e[0] != g], xb)
            reps["new_gene"].append(rep); logs["new_gene"].append({"held_out": g, **rep.log})
            own = rd & (((A == 1) & (tg == g)) | (grp == kg))
            for k in KINDS:
                q = rep.output(k, x[rd], o[rd])
                P["new_gene"][k][g][rd] = q
                pool["new_gene"][k][own] = rep.output(k, x[own], o[own])
    out = {}
    for s in P:
        P[s]["as_trained"], pool[s]["as_trained"] = {g: o for g in genes}, o
        out[s] = setting(y, A, tg, genes, P[s], pool[s])
        if bm is not None:
            out[s]["bone_marrow"] = guard(reps[s], *bm)
            out[s]["outcome"] = judge(out[s], out[s]["bone_marrow"])
        out[s]["fits"] = logs[s]
    out["gate"] = gates
    out["n"] = {"knockout": int((A == 1).sum()), "control": int((A == 0).sum()), "by_gene": {g: int(((A == 1) & (tg == g)).sum()) for g in genes}}
    return out


# --- another condition: the repair computed where the table was written, read on the new cells ----------------------------------
def across(src, new, genes, xb, bm):
    """src: the cells of the table (x, y, o, A, tg, C: the control cells of the fit). new: the cells to read (x, y, o, A, tg)."""
    x, y, o, A, tg, C = src
    xn, yn, on, An, tn = new
    rng = np.random.default_rng(SEED)
    grp = np.full(len(An), -1)
    ci_ = rng.permutation(np.where(An == 0)[0])
    grp[ci_] = np.arange(len(ci_)) % len(genes)
    ev = [(g, np.where((A == 1) & (tg == g))[0]) for g in genes]
    full = Repair(x, y, o, C, ev, xb)
    nan = lambda: np.full(len(An), np.nan)
    P = {"same_genes": {k: {g: full.output(k, xn, on) for g in genes} for k in KINDS}, "new_gene": {k: {g: nan() for g in genes} for k in KINDS}}
    pool = {"same_genes": {k: full.output(k, xn, on) for k in KINDS}, "new_gene": {k: nan() for k in KINDS}}
    reps, logs = {"same_genes": [full], "new_gene": []}, {"same_genes": [full.log], "new_gene": []}
    for kg, g in enumerate(genes):
        rep = Repair(x, y, o, C, [e for e in ev if e[0] != g], xb)
        reps["new_gene"].append(rep); logs["new_gene"].append({"held_out": g, **rep.log})
        own = ((An == 1) & (tn == g)) | (grp == kg)
        for k in KINDS:
            q = rep.output(k, xn, on)
            P["new_gene"][k][g] = q
            pool["new_gene"][k][own] = q[own]
    out = {}
    for s in P:
        P[s]["as_trained"], pool[s]["as_trained"] = {g: on for g in genes}, on
        out[s] = setting(yn, An, tn, genes, P[s], pool[s])
        out[s]["bone_marrow"] = guard(reps[s], *bm)
        out[s]["outcome"] = judge(out[s], out[s]["bone_marrow"])
        out[s]["fits"] = logs[s]
    out["gate"] = {"controls_read": float(full.gate(xn[An == 0]).mean()), "knockout_read": float(full.gate(xn[An == 1]).mean()), "knockout_read_lowest": float(full.gate(xn[An == 1]).min())}
    return out, full


def other_splits(x, y, o, A, tg, genes, xb):
    """The same reading with four other splits into halves: pooled values only. The rule reads the first split."""
    out = {}
    for seed in range(1, 5):
        r = within(x, y, o, A, tg, genes, xb, None, seed)
        out[str(seed)] = {s: {k: {q: r[s][k][q] for q in ("u", "ci", "reading", "ci_by_cell", "reading_by_cell")} for k in KINDS} for s in ("same_events", "new_gene")}
        print(f"  split {seed}: " + " | ".join(f"{s} " + " ".join(f"{k} {r[s][k]['u']:.3f} {r[s][k]['reading'][:5]}" for k in KINDS) for s in ("same_events", "new_gene")))
    return out


def show(name, table, genes):
    f = lambda v: f"{v:6.3f}"
    print(f"\n{name}" + (f"   -> {table['outcome']}" if "outcome" in table else ""))
    for kind in ["as_trained"] + KINDS:
        r = table[kind]
        print(f"  {kind:10s} {f(r['u'])} [{f(r['ci'][0])}, {f(r['ci'][1])}] {r['reading']:10s} by cell [{f(r['ci_by_cell'][0])}, {f(r['ci_by_cell'][1])}] {r['reading_by_cell']:10s} slope {r['slope']:6.2f} {'ok' if r['valid'] else 'NO'} | "
              + " ".join(f"{g}:{r['events'][g]['u']:.2f}" for g in genes)
              + (f" | bone marrow PD-L1 {table['bone_marrow'][kind]['pdl1_as_trained']:.3f} -> {table['bone_marrow'][kind]['pdl1_now']:.3f}" if kind in table.get("bone_marrow", {}) else ""))
    for b in ("controls", "cells"):
        d = table["rows"][f"abs_u_minus_{b}"]
        print(f"  rows - {b:8s} {f(d['difference'])} [{f(d['ci'][0])}, {f(d['ci'][1])}]")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--development", action="store_true")
    a = ap.parse_args()
    P0, names, j = weights()
    b, t_ = np.load(BM_FIT), np.load(BM_TEST)
    xb, xt, yt = unit(b["z_train"]).astype(np.float64), unit(t_["z"]).astype(np.float64), t_["adt"].astype(np.float64)
    bm = (xt, forward_np(P0, xt)[1], j, yt)
    head = lambda x: forward_np(P0, x)[1][:, j]
    RES.mkdir(exist_ok=True)
    if a.development:
        res = {"note": "development: the code of the test on cells already read; nothing here is a test", "margin": MARGIN, "grid": GRID, "data_sets": {}}
        for ds, genes in TARGETS.items():
            d = load(ds)
            x = d["x"].astype(np.float64)
            res["data_sets"][ds] = within(x, d["y"], head(x), d["A"], d["tg"], genes, xb, bm)
            for s in ("same_events", "new_gene"):
                show(f"{CELLS[ds]} (already read), {s}", res["data_sets"][ds][s], genes)
            res["data_sets"][ds]["other_splits"] = other_splits(x, d["y"], head(x), d["A"], d["tg"], genes, xb)
        (RES / "development.json").write_text(json.dumps(res, indent=1))
        return
    target = RES / "reading.json"
    if target.exists():
        raise SystemExit("results/reading.json exists: the reading was made once and is not made again")
    genes = TARGETS["frangieh2021"]
    cells = Path(os.environ["ANM_TEST_CELLS"])
    co = np.load(cells / "coculture.npz")
    xn, yn, An, tn = unit(co["Z"]).astype(np.float64), co["readout__PDL1"].astype(np.float64), co["A"], co["target"].astype(str)
    d = load("frangieh2021")
    assert not np.isin(co["ids"].astype(str), d["ids"]).any() and set(np.unique(tn[An == 1])) == set(genes)
    x = d["x"].astype(np.float64)
    res = {"margin": MARGIN, "guard": GUARD, "grid": GRID, "resamples": 1000,
           "cells": {"condition_of_the_table": "melanoma, IFN-gamma (read by earlier runs)", "condition_read_here": "melanoma, co-culture (not read before)"}}
    # 1. within the co-culture condition
    res["coculture"] = within(xn, yn, head(xn), An, tn, genes, xb, bm)
    res["coculture"]["other_splits"] = other_splits(xn, yn, head(xn), An, tn, genes, xb)
    # 2. the repair computed in the IFN-gamma condition, read on the co-culture cells
    src = (x, d["y"], head(x), d["A"], d["tg"], np.where(d["train_ctl"])[0])
    res["from_ifn"], full = across(src, (xn, yn, head(xn), An, tn), genes, xb, bm)
    # 3. the library of the IFN-gamma condition: a gene outside the pathway that moves PD-L1 (the probe), and the genes that do not
    lib = np.load(cells / "library.npz")
    newc = ~np.isin(lib["ids"].astype(str), d["ids"])
    xl, yl, tl = unit(lib["Z"][newc]).astype(np.float64), lib["readout__PDL1"][newc].astype(np.float64), lib["target"][newc].astype(str)
    ol = head(xl)
    ac = np.where(d["audit_ctl"])[0]                                     # control cells of the IFN-gamma condition that no fit has seen
    outs = {"as_trained": (ol, head(x)[ac])} | {k: (full.output(k, xl, ol), full.output(k, x[ac], head(x)[ac])) for k in KINDS}
    T = tl == PROBE
    res["probe"] = {"gene": PROBE, "arms": {k: one_event(np.r_[yl[T], d["y"][ac]], np.r_[pk[T], pc_], np.r_[np.ones(T.sum(), int), np.zeros(len(ac), int)]) for k, (pk, pc_) in outs.items()}}
    others = sorted(g for g in np.unique(tl) if g != PROBE)
    yc, sdc = d["y"][ac].mean(), d["y"][ac].std(ddof=1)
    meas = np.array([yl[tl == g].mean() - yc for g in others])
    spec = {"n_genes": len(others), "n_cells": int(np.isin(tl, others).sum()), "control_sd": float(sdc), "measured": {"median_abs": float(np.median(np.abs(meas))), "largest_abs": float(np.abs(meas).max())}, "arms": {}}
    for k, (pk, pc_) in outs.items():
        slope = np.cov(pc_, d["y"][ac])[0, 1] / pc_.var(ddof=1)
        pr = np.array([slope * (pk[tl == g].mean() - pc_.mean()) for g in others])
        spec["arms"][k] = {"predicted_median_abs": float(np.median(np.abs(pr))), "predicted_largest_abs": float(np.abs(pr).max()), "largest_gene": others[int(np.abs(pr).argmax())],
                           "rms_of_predicted_minus_measured": float(np.sqrt(((pr - meas) ** 2).mean())), "rms_of_measured": float(np.sqrt((meas ** 2).mean())),
                           "correlation_with_measured": float(np.corrcoef(pr, meas)[0, 1]), "genes_predicted_beyond_0.2": int((np.abs(pr) > 0.2).sum()),
                           "genes_predicted_beyond_0.5": int((np.abs(pr) > 0.5).sum())}
    res["other_genes"] = spec
    res["outcomes"] = {"same_events_new_cells": res["coculture"]["same_events"]["outcome"], "new_gene": res["coculture"]["new_gene"]["outcome"],
                       "new_condition": res["from_ifn"]["same_genes"]["outcome"], "new_gene_and_new_condition": res["from_ifn"]["new_gene"]["outcome"]}
    target.write_text(json.dumps(res, indent=1))
    show("co-culture, same events on new cells", res["coculture"]["same_events"], genes)
    show("co-culture, a new gene", res["coculture"]["new_gene"], genes)
    show("from IFN-gamma to co-culture, same genes", res["from_ifn"]["same_genes"], genes)
    show("from IFN-gamma to co-culture, a new gene", res["from_ifn"]["new_gene"], genes)
    print("\ngate:", json.dumps(res["coculture"]["gate"]), json.dumps(res["from_ifn"]["gate"]))
    print("probe", PROBE, {k: (round(v["u"], 3), [round(c, 3) for c in v["ci"]], v["reading"]) for k, v in res["probe"]["arms"].items()})
    print("other genes:", json.dumps(spec, indent=1))
    print("\noutcomes:", json.dumps(res["outcomes"]))


if __name__ == "__main__":
    main()
