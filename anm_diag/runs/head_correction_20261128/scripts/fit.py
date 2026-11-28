#!/usr/bin/env python3
"""The arms: three rows of the head's last layer, fitted on the teach control cells, with or without the loop's rows.
No cell of a drawn unit and no judge control cell is touched here.

  python fit.py --cells CELLS.tsv --units units.json --states STATES... --packs PACK... --out rows.npz
                [--steps 4000] [--selftest]

H0 the head as shipped; H1 the rows fitted on the teach control cells; H2 the same with the loop's rows (every event
recorded in these two conditions: the shift the output shows must equal the shift the loop measured); H3 as H2 with
the measured shifts permuted among the events of a readout (seed 20261128).

The fit is by Adam on the head's own output (log of softplus), toward the measured readout, with a ridge term that
pulls the row toward the one the head was shipped with. Its weight is chosen inside the teach control cells alone:
they are split in two by a hash of the cell's name, every weight of a fixed grid is fitted on one half and measured
on the other, and the one with the smallest error is used for every arm of that readout. The weight of the loop's
rows is chosen the same way, by splitting the recorded events in two: the fit is made with the rows of one half and
judged by how well it gives the shifts of the other. No judge control cell and no drawn unit takes part. The fit is carried out with the activation centred and scaled by the teach control cells
and the row written back into the head's own numbers afterwards; that only conditions the search and changes no
result. --selftest fits a made-up world and checks the code alone.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _hc as K  # noqa: E402

LR = 0.05
LAMBDAS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
RHOS = (0.03, 0.1, 0.3, 1.0, 3.0)
inner = lambda n: int(__import__("hashlib").sha256(("lam|" + str(n)).encode()).hexdigest(), 16) % 2 == 0


def adam(grad, w0, steps, lr=LR):
    w = w0.copy()
    m, v = np.zeros_like(w), np.zeros_like(w)
    for t in range(1, steps + 1):
        g = grad(w)
        m, v = 0.9 * m + 0.1 * g, 0.999 * v + 0.001 * g * g
        w -= lr * (m / (1 - 0.9 ** t)) / (np.sqrt(v / (1 - 0.999 ** t)) + 1e-8)
    return w


def fit_one(a_ctl, y_ctl, rows, w0, b0, lam, rho, steps):
    """One readout. rows: list of (activations of the event's cells, the measured shift). Empty for H1."""
    scale = float(np.linalg.norm(w0) ** 2 + b0 ** 2) + 1e-12
    mc = a_ctl.mean(0)

    def value(p):
        w, b = p[:-1], p[-1]
        o = K.out_of(a_ctl, w, b)
        loss = float(((o - y_ctl) ** 2).mean()) + lam * float(((p - np.r_[w0, b0]) ** 2).sum()) / scale
        for A, s in rows:
            loss += rho / len(rows) * (float(K.out_of(A, w, b).mean() - K.out_of(mc[None, :], w, b)[0]) - s) ** 2
        return loss

    def grad(p):
        w, b, g = p[:-1], p[-1], np.zeros_like(p)
        z = a_ctl @ w + b
        sp = np.logaddexp(0.0, z)
        d = (np.log(np.maximum(sp, 1e-6)) - y_ctl) * 2.0 / len(z) / np.maximum(sp, 1e-6) / (1.0 + np.exp(-z))
        g[:-1] = a_ctl.T @ d
        g[-1] = d.sum()
        g += 2.0 * lam * (p - np.r_[w0, b0]) / scale
        for A, s in rows:
            zc, za = float(mc @ w + b), A @ w + b
            oc = np.log(max(np.logaddexp(0.0, zc), 1e-6))
            oa = np.log(np.maximum(np.logaddexp(0.0, za), 1e-6))
            e = 2.0 * rho / len(rows) * (oa.mean() - oc - s)
            da = 1.0 / np.maximum(np.logaddexp(0.0, za), 1e-6) / (1.0 + np.exp(-za)) / len(za)
            dc = 1.0 / max(np.logaddexp(0.0, zc), 1e-6) / (1.0 + np.exp(-zc))
            g[:-1] += e * (A.T @ da - dc * mc)
            g[-1] += e * (da.sum() - dc)
        return g

    p = adam(grad, np.r_[w0, b0], steps)
    return p[:-1], float(p[-1]), {"loss_at_the_start": value(np.r_[w0, b0]), "loss_at_the_end": value(p)}


def selftest():
    rng = np.random.default_rng(0)
    n, p = 3000, 512
    a = rng.normal(size=(n, p)) / np.sqrt(p)
    w = rng.normal(size=p); b = 1.0
    y = K.out_of(a, w, b) + rng.normal(size=n) * 0.05
    w0 = rng.normal(size=p) * 0.1
    wf, bf, info = fit_one(a, y, [], w0, 0.0, 1e-6, 0.0, 3000)
    e0, e1 = float(((K.out_of(a, w0, 0.0) - y) ** 2).mean()), float(((K.out_of(a, wf, bf) - y) ** 2).mean())
    A = rng.normal(size=(200, p)) / np.sqrt(p) + 0.3 * w / np.linalg.norm(w)
    s = float(K.out_of(A, w, b).mean() - K.out_of(a.mean(0)[None, :], w, b)[0])
    w2, b2, _ = fit_one(a, y, [(A, s)], w0, 0.0, 1e-6, 10.0, 3000)
    got = float(K.out_of(A, w2, b2).mean() - K.out_of(a.mean(0)[None, :], w2, b2)[0])
    print(f"selftest: error on the cells {e0:.3f} -> {e1:.4f}; the row asked for {s:+.3f}, the fit gives {got:+.3f}")
    assert e1 < 0.01 and abs(got - s) < 0.05, "the fit does not reach a world it could reach"
    return {"error_before": e0, "error_after": e1, "row_asked": s, "row_given": got}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path)
    ap.add_argument("--units", type=Path)
    ap.add_argument("--states", type=Path, nargs="+")
    ap.add_argument("--packs", type=Path, nargs="+")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--steps", type=int, default=4000)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        print(json.dumps(selftest(), indent=1))
        return
    if a.out.exists():
        raise SystemExit(f"{a.out.name} exists: the arms are fitted once")
    C = K.cells(a.cells)
    rec = K.recorded_units(C)
    St, Pk = K.B3.States(a.states), K.B3.Packs(a.packs)
    rng = np.random.default_rng(K.SEED)
    out, report = {}, {"lambdas_tried": list(LAMBDAS), "rhos_tried": list(RHOS), "steps": a.steps, "lr": LR, "conditions": {}}
    act, meas = {}, {}
    for c in K.CONDITIONS:
        g = C[c]
        keep = ~g["drawn"]                                                   # control cells and the recorded events; no drawn cell
        ids = g["ids"][keep]
        A, rows0 = K.head(St.rows(ids), ids)
        clr, _log, _bg = Pk.readouts(K.STUDY, ids, list(K.READOUTS))
        act[c] = {"A": A, "unit": g["unit"][keep], "teach": g["teach"][keep], "control": g["control"][keep], "y": clr, "inner": np.array([inner(i) for i in ids])}
        report["conditions"][c] = {"cells": int(keep.sum()), "teach_control_cells": int((g["control"] & g["teach"])[keep].sum()), "recorded_events": len(rec[c])}
        meas[c] = rows0
    rows0 = meas[K.CONDITIONS[0]]
    A_all = np.vstack([act[c]["A"][act[c]["control"] & act[c]["teach"]] for c in K.CONDITIONS])
    mu, sd = A_all.mean(0), A_all.std(0) + 1e-9                          # the activation as the teach control cells see it; the search is made there
    for c in K.CONDITIONS:
        act[c]["A"] = (act[c]["A"] - mu) / sd
    report["conditioning"] = {"the_activation_is_centred_and_scaled_by_the_teach_control_cells": True}
    for r in K.READOUTS:
        A_ctl = np.vstack([act[c]["A"][act[c]["control"] & act[c]["teach"]] for c in K.CONDITIONS])
        y_ctl = np.concatenate([act[c]["y"][r][act[c]["control"] & act[c]["teach"]] for c in K.CONDITIONS])
        true_rows, names = [], []
        for c in K.CONDITIONS:
            d, mc = act[c], act[c]["control"] & act[c]["teach"]
            base = float(d["y"][r][mc].mean())
            for u in rec[c]:
                m = d["unit"] == u
                true_rows.append((d["A"][m], float(d["y"][r][m].mean() - base)))
                names.append(f"{u}@{c}")
        shifts = np.array([s for _A, s in true_rows])
        perm = rng.permutation(len(shifts))
        w0, b0 = rows0[r]
        out[f"H0__{r}__w"], out[f"H0__{r}__b"] = w0, b0
        wt0, bt0 = sd * w0, b0 + float(mu @ w0)                          # the shipped row, written in the conditioned space
        inner_a = np.concatenate([act[c]["inner"][act[c]["control"] & act[c]["teach"]] for c in K.CONDITIONS])
        tried = {}
        for lam in LAMBDAS:                                              # the weight of the ridge, chosen inside the teach control cells
            wi, bi, _ = fit_one(A_ctl[inner_a], y_ctl[inner_a], [], wt0, bt0, lam, 0.0, a.steps)
            tried[lam] = float(((K.out_of(A_ctl[~inner_a], wi, bi) - y_ctl[~inner_a]) ** 2).mean())
        lam = min(tried, key=tried.get)
        report.setdefault("lambda", {})[r] = {"chosen": lam, "error_on_the_other_half": {str(k): v for k, v in tried.items()}}
        print(f"{r}: ridge weight {lam} (error on the other half " + ", ".join(f"{k:g}:{v:.3f}" for k, v in tried.items()) + ")", flush=True)
        inner_r = np.array([inner(n) for n in names])
        got = {}
        for rho in RHOS:                                                 # the weight of the rows, chosen by the events the fit was not given
            wi, bi, _ = fit_one(A_ctl, y_ctl, [true_rows[i] for i in np.where(inner_r)[0]], wt0, bt0, lam, rho, a.steps)
            mc = A_ctl.mean(0)
            got[rho] = float(np.mean([(float(K.out_of(true_rows[i][0], wi, bi).mean() - K.out_of(mc[None, :], wi, bi)[0]) - true_rows[i][1]) ** 2 for i in np.where(~inner_r)[0]]))
        rho = min(got, key=got.get)
        report.setdefault("rho", {})[r] = {"chosen": rho, "error_on_the_other_events": {str(k): v for k, v in got.items()}}
        print(f"{r}: weight of the rows {rho} (error on the other events " + ", ".join(f"{k:g}:{v:.3f}" for k, v in got.items()) + ")", flush=True)
        fits = {"H1": [], "H2": true_rows, "H3": [(A, float(shifts[perm[i]])) for i, (A, _s) in enumerate(true_rows)]}
        for arm, rws in fits.items():
            wt, bt, info = fit_one(A_ctl, y_ctl, rws, wt0, bt0, lam, rho if rws else 0.0, a.steps)
            w, b = wt / sd, float(bt - (mu / sd) @ wt)                   # written back into the head's own numbers
            out[f"{arm}__{r}__w"], out[f"{arm}__{r}__b"] = w, b
            report.setdefault("fits", {}).setdefault(r, {})[arm] = info
            print(f"{r} {arm}: loss {info['loss_at_the_start']:.4f} -> {info['loss_at_the_end']:.4f}", flush=True)
        report.setdefault("rows", {})[r] = {"n": len(true_rows), "measured_shift": {"mean": float(shifts.mean()), "spread": float(shifts.std()), "largest": float(np.abs(shifts).max())},
                                            "correlation_of_the_shuffled_rows_with_the_true_ones": float(np.corrcoef(shifts, shifts[perm])[0, 1]), "events": names}
    np.savez_compressed(a.out, **out)
    (a.out.parent / (a.out.stem + "_report.json")).write_text(json.dumps(report, indent=1) + "\n")
    print("rows written:", a.out.name, "| sha256", K.sha(a.out)[:16])


if __name__ == "__main__":
    main()
