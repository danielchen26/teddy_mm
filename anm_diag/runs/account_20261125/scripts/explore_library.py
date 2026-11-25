#!/usr/bin/env python3
"""Section 3 (explored, not counted): the 206 knockout events of the melanoma library under IFN-gamma, event by event.
results/explored_library.json.

  python explore_library.py --library LIBRARY.npz --own-gene OWN_GENE.json --site-counts SITE_COUNTS.npz

LIBRARY as for library.py. OWN_GENE and SITE_COUNTS are written by own_gene.py and site_counts.py from the RNA counts of
the screen (labels and counts only).
For every event: how far its shift of the mean state stands out, its change of every readout, its own transcript.
For the events that change a readout most: what a readout of control cells reaches from the state, with the count of the
readout's own genes beside it, and with the count of the event's own site beside it; what is left with the event's own
record under the rule of run 40; and the forecast of the direct part by an edge known beforehand, one pair left out.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from _read import RUN, RUNS, unit
from repair import cross, cycle, u_of

OWN = {"PDL1": ["CD274"], "MHC1": ["HLA-A", "HLA-B", "HLA-C", "B2M"], "HLAE": ["HLA-E", "B2M"]}
LAMS = 10.0 ** np.arange(-3, 3.01, 0.5)


def fit(F, y, C, pen):
    """A ridge readout of control cells; the weight by five folds; pen: one for a penalised column."""
    nc = len(C)
    Fc, yc, idx, best = F[C], y[C] - y[C].mean(), np.arange(nc) % 5, None
    for lam in LAMS:
        err = 0.0
        for f in range(5):
            tr = idx != f
            th = np.linalg.solve(Fc[tr].T @ Fc[tr] / tr.sum() + lam * np.diag(pen), Fc[tr].T @ yc[tr] / tr.sum())
            err += ((yc[~tr] - Fc[~tr] @ th) ** 2).sum()
        if best is None or err < best[0]:
            best = (err, lam)
    return np.linalg.solve(Fc.T @ Fc / nc + best[1] * np.diag(pen), Fc.T @ yc / nc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--own-gene", type=Path, required=True)
    ap.add_argument("--site-counts", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    d, own, sc = np.load(a.library), json.loads(a.own_gene.read_text())["genes"], np.load(a.site_counts)
    assert (sc["ids"].astype(str) == d["ids"].astype(str)).all()
    Z = unit(d["Z"]).astype(np.float64)
    tg, Cm = d["target"].astype(str), d["A"] == 0
    C, nc, p = np.where(Cm)[0], int(Cm.sum()), Z.shape[1]
    X = (Z - Z[C].mean(0)) / Z[C].std(0)
    ev = [t for t in sorted(set(tg)) if t != "control"]
    Y = {k: d["readout__" + k].astype(np.float64) for k in OWN}
    genes, cnt, tot = list(sc["genes"].astype(str)), sc["counts"].astype(np.float64), sc["total_count"].astype(np.float64)
    lc = np.log1p(cnt / tot[:, None] * 1e4)
    lc = (lc - lc[C].mean(0)) / np.where(lc[C].std(0) > 1e-9, lc[C].std(0), 1.0)
    pres = (cnt > 0).astype(float)
    pres = (pres - pres[C].mean(0)) / np.where(pres[C].std(0) > 1e-9, pres[C].std(0), 1.0)
    events = {}
    for t in ev:
        T = tg == t
        n = int(T.sum())
        r = {"cells": n, "shift_of_the_state_times_its_noise": float((X[T].mean(0) ** 2).sum() / (p * (1 / n + 1 / nc)))}
        for k, y in Y.items():
            dy, se = y[T].mean() - y[C].mean(), np.sqrt(y[T].var(ddof=1) / n + y[C].var(ddof=1) / nc)
            r[k] = {"change": float(dy), "errors": float(dy / se)}
        q = own.get(t)
        if q and "cells" in q and q["mean_in_control_cells"] * nc >= 20:
            ratio = q["mean_in_its_cells"] / q["mean_in_control_cells"]
            se = ratio * np.sqrt(1 / max(q["mean_in_its_cells"] * n, 1) + 1 / (q["mean_in_control_cells"] * nc))      # as for counts; a lower bound of the error
            r["own_transcript"] = {"per_cell_in_its_cells": q["mean_in_its_cells"], "per_cell_in_control_cells": q["mean_in_control_cells"], "ratio": ratio, "errors_from_one": float((ratio - 1) / max(se, 1e-9)),
                                   "detected_in_its_cells": q["detected_in_its_cells"], "detected_in_control_cells": q["detected_in_control_cells"]}
        events[t] = r
    have = [r for r in events.values() if "own_transcript" in r]
    ratio = np.array([r["own_transcript"]["ratio"] for r in have])
    lower = (ratio < 0.8) & (np.array([r["own_transcript"]["errors_from_one"] for r in have]) < -3)
    moved = np.array([r["shift_of_the_state_times_its_noise"] for r in have])
    res = {"note": "development units; explored: reported, not counted", "cells": int(len(tg)), "control_cells": nc, "events": len(ev),
           "events_whose_shift_of_the_state_is_beyond": {"2x_its_noise": int(sum(r["shift_of_the_state_times_its_noise"] > 2 for r in events.values())),
                                                         "5x_its_noise": int(sum(r["shift_of_the_state_times_its_noise"] > 5 for r in events.values()))},
           "own_transcript": {"events_with_enough_counts": len(have), "median_ratio": float(np.median(ratio)), "lowered_below_0.8_and_beyond_3_errors": int(lower.sum()),
                              "of_those_the_state_moves_beyond_2x_noise": int((lower & (moved > 2)).sum()), "state_moves_beyond_2x_noise": int((moved > 2).sum())},
           "by_event": events, "the_events_that_change_a_readout_most": {}}
    rng = np.random.default_rng(0)
    for k in ("PDL1", "MHC1"):
        y = Y[k]
        th = fit(X, y, C, np.ones(p))
        gi = [genes.index(g) for g in OWN[k]]
        Fo = np.hstack([X, lc[:, gi], pres[:, gi]])
        tho = fit(Fo, y, C, np.r_[np.ones(p), np.zeros(2 * len(gi))])
        rows = {}
        for t in ev:
            if abs(events[t][k]["errors"]) < 4:
                continue
            T, dy = tg == t, events[t][k]["change"]
            r = {"change": dy, "reached_from_the_state": float(X[T].mean(0) @ th / dy), "with_the_readouts_own_genes": float(Fo[T].mean(0) @ tho / dy)}
            if t in genes:
                j = genes.index(t)
                Fs = np.hstack([X, lc[:, [j]], pres[:, [j]]])
                r["with_the_events_own_site"] = float(Fs[T].mean(0) @ fit(Fs, y, C, np.r_[np.ones(p), 0, 0]) / dy)
            m = Cm | T                                           # its own record, the rule of run 40, fitted on one half and read on the other
            first = rng.permutation(int(m.sum())) < int(m.sum()) // 2
            o, _inside, fits = cross(Z[m], y[m], y[m], T[m], np.where(T[m], t, "control"), first, [t])
            r["u_with_no_record"] = float(u_of(y[m], o["controls"], T[m])[0][0])
            u = u_of(y[m], o["rows"], T[m])[0]
            r["u_with_its_own_record_run_40"] = {"u": float(u[0]), "ci": cycle.interval(u[1:]), "row_closed_on_the_fitted_cells": float(np.mean([f["rows"][0]["share_of_the_residual_closed"] for f in fits]))}
            rows[t] = r
            print(f"{k:5s} {t:7s} change {dy:+.2f} | state {r['reached_from_the_state']:4.0%} | + own genes {r['with_the_readouts_own_genes']:4.0%} | + own site {r.get('with_the_events_own_site', float('nan')):4.0%} | u no record {r['u_with_no_record']:+.2f} | own record {u[0]:+.2f}", flush=True)
        res["the_events_that_change_a_readout_most"][k] = rows
    # the direct part of every pair, in spreads of the readout; an edge known beforehand; one pair left out
    rec = []
    for k, y in Y.items():
        th, sd = fit(X, y, C, np.ones(p)), y[C].std()
        for t in ev:
            T = tg == t
            rec.append((k, t, float((events[t][k]["change"] - X[T].mean(0) @ th) / sd), float(events[t][k]["change"] / events[t][k]["errors"] / sd), t in OWN[k]))
    dr, se, ed = np.array([r[2] for r in rec]), np.array([r[3] for r in rec]), np.array([r[4] for r in rec])
    nlp = lambda x, m, s: 0.5 * np.log(2 * np.pi * s * s) + 0.5 * ((x - m) / s) ** 2
    pairs, t0, t1 = [], 0.0, 0.0
    for i in np.where(ed)[0]:
        o = np.arange(len(rec)) != i
        m0, s0 = np.median(dr[o]), np.sqrt(np.mean((dr[o] - np.median(dr[o])) ** 2) + se[i] ** 2)
        oe = o & ed
        m1, s1 = dr[oe].mean(), np.sqrt(dr[oe].var(ddof=1) * (1 + 1 / oe.sum()) + se[i] ** 2)
        t0, t1 = t0 + nlp(dr[i], m0, s0), t1 + nlp(dr[i], m1, s1)
        pairs.append({"event": rec[i][1], "readout": rec[i][0], "direct": float(dr[i]), "without_edges": {"forecast": float(m0), "sd": float(s0), "errors": float((dr[i] - m0) / s0)},
                      "with_the_edge": {"forecast": float(m1), "sd": float(s1), "errors": float((dr[i] - m1) / s1)}})
    res["an_edge_known_beforehand"] = {"pairs": int(len(rec)), "with_an_edge": int(ed.sum()), "direct_without_an_edge": {"median_size": float(np.median(np.abs(dr[~ed]))), "95th_percentile": float(np.percentile(np.abs(dr[~ed]), 95))},
                                       "one_left_out": pairs, "sum_of_minus_log_density": {"without_edges": float(t0), "with_edges": float(t1)},
                                       "without_an_edge_and_beyond_5_errors": sorted([{"event": r[1], "readout": r[0], "direct": r[2], "errors": r[2] / r[3]} for r in rec if not r[4] and abs(r[2] / r[3]) > 5], key=lambda q: -abs(q["errors"]))}
    (a.results / "explored_library.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps({k: res[k] for k in ("events_whose_shift_of_the_state_is_beyond", "own_transcript")}), json.dumps(res["an_edge_known_beforehand"]["sum_of_minus_log_density"]))


if __name__ == "__main__":
    main()
