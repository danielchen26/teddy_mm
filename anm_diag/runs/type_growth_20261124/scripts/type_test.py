#!/usr/bin/env python3
"""Can a finer type cover a donor that has no record of its own? The reading, made once.

  python type_test.py --work WORK [--results RESULTS] [--made-up]

WORK: the work folder of runs/unit_events_20261123 (donors.npz: the cells outside the selection third, with their
donor, batch, state, measured proteins and the head's own outputs). The pairs of a donor and a protein are those the
round of 23 November admitted (results/selection.json there), and its readings are the anchors (results/round.json).

A type here is a grouping of the eight donors into two groups of four. For a pair (donor U, protein), the record of
U's type is what its three type-mates measured. DESIGN.md, Sections 3 to 5.

  part 1, the record alone   u_type = 1 - (the mates' change) / (U's change), both against the cells of the seven
                             other donors. No state and no readout enter.
  part 2, with the repair    the arm `typed`: the round's repair (anm_diag.cycle.UnitRows) fitted with no protein of
                             U, in which the mean of U is set to the mean of its type-mates. Read as the round reads
                             an output, beside `unseen` (no record of U, no type) and `rows` (U's own record).

--made-up: for a trial of the script on made-up values only; it skips the check against the recorded round.
"""
import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN, RUNS = HERE.parent, HERE.parents[1]
sys.path.insert(0, str(RUNS.parent))
from anm_diag import cycle  # noqa: E402

SEED, DRAWS, STEP, MARGIN = 0, 1000, cycle.STEP, cycle.MARGIN
unit_length = lambda z: z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)
iv = lambda d: [float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))]
verdict = lambda c: "support" if (c[0] > -MARGIN and c[1] < MARGIN) else ("reject" if (c[0] > MARGIN or c[1] < -MARGIN) else "unresolved")
step = lambda c: "a step" if c[0] > STEP else ("no step" if c[1] < STEP else "unresolved")


def groupings(donors):
    """Every way to put eight donors into two groups of four: 35."""
    d0, out = donors[0], []
    for rest in itertools.combinations(donors[1:], 3):
        a = (d0,) + rest
        out.append((a, tuple(d for d in donors if d not in a)))
    return out


def mates(g, u):
    side = g[0] if u in g[0] else g[1]
    return [d for d in side if d != u]


class Typed(cycle.UnitRows):
    """cycle.UnitRows with one more arm for the units named in `plan` ({unit: {name: its type-mates}}):
    typed|unit|name: a fit with no protein of the unit. Its readout of the variation among cells is that of the other
    units (as in `unseen`); the rows are those of every unit against the others, with the mean of the unit set to the
    mean of its type-mates. The arm `unseen` is computed again here and must equal the package's."""
    plan = {}

    def __init__(self, x, y, o, C, events, other=None):
        super().__init__(x, y, o, C, events, other)
        events = [(str(n), np.asarray(T)) for n, T in events if len(T)]
        names = [n for n, _T in events]
        A = np.concatenate([T for _n, T in events])
        scale = o[A].std() / y[A].std()
        t = o[A].mean() + scale * (y - y[A].mean())
        lam, cv = self.log["ridge_weight"], self.log["cv_error"]
        xw = np.vstack([x[T] - x[T].mean(0) for _n, T in events])
        pen = np.ones(x.shape[1]) * xw.var(0).mean()
        st = {}
        for n, T in events:
            xm, tm = x[T].mean(0), t[T].mean()
            st[n] = dict(n=len(T), sx=x[T].sum(0), st=t[T].sum(), sy=y[T].sum(), syy=(y[T] ** 2).sum(), W=(x[T] - xm).T @ (x[T] - xm), b=(x[T] - xm).T @ (t[T] - tm))
        var = lambda s1, s2, m: (s2 - s1 ** 2 / m) / np.maximum(m - 1, 1)
        for u, types in self.plan.items():
            if u not in st:
                continue
            keep = [n for n in names if n != u]
            N = float(sum(st[k]["n"] for k in keep))
            Minv = np.linalg.inv(sum(st[k]["W"] for k in keep) / N + lam * np.diag(pen))
            w = Minv @ (sum(st[k]["b"] for k in keep) / N)

            def rows(units, sy, syy, xc, tbar):
                n = np.array([st[k]["n"] for k in units], np.float64)
                sx = np.array([st[k]["sx"] for k in units])
                D = sx / n[:, None] - (sx.sum(0) - sx) / (n.sum() - n)[:, None]
                dy = sy / n - (sy.sum() - sy) / (n.sum() - n)
                se2 = var(sy, syy, n) / n + var(sy.sum() - sy, syy.sum() - syy, n.sum() - n) / (n.sum() - n)
                K = Minv @ D.T @ np.linalg.inv(D @ Minv @ D.T + np.diag(N * scale ** 2 * se2 / cv))
                w2 = w + K @ (scale * dy - D @ w)
                return np.r_[w2, tbar + (self.xc - xc) @ w2]

            sy, syy = np.array([st[k]["sy"] for k in keep]), np.array([st[k]["syy"] for k in keep])
            again = rows(keep, sy, syy, sum(st[k]["sx"] for k in keep) / N, sum(st[k]["st"] for k in keep) / N)
            assert np.allclose(again, self.th[f"unseen|{u}"], rtol=1e-6, atol=1e-8), "the arm `unseen` computed here is not the package's"
            for name, ms in types.items():                       # the mean of u is set to the mean of its type-mates; its spread to theirs
                nm = float(sum(st[k]["n"] for k in ms))
                ym = sum(st[k]["sy"] for k in ms) / nm
                vm = sum(st[k]["syy"] - st[k]["sy"] ** 2 / st[k]["n"] for k in ms) / (nm - len(ms))
                nu = float(st[u]["n"])
                sy8 = np.array([st[k]["sy"] if k != u else nu * ym for k in names])
                syy8 = np.array([st[k]["syy"] if k != u else (nu - 1) * vm + nu * ym ** 2 for k in names])
                N8 = float(sum(st[k]["n"] for k in names))
                tbar = o[A].mean() + scale * (sy8.sum() / N8 - y[A].mean())
                self.th[f"typed|{u}|{name}"] = rows(names, sy8, syy8, sum(st[k]["sx"] for k in names) / N8, tbar)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    ap.add_argument("--made-up", action="store_true")
    ap.add_argument("--selection", type=Path, default=RUNS / "unit_events_20261123" / "results" / "selection.json")
    ap.add_argument("--round", type=Path, default=RUNS / "unit_events_20261123" / "results" / "round.json")
    a = ap.parse_args()
    out_file = a.results / "types.json"
    if out_file.exists():
        raise SystemExit("types.json exists: this reading was made once and is not made again")
    d = np.load(a.work / "donors.npz")
    unit, batch = d["unit"].astype(str), d["batch"].astype(str)
    donors = sorted(set(unit))
    assert len(donors) == 8
    x, Y, O = unit_length(d["Z"].astype(np.float64)), d["Y"].astype(np.float64), d["frozen"].astype(np.float64)
    sel = json.loads(a.selection.read_text())["proteins"]
    rec = json.loads(a.round.read_text())["targets"]
    pairs = [(n, s) for n, q in sel.items() for s in q["admitted"]]
    A3 = {(n, s): bool(rec[n]["events"][s]["A3"]["holds"]) for n, s in pairs}
    need = {(n, s): rec[n]["events"][s].get("rows_minus_unseen", {}).get("reading") == "a step" for n, s in pairs}
    sets = {"own_row_needed": [p for p in pairs if A3[p] and need[p]], "a_step_to_remove": [p for p in pairs if A3[p]], "admitted": pairs}
    G = groupings(donors)
    key = lambda g: "|".join("+".join(x_.split("_")[-1] for x_ in side) for side in g)
    # --- the candidates
    b_of = {u: sorted(set(batch[unit == u])) for u in donors}
    assert all(len(v) == 1 for v in b_of.values())
    by_batch = sorted({v[0] for v in b_of.values()})
    g_batch = next(g for g in G if {b_of[u][0] for u in g[0]} == {b_of[donors[0]][0]}) if len(by_batch) == 2 and sum(b_of[u][0] == by_batch[0] for u in donors) == 4 else None
    M = np.array([x[unit == u].mean(0) for u in donors])          # the donors' mean states: no protein enters
    def separation(g):
        i0, i1 = [donors.index(u) for u in g[0]], [donors.index(u) for u in g[1]]
        within = np.mean([((M[i] - M[i].mean(0)) ** 2).sum(1).mean() for i in (i0, i1)])
        return float(((M[i0].mean(0) - M[i1].mean(0)) ** 2).sum() / within)
    sep = {key(g): separation(g) for g in G}
    g_state = max(G, key=separation)
    cand = {"batch": g_batch, "state": g_state}
    res = {"margin": MARGIN, "step_margin": STEP, "note": "exploratory; development donors; cells outside the selection third; a type is a grouping of the eight donors into two groups of four",
           "donors": donors, "batch_of_donor": {u: v[0] for u, v in b_of.items()}, "pairs": {k: len(v) for k, v in sets.items()},
           "candidates": {"batch": None if g_batch is None else key(g_batch), "state": key(g_state), "the_same_grouping": bool(g_batch is not None and key(g_batch) == key(g_state)),
                          "separation_of_the_donors_mean_states": {"batch": None if g_batch is None else sep[key(g_batch)], "state": sep[key(g_state)], "median_of_all_groupings": float(np.median(list(sep.values())))}}}
    if a.made_up:
        res["MADE_UP"] = "a trial of the script on made-up values; not a result"
    # --- part 1: the record alone. Cells are resampled within a donor, the same resamples for every protein
    rng = np.random.default_rng(SEED + 1)
    idx = {u: np.where(unit == u)[0] for u in donors}
    n = {u: len(idx[u]) for u in donors}
    S = {u: np.vstack([Y[idx[u]].sum(0)[None], rng.multinomial(n[u], np.ones(n[u]) / n[u], size=DRAWS) @ Y[idx[u]]]) for u in donors}     # sums of every protein, row 0 the data
    def u_type(name, s, ms):
        j = sel[name]["index"]
        others = [u for u in donors if u != s]
        mo = sum(S[u][:, j] for u in others) / sum(n[u] for u in others)
        return 1 - (sum(S[u][:, j] for u in ms) / sum(n[u] for u in ms) - mo) / (S[s][:, j] / n[s] - mo)
    P1 = {key(g): {p: u_type(p[0], p[1], mates(g, p[1])) for p in pairs} for g in G}
    med = lambda k, ps: np.median(np.abs(np.array([P1[k][p] for p in ps])), axis=0)       # the median over pairs, on every resample
    part1 = {"by_grouping": {k: {s_: float(med(k, ps)[0]) for s_, ps in sets.items()} for k in P1}, "candidates": {}}
    for c, g in cand.items():
        if g is None:
            continue
        k = key(g)
        e = {"grouping": k, "by_pair": {f"{p[0]}|{p[1]}": {"u": float(P1[k][p][0]), "ci": iv(P1[k][p][1:]), "verdict": verdict(iv(P1[k][p][1:])), "mates": mates(g, p[1]),
                                                           "own_row_needed": need[p], "a_step_to_remove": A3[p]} for p in pairs}, "sets": {}}
        for s_, ps in sets.items():
            m = med(k, ps)
            v_ = [e["by_pair"][f"{p[0]}|{p[1]}"]["verdict"] for p in ps]
            order = sorted(part1["by_grouping"], key=lambda q: part1["by_grouping"][q][s_])
            e["sets"][s_] = {"pairs": len(ps), "median_abs_u": float(m[0]), "ci": iv(m[1:]), "supported": v_.count("support"), "unresolved": v_.count("unresolved"), "rejected": v_.count("reject"),
                             "rank_among_the_35_groupings": order.index(k) + 1, "best_grouping": order[0], "median_abs_u_of_the_best": part1["by_grouping"][order[0]][s_],
                             "median_over_the_35_groupings": float(np.median([part1["by_grouping"][q][s_] for q in order]))}
        part1["candidates"][c] = e
    res["part1_the_record_alone"] = part1
    # --- part 2: with the repair
    Typed.plan = {}
    part2, worst, anchors = {"by_pair": {}, "candidates": {}}, 0.0, {}
    first = cycle.halves(unit, SEED)
    for name in sorted({p[0] for p in pairs}):
        j = sel[name]["index"]
        y, o = Y[:, j], O[:, j]
        Typed.plan = {s: {key(g): mates(g, s) for g in G} for s in sel[name]["admitted"]}
        out, _inside, _fits = cycle.cross_fit(x, y, o, np.ones(len(y), bool), unit, donors, first, Typed)
        for s in sel[name]["admitted"]:
            U = {arm: cycle.output_u(y, out[k_], unit == s, None, SEED, unit)[0] for arm, k_ in (("rows", "rows"), ("unseen", f"unseen|{s}"))}
            T = {key(g): cycle.output_u(y, out[f"typed|{s}|{key(g)}"], unit == s, None, SEED, unit)[0] for g in G}
            if not a.made_up:
                r = rec[name]["events"][s]["arms"]
                worst = max(worst, abs(float(U["rows"][0]) - r["rows"]["u"]), abs(float(U["unseen"][0]) - r["unseen"]["u"]))
            anchors[(name, s)] = (U, T)
            print(f"{name:12s} {s[-2:]}: rows {float(U['rows'][0]):+.3f} unseen {float(U['unseen'][0]):+.3f} | typed, batch {float(T[key(g_batch)][0]) if g_batch else float('nan'):+.3f} state {float(T[key(g_state)][0]):+.3f}", flush=True)
    res["checks"] = {"rows_and_unseen_repeat_the_recorded_round": worst, "limit": 0.005, "hold": bool(a.made_up or worst <= 0.005)}
    if not res["checks"]["hold"]:
        raise SystemExit(f"the arms rows and unseen do not repeat the recorded round ({worst:.4f}): stopped, nothing written")
    absu = lambda v_: np.abs(np.asarray(v_, np.float64))
    def gain(a_, b_):                                           # |a| - |b| on paired resamples: what b removes that a leaves
        d_ = absu(a_) - absu(b_)
        c = iv(d_[1:])
        return {"gain": float(d_[0]), "ci": c, "reading": step(c)}
    medT = lambda k, ps: np.median(np.array([absu(anchors[p][1][k]) for p in ps]), axis=0)
    medU = lambda arm, ps: np.median(np.array([absu(anchors[p][0][arm]) for p in ps]), axis=0)
    part2["by_grouping"] = {key(g): {s_: float(medT(key(g), ps)[0]) for s_, ps in sets.items()} for g in G}
    part2["anchors"] = {s_: {arm: {"median_abs_u": float(medU(arm, ps)[0]), "ci": iv(medU(arm, ps)[1:])} for arm in ("unseen", "rows")} for s_, ps in sets.items()}
    for c, g in cand.items():
        if g is None:
            continue
        k = key(g)
        e = {"grouping": k, "by_pair": {}, "sets": {}}
        for p in pairs:
            U, T = anchors[p]
            e["by_pair"][f"{p[0]}|{p[1]}"] = {"typed": {"u": float(T[k][0]), "ci": iv(T[k][1:])}, "unseen": float(U["unseen"][0]), "rows": float(U["rows"][0]),
                                             "gain_over_unseen": gain(U["unseen"], T[k]), "own_row_still_adds": gain(T[k], U["rows"])}
        for s_, ps in sets.items():
            m = medT(k, ps)
            g_ = [e["by_pair"][f"{p[0]}|{p[1]}"]["gain_over_unseen"]["reading"] for p in ps]
            r_ = [e["by_pair"][f"{p[0]}|{p[1]}"]["own_row_still_adds"]["reading"] for p in ps]
            dm = medU("unseen", ps) - m
            order = sorted(part2["by_grouping"], key=lambda q: part2["by_grouping"][q][s_])
            e["sets"][s_] = {"pairs": len(ps), "median_abs_u": float(m[0]), "ci": iv(m[1:]), "gain_of_the_median_over_unseen": {"gain": float(dm[0]), "ci": iv(dm[1:]), "reading": step(iv(dm[1:]))},
                             "pairs_that_gain_a_step_over_unseen": g_.count("a step"), "pairs_with_no_step": g_.count("no step"), "pairs_unresolved": g_.count("unresolved"),
                             "pairs_where_the_own_row_still_adds_a_step": r_.count("a step"), "pairs_where_it_adds_no_step": r_.count("no step"),
                             "rank_among_the_35_groupings": order.index(k) + 1, "best_grouping": order[0], "median_abs_u_of_the_best": part2["by_grouping"][order[0]][s_],
                             "median_over_the_35_groupings": float(np.median([part2["by_grouping"][q][s_] for q in order]))}
        part2["candidates"][c] = e
    res["part2_with_the_repair"] = part2
    # --- the statements
    def h1(c):
        q = part1["candidates"][c]["sets"]["own_row_needed"]
        return "confirmed" if q["supported"] * 2 >= q["pairs"] else ("refuted" if q["rejected"] * 2 > q["pairs"] else "unresolved")
    def h2(c):
        q = part2["candidates"][c]["sets"]["own_row_needed"]
        return "confirmed" if q["pairs_that_gain_a_step_over_unseen"] * 2 >= q["pairs"] else ("refuted" if q["pairs_with_no_step"] * 2 > q["pairs"] else "unresolved")
    res["statements"] = {c: {"H1_the_types_record_covers_the_member": h1(c), "H2_the_repair_with_the_types_record_gains_a_step": h2(c),
                             "not_a_random_grouping_by_the_record": bool(part1["candidates"][c]["sets"]["own_row_needed"]["rank_among_the_35_groupings"] <= 2),
                             "not_a_random_grouping_with_the_repair": bool(part2["candidates"][c]["sets"]["own_row_needed"]["rank_among_the_35_groupings"] <= 2)} for c in part1["candidates"]}
    a.results.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print(json.dumps(res["candidates"], indent=1))
    for c in part1["candidates"]:
        print(f"\n== {c}: {part1['candidates'][c]['grouping']}")
        for s_ in sets:
            q1, q2 = part1["candidates"][c]["sets"][s_], part2["candidates"][c]["sets"][s_]
            print(f"  {s_:18s} ({q1['pairs']:2d}) record alone: median |u| {q1['median_abs_u']:.3f} [{q1['ci'][0]:.3f}, {q1['ci'][1]:.3f}], supported {q1['supported']}, rejected {q1['rejected']}, rank {q1['rank_among_the_35_groupings']} of 35 (best {q1['median_abs_u_of_the_best']:.3f})")
            print(f"  {'':23s} with the repair: median |u| {q2['median_abs_u']:.3f} (unseen {part2['anchors'][s_]['unseen']['median_abs_u']:.3f}, rows {part2['anchors'][s_]['rows']['median_abs_u']:.3f}); a step over unseen in {q2['pairs_that_gain_a_step_over_unseen']}, no step in {q2['pairs_with_no_step']}; rank {q2['rank_among_the_35_groupings']} of 35")
    print("\nstatements:", json.dumps(res["statements"]))


if __name__ == "__main__":
    main()
