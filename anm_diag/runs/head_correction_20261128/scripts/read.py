#!/usr/bin/env python3
"""The judge, once: on the 93 drawn units, how much of each change every arm's output carries. `account.reading`,
unchanged. The arms never saw a drawn unit or a judge control cell.

  python read.py --cells CELLS.tsv --units units.json --rows rows.npz --states STATES... --packs PACK...
                 [--results DIR] [--rehearsal]

-> results/reading.json, results/table.tsv, results/outcome.json. Refuses to run if reading.json exists.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _hc as K  # noqa: E402
from anm_diag import account as AC  # noqa: E402

half = lambda n: int(hashlib.sha256(("fit52|" + str(n)).encode()).hexdigest(), 16) % 2 == 0


def within(y, u):
    """y with the mean of its unit taken out."""
    out = np.array(y, float)
    for g in set(u.tolist()):
        m = u == g
        out[m] -= out[m].mean()
    return out


def pooled(pairs, key, draws, seed):
    """The share the place carries, over pairs; the interval by resampling the drawn units."""
    use = [p for p in pairs if p["arms"][key]["share"] is not None]
    if not use:
        return None
    by = {}
    for p in use:
        s = by.setdefault(p["unit"], [0.0, 0.0])
        s[0] += np.sign(p["change"]) * p["arms"][key]["share"] * p["change"]
        s[1] += abs(p["change"])
    names = sorted(by)
    rng = np.random.default_rng(seed)
    val = lambda d: sum(by[g][0] for g in d) / max(sum(by[g][1] for g in d), 1e-300)
    boot = [val([names[i] for i in rng.integers(0, len(names), len(names))]) for _ in range(draws)]
    return {"share": float(val(names)), "ci": [float(x) for x in np.percentile(boot, [2.5, 97.5])], "pairs": len(use), "units": len(names)}


def difference(pairs, a, b, draws, seed):
    """Arm a minus arm b, paired by unit, the interval by resampling the drawn units."""
    use = [p for p in pairs if p["arms"][a]["share"] is not None and p["arms"][b]["share"] is not None]
    by = {}
    for p in use:
        s = by.setdefault(p["unit"], [0.0, 0.0, 0.0])
        s[0] += np.sign(p["change"]) * p["arms"][a]["share"] * p["change"]
        s[1] += np.sign(p["change"]) * p["arms"][b]["share"] * p["change"]
        s[2] += abs(p["change"])
    names = sorted(by)
    rng = np.random.default_rng(seed)
    val = lambda d: (sum(by[g][0] for g in d) - sum(by[g][1] for g in d)) / max(sum(by[g][2] for g in d), 1e-300)
    boot = [val([names[i] for i in rng.integers(0, len(names), len(names))]) for _ in range(draws)]
    return {"difference": float(val(names)), "ci": [float(x) for x in np.percentile(boot, [2.5, 97.5])], "pairs": len(use), "units": len(names)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--units", type=Path, required=True)
    ap.add_argument("--rows", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--results", type=Path, default=K.RUN / "results")
    ap.add_argument("--rehearsal", action="store_true")
    ap.add_argument("--draws", type=int, default=K.RESAMPLES)
    a = ap.parse_args()
    assert not a.rehearsal or a.results.resolve() != (K.RUN / "results").resolve(), "a rehearsal does not write into the results of the run"
    out_file = a.results / "reading.json"
    if out_file.exists():
        raise SystemExit("reading.json exists: the judge was read once and is not read again")
    R = dict(np.load(a.rows))
    C = K.cells(a.cells)
    drawn = K.drawn_units(C)
    assert sum(len(v) for v in drawn.values()) == json.loads(a.units.read_text())["drawn"]["n"]
    St, Pk = K.B3.States(a.states), K.B3.Packs(a.packs)
    res = {"rehearsal_on_made_up_values": bool(a.rehearsal), "rows_sha256": K.sha(a.rows), "min_errors": AC.MIN_SNR, "arms": list(K.ARMS), "conditions": {}, "pairs": [], "quiet": [], "control": {}}
    for c in K.CONDITIONS:
        g = C[c]
        keep = g["drawn"] | (g["control"] & ~g["teach"])                       # the drawn units and the judge control cells only
        ids = g["ids"][keep]
        Z = K.B3.V.unit_length(St.rows(ids))
        A, _rows0 = K.head(St.rows(ids), ids)
        clr, _log, _bg = Pk.readouts(K.STUDY, ids, list(K.READOUTS))
        unit = np.where(g["control"][keep], "control", g["unit"][keep])
        first = np.array([half(i) for i in ids])
        o = {arm: {r: K.out_of(A, R[f"{arm}__{r}__w"], float(R[f"{arm}__{r}__b"])) for r in K.READOUTS} for arm in K.ARMS}
        res["conditions"][c] = {"drawn_units": len(drawn[c]), "cells_of_the_drawn_units": int(g["drawn"][keep].sum()), "judge_control_cells": int((unit == "control").sum())}
        ctl = unit == "control"
        res["control"][c] = {arm: {r: float(np.std(within(clr[r][ctl], np.array(["control"] * int(ctl.sum()))) - np.polyfit(within(o[arm][r][ctl], np.array(["c"] * int(ctl.sum()))), within(clr[r][ctl], np.array(["c"] * int(ctl.sum()))), 1)[0] * within(o[arm][r][ctl], np.array(["c"] * int(ctl.sum()))))) for r in K.READOUTS} for arm in K.ARMS}
        for u in drawn[c]:
            for r in K.READOUTS:
                m = ctl | (unit == u)
                base = AC.reading(Z[m], clr[r][m], unit[m], first[m], u, o=o["H0"][r][m], seed=K.SEED)
                row = {"unit": f"{u}@{c}", "readout": r, "change": base["change"], "errors": base["snr"], "cells": base["cells"], "carried": base["shares"]["carried"]["share"], "arms": {}}
                if base["snr"] >= AC.MIN_SNR:
                    for arm in K.ARMS:
                        s = base if arm == "H0" else AC.reading(Z[m], clr[r][m], unit[m], first[m], u, o=o[arm][r][m], seed=K.SEED)
                        row["arms"][arm] = {"share": s["shares"]["as_trained"]["share"], "ci": s["shares"]["as_trained"]["ci"]}
                        if arm != "H0":
                            row.setdefault("carried_every_arm", []).append(s["shares"]["carried"]["share"])
                    res["pairs"].append(row)
                    print(f"{c:11s} {u:9s} {r:5s} change {row['change']:+.3f} ({row['errors']:4.1f}) " + " ".join(f"{k} {v['share']:+.2f}" for k, v in row["arms"].items()), flush=True)
                else:
                    sp = {arm: float((o[arm][r][unit == u].mean() - o[arm][r][ctl].mean()) / (o[arm][r][ctl].std() + 1e-12)) for arm in K.ARMS}
                    res["quiet"].append({"unit": f"{u}@{c}", "readout": r, "change": base["change"], "errors": base["snr"], "the_output_moves_in_spreads": sp})
    P = res["pairs"]
    res["pooled"] = {arm: pooled(P, arm, a.draws, K.SEED) for arm in K.ARMS}
    res["differences"] = {f"{x}-{y}": difference(P, x, y, a.draws, K.SEED) for x, y in (("H2", "H1"), ("H1", "H0"), ("H2", "H0"), ("H3", "H1"))}
    best = max(("H1", "H2"), key=lambda k: res["pooled"][k]["share"])
    c1 = {"says": "a correction at the head closes the loss, on units never read", "best_arm": best, "best": res["pooled"][best]["share"], "H0": res["pooled"]["H0"]["share"]}
    c1["holds"] = bool(c1["best"] >= 0.60 and c1["H0"] <= 0.25)
    c1["refuted"] = bool(c1["best"] < 0.40 or c1["H0"] > 0.40)
    d = res["differences"]["H2-H1"]
    c2 = {"says": "the loop's rows teach something a plain refit does not", **d}
    c2["holds"] = bool(d["ci"][0] > 0)
    c2["refuted"] = bool(d["ci"][1] < 0)
    g1 = {"says": "the state did not move", "largest_difference": max((abs(x - p["carried"]) for p in P for x in p.get("carried_every_arm", [])), default=0.0)}
    g1["holds"] = bool(g1["largest_difference"] <= 0.02)
    worst = max((max(abs(v) for v in q["the_output_moves_in_spreads"].values()) for q in res["quiet"]), default=0.0)
    g2 = {"says": "nothing got worse where nothing was wrong", "largest_move_of_an_output_on_a_quiet_pair_in_spreads": worst, "quiet_pairs": len(res["quiet"]),
          "error_on_the_judge_control_cells": {arm: float(np.mean([res["control"][c][arm][r] for c in K.CONDITIONS for r in K.READOUTS])) for arm in K.ARMS}}
    g2["holds"] = bool(worst <= 1.5 and all(g2["error_on_the_judge_control_cells"][arm] <= g2["error_on_the_judge_control_cells"]["H0"] + 1e-9 for arm in K.ARMS))
    res["lines"] = {"C1": c1, "C2": c2, "G1": g1, "G2": g2}
    res["passed"] = bool(c1["holds"])
    a.results.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=1, ensure_ascii=False, default=float) + "\n")
    head = ["unit", "readout", "change", "errors"] + list(K.ARMS) + ["the state carries", "cells of the unit", "judge control cells"]
    lines = ["\t".join(head)] + ["\t".join([p["unit"], p["readout"], f"{p['change']:+.3f}", f"{p['errors']:.1f}"] + [f"{p['arms'][k]['share']:+.2f}" for k in K.ARMS] + [f"{p['carried']:.2f}", str(p["cells"]["event"]), str(p["cells"]["control"])]) for p in sorted(P, key=lambda p: -abs(p["change"]))]
    (a.results / "table.tsv").write_text("\n".join(lines) + "\n")
    (a.results / "outcome.json").write_text(json.dumps({"rehearsal_on_made_up_values": bool(a.rehearsal), "passed": res["passed"], "lines": res["lines"], "pooled": res["pooled"], "differences": res["differences"],
                                                        "pairs": len(P), "quiet_pairs": len(res["quiet"]), "units": res["conditions"]}, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({"pooled": {k: (None if v is None else round(v["share"], 3)) for k, v in res["pooled"].items()}, "H2-H1": [round(x, 3) for x in d["ci"]],
                      "lines": {k: ("holds" if v.get("holds") else ("refuted" if v.get("refuted") else "unresolved")) for k, v in res["lines"].items()}, "passed": res["passed"], "pairs": len(P)}, indent=1))


if __name__ == "__main__":
    main()
