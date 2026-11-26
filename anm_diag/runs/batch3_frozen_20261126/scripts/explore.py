#!/usr/bin/env python3
"""Explored after the reading: why the intervals of version 5 did not hold. results/explored.json.

  python explore.py

Everything here was looked at after stage B was read. It is no test and changes no verdict of score.py. It uses
reading.json and forecasts.json alone.

1. The arithmetic of the readout. A readout of the melanoma screen is the log count of one antibody minus the mean over
   the 24 columns of the cell. If a unit changes the log count of the columns in A and of no other column, every readout
   outside A moves by

       - (sum of the changes of the readouts in A, in log units) / (24 - |A|)

   with no constant to fit. A: the readouts the unit has an edge to (declared, or admitted in run 43). Two variants:
   `at forecast time` uses what was known before stage B (the recorded change of a pair read before, else the forecast
   of S2 for the pair); `measured` uses the measured changes of the readouts in A.
2. Where the misses of S2 are: by unit, before and after that term; the pairs still far outside.
3. The spread of a pair with no edge and no type that the rows of this batch give by rule 5a (`account._calibrated_spread`,
   the median of z squared), per stratum, against the spread S2 used; and how many pairs lie inside with it. This is
   in-sample: the same rows set the spread and are counted.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _b3 as B  # noqa: E402
from anm_diag import account as AC  # noqa: E402

COLUMNS = 24
MEL = ("mel_co", "mel_ctl")
rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None


def stats(err, sd):
    z2 = (err / sd) ** 2
    return {"pairs": int(len(err)), "error": rms(err), "mean_error": float(err.mean()), "inside": int((np.abs(err) <= 1.96 * sd).sum()), "share_inside": float((np.abs(err) <= 1.96 * sd).mean()),
            "mean_z2": float(z2.mean()), "median_z2": float(np.median(z2)), "beyond_3": int((z2 > 9).sum()), "beyond_5": int((z2 > 25).sum())}


def main():
    res_dir = B.RUN / "results"
    R, fc = json.loads((res_dir / "reading.json").read_text()), json.loads((res_dir / "forecasts.json").read_text())
    r43 = json.loads((B.RUN.parent / "one_loop_20261125" / "results" / "retro.json").read_text())["edges"]
    rows = R["rows"]
    by = {(r["stratum"], r["event"], r["readout"]): r for r in rows}
    edges = {g: {} for g in MEL}
    for e in fc["edges"]["declared"]:
        for g in MEL:
            edges[g].setdefault(e["site"], set()).add(e["readout"])
    for e in r43:
        if e["status"] == "admitted":
            edges[e["stratum"]].setdefault(e["site"], set()).add(e["readout"])

    def term(r, how):
        g, e, j = r["stratum"], r["event"], r["readout"]
        if g not in MEL:
            return 0.0
        A = sorted(i for i in edges[g].get(e, ()) if i != j)
        if not A:
            return 0.0
        sd, tot = R["spread"][g], 0.0
        for i in A:
            f = by.get((g, e, i))
            known = R["records"][g][e]["changes"][i][0] if (how == "measured" or f is None) else f["S2"]["forecast"]
            tot += known * sd[i]
        return -tot / ((COLUMNS - len(A)) * sd[j])

    out = {"note": "explored after the reading; no test", "the_arithmetic_of_the_readout": {}}
    no_edge = [r for r in rows if r["edge"] is None]
    err, sd = np.array([r["S2"]["error"] for r in no_edge]), np.array([r["S2"]["sd"] for r in no_edge])
    mel = np.array([r["stratum"] in MEL for r in no_edge])
    for how in ("at forecast time", "measured"):
        c = np.array([term(r, how) for r in no_edge])
        has = c != 0
        o = {"pairs_with_no_edge_in_the_melanoma_screen": {"as_scored": stats(err[mel], sd[mel]), "with_the_term": stats((err - c)[mel], sd[mel])},
             "of_units_with_an_edge_to_a_readout_of_the_panel": {"as_scored": stats(err[has], sd[has]), "with_the_term": stats((err - c)[has], sd[has])},
             "of_the_other_units": stats(err[mel & ~has], sd[mel & ~has]),
             "all_pairs_with_no_edge": {"as_scored": stats(err, sd), "with_the_term": stats(err - c, sd)}, "units": []}
        for g in MEL:
            for u in sorted(edges[g]):
                m = np.array([r["stratum"] == g and r["event"] == u for r in no_edge])
                if m.sum():
                    o["units"].append({"stratum": g, "unit": u, "edges_to": sorted(edges[g][u]), "pairs": int(m.sum()), "mean_error_of_S2": float(err[m].mean()), "the_term": float(c[m].mean()),
                                       "outside": int((np.abs(err[m]) > 1.96 * sd[m]).sum()), "outside_with_the_term": int((np.abs((err - c)[m]) > 1.96 * sd[m]).sum())})
        out["the_arithmetic_of_the_readout"][how] = o
        if how == "at forecast time":
            c_fc = c
    # every pair, as H1 counts them, with the term known at forecast time
    c_all = {id(r): v for r, v in zip(no_edge, c_fc)}
    e_all = np.array([r["S2"]["error"] - c_all.get(id(r), 0.0) for r in rows])
    s_all = np.array([r["S2"]["sd"] for r in rows])
    out["all_698_pairs"] = {"S1_as_scored": stats(np.array([r["S1"]["error"] for r in rows]), np.array([r["S1"]["sd"] for r in rows])),
                            "S2_as_scored": stats(np.array([r["S2"]["error"] for r in rows]), s_all), "S2_with_the_term_at_forecast_time": stats(e_all, s_all)}
    # what is still far outside with the term: edges that are not on the graph
    tested = {(e["site"], e["readout"], e["stratum"]): e["status"] for e in R["edges"]["tested"]}
    far = sorted(((abs(e / s), r, e) for r, e, s in zip(rows, e_all, s_all) if abs(e / s) > 3), key=lambda t: -t[0])
    out["still_beyond_three_with_the_term"] = [{"stratum": r["stratum"], "event": r["event"], "readout": r["readout"], "edge": r["edge"], "type": r["type"], "change": r["change"], "se": r["se"], "forecast_S2": r["S2"]["forecast"],
                                                "miss_with_the_term": float(e), "z": float(np.sign(e) * z), "proposed_before_the_batch": tested.get((r["event"], r["readout"], r["stratum"]))} for z, r, e in far]
    # the spread with no edge and no type that the rows of this batch give, by rule 5a
    out["spread_with_no_edge_from_the_rows_of_the_batch"] = {}
    plain = [(r, e) for r, e in zip(rows, e_all) if r["edge"] is None and r["type"] is None]
    parts = {(q["stratum"], q["event"], q["readout"]): q["S2"]["parts_of_the_spread"] for q in fc["rows"]}
    for g in B.STRATA:
        for grp in ("read before", "never read"):
            x = [(r, e) for r, e in plain if r["stratum"] == g and r["group"] == grp]
            if len(x) < AC.MIN_ROWS:
                continue
            e = np.array([v for _r, v in x])
            f = [parts[(r["stratum"], r["event"], r["readout"])] for r, _v in x]
            base = np.array([q["readout"] ** 2 + q["row"] ** 2 for q in f])
            used = float(np.median([q["direct"] for q in f]))
            s = float(AC._calibrated_spread(e ** 2, base))
            sd_new = np.sqrt(base + s ** 2)
            out["spread_with_no_edge_from_the_rows_of_the_batch"][f"{g} | {grp}"] = {"pairs": len(x), "S2_used": used, "these_rows_give": s, "with_the_spread_S2_used": stats(e, np.sqrt(base + used ** 2)), "with_the_spread_these_rows_give": stats(e, sd_new)}
    (res_dir / "explored.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for how, o in out["the_arithmetic_of_the_readout"].items():
        a, b = o["pairs_with_no_edge_in_the_melanoma_screen"]["as_scored"], o["pairs_with_no_edge_in_the_melanoma_screen"]["with_the_term"]
        print(f"{how:17s} melanoma pairs with no edge {a['pairs']}: inside {a['inside']} -> {b['inside']} | error {a['error']:.3f} -> {b['error']:.3f} | mean z2 {a['mean_z2']:.2f} -> {b['mean_z2']:.2f} | median {a['median_z2']:.2f} -> {b['median_z2']:.2f}")
    for k, v in out["all_698_pairs"].items():
        print(f"  {k:34s} inside {v['inside']} ({v['share_inside']:.3f}) mean z2 {v['mean_z2']:.2f} median {v['median_z2']:.2f} beyond 3: {v['beyond_3']} beyond 5: {v['beyond_5']}")
    print("still beyond three with the term:", len(out["still_beyond_three_with_the_term"]), [(r["event"], r["readout"], r["stratum"], round(r["z"], 1), r["proposed_before_the_batch"]) for r in out["still_beyond_three_with_the_term"][:14]])
    for k, v in out["spread_with_no_edge_from_the_rows_of_the_batch"].items():
        a, b = v["with_the_spread_S2_used"], v["with_the_spread_these_rows_give"]
        print(f"  {k:22s} n={v['pairs']:3d} spread S2 used {v['S2_used']:.3f} | these rows give {v['these_rows_give']:.3f} | inside {a['inside']} -> {b['inside']} ({b['share_inside']:.3f}) | mean z2 {a['mean_z2']:.2f} -> {b['mean_z2']:.2f}")


if __name__ == "__main__":
    main()
