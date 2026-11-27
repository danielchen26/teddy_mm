#!/usr/bin/env python3
"""The hypotheses of DESIGN.md, Section 4, against the reading. results/outcome.json. Nothing here is chosen by hand.

  python score.py [--results DIR]

The sets the hypotheses are scored on were written into forecasts.json before the final cells were exported. What
the design left open is fixed here, before the reading:

- F2. The units of the set are resampled 2,000 times with seed 20261127. For a pair of systems the difference of the
  summed squared errors is taken over the resampled units; its 95% interval is between the 2.5th and the 97.5th
  percentile. Holds: the interval of S0 minus S3 lies above zero, and the interval of S2 minus S3 does not lie below
  zero. Refuted: the error of S3 on the set is not below that of S0, or the interval of S2 minus S3 lies below zero.
- F3 and F4 have no middle. F3 is not tested with fewer than four pairs.
- F5 counts by test (a proposal asked of one stratum's row) and by proposal (a pair is admitted if a test admits it).
  It holds if at least half are admitted in both counts, does not hold if fewer than half in both, and is unresolved
  if the counts disagree; not tested with fewer than four tests.
- The audit: passed if F1 and F2 hold; not passed if either is refuted; else unresolved.
- `explored_against_the_background` scores the melanoma forecasts against the change measured with the isotype
  antibodies as the reference. Declared here; it enters no hypothesis.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fa as F  # noqa: E402

EXPECTED = {"F1_holds": 0.8, "F1_unresolved": 0.12, "F1_refuted": 0.08, "F2_holds": 0.7, "the_audit_is_passed": 0.6, "F3_holds": 0.5, "F4_holds": 0.65, "F5_holds": 0.65}
INSIDE, INSIDE_REFUTED, MEDIAN, MEDIAN_REFUTED = 0.90, 0.85, (0.23, 0.91), (0.15, 1.4)
DRAWS, SEED = 2000, 20261127
rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None
key = lambda r: f"{r['stratum']}|{r['event']}|{r['readout']}"


def errors(rows):
    return {"pairs": len(rows), "no_change": rms([r["change"] for r in rows]), **{s: rms([r[s]["error"] for r in rows]) for s in ("S0", "S2", "S3")}}


def intervals(rows, s="S3"):
    if not rows:
        return {"pairs": 0}
    z2 = [r[s]["z"] ** 2 for r in rows]
    n = sum(r[s]["inside"] for r in rows)
    return {"pairs": len(rows), "inside": n, "share_inside": n / len(rows), "mean_z2": float(np.mean(z2)), "median_z2": float(np.median(z2)),
            "median_half_width": float(np.median([(r[s]["interval"][1] - r[s]["interval"][0]) / 2 for r in rows]))}


def verdict(c):
    if c["share_inside"] >= INSIDE and MEDIAN[0] <= c["median_z2"] <= MEDIAN[1]:
        return "holds"
    if c["share_inside"] < INSIDE_REFUTED or not (MEDIAN_REFUTED[0] <= c["median_z2"] <= MEDIAN_REFUTED[1]):
        return "refuted"
    return "unresolved"


def both(x):
    return {"errors": errors(x), "S2": intervals(x, "S2"), "S3": intervals(x, "S3")}


def brief(r):
    return {"stratum": r["stratum"], "event": r["event"], "readout": r["readout"], "edge": r["edge"], "type": r["type"], "cells": r["cells"], "change": r["change"], "se": r["se"], "snr": r["snr"],
            "S0": r["S0"]["forecast"], "S2": r["S2"]["forecast"], "S3": r["S3"]["forecast"], "interval_S3": r["S3"]["interval"], "inside_S2": r["S2"]["inside"], "inside_S3": r["S3"]["inside"], "z_S3": r["S3"]["z"]}


def resampled(rows, a, b):
    """The difference of the summed squared errors of systems a and b over resampled units: the value and its interval."""
    units = sorted({(r["stratum"], r["event"]) for r in rows})
    d = np.array([sum(r[a]["error"] ** 2 - r[b]["error"] ** 2 for r in rows if (r["stratum"], r["event"]) == u) for u in units])
    rng = np.random.default_rng(SEED)
    draws = d[rng.integers(0, len(d), (DRAWS, len(d)))].sum(1)
    return {"units": len(units), "difference": float(d.sum()), "interval": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))], "units_where_the_second_is_closer": int((d > 0).sum()), "units_where_it_is_further": int((d < 0).sum())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=F.RUN / "results")
    a = ap.parse_args()
    fc = json.loads((F.RUN / "results" / "forecasts.json").read_text())
    R = json.loads((a.results / "reading.json").read_text())
    rows = R["rows"]
    by = {key(r): r for r in rows}
    assert len(by) == len(rows) == fc["pairs"]
    sets = {k: [by[i] for i in v] for k, v in fc["sets"].items()}
    D = sets["the_systems_can_differ"]
    other = [r for r in rows if key(r) not in set(fc["sets"]["the_systems_can_differ"])]
    out = {"pairs": len(rows), "rehearsal_on_made_up_values": R["rehearsal_on_made_up_values"], **both(rows),
           "by_stratum": {g: both([r for r in rows if r["stratum"] == g]) for g in F.STRATA},
           "pairs_where_the_systems_can_differ": both(D), "all_other_pairs": both(other),
           "pairs_with_an_edge": {**both(sets["with_an_edge"]), "rows": [brief(r) for r in sets["with_an_edge"]]},
           "pairs_with_a_change": {**both(x), "rows": [brief(r) for r in sorted(x, key=lambda r: -r["snr"])]} if (x := [r for r in rows if r["snr"] >= R["min_snr"]]) else {"errors": errors([])},
           "outside_S3": [brief(r) for r in sorted(rows, key=lambda r: -abs(r["S3"]["z"])) if not r["S3"]["inside"]], "outside_S3_by_unit": {}, "outside_S3_by_readout": {}, "hypotheses": {}}
    for r in rows:
        if not r["S3"]["inside"]:
            u, q = f"{r['stratum']}|{r['event']}", f"{r['stratum']}|{r['readout']}"
            out["outside_S3_by_unit"][u] = out["outside_S3_by_unit"].get(u, 0) + 1
            out["outside_S3_by_readout"][q] = out["outside_S3_by_readout"].get(q, 0) + 1
    H = out["hypotheses"]
    # --- F1
    c = out["S3"]
    H["F1_the_intervals_hold"] = dict(c, asked={"share_inside_at_least": INSIDE, "median_z2_between": MEDIAN, "refuted_below": INSIDE_REFUTED, "or_median_z2_outside": MEDIAN_REFUTED}, outcome=verdict(c),
                                      S2_for_comparison=dict(out["S2"], as_if_scored=verdict(out["S2"])), by_stratum={g: dict(v["S3"], as_if_scored=verdict(v["S3"])) for g, v in out["by_stratum"].items() if v["S3"]["pairs"]})
    # --- F2
    e = errors(D)
    over_s0, over_s2 = resampled(D, "S0", "S3"), resampled(D, "S2", "S3")
    holds = over_s0["interval"][0] > 0 and not over_s2["interval"][1] < 0
    refuted = not e["S3"] < e["S0"] or over_s2["interval"][1] < 0
    H["F2_the_loop_corrects"] = {"errors": e, "S0_minus_S3": over_s0, "S2_minus_S3": over_s2, "on_all_pairs": out["errors"], "outcome": "holds" if holds else "refuted" if refuted else "unresolved"}
    # --- F3
    x = sets["F3"]
    e3, c3 = errors(x), intervals(x)
    H["F3_an_edge_has_the_size_of_its_pair"] = {"errors": e3, "intervals": c3, "rows": [brief(r) for r in x],
                                                "outcome": "not tested: fewer than four such pairs" if len(x) < 4 else "holds" if e3["S3"] < e3["S0"] and c3["share_inside"] >= 0.75 else "does not hold"}
    # --- F4
    x = sets["F4"]
    H["F4_the_arithmetic_of_the_readout"] = {"errors": errors(x), "S3": intervals(x, "S3"), "S2": intervals(x, "S2"), "by_unit": {u: both([r for r in x if f"{r['stratum']}|{r['event']}" == u])["S3"] for u in sorted({f"{r['stratum']}|{r['event']}" for r in x})},
                                             "outcome": "not tested: no such pair" if not x else "holds" if intervals(x)["share_inside"] >= 0.90 else "does not hold"}
    # --- F5
    t = R["edges"]["tested"]
    pairs = sorted({(x["site"], x["readout"]) for x in t})
    adm_t = sum(x["status"] == "admitted" for x in t)
    adm_p = sum(any(x["status"] == "admitted" for x in t if (x["site"], x["readout"]) == p) for p in pairs)
    half_t, half_p = (adm_t >= 0.5 * len(t)) if t else None, (adm_p >= 0.5 * len(pairs)) if pairs else None
    H["F5_an_edge_proposed_elsewhere_holds_here"] = {"tests": len(t), "admitted_tests": adm_t, "proposals_with_a_row": len(pairs), "admitted_proposals": adm_p, "with_no_row": [x["edge_id"] for x in R["edges"]["with_no_row"]],
                                                     "rows": [{"edge": x["edge_id"], "status": x["status"], "proposed_with_z": x["evidence"]["z"], **x["test"]} for x in t],
                                                     "outcome": "not tested: fewer than four tests" if len(t) < 4 else "holds" if half_t and half_p else "does not hold" if not half_t and not half_p else "unresolved"}
    f1, f2 = H["F1_the_intervals_hold"]["outcome"], H["F2_the_loop_corrects"]["outcome"]
    out["the_audit"] = {"F1": f1, "F2": f2, "verdict": "passed" if f1 == f2 == "holds" else "not passed" if "refuted" in (f1, f2) else "unresolved"}
    out["expected_before_any_cell_was_listed"] = EXPECTED
    # --- reported apart: THP-1, its two own-gene pairs, and the melanoma pairs against the background
    thp_own = [r for r in rows if r["stratum"] == "thp1" and (r["event"], r["readout"]) in (("CD86", "CD86"), ("PDCD1LG2", "PDL2"))]
    out["thp1_pairs_of_a_unit_that_removes_its_own_protein"] = [brief(r) for r in thp_own]

    def bg(rs):
        rs = [r for r in rs if "against_the_background" in r]
        if not rs:
            return None
        d = lambda r, s: r["against_the_background"]["change"] - r[s]["forecast"]
        wide = lambda r: np.sqrt(max(r["S3"]["sd"] ** 2 + r["against_the_background"]["se"] ** 2 - r["se"] ** 2, 1e-12)) * (r["S3"]["interval"][1] - r["S3"]["interval"][0]) / 2 / r["S3"]["sd"]
        return {"pairs": len(rs), "error": {"no_change": rms([r["against_the_background"]["change"] for r in rs]), **{s: rms([d(r, s) for r in rs]) for s in ("S0", "S2", "S3")}},
                "S3_inside_with_the_error_of_this_measure": int(sum(abs(d(r, "S3")) <= wide(r) for r in rs)), "S3_inside_as_scored": int(sum(r["S3"]["inside"] for r in rs))}
    out["explored_against_the_background"] = {"all_melanoma_pairs": bg(rows), "F4_pairs": bg(sets["F4"])}
    acc_file = a.results / "accounts.json"
    if acc_file.exists():
        A = json.loads(acc_file.read_text())
        says, largest = {}, {}
        for v in A["accounts"].values():
            says[v["chain"]["says"]] = says.get(v["chain"]["says"], 0) + 1
            part = v["chain"]["says"].split(" (")[0].split(";")[0].rsplit(" ", 1)[0]
            largest[part] = largest.get(part, 0) + 1
        states = {}
        for it in A["open_parts"]:
            states[f"{it.get('state')} | {it.get('lever')}"] = states.get(f"{it.get('state')} | {it.get('lever')}", 0) + 1
        out["accounts"] = {"pairs_with_a_change": A["pairs_with_a_change"], "accounts": len(A["accounts"]), "by_the_first_part_named": largest,
                           "rows": [{"pair": k, "change": v["change_in_spreads"], "snr": v["snr"], "edge": v["edge"], "carried": {s: v["stages"][s]["shares"]["carried"]["share"] for s in v["stages"]}, "says": v["chain"]["says"]} for k, v in A["accounts"].items()],
                           "open_parts": {"count": len(A["open_parts"]), "by_state_and_lever": states}}
    (a.results / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    e, c2 = out["errors"], out["S2"]
    print(f"pairs {out['pairs']} | error: no change {e['no_change']:.3f}, S0 {e['S0']:.3f}, S2 {e['S2']:.3f}, S3 {e['S3']:.3f}")
    print(f"  inside S3 {c['inside']} ({c['share_inside']:.3f}), median z2 {c['median_z2']:.2f}, mean z2 {c['mean_z2']:.2f} | S2 {c2['inside']} ({c2['share_inside']:.3f}), median z2 {c2['median_z2']:.2f}")
    for g, v in H["F1_the_intervals_hold"]["by_stratum"].items():
        print(f"    {g}: {v['inside']}/{v['pairs']} ({v['share_inside']:.3f}), median z2 {v['median_z2']:.2f}")
    h = H["F2_the_loop_corrects"]
    print(f"  on the {h['errors']['pairs']} pairs where the systems can differ: error S0 {h['errors']['S0']:.3f}, S2 {h['errors']['S2']:.3f}, S3 {h['errors']['S3']:.3f} | S0-S3 {h['S0_minus_S3']['difference']:.2f} {np.round(h['S0_minus_S3']['interval'], 2)} | S2-S3 {h['S2_minus_S3']['difference']:.3f} {np.round(h['S2_minus_S3']['interval'], 3)}")
    for k, v in H.items():
        print(f"  {k}: {v['outcome']}")
    print("  THE AUDIT:", out["the_audit"])


if __name__ == "__main__":
    main()
