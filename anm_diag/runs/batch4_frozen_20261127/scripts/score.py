#!/usr/bin/env python3
"""The gate of DESIGN.md, Section 4, against the reading. results/outcome.json. Nothing here is chosen by hand.

  python score.py [--results DIR]

What the design left open is fixed here, before the reading:

- H1 is scored on the 600 pairs with the interval of S3 as committed in forecasts.json. The 480 pairs of the units
  chosen by lot, each stratum, and the intervals of S2 are reported beside it and enter no verdict.
- "The errors of the three systems agree to two digits": the root mean square error of S0, S2 and S3 over the 600 pairs,
  each rounded to two digits, are the same number.
- The gate by the numbers: H1 holds, and the three parts of the turn by rule exist (reading.json, by_rule.json,
  accounts.json). Whether a fault of the method changes a verdict is said in the README.
- The readout is a centred log-ratio. `explored_against_the_background` scores the same forecasts against the change
  measured with the isotype antibodies as the reference. It is declared here, enters no hypothesis and changes no
  verdict.

By rule the batch leaves: its rows as records; the edges its rows admitted or refused; the edges its flags proposed,
with what the library's rows said of them; the open parts of its accounts, each with where it goes; and the pairs now
measured twice, which give the spread of a carried share.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _b4 as B  # noqa: E402
from anm_diag import events as EV  # noqa: E402

EXPECTED = {"H1_holds": 0.7, "H1_unresolved": 0.2, "H1_refuted": 0.1, "the_errors_of_the_three_systems_agree_to_two_digits": 0.9, "at_least_one_of_the_three_edge_tests_admits": 0.5}
INSIDE, INSIDE_REFUTED = 0.90, 0.85
MEDIAN, MEDIAN_REFUTED = (0.23, 0.91), (0.15, 1.4)
rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None
key = lambda r: f"{r['stratum']}|{r['event']}|{r['readout']}"


def errors(rows):
    return {"pairs": len(rows), "no_change": rms([r["change"] for r in rows]), **{s: rms([r[s]["error"] for r in rows]) for s in ("S0", "S2", "S3")}}


def intervals(rows, s="S3"):
    if not rows:
        return {"pairs": 0}
    z2 = [r[s]["z"] ** 2 for r in rows]
    n = sum(r[s]["inside"] for r in rows)
    hw = [(r[s]["interval"][1] - r[s]["interval"][0]) / 2 for r in rows]
    return {"pairs": len(rows), "inside": n, "share_inside": n / len(rows), "mean_z2": float(np.mean(z2)), "median_z2": float(np.median(z2)), "median_half_width": float(np.median(hw))}


def verdict(c):
    lo, hi = MEDIAN
    if c["share_inside"] >= INSIDE and lo <= c["median_z2"] <= hi:
        return "holds"
    if c["share_inside"] < INSIDE_REFUTED or not (MEDIAN_REFUTED[0] <= c["median_z2"] <= MEDIAN_REFUTED[1]):
        return "refuted"
    return "unresolved"


def brief(r):
    return {"stratum": r["stratum"], "event": r["event"], "readout": r["readout"], "chosen": r["chosen"], "edge": r["edge"], "cells": r["cells"], "change": r["change"], "se": r["se"], "snr": r["snr"],
            "S0": r["S0"]["forecast"], "S2": r["S2"]["forecast"], "S3": r["S3"]["forecast"], "interval_S3": r["S3"]["interval"], "inside_S2": r["S2"]["inside"], "inside_S3": r["S3"]["inside"], "z_S3": r["S3"]["z"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=B.RUN / "results")
    a = ap.parse_args()
    fc = json.loads((B.RUN / "results" / "forecasts.json").read_text())
    R = json.loads((a.results / "reading.json").read_text())
    rows = R["rows"]
    by = {key(r): r for r in rows}
    assert len(by) == len(rows) == fc["pairs"]
    sets = {k: [by[i] for i in v] for k, v in fc["sets"].items()}
    both = lambda x: {"errors": errors(x), "S2": intervals(x, "S2"), "S3": intervals(x, "S3")}
    out = {"pairs": len(rows), "rehearsal_on_made_up_values": R["rehearsal_on_made_up_values"], "errors": errors(rows), "intervals": {s: intervals(rows, s) for s in ("S2", "S3")},
           "by_how_the_unit_was_chosen": {k: both(sets[k]) for k in ("by lot", "by an item on record")},
           "by_stratum": {g: both([r for r in rows if r["stratum"] == g]) for g in B.STRATA},
           "pairs_with_an_edge": {**both(sets["with an edge"]), "rows": [brief(r) for r in sets["with an edge"]]},
           "pairs_with_a_change": {**both(x), "rows": [brief(r) for r in sorted(x, key=lambda r: -r["snr"])]} if (x := [r for r in rows if r["snr"] >= R["min_snr"]]) else {"errors": errors([])},
           "outside_S3": [brief(r) for r in sorted(rows, key=lambda r: -abs(r["S3"]["z"])) if not r["S3"]["inside"]], "outside_S3_by_unit": {}, "outside_S3_by_readout": {}, "hypotheses": {}}
    for r in rows:
        if not r["S3"]["inside"]:
            u, q = f"{r['stratum']}|{r['event']}", f"{r['stratum']}|{r['readout']}"
            out["outside_S3_by_unit"][u] = out["outside_S3_by_unit"].get(u, 0) + 1
            out["outside_S3_by_readout"][q] = out["outside_S3_by_readout"].get(q, 0) + 1
    # --- H1, the gate
    c = out["intervals"]["S3"]
    H = out["hypotheses"]
    H["H1_the_intervals_of_version_6_hold"] = dict(c, asked={"share_inside_at_least": INSIDE, "median_z2_between": MEDIAN, "refuted_below": INSIDE_REFUTED, "or_median_z2_outside": MEDIAN_REFUTED}, outcome=verdict(c),
                                                   S2_for_comparison=dict(out["intervals"]["S2"], as_if_scored=verdict(out["intervals"]["S2"])),
                                                   by_lot=dict(out["by_how_the_unit_was_chosen"]["by lot"]["S3"], as_if_scored=verdict(out["by_how_the_unit_was_chosen"]["by lot"]["S3"])),
                                                   by_stratum={g: dict(v["S3"], as_if_scored=verdict(v["S3"])) for g, v in out["by_stratum"].items()})
    # --- what the author expected beside it
    e = out["errors"]
    t = R["edges"]["tested"]
    out["expected_before_the_states"] = EXPECTED
    out["expectations"] = {"the_errors_of_the_three_systems_agree_to_two_digits": len({round(e[s], 2) for s in ("S0", "S2", "S3")}) == 1,
                           "at_least_one_of_the_three_edge_tests_admits": any(x["status"] == "admitted" for x in t)}
    out["edges"] = {"tests": len(t), "admitted": [x["edge_id"] for x in t if x["status"] == "admitted"], "refused": [x["edge_id"] for x in t if x["status"] == "refused"],
                    "undecided": [x["edge_id"] for x in t if x["status"] == "undecided"], "with_no_row_in_the_batch": [x["edge_id"] for x in R["edges"]["with_no_row_in_the_batch"]],
                    "rows": [{"edge": x["edge_id"], "status": x["status"], "proposed_with_z": x["evidence"]["z"], **x["test"]} for x in t]}
    # --- the turn by rule
    parts = {"reading.json": True, "by_rule.json": (a.results / "by_rule.json").exists(), "accounts.json": (a.results / "accounts.json").exists()}
    lib_share = {(s["event"], s["readout"]): s["share"] for s in fc["library"]["shares"]}
    own = []
    for r in sets["with an edge"]:
        if r["scale"] is not None and (r["event"], r["readout"]) in lib_share:
            own.append({"stratum": r["stratum"], "event": r["event"], "readout": r["readout"], "scale": r["scale"], "share_under_interferon": lib_share[(r["event"], r["readout"])],
                        "share_here": (r["change"] - (r["S3"]["forecast"] - r["S3"]["by_its_edge"])) / r["scale"]})
    admitted = []
    for x in t:
        if x["status"] != "admitted":
            continue
        s_here, s_lib = x["test"]["scale_of_the_readout_here"], fc["scale"][x["readout"]]["library"]
        admitted.append({"edge": x["edge_id"], "class": x["class"], "direct_part_here": x["test"]["direct_part"], "scale_here": s_here, "share_here": x["test"]["direct_part"] / s_here if s_here >= EV.MIN_SCALE else None,
                         "direct_part_under_interferon": x["evidence"]["direct"], "scale_under_interferon": s_lib, "share_under_interferon": x["evidence"]["direct"] / s_lib if s_lib >= EV.MIN_SCALE else None})
    twice = [(o["share_under_interferon"], o["share_here"]) for o in own] + [(x["share_under_interferon"], x["share_here"]) for x in admitted if x["share_here"] is not None and x["share_under_interferon"] is not None]
    leaves = {"records": {g: {"units": len(v), "pairs": sum(len(u["changes"]) for u in v.values())} for g, v in R["records"].items()},
              "edges_admitted": admitted, "edges_refused": out["edges"]["refused"], "edges_still_proposed": out["edges"]["with_no_row_in_the_batch"],
              "own_gene_pairs_now_measured_twice": own, "spread_of_a_carried_share": {"before": fc["pairs_measured_twice"], "pairs_added": len(twice), "of_the_pairs_added_alone": EV.carried_spread(twice)}}
    if parts["by_rule.json"]:
        Q = json.loads((a.results / "by_rule.json").read_text())
        lt = Q["edges"]["tested_on_the_library"]
        leaves["the_laws_fitted_again"] = {g: {k: v[k] for k in ("recorded_events", "spread_with_no_edge", "spread_from")} | {"width": v["width"]["no type"], "flags": len(v["flags"]), "flags_of_units_of_the_batch": sum(f["a_unit_of_the_batch"] for f in v["flags"]),
                                               "flags_that_ask_for_a_record_of_a_type": [f"{f['event']}->{f['readout']}" for f in v["flags"] if f["flag"] != "a missing edge"]} for g, v in Q["strata"].items()}
        leaves["edges_proposed_by_the_flags"] = [{"edge": x["edge_id"], "z_here": x["evidence"]["z"], "on_the_library's_row": x["status"], **({"z_there": x["test"]["z"]} if "z" in x.get("test", {}) else {"why": x.get("test", {}).get("why")})} for x in lt]
    if parts["accounts.json"]:
        A = json.loads((a.results / "accounts.json").read_text())
        says = {}
        for v in A["accounts"].values():
            says[v["chain"]["says"]] = says.get(v["chain"]["says"], 0) + 1
        out["accounts"] = {"pairs_with_a_change": A["pairs_with_a_change"], "accounts": len(A["accounts"]), "what_they_say": says,
                           "rows": [{"pair": k, "chosen": v["chosen"], "change": v["change_in_spreads"], "snr": v["snr"], "carried": {s: v["stages"][s]["shares"]["carried"]["share"] for s in v["stages"]}, "says": v["chain"]["says"]} for k, v in A["accounts"].items()],
                           "open_parts": [{k: it.get(k) for k in ("id", "event", "readout", "group", "reading", "lever", "question", "state", "note")} for it in A["open_parts"]]}
        states = {}
        for it in A["open_parts"]:
            states[f"{it.get('state')} | {it.get('lever')}"] = states.get(f"{it.get('state')} | {it.get('lever')}", 0) + 1
        leaves["open_parts"] = {"count": len(A["open_parts"]), "by_state_and_lever": states}
    out["by_rule_the_batch_leaves"] = leaves
    out["gate_by_the_numbers"] = {"H1": H["H1_the_intervals_of_version_6_hold"]["outcome"], "the_turn_by_rule": parts,
                                  "passed": H["H1_the_intervals_of_version_6_hold"]["outcome"] == "holds" and all(parts.values())}
    # --- explored, declared before the reading: the isotype antibodies as the reference
    def bg(rs):
        if not rs:
            return None
        d = lambda r, s: r["against_the_background"]["change"] - r[s]["forecast"]
        wide = lambda r: np.sqrt(max(r["S3"]["sd"] ** 2 + r["against_the_background"]["se"] ** 2 - r["se"] ** 2, 1e-12)) * (r["S3"]["interval"][1] - r["S3"]["interval"][0]) / 2 / r["S3"]["sd"]
        return {"pairs": len(rs), "error": {"no_change": rms([r["against_the_background"]["change"] for r in rs]), **{s: rms([d(r, s) for r in rs]) for s in ("S0", "S2", "S3")}},
                "S3_inside_with_the_error_of_this_measure": int(sum(abs(d(r, "S3")) <= wide(r) for r in rs)), "S3_inside_as_scored": int(sum(r["S3"]["inside"] for r in rs))}
    out["explored_against_the_background"] = {"all_pairs": bg(rows), "by_stratum": {g: bg([r for r in rows if r["stratum"] == g]) for g in B.STRATA}}
    (a.results / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    c2 = out["intervals"]["S2"]
    print(f"pairs {out['pairs']} | error: no change {e['no_change']:.3f}, S0 {e['S0']:.3f}, S2 {e['S2']:.3f}, S3 {e['S3']:.3f}")
    print(f"  inside S3 {c['inside']} ({c['share_inside']:.3f}), median z2 {c['median_z2']:.2f}, mean z2 {c['mean_z2']:.2f} | S2 {c2['inside']} ({c2['share_inside']:.3f}), median z2 {c2['median_z2']:.2f}")
    h = H["H1_the_intervals_of_version_6_hold"]
    print(f"  H1: {h['outcome']} | by lot {h['by_lot']['inside']}/{h['by_lot']['pairs']} median {h['by_lot']['median_z2']:.2f} | " + " | ".join(f"{g} {v['inside']}/{v['pairs']} median {v['median_z2']:.2f}" for g, v in h["by_stratum"].items()))
    print("  edges:", [(x["edge"], x["status"], round(x["z"], 1)) for x in out["edges"]["rows"]], "| expectations:", out["expectations"])
    print("  gate by the numbers:", out["gate_by_the_numbers"])


if __name__ == "__main__":
    main()
