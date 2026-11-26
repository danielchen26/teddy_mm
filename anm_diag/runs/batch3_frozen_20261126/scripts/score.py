#!/usr/bin/env python3
"""The hypotheses of DESIGN.md, Section 7, against the reading. results/outcome.json. Nothing here is chosen by hand.

  python score.py [--results DIR]

The sets the hypotheses are scored on were fixed in forecasts.json before stage B was opened. What the design left
open is fixed here, before the reading:

- H3, H4 and H5 have no middle: each holds or does not hold.
- H6 counts in two ways, because the same proposed pair can have a row in both conditions. By test: a proposal asked of
  one stratum's row. By proposal: a pair is admitted if a test admits it. H6 holds if at least half are admitted in
  both counts, does not hold if fewer than half in both, and is unresolved if the counts disagree.
- H7 (b): the state carries a step less than the input where the interval of carried(V) - carried(Z) lies above 0.10;
  the typed state brings back half if the gain of T over Z is at least half of that loss.
- The gate by the numbers: H1 holds and H2 holds. Whether an error of the method changes a verdict is said in the README.

Known before the reading and left as it is: the readout of the melanoma screen is a centred log-ratio, which couples
the 24 columns of a cell. A unit that removes a large protein raises every other readout a little, and no system here
forecasts that. `explored_against_the_background` scores the same forecasts against the change measured with the
isotype antibodies as the reference instead of the mean of the columns. It is declared here, enters no hypothesis and
changes no verdict.

By rule, the batch leaves: its rows as records; the admitted proposals as edges, with the share of the scale they
showed; and the pairs now measured twice, which give the spread of a carried share (rule 5d).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _b3 as B  # noqa: E402
from anm_diag import events as EV  # noqa: E402

EXPECTED = {"H1": 0.55, "H2": 0.75, "H3": 0.6, "H4": 0.5, "H5": 0.8, "H6": 0.4, "H7a": 0.8, "H7b_not_tested": 0.6, "gate": 0.45}
rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None
key = lambda r: f"{r['stratum']}|{r['event']}|{r['readout']}"


def errors(rows):
    return {"pairs": len(rows), "no_change": rms([r["change"] for r in rows]), "S0": rms([r["S0"]["error"] for r in rows]), "S1": rms([r["S1"]["error"] for r in rows]), "S2": rms([r["S2"]["error"] for r in rows])}


def intervals(rows, s="S2"):
    z2 = [r[s]["z"] ** 2 for r in rows]
    n = sum(r[s]["inside"] for r in rows)
    return {"pairs": len(rows), "inside": n, "share_inside": n / len(rows) if rows else None, "mean_z2": float(np.mean(z2)) if rows else None, "median_z2": float(np.median(z2)) if rows else None}


def brief(r):
    return {"stratum": r["stratum"], "event": r["event"], "readout": r["readout"], "kind": r["kind"], "group": r["group"], "edge": r["edge"], "type": r["type"], "change": r["change"], "se": r["se"],
            "S0": r["S0"]["forecast"], "S1": r["S1"]["forecast"], "S2": r["S2"]["forecast"], "interval_S2": r["S2"]["interval"], "inside_S1": r["S1"]["inside"], "inside_S2": r["S2"]["inside"], "z_S2": r["S2"]["z"]}


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
    out = {"pairs": len(rows), "errors": errors(rows), "intervals": {s: intervals(rows, s) for s in ("S1", "S2")},
           "by_stratum": {g: {"errors": errors(x), "S1": intervals(x, "S1"), "S2": intervals(x, "S2")} for g in B.STRATA for x in [[r for r in rows if r["stratum"] == g]]},
           "by_group": {f"{g} | {grp}": {"errors": errors(x), "S1": intervals(x, "S1"), "S2": intervals(x, "S2")} for g in B.STRATA for grp in ("read before", "never read")
                        for x in [[r for r in rows if r["stratum"] == g and r["group"] == grp]] if x},
           "pairs_with_a_change": {"errors": errors(x), "S1": intervals(x, "S1"), "S2": intervals(x, "S2"), "rows": [brief(r) for r in sorted(x, key=lambda r: -r["snr"])]} if (x := [r for r in rows if r["snr"] >= R["min_snr"]]) else {"errors": errors([])},
           "outside_S2": [brief(r) for r in sorted(rows, key=lambda r: -abs(r["S2"]["z"])) if not r["S2"]["inside"]],
           "outside_S2_by_unit": {}, "outside_S2_by_readout": {}, "hypotheses": {}}
    for r in rows:
        if not r["S2"]["inside"]:
            u, q = f"{r['stratum']}|{r['event']}", f"{r['stratum']}|{r['readout']}"
            out["outside_S2_by_unit"][u] = out["outside_S2_by_unit"].get(u, 0) + 1
            out["outside_S2_by_readout"][q] = out["outside_S2_by_readout"].get(q, 0) + 1
    H = out["hypotheses"]
    # --- H1
    c = out["intervals"]["S2"]
    ok = c["share_inside"] >= 0.90 and 0.5 <= c["mean_z2"] <= 2
    bad = c["share_inside"] < 0.85 or not (0.33 <= c["mean_z2"] <= 3)
    H["H1_the_intervals_of_version_5_hold"] = dict(c, outcome="holds" if ok else "refuted" if bad else "unresolved", S1_for_comparison=out["intervals"]["S1"])
    # --- H2
    e_all, e_dif = errors(rows), errors(sets["the_systems_can_differ"])
    least = lambda e: e["S2"] < e["S1"] and e["S2"] < e["S0"]
    H["H2_the_loop_corrects"] = {"all_pairs": e_all, "pairs_on_which_the_systems_can_differ": e_dif,
                                 "outcome": "holds" if least(e_all) and least(e_dif) else "refuted" if not least(e_dif) else "unresolved"}
    # --- H3
    x = sets["H3"]
    e, c = errors(x), intervals(x)
    H["H3_an_edge_has_the_size_of_its_pair"] = {"errors": e, "intervals": c, "rows": [brief(r) for r in x],
                                                "outcome": "not tested: fewer than four such pairs" if len(x) < 4 else "holds" if least(e) and c["share_inside"] >= 0.75 else "does not hold"}
    # --- H4
    x = sets["H4"]
    e, c = errors(x), intervals(x)
    H["H4_the_records_of_a_type_forecast_its_other_events"] = {"errors": e, "intervals": c, "by_stratum": {g: {"errors": errors(y), "intervals": intervals(y)} for g in B.STRATA for y in [[r for r in x if r["stratum"] == g]] if y},
                                                               "rows": [brief(r) for r in x], "outcome": "holds" if e["S2"] < e["S0"] and c["share_inside"] >= 0.80 else "does not hold"}
    # --- H5
    x = sets["H5"]
    e = errors(x)
    H["H5_they_move_no_event_of_another_type"] = {"errors": e, "intervals": {s: intervals(x, s) for s in ("S1", "S2")}, "outcome": "holds" if e["S2"] <= e["S1"] else "does not hold"}
    # --- H6
    t = R["edges"]["tested"]
    pairs = sorted({(e["site"], e["readout"]) for e in t})
    adm_t = sum(e["status"] == "admitted" for e in t)
    adm_p = sum(any(e["status"] == "admitted" for e in t if (e["site"], e["readout"]) == p) for p in pairs)
    half_t, half_p = (adm_t >= 0.5 * len(t)) if t else None, (adm_p >= 0.5 * len(pairs)) if pairs else None
    H["H6_a_proposed_edge_holds_elsewhere"] = {"tests": len(t), "admitted_tests": adm_t, "proposals_with_a_row": len(pairs), "admitted_proposals": adm_p,
                                               "with_no_row_in_stage_B": [e["edge_id"] for e in R["edges"]["with_no_row_in_stage_B"]],
                                               "rows": [{"edge": e["edge_id"], "status": e["status"], "proposed_with_z": e["evidence"]["z"], **e["test"]} for e in t],
                                               "outcome": "not tested: fewer than four tests" if len(t) < 4 else "holds" if half_t and half_p else "does not hold" if not half_t and not half_p else "unresolved"}
    # --- H7
    acc_file = a.results / "accounts.json"
    if acc_file.exists():
        A = json.loads(acc_file.read_text())
        x = []
        for k, v in A["accounts"].items():
            if v["stratum"] not in ("mel_co", "thp1") or v["edge"] is not None:
                continue
            tz, vz = v["gains"]["T over Z"]["carried"], v["gains"]["V over Z"]["carried"]
            lost = vz["ci"][0] > 0.10
            x.append({"pair": k, "kind": v["kind"], "type": v["type"], "change": v["change_in_spreads"], "carried": {s: v["stages"][s]["shares"]["carried"]["share"] for s in ("V", "Z", "T", "R")},
                      "T_over_Z": tz, "T_over_R": v["gains"]["T over R"]["carried"], "V_over_Z": vz, "the_typed_state_carries_less": bool(tz["ci"][1] < -0.10),
                      "the_state_loses_a_step": bool(lost), "brought_back": (tz["gain"] / vz["gain"]) if lost else None, "says": v["chain"]["says"]})
        lost = [r for r in x if r["the_state_loses_a_step"]]
        a_ok = not any(r["the_typed_state_carries_less"] for r in x)
        b = "not tested: no pair loses a step before the state" if not lost else "holds" if all(r["brought_back"] >= 0.5 for r in lost) else "does not hold"
        H["H7_the_state_kept_by_type"] = {"pairs": len(x), "a_no_pair_where_the_typed_state_carries_less": a_ok, "pairs_where_the_state_loses_a_step": len(lost), "b_the_typed_state_brings_back_half": b, "rows": x,
                                          "outcome": "not tested: no pair with a change and no edge" if not x else f"(a) {'holds' if a_ok else 'does not hold'}; (b) {b}"}
        out["accounts"] = {"pairs_with_a_change": A["pairs_with_a_change"], "accounts": len(A["accounts"]),
                           "says": {k: v["chain"]["says"] for k, v in A["accounts"].items()},
                           "open_parts": [{k: it.get(k) for k in ("id", "event", "readout", "group", "reading", "lever", "question", "state", "note")} for it in A["open_parts"]]}
    else:
        H["H7_the_state_kept_by_type"] = {"outcome": "pending: the accounts are not read yet"}
    out["gate_by_the_numbers"] = {"H1": H["H1_the_intervals_of_version_5_hold"]["outcome"], "H2": H["H2_the_loop_corrects"]["outcome"],
                                  "passed": H["H1_the_intervals_of_version_5_hold"]["outcome"] == "holds" and H["H2_the_loop_corrects"]["outcome"] == "holds"}
    out["expected_before_stage_A"] = EXPECTED
    # --- explored, declared before the reading: the isotype antibodies as the reference (melanoma)
    x = [r for r in rows if "against_the_background" in r]
    def bg(rs):
        if not rs:
            return None
        d = lambda r, s: r["against_the_background"]["change"] - r[s]["forecast"]
        wide = lambda r: np.sqrt(max(r["S2"]["sd"] ** 2 + r["against_the_background"]["se"] ** 2 - r["se"] ** 2, 1e-12))
        return {"pairs": len(rs), "error": {"no_change": rms([r["against_the_background"]["change"] for r in rs]), **{s: rms([d(r, s) for r in rs]) for s in ("S0", "S1", "S2")}},
                "error_as_scored": {s: rms([r[s]["error"] for r in rs]) for s in ("S0", "S1", "S2")},
                "S2_inside_with_the_error_of_this_measure": int(sum(abs(d(r, "S2")) <= 1.96 * wide(r) for r in rs)), "S2_inside_as_scored": int(sum(r["S2"]["inside"] for r in rs))}
    out["explored_against_the_background"] = {"all_melanoma_pairs": bg(x), "pairs_with_no_edge_on_the_readouts_never_read": bg([r for r in x if r["edge"] is None and r["group"] == "never read"]),
                                              "by_stratum": {g: bg([r for r in x if r["stratum"] == g]) for g in ("mel_co", "mel_ctl")}, "units": []}
    for g in ("mel_co", "mel_ctl"):
        for u in sorted({r["event"] for r in x if r["stratum"] == g}):
            rs = [r for r in x if r["stratum"] == g and r["event"] == u and r["edge"] is None and r["group"] == "never read"]
            if rs:
                out["explored_against_the_background"]["units"].append({"stratum": g, "unit": u, "pairs": len(rs), "mean_error_S2": float(np.mean([r["S2"]["error"] for r in rs])),
                                                                        "mean_error_S2_against_the_background": float(np.mean([r["against_the_background"]["change"] - r["S2"]["forecast"] for r in rs])),
                                                                        "outside_S2": int(sum(not r["S2"]["inside"] for r in rs))})
    # --- by rule: what the batch leaves
    twice_before = fc["pairs_measured_twice"]
    lib_share = {(s["event"], s["readout"]): s["share"] for s in fc["library"]["shares"]}
    own = []
    for r in sets["H3"]:
        here = (r["change"] - (r["S2"]["forecast"] - r["S2"]["by_its_edge"])) / r["scale"]
        own.append({"stratum": r["stratum"], "event": r["event"], "readout": r["readout"], "scale": r["scale"], "share_under_interferon": lib_share[(r["event"], r["readout"])], "share_here": here})
    admitted = []
    for e in R["edges"]["tested"]:
        if e["status"] != "admitted":
            continue
        s_here, s_lib = e["test"]["scale_of_the_readout_here"], fc["scale"][e["readout"]]["library"]
        admitted.append({"edge": e["edge_id"], "class": e["class"], "direct_part_here": e["test"]["direct_part"], "scale_here": s_here, "share_here": e["test"]["direct_part"] / s_here if s_here >= EV.MIN_SCALE else None,
                         "direct_part_under_interferon": e["evidence"]["direct"], "scale_under_interferon": s_lib, "share_under_interferon": e["evidence"]["direct"] / s_lib if s_lib >= EV.MIN_SCALE else None})
    new_twice = [(o["share_under_interferon"], o["share_here"]) for o in own] + [(e["share_under_interferon"], e["share_here"]) for e in admitted if e["share_here"] is not None and e["share_under_interferon"] is not None]
    out["by_rule_the_batch_leaves"] = {"records": {g: {"units": len(v), "pairs": sum(len(u["changes"]) for u in v.values())} for g, v in R["records"].items()},
                                       "edges_admitted": admitted, "edges_refused": [e["edge_id"] for e in R["edges"]["tested"] if e["status"] == "refused"],
                                       "edges_undecided": [e["edge_id"] for e in R["edges"]["tested"] if e["status"] == "undecided"],
                                       "edges_still_proposed": [e["edge_id"] for e in R["edges"]["with_no_row_in_stage_B"]],
                                       "own_gene_pairs_now_measured_twice": own,
                                       "spread_of_a_carried_share": {"before": twice_before, "pairs_added": len(new_twice), "of_the_pairs_added_alone": EV.carried_spread(new_twice)}}
    (a.results / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    e, c = out["errors"], out["intervals"]
    print(f"pairs {out['pairs']} | error: no change {e['no_change']:.3f}, S0 {e['S0']:.3f}, S1 {e['S1']:.3f}, S2 {e['S2']:.3f} | inside: S1 {c['S1']['inside']}, S2 {c['S2']['inside']} | mean z2: S1 {c['S1']['mean_z2']:.2f}, S2 {c['S2']['mean_z2']:.2f}")
    for k, v in H.items():
        print(f"  {k}: {v['outcome']}")
    print("  gate by the numbers:", out["gate_by_the_numbers"])


if __name__ == "__main__":
    main()
