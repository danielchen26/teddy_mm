#!/usr/bin/env python3
"""The forecasts of part A against results/round_a.json: results/outcome_a.json.

  python score_a.py [--results RESULTS]
"""
import argparse
import json
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    dec = json.loads((RUN / "results" / "declared.json").read_text())["part_A"]
    R = json.loads((a.results / "round_a.json").read_text())
    T, H = R["targets"], R["hypotheses"]
    new = lambda t, st: T[t]["stages"][st]["with_the_rows_applied"]
    ev = lambda t, arm: [abs(q["Z"]["arms"][arm]["u"]) for q in T[t]["events_one_by_one"].values()]
    ratios = [r for r in H["H7_the_rule_on_the_records"]["derived_over_measured"].values()]
    got = {**{h[:2]: v["reading"] for h, v in H.items()}, "grade-T2": new("T2", "Z")["closure"]["grade"],
           "type-or-derived": "type" if np.median(ev("T2", "type")) <= np.median(ev("T2", "derived")) else "derived",
           "u-T2-Z": new("T2", "Z")["arms"]["rows"]["u"], "u-T2-ZV": new("T2", "ZV")["arms"]["rows"]["u"], "gain-T2-Z": T["T2"]["stages"]["Z"]["new_rows_over_old_rows"]["gain"],
           "u-C2-Z": new("C2", "Z")["arms"]["rows"]["u"], "u-C3-Z": new("C3", "Z")["arms"]["rows"]["u"], "derived-MHC1": float(np.median(ev("T2", "derived"))),
           "type-MHC1": float(np.median(ev("T2", "type"))), "ratio-15": float(np.median(ratios))}
    rows, hit, n_i, top, n_o = [], 0, 0, 0, 0
    for f in dec:
        v = got[f["id"]]
        if f["kind"] == "interval":
            ok = bool(f["interval"][0] <= v <= f["interval"][1]); hit += ok; n_i += 1
            rows.append({**f, "value": v, "inside": ok})
        else:
            best = max(f["probabilities"], key=f["probabilities"].get); ok = bool(best == v); top += ok; n_o += 1
            rows.append({**f, "outcome": v, "probability_given": f["probabilities"].get(v, 0.0), "the_most_probable_class_happened": ok})
    out = {"intervals": {"held": hit, "of": n_i}, "outcomes": {"the_most_probable_class_happened": top, "of": n_o}, "rows": rows}
    (a.results / "outcome_a.json").write_text(json.dumps(out, indent=1) + "\n")
    print(f"intervals: {hit} of {n_i} held the value | outcomes: the class given the highest probability happened in {top} of {n_o}")
    for r in rows:
        print(f"  {r['id']:16s} {str(r.get('outcome', round(r.get('value', 0), 3))):26s} forecast {json.dumps(r.get('probabilities', r.get('interval')))}")


if __name__ == "__main__":
    main()
