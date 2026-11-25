#!/usr/bin/env python3
"""The forecasts of parts B, D and E against results/round_b.json and round_e.json: results/outcome_b.json.

  python score_b.py [--results RESULTS]
"""
import argparse
import json
from pathlib import Path

RUN = Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    dec = json.loads((RUN / "results" / "declared.json").read_text())["parts_B_D_E"]
    B, E = json.loads((a.results / "round_b.json").read_text()), json.loads((a.results / "round_e.json").read_text())
    T, HB, T1e = B["targets"], B["hypotheses"], E["targets"]["T1"]["stages"]
    got = {**{h.split("_")[0]: v["reading"] for h, v in HB.items()}, "HE": E["hypotheses"]["HE_a_trained_last_layer_adds_a_step_under_the_new_repair"]["reading"],
           "u-T2-T": T["T2"]["stages"]["T"]["arms"]["rows"]["u"], "u-T2-P": T["T2"]["stages"]["P"]["arms"]["rows"]["u"], "u-IRF1-T": T["T1"]["stages"]["T"]["arms"]["rows"]["u"],
           "u-IRF1-T-derived": T["T1"]["events_one_by_one"]["IRF1"]["T"]["derived"]["u"], "u-IRF1-rows5": T1e["rows5"]["arms"]["rows"]["u"], "u-IRF1-rows5-none": T1e["rows5"]["arms"]["controls"]["u"]}
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
    (a.results / "outcome_b.json").write_text(json.dumps(out, indent=1) + "\n")
    print(f"intervals: {hit} of {n_i} held the value | outcomes: the class given the highest probability happened in {top} of {n_o}")
    for r in rows:
        print(f"  {r['id']:18s} {str(r.get('outcome', round(r.get('value', 0), 3))):12s} forecast {json.dumps(r.get('probabilities', r.get('interval')))}")


if __name__ == "__main__":
    main()
