#!/usr/bin/env python3
"""The forecasts of DESIGN.md against what was read: results/outcome.json and a table.

  python score.py [--results RESULTS]

An interval forecast is scored inside or outside; an outcome forecast by whether the class given the highest probability happened.
"""
import argparse
import json
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent


def observed(i, R, M):
    T = R["targets"]
    if i == "T1-ref-C-u":
        return T["T1"]["C"]["reference"]["u"]
    if i == "T1-ref-ZC-u":
        return T["T1"]["ZC"]["reference"]["u"]
    if i in ("T1-rows-ZC-u", "T1-cells-ZC-u", "T1-controls-ZC-u", "T2-rows-ZC-u", "T1-rows-C-u"):
        t, arm, st, _u = i.split("-")
        return T[t][st]["arms"][arm]["u"]
    if i == "T1-rows-gain-ZC":
        return T["T1"]["gain"]["ZC"]["rows"]["gain"]
    if i in ("T1-rows-gain-ZC-reading", "T2-rows-gain-ZC-reading"):
        return T[i[:2]]["gain"]["ZC"]["rows"]["reading"]
    if i in ("H1", "H2", "C1", "C2", "C3"):
        return R["hypotheses"][i]
    if i == "T1-ZC-against-ZV":
        return T["T1"]["annotation_position_against_input_view"]["rows"]["reading"]
    if i == "grades-kept":
        return "yes" if all(R["hypotheses"]["grades_kept"].values()) else "no"
    if i == "records-repeat":
        c = R["checks_against_the_records"]
        return "yes" if max(c["largest_difference_of_an_arm_at_Z"], c["largest_difference_of_an_arm_at_ZV"]) <= c["limit"] else "no"
    if i == "mel-says-disease":
        s = M["strata"]["mel_co"]["the_model_says_at_the_placeholders"]["disease"]["ranked_first_share_of_cells"]
        w = max(s, key=s.get)
        return w if w in ("cancer_disease", "healthy") else "another"
    raise KeyError(i)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    D = json.loads((RUN / "results" / "declared.json").read_text())["forecasts"]
    R = json.loads((a.results / "round.json").read_text())
    M = json.loads((a.results / "stages_manifest.json").read_text())
    rows, n_in, n_iv, n_top, n_out = [], 0, 0, 0, 0
    for f in D:
        x = observed(f["id"], R, M)
        if f["kind"] == "interval":
            ok = bool(f["interval"][0] <= x <= f["interval"][1])
            n_iv += 1
            n_in += ok
            rows.append({"id": f["id"], "what": f["what"], "forecast": f["interval"], "observed": x, "inside": ok})
        else:
            top = max(f["probabilities"], key=f["probabilities"].get)
            ok = bool(x == top)
            n_out += 1
            n_top += ok
            rows.append({"id": f["id"], "what": f["what"], "forecast": f["probabilities"], "observed": x, "probability_given": f["probabilities"].get(x, 0.0), "most_likely_happened": ok})
    out = {"intervals": {"inside": n_in, "of": n_iv}, "outcomes": {"most_likely_happened": n_top, "of": n_out}, "rows": rows}
    (a.results / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(f"intervals: {n_in} of {n_iv} held the value | outcomes: the class given the highest probability happened in {n_top} of {n_out}")
    for r in rows:
        x = r["observed"]
        print(f"  {r['id']:26s} {('%.3f' % x) if isinstance(x, float) else x!s:20s} forecast {json.dumps(r['forecast'])} -> {r.get('inside', r.get('most_likely_happened'))}")


if __name__ == "__main__":
    main()
