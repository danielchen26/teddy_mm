#!/usr/bin/env python3
"""The forecasts of DESIGN.md against what was read: results/outcome.json and a table.

  python score.py

An interval forecast is scored inside or outside; an outcome forecast by whether the class given the highest probability happened.
"""
import json
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent
D = json.loads((RUN / "results" / "declared.json").read_text())["forecasts"]
R = json.loads((RUN / "results" / "round.json").read_text())
T = R["targets"]


def observed(i):
    if i == "CHECK-Z":
        return R["check_against_the_record"]["largest_difference_of_an_arm_at_Z"]
    if i in R["hypotheses"]:
        return R["hypotheses"][i]
    if i == "C3-support":
        return "yes" if T["C3"]["ZV"]["arms"]["rows"]["verdict"] == "support" else "no"
    t, what, *more = i.split("-")
    x = T[t]
    if what == "ref":
        return x["ZV"]["reference"]["u"]
    if what in ("rows", "cells"):
        return x["ZV"]["arms"][what]["u"]
    if what == "gain" and more[-1] == "reading":
        return x["gain"][more[0]]["reading"]
    if what == "gain":
        return x["gain"][more[0]]["gain"]
    if what == "grade":
        g = x["ZV"]["closure"]["grade"]
        return {"relative": "relative or absolute", "absolute": "relative or absolute", "not a target for a repair": "not a target", "nothing to close": "not a target"}.get(g, g)
    raise KeyError(i)


def main():
    rows, n_in, n_iv, n_top, n_out = [], 0, 0, 0, 0
    for f in D:
        v = observed(f["id"])
        if f["kind"] == "interval":
            ok = bool(f["interval"][0] <= v <= f["interval"][1])
            n_iv += 1
            n_in += ok
            rows.append({"id": f["id"], "what": f["what"], "forecast": f["interval"], "observed": v, "inside": ok})
        else:
            top = max(f["probabilities"], key=f["probabilities"].get)
            ok = bool(v == top)
            n_out += 1
            n_top += ok
            rows.append({"id": f["id"], "what": f["what"], "forecast": f["probabilities"], "observed": v, "probability_given": f["probabilities"].get(v, 0.0), "most_likely_happened": ok})
    out = {"intervals": {"inside": n_in, "of": n_iv}, "outcomes": {"most_likely_happened": n_top, "of": n_out}, "rows": rows}
    (RUN / "results" / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(f"intervals: {n_in} of {n_iv} held the value | outcomes: the class given the highest probability happened in {n_top} of {n_out}")
    for r in rows:
        v = r["observed"]
        print(f"  {r['id']:22s} {('%.3f' % v) if isinstance(v, float) else v!s:28s} forecast {json.dumps(r['forecast'])} -> {r.get('inside', r.get('most_likely_happened'))}")


if __name__ == "__main__":
    main()
