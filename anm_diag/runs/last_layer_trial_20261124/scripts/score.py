#!/usr/bin/env python3
"""The forecasts of DESIGN.md against what was read: results/outcome.json and a table.

  python score.py [--part A|B] [--results RESULTS]

Every forecast names its part (A: results/trial.json, B: results/trial_B.json) and where its value stands there
(`at`: the keys, in order). The forecasts of one part are scored into results/outcome_<part>.json. An interval forecast is
scored inside or outside; an outcome forecast by whether the class given the highest probability happened.
"""
import argparse
import json
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent


def at(R, path):
    for k in path:
        R = R[k]
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=RUN / "results")
    ap.add_argument("--part", choices=("A", "B"), default="A")
    a = ap.parse_args()
    D = [f for f in json.loads((RUN / "results" / "declared.json").read_text())["forecasts"] if f["part"] == a.part]
    R = json.loads((a.results / ("trial.json" if a.part == "A" else "trial_B.json")).read_text())
    rows, n_in, n_iv, n_top, n_out = [], 0, 0, 0, 0
    for f in D:
        x = at(R, f["at"])
        if isinstance(x, bool):
            x = "yes" if x else "no"
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
    (a.results / f"outcome_{a.part}.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(f"intervals: {n_in} of {n_iv} held the value | outcomes: the class given the highest probability happened in {n_top} of {n_out}")
    for r in rows:
        x = r["observed"]
        print(f"  {r['id']:26s} {('%.3f' % x) if isinstance(x, float) else x!s:22s} forecast {json.dumps(r['forecast'])} -> {r.get('inside', r.get('most_likely_happened'))}")


if __name__ == "__main__":
    main()
