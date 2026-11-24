#!/usr/bin/env python3
"""The forecasts of DESIGN.md against what was read: results/outcome.json and a table.

  python score.py

An interval forecast is scored inside or outside; an outcome forecast by whether the class given the highest probability happened.
"""
import json
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent
D = json.loads((RUN / "results" / "declared.json").read_text())["forecasts"]
R = json.loads((RUN / "results" / "stages.json").read_text())
S = R["IRF1"]["stages"]


def observed(i):
    if i.endswith("-u"):
        return S[i[:-2]]["u"]
    if i.endswith("-from-input"):
        return S[i.split("-")[0]]["from_the_input"]["reading"]
    if i == "first-loss":
        f = R["located"]["first_stage_that_is_a_loss"]
        return "none" if f is None else ("a layer" if f.startswith("L") else f)
    if i == "layers-bring-back":
        return R["layers"]["L12_against_E0"]["reading"]
    if i == "mean-at-the-end":
        return R["pooling"]["at_the_end"]["the_single_mean_loses_against_eight_windows"]["reading"]
    if i == "control":
        return "some" if R["control"]["stages_that_are_a_loss_from_the_input"] else "none"
    if i == "PERM":
        return "yes" if R["checks"]["control_view_reads_nothing"] else "no"
    raise KeyError(i)


def main():
    rows, n_in, n_iv, n_top, n_out = [], 0, 0, 0, 0
    for f in D:
        x = observed(f["id"])
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
    (RUN / "results" / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(f"intervals: {n_in} of {n_iv} held the value | outcomes: the class given the highest probability happened in {n_top} of {n_out}")
    for r in rows:
        x = r["observed"]
        print(f"  {r['id']:18s} {('%.3f' % x) if isinstance(x, float) else x!s:16s} forecast {json.dumps(r['forecast'])} -> {r.get('inside', r.get('most_likely_happened'))}")


if __name__ == "__main__":
    main()
