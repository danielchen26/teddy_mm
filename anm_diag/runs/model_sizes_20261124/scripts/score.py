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
    S = R["IRF1"]["stages"]
    if i.endswith("-u"):
        return S[i[:-2]]["u"]
    if i.endswith("-from-input"):
        return S[i.split("-")[0]]["from_the_input"]["reading"]
    m = i.split("-")[-1]
    if i.startswith("Z") and i.endswith("-against-small"):
        return R["sizes"][i[1:4]]["against_the_small_model"]["reading"]
    if i.startswith("statement-"):
        return R["sizes"][m]["statement"]
    if i.startswith("last-layer-"):
        return R["sizes"][m]["last_layer"]["reading"]
    if i.startswith("C") and i.endswith("-against-state"):
        return R["sizes"][i[1:4]]["annotation_position_against_the_state"]["reading"]
    if i.startswith("placeholders-"):
        return R["sizes"][m]["placeholders_against_none"]["reading"]
    if i == "control":
        return "holds" if R["control"]["holds"] else "fails"
    if i == "records-repeat":
        return "yes" if R["checks"]["hold"] else "no"
    if i.startswith("says-disease-"):
        s = M["sizes"][m]["the_model_says_at_the_placeholders"]["disease"]["ranked_first_share_of_cells"]
        w = max(s, key=s.get)
        return w if w in ("infectious_disease", "cancer_disease") else "another"
    raise KeyError(i)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    D = json.loads((RUN / "results" / "declared.json").read_text())["forecasts"]
    R = json.loads((a.results / "sizes.json").read_text())
    M = json.loads((a.results / "sizes_manifest.json").read_text())
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
        print(f"  {r['id']:22s} {('%.3f' % x) if isinstance(x, float) else x!s:34s} forecast {json.dumps(r['forecast'])} -> {r.get('inside', r.get('most_likely_happened'))}")


if __name__ == "__main__":
    main()
