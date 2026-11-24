#!/usr/bin/env python3
"""The forecasts of DESIGN.md against what was read: results/outcome.json and a table.

  python score.py [--results RESULTS]

An interval forecast is scored inside or outside; an outcome forecast by whether the class given the highest probability happened.
"""
import argparse
import json
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent
ORDER = ["disease", "tissue", "cell_type", "sex"]


def observed(i, R, M):
    S = R["IRF1"]["stages"]
    said = M["the_model_says_at_the_placeholders"]
    first = lambda name: max(said[name]["ranked_first_share_of_cells"], key=said[name]["ranked_first_share_of_cells"].get)
    if i.endswith("-u"):
        return S[i[:-2]]["u"]
    if i.endswith("-from-input"):
        return S[i.split("-")[0]]["from_the_input"]["reading"]
    if i.endswith("-against-none"):
        return S[i.split("-")[0]]["against_no_annotations"]["reading"]
    if i.endswith("-moves"):
        return M["the_state_moves"][i.split("-")[0]]["common_shift_taken_out"]["median"]
    if i == "statement":
        w = R["annotations_in_the_input"]["statement"]
        return "brings back" if w.startswith("brings") else w
    if i == "content":
        return R["content"]["labels_of_the_line_against_wrong_labels"]["reading"]
    if i == "control":
        return "holds" if R["control"]["holds"] else "fails"
    if i == "PERM":
        return "yes" if R["checks"]["control_view_reads_nothing"] else "no"
    if i == "own-category":
        groups = {"disease": lambda w: w.endswith("_disease") or w == "healthy", "tissue": lambda w: w.endswith("_tissue"), "cell_type": lambda w: w.endswith("_cell") or w == "skeletal_muscle",
                  "sex": lambda w: w in ("female", "male")}
        ok = all(sum(s for w, s in said[name]["ranked_first_share_of_cells"].items() if groups[name](w)) >= 0.9 for name in ORDER)
        return "yes" if ok else "no"
    if i.startswith("says-"):
        name = i[5:]
        w = first(name)
        known = {"disease": ("cancer_disease", "healthy"), "tissue": ("hematopoietic_tissue", "immune_tissue"), "cell_type": ("precursor_cell", "immune_cell", "hematopoietic_cell"), "sex": ("male", "female")}[name]
        return w if w in known else "another"
    if i == "labels-untrained":
        return "yes" if M["input_embeddings"]["label_tokens_look_untrained"] else "no"
    raise KeyError(i)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    D = json.loads((RUN / "results" / "declared.json").read_text())["forecasts"]
    R = json.loads((a.results / "modes.json").read_text())
    M = json.loads((a.results / "modes_manifest.json").read_text())
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
        print(f"  {r['id']:20s} {('%.3f' % x) if isinstance(x, float) else x!s:22s} forecast {json.dumps(r['forecast'])} -> {r.get('inside', r.get('most_likely_happened'))}")


if __name__ == "__main__":
    main()
