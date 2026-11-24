#!/usr/bin/env python3
"""The forecasts of DESIGN.md (Section 6) against what was read: results/outcome.json and a table.

  python score.py

An interval forecast is scored inside or outside; an outcome forecast by whether the class given the highest probability happened.
"""
import json
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent
D = json.loads((RUN / "results" / "declared.json").read_text())["forecasts"]
S = json.loads((RUN / "results" / "round.json").read_text())["summary"]
A = S["among_A3"]
OBS = {"N-pairs": S["pairs"], "N-proteins": S["proteins"], "T-rejected": S["as_trained"]["rejected"], "T-median": S["as_trained"]["median_abs_u"], "Z-supported": S["reference"]["support"],
       "Z-rejected": S["reference"]["reject"], "A3-share": S["A3"]["holds"], "R-closed": A["closed"], "R-absolute": A["absolute"], "R-fall": A["median_fall_to_rows"],
       "U-needed": A["own_row_needed"], "U-not-needed": A["own_row_not_needed"], "U-median": A["median_unseen_minus_rows"], "C-rows": A["rows_add_a_step_over_controls"],
       "P-unreadable": S["among_reference_rejected"]["unreadable"], "P-upstream": S["among_reference_rejected"]["upstream"],
       "K-control": "yes" if S["control_view_reads_nothing_for_every_pair"] else "no", "K-positive": "yes" if S["positive_control_holds_for_some_protein"] else "no"}


def main():
    rows, n_in, n_iv, n_top, n_out, skipped = [], 0, 0, 0, 0, 0
    for f in D:
        v = OBS[f["id"]]
        if v is None:                                           # its condition did not occur (no pair with A3, or none with a rejected reference)
            skipped += 1
            rows.append({"id": f["id"], "what": f["what"], "forecast": f.get("interval", f.get("probabilities")), "observed": None, "scored": False})
        elif f["kind"] == "interval":
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
    out = {"intervals": {"inside": n_in, "of": n_iv}, "outcomes": {"most_likely_happened": n_top, "of": n_out}, "not_scored": skipped, "rows": rows}
    (RUN / "results" / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print(f"intervals: {n_in} of {n_iv} held the value | outcomes: the class given the highest probability happened in {n_top} of {n_out} | not scored: {skipped}")
    for r in rows:
        v = r["observed"]
        print(f"  {r['id']:14s} {('%.3f' % v) if isinstance(v, float) else v!s:10s} forecast {json.dumps(r['forecast'])} -> {r.get('inside', r.get('most_likely_happened'))}")


if __name__ == "__main__":
    main()
