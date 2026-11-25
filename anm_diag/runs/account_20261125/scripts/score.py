#!/usr/bin/env python3
"""The hypotheses of DESIGN.md, Section 8, as declared, and the forecasts. results/outcome.json.

  python score.py

Reads results/melanoma_as_first_read.json (the reading of Section 7.4 with the chain as it was declared: its first part
was named `direct` and the input was the view V), results/planted.json and results/declared.json. It also checks that
results/melanoma.json, the same reading put together by the module as amended after the reading (README, "What was
not as written"), holds the same readings stage by stage.
"""
import json
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
R = RUN / "results"
MHC, PDL1 = ["MHC1_coculture"], ["PDL1_coculture", "PDL1_ifn"]
three = lambda ok, bad: "confirmed" if ok else ("refuted" if bad else "unresolved")


def main():
    first, amended, dec = (json.loads((R / f).read_text()) for f in ("melanoma_as_first_read.json", "melanoma.json", "declared.json"))
    planted = json.loads((R / "planted.json").read_text())
    worst = 0.0
    for t, T in first["targets"].items():
        for e, q in T["events"].items():
            for st, r in q["stages"].items():
                for k, v in r["shares"].items():
                    worst = max(worst, abs(v["share"] - amended["targets"][t]["events"][e]["stages"][st]["shares"][k]["share"]))
    acc = {(t, e): q for t, T in first["targets"].items() for e, q in T["events"].items()}
    made = {k: q for k, q in acc.items() if q["chain"]["snr"] >= 4.0}
    left_out = [f"{t} {e}" for (t, e), q in acc.items() if q["chain"]["snr"] < 4.0]
    mhc = {e: q for (t, e), q in made.items() if t in MHC}
    pdl = {f"{t} {e}": q for (t, e), q in made.items() if t in PDL1}
    z = lambda q: q["stages"]["Z"]
    h = {}
    fa = {e: q["found_again_over_carried_at_the_state"] for e, q in mhc.items()}
    n1 = sum(v["reading"] == "a step" for v in fa.values())
    h["HA1_the_mean_does_not_see"] = {"reading": three(n1 >= 4, n1 < 2), "MHC1_events_where_found_again_exceeds_carried_by_a_step": n1, "of": len(mhc), "by_event": fa}
    tz = {e: q["gains"]["T over Z"]["carried"] for e, q in mhc.items()}
    n2, med2 = sum(v["reading"] == "a step" for v in tz.values()), float(np.median([v["gain"] for v in tz.values()]))
    h["HA2_the_typed_state_carries_more"] = {"reading": three(n2 >= 3, n2 == 0 or med2 < 0.03), "MHC1_events_with_a_step": n2, "median_gain": med2, "by_event": tz}
    tr = {e: q["gains"]["T over R"]["carried"] for e, q in mhc.items()}
    n3, med3 = sum(v["reading"] == "a step" for v in tr.values()), float(np.median([v["gain"] for v in tr.values()]))
    h["HA3_it_is_the_type"] = {"reading": three(n3 >= 3, n3 == 0), "MHC1_events_with_a_step": n3, "median_gain": med3, "by_event": tr}
    lost = {e: q["chain"]["parts"]["lost_before_the_state"] for e, q in mhc.items()}
    n4, med4 = sum(v["reading"] == "a step" for v in lost.values()), float(np.median([v["share"] for v in lost.values()]))
    h["HA4_where_MHC1_is_lost"] = {"reading": three(n4 >= 3, n4 == 0), "MHC1_events_where_lost_before_the_state_is_open": n4, "median_share": med4, "by_event": lost}
    nothing = {k: (not q["chain"]["open"]) for k, q in pdl.items()}
    n5 = sum(nothing.values())
    h["HA5_PDL1_has_nothing_open"] = {"reading": three(n5 >= 7, n5 < 4), "accounts_with_no_open_part": n5, "of": len(pdl), "what_is_open": {k: [p["part"] for p in q["chain"]["open"]] for k, q in pdl.items()}}
    perm = {f"{t} {e}": q["must_fail"]["permuted_stage"]["carried"]["share"] for (t, e), q in acc.items()}
    fake = {f"{t} {e}": q["must_fail"]["made_up_event"]["no_change_to_account_for"] for (t, e), q in acc.items()}
    h["HA6_the_controls"] = {"reading": "confirmed" if max(abs(v) for v in perm.values()) < 0.15 and all(fake.values()) else "refuted", "largest_carried_by_a_permuted_stage": float(max(abs(v) for v in perm.values())),
                             "made_up_events_with_no_change_to_account_for": int(sum(fake.values())), "of": len(fake), "permuted_stage": perm}
    ov = {f"{t} {e}": z(q)["overlap"] for (t, e), q in acc.items()}
    n7 = sum(v >= 0.10 for v in ov.values())
    h["HA7_the_cells_can_be_read"] = {"reading": three(n7 == len(ov), n7 < 12), "pairs_with_overlap_at_least_0.10": n7, "of": len(ov), "smallest": float(min(ov.values())), "overlap": ov}
    oth = {f"{t} {e}": z(q)["shares"]["others"]["share"] for (t, e), q in made.items()}
    n8 = sum(v >= 0.75 for v in oth.values())
    h["HA8_transfer_inside_the_type"] = {"reading": three(n8 >= 10, n8 < 6), "pairs_with_others_at_least_0.75": n8, "of": len(oth), "others_at_Z": oth}
    ap = planted["as planted"]
    h["HP_planted_suite"] = {"reading": three(ap["wrong"] / ap["checks"] <= 0.02, ap["wrong"] / ap["checks"] > 0.05), "wrong": ap["wrong"], "checks": ap["checks"], "not_counted": dec["forecasts"]["HP_planted_suite"]["not_counted"]}
    value = {"HA1_the_mean_does_not_see": n1, "HA2_the_typed_state_carries_more": med2, "HA3_it_is_the_type": med3, "HA4_where_MHC1_is_lost": med4, "HA5_PDL1_has_nothing_open": n5,
             "HA7_the_cells_can_be_read": float(min(ov.values())), "HA8_transfer_inside_the_type": n8}
    fc, inside, n_iv, right, n_out = {}, 0, 0, 0, 0
    for k, f in dec["forecasts"].items():
        if "not_counted" in f:
            continue
        r = {"probability_confirmed": f["probability_confirmed"], "reading": h[k]["reading"]}
        said = "confirmed" if f["probability_confirmed"] >= 0.5 else "not confirmed"
        r["the_outcome_given_the_higher_probability"] = said
        r["happened"] = bool((h[k]["reading"] == "confirmed") == (said == "confirmed"))
        right, n_out = right + r["happened"], n_out + 1
        if "interval" in f:
            r.update(interval=f["interval"], value=value[k], inside=bool(f["interval"][0] <= value[k] <= f["interval"][1]))
            inside, n_iv = inside + r["inside"], n_iv + 1
        fc[k] = r
    out = {"readings_of_the_two_files_differ_by_at_most": worst, "pairs_with_no_change_to_account_for": left_out, "hypotheses": h,
           "counts": {k: sum(v["reading"] == k for kk, v in h.items() if kk != "HP_planted_suite") for k in ("confirmed", "refuted", "unresolved")},
           "forecasts": {"by_hypothesis": fc, "intervals_that_held_the_value": [inside, n_iv], "outcomes_given_the_higher_probability_that_happened": [right, n_out]}}
    (R / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print("readings differ by at most", worst, "| left out:", left_out)
    for k, v in h.items():
        print(f"{k:36s} {v['reading']:10s} " + json.dumps({a: b for a, b in v.items() if not isinstance(b, dict) and a != "reading"}))
    print(out["counts"], out["forecasts"]["intervals_that_held_the_value"], out["forecasts"]["outcomes_given_the_higher_probability_that_happened"])


if __name__ == "__main__":
    main()
