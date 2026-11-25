#!/usr/bin/env python3
"""The hypotheses of DESIGN.md, Section 7, against the reading. results/outcome.json.

  python score.py
"""
import json

import numpy as np

from _b2 import RUN

OWN_PAIRS = 8
FLAGGED = [("mel_co", "STAT3", "MHC1"), ("mel_ctl", "SSR2", "MHC1"), ("mel_ctl", "CD44", "MHC1")]
rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None


def main():
    R = json.loads((RUN / "results" / "reading.json").read_text())
    D = json.loads((RUN / "results" / "declared.json").read_text())
    pairs = [dict(p, stratum=g) for g, s in R["strata"].items() for p in s["pairs"]]
    assert len(pairs) == 79
    err = lambda ps: {"no_change": rms([p["change"] for p in ps]), "control_readout": rms([p["change"] - p["read"] for p in ps]), "forecast": rms([p["change"] - p["forecast"] for p in ps])}
    out = {"pairs": len(pairs), "inside": sum(p["inside"] for p in pairs), "error": err(pairs),
           "by_stratum": {g: {"pairs": len(s["pairs"]), "inside": sum(p["inside"] for p in s["pairs"]), "error": err(s["pairs"]), "pairs_with_a_change": sum(p["snr"] >= R["min_snr"] for p in s["pairs"])}
                          for g, s in R["strata"].items()},
           "pairs_with_a_change": {"pairs": sum(p["snr"] >= R["min_snr"] for p in pairs), "inside": sum(p["inside"] for p in pairs if p["snr"] >= R["min_snr"]),
                                   "error": err([p for p in pairs if p["snr"] >= R["min_snr"]])},
           "hypotheses": {}}
    H = out["hypotheses"]
    n = out["inside"]
    H["HB1_the_intervals_hold"] = {"inside": n, "of": 79, "outcome": "confirmed" if n >= 72 else "refuted" if n <= 67 else "unresolved"}
    e = out["error"]
    better = (e["forecast"] < e["no_change"]) + (e["forecast"] < e["control_readout"])
    H["HB2_the_forecast_uses_what_it_adds"] = dict(e, outcome="confirmed" if better == 2 else "refuted" if better == 0 else "unresolved")
    own = [p for p in pairs if p["edge"] == "own"]
    assert len(own) == OWN_PAIRS
    eo = {"with_the_edge": rms([p["change"] - p["forecast"] for p in own]), "without_it": rms([p["change"] - p["without_its_edge"] for p in own])}
    ins = sum(p["inside"] for p in own)
    H["HB3_an_edge_in_another_condition"] = {"inside": ins, "of": OWN_PAIRS, "error": eo, "outcome": "confirmed" if ins >= 6 and eo["with_the_edge"] < eo["without_it"] else "refuted",
                                             "pairs": [{k: p[k] for k in ("stratum", "event", "readout", "change", "se", "forecast", "interval", "inside", "z", "by_its_edge", "read")} for p in own]}
    rows = []
    for g, s in R["strata"].items():
        for key, a in s["accounts"].items():
            if not (a["own_gene"] and a["of_the_batch"]):
                continue
            z, ch = a["stages"]["Z"]["shares"], a["chain"]
            largest = max(ch["parts"], key=lambda k: ch["parts"][k]["share"])
            ok = z["read"]["ci"][1] < 0.5 and largest == "not_in_the_input" and z["others"]["ci"][1] < 0.5
            rows.append({"stratum": g, "pair": key, "read_at_the_state": z["read"], "others_at_the_state": z["others"], "largest_part": largest,
                         "parts": {k: v["share"] for k, v in ch["parts"].items()}, "says": ch["says"], "as_forecast": bool(ok)})
    k = sum(r["as_forecast"] for r in rows)
    H["HB4_an_own_gene_change_is_outside_the_input"] = {"pairs_with_a_change": len(rows), "as_forecast": k, "pairs": rows,
                                                         "outcome": "not tested: fewer than four such pairs" if len(rows) < 4 else "confirmed" if k >= 0.75 * len(rows) else "refuted" if k < 0.5 * len(rows) else "unresolved"}
    fl = [next(p for p in pairs if (p["stratum"], p["event"], p["readout"]) == f) for f in FLAGGED]
    k = sum(not p["inside"] for p in fl)
    H["HB5_a_flag_in_another_condition"] = {"outside": k, "of": 3, "outcome": "confirmed" if k >= 2 else "refuted" if k == 0 else "unresolved",
                                            "pairs": [{q: p[q] for q in ("stratum", "event", "readout", "change", "se", "forecast", "interval", "inside", "z")} for p in fl]}
    rows = []
    for g in ("mel_co", "thp1"):
        for key, a in R["strata"].get(g, {}).get("accounts", {}).items():
            if a["own_gene"]:
                continue
            tz, vz = a["gains"]["T over Z"]["carried"], a["gains"]["V over Z"]["carried"]
            rows.append({"stratum": g, "pair": key, "of_the_batch": a["of_the_batch"], "carried": {k: a["stages"][k]["shares"]["carried"]["share"] for k in ("V", "Z", "T", "R")},
                         "T_over_Z": tz, "T_over_R": a["gains"]["T over R"]["carried"], "V_over_Z": vz, "lost_before_the_state": a["chain"]["parts"]["lost_before_the_state"],
                         "the_typed_state_carries_less": bool(tz["ci"][1] < -0.10), "the_state_loses_a_step": bool(vz["ci"][0] > 0.10),
                         "brought_back": (tz["gain"] / vz["gain"]) if vz["ci"][0] > 0.10 else None})
    lost = [r for r in rows if r["the_state_loses_a_step"]]
    a_ok = not any(r["the_typed_state_carries_less"] for r in rows)
    b = "not tested: no pair loses a step before the state" if not lost else "holds" if all(r["brought_back"] >= 0.5 for r in lost) else "does not hold"
    H["HB6_the_typed_state_where_the_mean_state_loses"] = {"pairs": len(rows), "a_no_pair_where_the_typed_state_carries_less": a_ok, "pairs_where_the_state_loses_a_step": len(lost), "b_brings_back_at_least_half": b,
                                                           "outcome": ("confirmed" if a_ok and b == "holds" else "refuted" if (not a_ok) or b == "does not hold" else "(a) holds; (b) not tested"), "rows": rows}
    out["forecasts_of_the_hypotheses"] = D["forecasts"]
    (RUN / "results" / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    print("inside", out["inside"], "of 79 |", out["error"], "| by stratum", {g: (v["inside"], v["pairs"]) for g, v in out["by_stratum"].items()})
    for k, v in H.items():
        print(k, "->", v["outcome"])


if __name__ == "__main__":
    main()
