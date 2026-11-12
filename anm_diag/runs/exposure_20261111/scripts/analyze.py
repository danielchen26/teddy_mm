"""Read the six records and apply the rule of DESIGN.md. Writes results/exposure.json."""
import json

from _paths import ARMS, RES, RUN, TARGETS

out = {"runs": {}}
for g in TARGETS:
    for arm in ARMS:
        r = json.load(open(RUN / f"{arm}_{g}" / "anm_record.json"))
        ro = r["readouts"]["PDL1"]
        out["runs"][f"{arm}_{g}"] = {
            "arm": arm, "target": g, "n_cells": r["n_cells"], "n_treated": r["n_treated"], "n_guides": r["n_treated_clusters"],
            "fold_seed": r.get("fold_seed"), "delta": ro["delta"], "delta_ci": ro["delta_ci"], "runtime_s": r["runtime_s"],
            "problems": [p["code"] + ": " + p["message"] for p in r["problems"]],
            "stages": {s: {"u": v["u"], "ci": v["ci"], "verdict": v["verdict"], "best_class": v["best_class"]}
                       for s, v in ro["states"].items() if v.get("role") == "chain"}}

head_out = lambda arm, g="IRF1": out["runs"][f"{arm}_{g}"]["stages"]["HEAD_OUT"]
u0 = head_out("as_trained")
mark = 0.5 * u0["u"]
meets = lambda arm: head_out(arm)["verdict"] != "reject" and head_out(arm)["u"] <= mark
still = lambda arm: head_out(arm)["verdict"] == "reject" and head_out(arm)["u"] > mark
if u0["verdict"] != "reject":
    outcome = "cannot be put on these cells"
elif meets("without_IRF1"):
    outcome = "not supported"
elif meets("with_IRF1") and still("without_IRF1"):
    outcome = "supported"
elif still("with_IRF1"):
    outcome = "not answered"
else:
    outcome = "between the marks"
out["rule"] = {"u0": u0["u"], "half_of_u0": mark, "meets_mark": {a: bool(meets(a)) for a in ARMS[1:]},
               "still_rejected_above_mark": {a: bool(still(a)) for a in ARMS[1:]}, "outcome": outcome}
json.dump(out, open(RES / "exposure.json", "w"), indent=1)
for g in TARGETS:
    for arm in ARMS:
        x = out["runs"][f"{arm}_{g}"]
        print(f"{g:5s} {arm:13s} n {x['n_treated']:4d}/{x['n_cells']:4d} delta {x['delta']:+.2f} | " + " | ".join(
            f"{s} {v['u']:.2f} [{v['ci'][0]:.2f}, {v['ci'][1]:.2f}] {v['verdict']}" for s, v in x["stages"].items()))
print("rule:", out["rule"])
