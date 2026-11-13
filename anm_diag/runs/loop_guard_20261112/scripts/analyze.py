"""Read the records and results/bone_marrow.json and apply the rule of DESIGN.md. Writes results/loop_guard.json.
The head as trained and the unguarded repair on the same cells are read from runs/loop_20261112/results/loop.json."""
import json

import numpy as np

from _paths import ARM, LOOP, RES, RUN, TARGETS, audits

MARGIN, GUARD = 0.10, 0.02


def reading(ro):
    """One readout of one audit, read as in runs/loop_20261112."""
    st = {s: v for s, v in ro["states"].items() if v.get("role") == "chain"}
    rc = ro["root_cause"]
    added = st["HEAD_OUT"]["u"] - st["Z"]["u"]
    if rc["cause"] in ("no_effect", "control_failed"):
        status = "no reading"
    elif rc["cause"] in ("loss_at", "relative_loss_at"):
        status = "located in the head"
    elif added <= MARGIN:
        status = "no loss in the head"
    else:
        status = "not resolved"
    return {"delta": ro["delta"], "delta_ci": ro["delta_ci"], "stages": {s: {"u": v["u"], "ci": v["ci"], "verdict": v["verdict"], "best_class": v["best_class"]} for s, v in st.items()},
            "root_cause": {k: rc[k] for k in ("cause", "from", "to") if k in rc}, "added_by_head": added, "status": status,
            "below_floor": bool(st["HEAD_OUT"]["u"] < st["Z"]["u"] - MARGIN)}


first = json.load(open(LOOP / "results" / "loop.json"))["audits"]
bm = json.load(open(RES / "bone_marrow.json"))
out = {"margin": MARGIN, "guard": GUARD, "audits": {}, "missing": [], "rule": {}}
for name, (kind, ds, g) in audits().items():
    p = RUN / name / "anm_record.json"
    if not p.exists():
        out["missing"].append(name); continue
    r = json.load(open(p))
    out["audits"][name] = {"kind": kind, "data": ds, "arm": ARM, "target": g, "n_cells": r["n_cells"], "n_treated": r["n_treated"], "n_clusters": r["n_treated_clusters"],
                           "fold_seed": r.get("fold_seed"), "runtime_s": r["runtime_s"], "problems": sorted({p["code"] for p in r["problems"]}),
                           "readouts": {rid: reading(ro) for rid, ro in r["readouts"].items()}}
A = out["audits"]
for ds, genes in TARGETS.items():
    if f"pooled_{ds}_{ARM}" not in A:
        continue
    g = A[f"pooled_{ds}_{ARM}"]["readouts"]; b = first[f"pooled_{ds}_as_trained"]["readouts"]["PDL1"]
    scores = [bm["heads"][f"{ds}|{t}"]["mean_pearson"] for t in genes]
    fall = bm["as_trained"]["mean_pearson"] - float(np.mean(scores))
    closed, kept = g["PDL1"]["status"] == "no loss in the head", fall <= GUARD
    if b["status"] != "located in the head":
        outcome = "cannot be put on these cells"
    else:
        outcome = {(True, True): "closed within the guard", (True, False): "closed, guard broken", (False, True): "guard kept, not closed", (False, False): "neither"}[closed, kept]
    single = [A[f"target_{ds}_{t}_{ARM}"]["readouts"]["PDL1"]["status"] for t in genes if f"target_{ds}_{t}_{ARM}" in A]
    out["rule"][ds] = {"status": {"as_trained": b["status"], "mlp_only": first[f"pooled_{ds}_mlp_only"]["readouts"]["PDL1"]["status"], ARM: g["PDL1"]["status"]},
                       "added_by_head": {"as_trained": b["added_by_head"], "mlp_only": first[f"pooled_{ds}_mlp_only"]["readouts"]["PDL1"]["added_by_head"], ARM: g["PDL1"]["added_by_head"]},
                       "bone_marrow": {"as_trained": bm["as_trained"]["mean_pearson"], "guarded_heads": scores, "mean_fall": fall, "largest_fall": bm["as_trained"]["mean_pearson"] - min(scores)},
                       "outcome": outcome, "single_head_audits": {s: single.count(s) for s in ("located in the head", "no loss in the head", "not resolved", "no reading")},
                       "other_readouts": {rid: r["status"] for rid, r in g.items() if rid != "PDL1"}}
out["below_floor"] = [f"{n}/{rid}" for n, x in A.items() for rid, r in x["readouts"].items() if r["below_floor"]]
json.dump(out, open(RES / "loop_guard.json", "w"), indent=1)
for name, x in A.items():
    for rid, r in x["readouts"].items():
        print(f"{name:44s} {rid:5s} d {r['delta']:+.2f} | " + " ".join(f"{s} {v['u']:.2f} {v['verdict'][:3]}" for s, v in r["stages"].items()) + f" | head adds {r['added_by_head']:+.2f} | {r['root_cause']['cause']} | {r['status']}")
print(json.dumps(out["rule"], indent=1)); print("below floor:", out["below_floor"], "| missing:", len(out["missing"]))
