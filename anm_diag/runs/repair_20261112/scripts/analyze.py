"""Read the records and results/bone_marrow.json and apply the rule of DESIGN.md. Writes results/repair.json, which also
carries the earlier rounds (runs/loop_20261112, runs/loop_guard_20261112) for the same cells, so the whole loop is in one file."""
import json

import numpy as np

from _paths import ARM, GENES, GUARD, LOOP, RES, RUN, SEEDS, TARGETS, audits

MARGIN, GUARD_FALL = 0.10, 0.02


def reading(ro):
    """One readout of one audit, read as in the earlier rounds."""
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


r1 = json.load(open(LOOP / "results" / "loop.json"))["audits"]; r2 = json.load(open(GUARD / "results" / "loop_guard.json"))["audits"]
bm = json.load(open(RES / "bone_marrow.json")); base = bm["as_trained"]["mean_pearson"]
out = {"margin": MARGIN, "guard": GUARD_FALL, "audits": {}, "missing": [], "rule": {}, "rounds": {}, "replication": {}}
for name, (kind, ds, g, seed) in audits().items():
    p = RUN / name / "anm_record.json"
    if not p.exists():
        out["missing"].append(name); continue
    r = json.load(open(p))
    out["audits"][name] = {"kind": kind, "data": ds, "arm": ARM, "target": g, "seed": seed, "n_cells": r["n_cells"], "n_treated": r["n_treated"], "n_clusters": r["n_treated_clusters"],
                           "fold_seed": r.get("fold_seed"), "runtime_s": r["runtime_s"], "problems": sorted({p["code"] for p in r["problems"]}),
                           "readouts": {rid: reading(ro) for rid, ro in r["readouts"].items()}}
A = out["audits"]
fall = {s: (lambda v: {"heads": v, "mean_fall": base - float(np.mean(v)), "largest_fall": base - min(v)} if v else None)([bm["heads"][f"s{s}__{g}"]["mean_pearson"] for g in GENES if f"s{s}__{g}" in bm["heads"]]) for s in SEEDS}
out["bone_marrow"] = {"as_trained": base, "by_seed": fall, "final_head": bm["heads"].get("s0__all")}


def outcome(ds, s):
    g = A[f"pooled_{ds}_{ARM}_s{s}"]["readouts"]
    kept = g["PDL1"]["status"] == "no loss in the head" and fall[s]["mean_fall"] <= GUARD_FALL
    located = sorted(rid for rid, r in g.items() if r["status"] == "located in the head")
    return {"pdl1": g["PDL1"]["status"], "guard_kept": bool(fall[s]["mean_fall"] <= GUARD_FALL), "located_in_the_head": located,
            "outcome": "regressed" if not kept else ("kept; nothing located in the head" if not located else "kept; still located in the head: " + ", ".join(located))}


for ds, genes in TARGETS.items():
    if f"pooled_{ds}_{ARM}_s0" not in A:
        continue
    out["rule"][ds] = outcome(ds, 0)
    single = [A[f"target_{ds}_{t}_{ARM}_s0"]["readouts"]["PDL1"] for t in genes if f"target_{ds}_{t}_{ARM}_s0" in A]
    out["rule"][ds]["single_head_audits"] = {s: [x["status"] for x in single].count(s) for s in ("located in the head", "no loss in the head", "not resolved", "no reading")}
    out["rule"][ds]["left_outside_the_head"] = [{"target": t, "u_Z": x["stages"]["Z"]["u"], "verdict_Z": x["stages"]["Z"]["verdict"], "u_output": x["stages"]["HEAD_OUT"]["u"], "status": x["status"]}
                                                for t, x in zip(genes, single) if x["status"] != "no loss in the head"]
    out["replication"][ds] = {f"s{s}": outcome(ds, s) for s in SEEDS[1:] if f"pooled_{ds}_{ARM}_s{s}" in A and fall[s]}
    rows = {}
    for rid in A[f"pooled_{ds}_{ARM}_s0"]["readouts"]:
        rows[rid] = {"round 0, as trained": r1[f"pooled_{ds}_as_trained"]["readouts"][rid], "round 1, MLP repaired on PD-L1": r1[f"pooled_{ds}_mlp_only"]["readouts"][rid],
                     "round 2, with the bone-marrow outputs held": r2[f"pooled_{ds}_mlp_guarded"]["readouts"][rid], "round 3, one head for both cell types and all measured proteins": A[f"pooled_{ds}_{ARM}_s0"]["readouts"][rid]}
    out["rounds"][ds] = rows
out["one_head_for_both"] = bool(out["rule"]) and all(v["outcome"] != "regressed" for v in out["rule"].values())
out["below_floor"] = [f"{n}/{rid}" for n, x in A.items() for rid, r in x["readouts"].items() if r["below_floor"]]
json.dump(out, open(RES / "repair.json", "w"), indent=1)
for name, x in A.items():
    for rid, r in x["readouts"].items():
        print(f"{name:40s} {rid:5s} d {r['delta']:+.2f} | " + " ".join(f"{s} {v['u']:.2f} {v['verdict'][:3]}" for s, v in r["stages"].items()) + f" | head adds {r['added_by_head']:+.2f} | {r['root_cause']['cause']} | {r['status']}")
print(json.dumps({"rule": out["rule"], "replication": out["replication"], "bone_marrow": {k: (v and {"mean_fall": round(v["mean_fall"], 4), "largest_fall": round(v["largest_fall"], 4)}) for k, v in fall.items()}, "one_head_for_both": out["one_head_for_both"]}, indent=1))
print("below floor:", out["below_floor"], "| missing:", len(out["missing"]))
