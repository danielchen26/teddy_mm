"""Read the records and apply the rule of DESIGN.md. Writes results/loop.json."""
import json

from _paths import ARMS, RES, RUN, TARGETS, audits

MARGIN = 0.10


def reading(ro):
    """One readout of one audit, read the way DESIGN.md says."""
    st = {s: v for s, v in ro["states"].items() if v.get("role") == "chain"}
    rc = ro["root_cause"]
    added = st["HEAD_OUT"]["u"] - st["Z"]["u"]
    named = rc["cause"] in ("loss_at", "relative_loss_at")
    if rc["cause"] in ("no_effect", "control_failed"):
        status = "no reading"
    elif named:
        status = "located in the head"
    elif added <= MARGIN:
        status = "no loss in the head"
    else:
        status = "not resolved"
    return {"delta": ro["delta"], "delta_ci": ro["delta_ci"], "stages": {s: {"u": v["u"], "ci": v["ci"], "verdict": v["verdict"], "best_class": v["best_class"]} for s, v in st.items()},
            "root_cause": {k: rc[k] for k in ("cause", "from", "to") if k in rc}, "added_by_head": added, "status": status,
            "below_floor": bool(st["HEAD_OUT"]["u"] < st["Z"]["u"] - MARGIN)}


out = {"margin": MARGIN, "audits": {}, "missing": []}
for name, (kind, ds, arm, g) in audits().items():
    p = RUN / name / "anm_record.json"
    if not p.exists():
        out["missing"].append(name); continue
    r = json.load(open(p))
    out["audits"][name] = {"kind": kind, "data": ds, "arm": arm, "target": g, "n_cells": r["n_cells"], "n_treated": r["n_treated"], "n_clusters": r["n_treated_clusters"],
                           "fold_seed": r.get("fold_seed"), "runtime_s": r["runtime_s"], "problems": sorted({p["code"] for p in r["problems"]}),
                           "readouts": {rid: reading(ro) for rid, ro in r["readouts"].items()}}
A = out["audits"]
get = lambda name, rid="PDL1": A[name]["readouts"][rid] if name in A else None
out["rule"] = {}
for ds, genes in TARGETS.items():
    b, m, d = (get(f"pooled_{ds}_{a}") for a in ARMS)
    if not (b and m and d):
        continue
    fixed = lambda x: x["status"] == "no loss in the head"
    if b["status"] != "located in the head":
        outcome = "cannot be put on these cells"
    elif fixed(m) and fixed(d):
        outcome = "closed, but not specific to the located step"
    elif fixed(m):
        outcome = "closed at the located step"
    elif fixed(d):
        outcome = "closed by the other part: the location was wrong"
    else:
        outcome = "not closed"
    per = {a: {s: sum(1 for g in genes if get(f"target_{ds}_{g}_{a}") and get(f"target_{ds}_{g}_{a}")["status"] == s) for s in ("located in the head", "no loss in the head", "not resolved", "no reading")} for a in ARMS}
    out["rule"][ds] = {"located_step": b["root_cause"], "status": {a: x["status"] for a, x in zip(ARMS, (b, m, d))}, "added_by_head": {a: x["added_by_head"] for a, x in zip(ARMS, (b, m, d))},
                       "outcome": outcome, "single_head_audits": per,
                       "other_readouts": {rid: {a: A[f"pooled_{ds}_{a}"]["readouts"][rid]["status"] for a in ARMS} for rid in A[f"pooled_{ds}_as_trained"]["readouts"] if rid != "PDL1"},
                       "transfer": {a: (get(f"transfer_{ds}_{a}") or {}).get("status") for a in ARMS[1:]}}
out["below_floor"] = [f"{n}/{rid}" for n, x in A.items() for rid, r in x["readouts"].items() if r["below_floor"]]
json.dump(out, open(RES / "loop.json", "w"), indent=1)
for name, x in A.items():
    for rid, r in x["readouts"].items():
        print(f"{name:46s} {rid:5s} d {r['delta']:+.2f} | " + " ".join(f"{s} {v['u']:.2f} {v['verdict'][:3]}" for s, v in r["stages"].items()) + f" | head adds {r['added_by_head']:+.2f} | {r['root_cause']['cause']} | {r['status']}")
print(json.dumps(out["rule"], indent=1)); print("below floor:", out["below_floor"], "| missing:", len(out["missing"]))
