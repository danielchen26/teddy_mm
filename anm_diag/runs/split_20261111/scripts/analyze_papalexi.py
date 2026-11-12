"""THP-1 (Papalexi 2021): does the pooled verdict hold per target and per batch? Sealed pooled record vs one record per
target (same controls, all batches) and one per batch (all five targets, that batch's own controls)."""
import json, os
from _paths import FEATURES, RES, RUNS, SPLIT
from anm_diag import events as ev
O = str(SPLIT)
G = ["IFNGR1", "IFNGR2", "JAK2", "STAT1", "IRF1"]; B = ["rep1", "rep3", "rep4"]; RO = ["PDL1", "PDL2"]
sealed = json.load(open(RUNS / "papalexi2021_headfix" / "anm_record.json"))
load = lambda n: json.load(open(f"{O}/{n}/anm_record.json")) if os.path.exists(f"{O}/{n}/anm_record.json") else None
per = {g: r for g in G if (r := load(f"papalexi2021_target_{g}"))}; bat = {b: r for b in B if (r := load(f"papalexi2021_batch_{b}"))}
q = load("papalexi2021_pooled_quick")
meta = lambda d: {k: dict(n_cells=r["n_cells"], n_treated=r["n_treated"], n_clusters=r["n_treated_clusters"], runtime_s=round(r["runtime_s"])) for k, r in d.items()}
res = {"targets_done": list(per), "batches_done": list(bat), "n_targets": meta(per), "n_batches": meta(bat), "repro": None, "targets": {}, "batches": {},
       "pooled": dict(n_cells=sealed["n_cells"], n_treated=sealed["n_treated"]),
       "admission": json.load(open(RES / "papalexi_admission.json"))}
if q:
    rows = {rid: {s: dict(here=q["readouts"][rid]["states"][s]["u"], sealed=sealed["readouts"][rid]["states"][s]["u"], same=q["readouts"][rid]["states"][s]["verdict"] == sealed["readouts"][rid]["states"][s]["verdict"]) for s in ("Z", "HEAD_LAT", "HEAD_OUT")} for rid in RO}
    res["repro"] = dict(rows=rows, max_abs_diff=max(abs(v["here"] - v["sealed"]) for r in rows.values() for v in r.values()), verdicts_equal=all(v["same"] for r in rows.values() for v in r.values()))
for rid in RO:
    if per: res["targets"][rid] = ev.pooling_evidence(sealed, per, rid)
    if bat: res["batches"][rid] = ev.pooling_evidence(sealed, bat, rid)
A4 = ("RAW", "Z", "HEAD_LAT", "HEAD_OUT")
res["auc"] = {"pool": {k: sealed["forward"]["auc"][k] for k in A4}, **{g: {k: r["forward"]["auc"][k] for k in A4} for g, r in per.items()}, **{b: {k: r["forward"]["auc"][k] for k in A4} for b, r in bat.items()}}
res["notes"] = dict(batch_rep4_fold_seed=bat["rep4"].get("fold_seed") if "rep4" in bat else None, batch_rep4_reason="with the declared seed 0 the fold assignment puts all 5 treated target clusters into one fold; the runner falls back to the next seed for the folds", probes_dropped=["SET"], batch_not_run=["rep2 (73 cells)"])
json.dump(res, open(RES / "papalexi_split.json", "w"), indent=1)
f = lambda v: f"{v:+.3f}"
print("targets:", res["targets_done"], "| batches:", res["batches_done"], "| repro:", res["repro"] and (round(res["repro"]["max_abs_diff"], 5), res["repro"]["verdicts_equal"]))
for kind in ("targets", "batches"):
    for rid, e in res[kind].items():
        print(f"== {kind} {rid}: pooled {f(e['pooled_delta'])} | label {e['label']} | no response {e['no_response']} | opposite {e['opposite_direction']} | conflicts {e['verdict_conflicts']} | delta het p={e['delta_heterogeneity']['p']:.3g}")
        print("   delta:", {s: (f(d["v"]), [round(c, 2) for c in d["ci"]], d["n_treated"]) for s, d in e["delta"].items()})
        for s, v in e["stages"].items():
            print(f"   {s:9s} pooled {f(v['pooled']['u'])} {v['pooled']['verdict'][:3]} |", {g: (f(d["u"]), d["verdict"][:3]) for g, d in v["sites"].items()}, "| het p=%.3g" % v["heterogeneity"]["p"])
