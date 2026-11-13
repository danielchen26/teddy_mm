"""Melanoma (Frangieh 2021): pooling evidence on the measured side. The sealed pooled record against one record per
target (same controls), the quick pooled reruns as a check of the environment, and the admission trail of the real knockouts."""
import json, os
from _paths import RES, RUNS, SPLIT
from anm_diag import events as ev
out, sealed_p = str(SPLIT), str(RUNS / "frangieh2021_headfix" / "anm_record.json")
G = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]; RO = ["PDL1", "MHC1", "HLAE"]
sealed = json.load(open(sealed_p))
per = {g: json.load(open(f"{out}/frangieh2021_target_{g}/anm_record.json")) for g in G if os.path.exists(f"{out}/frangieh2021_target_{g}/anm_record.json")}
quick = {k: json.load(open(f"{out}/{k}/anm_record.json")) for k in ("frangieh2021_pooled_quick", "frangieh2021_pooled_quick_byguide") if os.path.exists(f"{out}/{k}/anm_record.json")}
res = {"targets_done": list(per), "n": {g: dict(n_cells=r["n_cells"], n_treated=r["n_treated"], n_guide_clusters=r["n_treated_clusters"], runtime_s=round(r["runtime_s"])) for g, r in per.items()},
       "repro": {}, "readouts": {}}
for k, q in quick.items():
    rows = {rid: {s: dict(here=q["readouts"][rid]["states"][s]["u"], sealed=sealed["readouts"][rid]["states"][s]["u"], verdict_here=q["readouts"][rid]["states"][s]["verdict"], verdict_sealed=sealed["readouts"][rid]["states"][s]["verdict"])
                  for s in ("Z", "HEAD_LAT", "HEAD_OUT")} for rid in RO}
    res["repro"][k] = dict(rows=rows, max_abs_diff=max(abs(v["here"] - v["sealed"]) for r in rows.values() for v in r.values()), verdicts_equal=all(v["verdict_here"] == v["verdict_sealed"] for r in rows.values() for v in r.values()))
for rid in RO:
    # a target whose own response is not resolved gets no stage verdicts from the runner ("no_effect"): keep it for the delta only
    full = {g: r for g, r in per.items() if r["readouts"][rid]["states"]}
    e = ev.pooling_evidence(sealed, per, rid)
    res["readouts"][rid] = dict(delta={g: dict(v=r["readouts"][rid]["delta"], ci=r["readouts"][rid]["delta_ci"], resolved=bool(r["readouts"][rid]["delta_ci"][1] < 0 or r["readouts"][rid]["delta_ci"][0] > 0)) for g, r in per.items()},
                                with_stage_verdicts=list(full), evidence=e)
res["admission_real"] = json.load(open(RES / "frangieh_admission.json"))
json.dump(res, open(RES / "frangieh_split.json", "w"), indent=1, ensure_ascii=False)
f = lambda v: f"{v:+.3f}"
print("targets:", res["targets_done"], "| repro:", {k: (round(v["max_abs_diff"], 4), v["verdicts_equal"]) for k, v in res["repro"].items()})
for rid, x in res["readouts"].items():
    print(f"== {rid} pooled {f(sealed['readouts'][rid]['delta'])} | per target:", {g: (f(d["v"]), [round(c, 2) for c in d["ci"]], "res" if d["resolved"] else "n.s.") for g, d in x["delta"].items()})
    e = x["evidence"]
    if e:
        print("   label:", e["label"], "| no response:", e["no_response"], "| conflicts:", e["verdict_conflicts"], "| opposite:", e["opposite_direction"], "| delta Q p=%.3g" % e["delta_heterogeneity"]["p"])
        for s, v in e["stages"].items():
            print(f"   {s:9s} pooled {f(v['pooled']['u'])} {v['pooled']['verdict'][:3]} |", {g: (f(d["u"]), d["verdict"][:3]) for g, d in v["sites"].items()}, "| het p=%.3g" % v["heterogeneity"]["p"])
