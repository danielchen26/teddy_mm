"""Run the event checks (closure, admission, join of the forward and the inverse half) on every finished record under runs/."""
import json, os
from _paths import RES, ROOT
root, out = str(ROOT), str(RES / "join_closure_admission.json")
from anm_diag import events as ev
from anm_diag.spec import load_declaration
fwd = json.load(open(os.path.join(root, "runs/forward_typed_events/forward_record.json")))
res = {"datasets": {}, "ks2": {}}
for ds in ("frangieh2021", "papalexi2021"):
    decl = load_declaration(os.path.join(root, f"declarations/teddy_mm_{ds}_headfix.sealed.yaml"))
    inv = json.load(open(os.path.join(root, f"runs/{ds}_headfix/anm_record.json")))
    spl = json.load(open(os.path.join(root, f"runs/{ds}_headsplit/anm_record.json")))
    sites = decl["source"]["treated"]; evs = ev.source_events(decl, inv) + ev.forward_events(fwd, ds)
    r = {"n_events": len(evs), "closure": ev.closure(evs)[f"{decl['name']}:A=1"], "problems": ev.check_source(decl), "sites": sites,
         "n_treated": inv["n_treated"], "n_cells": inv["n_cells"], "admission": ev.admission(fwd, ds, sites), "readouts": {}}
    for rid in inv["readouts"]:
        j = ev.join(fwd, inv, ds, rid)
        if j is None:
            r["readouts"][rid] = None; continue
        rows = ev.site_responses(fwd, ds, sites, j["forward_readout"])
        count = [x for x in rows if x["channel"] == "count" and x.get("n")]
        j["sites"] = rows; j["distinguishable_pairs"] = ev.distinguishable(rows)
        j["n_resolved"] = sum(1 for x in rows if x.get("resolved")); j["n_applied"] = sum(1 for x in rows if x.get("n"))
        j["insilico_ko_max"] = max((abs(2 * x["slope"]) for x in count), default=None)
        j["insilico_share_of_measured"] = j["insilico_ko_max"] / abs(j["measured_delta"]) if count else None
        js = ev.join(fwd, spl, ds, rid)
        j["headsplit_inverse"] = js["inverse"] if js else None
        r["readouts"][rid] = j
    res["datasets"][ds] = r
    ev.write_jsonl(evs, os.path.join(os.path.dirname(out), f"events_{ds}.jsonl"))
for k in ("ks2_teddy_mitoribo", "ks2_teddy_ribo"):
    d = load_declaration(os.path.join(root, f"declarations/{k}.yaml")); d.setdefault("name", k)
    res["ks2"][k] = {"closure": ev.closure(ev.source_events(d))[f"{d['name']}:A=1"], "problems": [p["message"][:120] for p in ev.check_source(d)]}
json.dump(res, open(out, "w"), indent=1, ensure_ascii=False)
f3 = lambda x: None if x is None else round(x, 3)
for ds, r in res["datasets"].items():
    print(f"== {ds}: {r['n_events']} events; closure {r['closure']}; token share " + ", ".join(f"{s} {round(100*r['admission'][s]['embed'])}%" for s in r["sites"]) + "; detected " + ", ".join(f"{s} {round(100*r['admission'][s]['count'])}%" for s in r["sites"]))
    for rid, j in r["readouts"].items():
        if j is None: print(f"  {rid}: not in the forward run"); continue
        print(f"  {rid}: measured {f3(j['measured_delta'])} | KO-dir {f3(j['head_response'])} {[f3(v) for v in j['head_response_ci']]} share {f3(j['share_of_measured'])} same_sign {j['same_sign']} step {f3(j['step_difference'])} pct {j['percentile_among_random']} rnd|med| {f3(j['random_abs_median'])} max {f3(j['random_abs_max'])}")
        print(f"       in-silico max {f3(j['insilico_ko_max'])} share {f3(j['insilico_share_of_measured'])} | resolved {j['n_resolved']}/{j['n_applied']} | disjoint pairs {j['distinguishable_pairs']}")
        print("       inverse headfix:", {s: (f3(v['u']), v['verdict'][:3]) for s, v in j["inverse"].items()})
        if rid == "PDL1":
            print("       inverse headsplit:", {s: (f3(v['u']), v['verdict'][:3]) for s, v in j["headsplit_inverse"].items()})
            print("       state profile KO/random:", {s: (f3(v['ko']), f3(v['random_median'])) for s, v in j["state_profile"].items()})
print("ks2:", res["ks2"])
