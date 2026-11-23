#!/usr/bin/env python3
"""The forecasts of DESIGN.md against what was read: results/outcome.json and a table.

  python score.py

An interval forecast is scored inside or outside; an outcome forecast by whether the class given the highest probability
happened. A forecast whose condition did not occur (the event was not admitted, A3 does not hold) is listed as not scored.
"""
import json
from pathlib import Path

RUN = Path(__file__).resolve().parent.parent
D = json.loads((RUN / "results" / "declared.json").read_text())["forecasts"]
SEL = json.loads((RUN / "results" / "selection.json").read_text())
R = json.loads((RUN / "results" / "round.json").read_text())
MAN = json.loads((RUN / "results" / "states_manifest.json").read_text())
EV = {"BRD4": ("PDL1", "thp1"), "CMTM6": ("PDL1", "thp1"), "TAPBP": ("MHC1", "mel_co"), "CD44": ("MHC1", "mel_co"), "B2M": ("MHC1", "mel_ctl")}
CAL = {"POU2F2": "thp1", "CDH19": "mel_co", "KLF4": "mel_ctl"}
adm = lambda g, r: SEL["strata"][g]["readouts"].get(r, {}).get("admitted", [])
ev = lambda s: R["targets"].get(f"{EV[s][0]} {EV[s][1]}", {}).get("events", {}).get(s)
BEFORE = ("open before the model", "undecided before the model")
CLOSED = ("closed (absolute)", "the model's part closed (relative)")


def observed(i):
    kind, _, s = i.partition("-")
    if i == "C-STATES":
        return MAN["check_largest_difference"]
    if i == "A-CAL":
        n = sum(u in adm(g, r) for u, g in CAL.items() for r in SEL["strata"][g]["readouts"])
        return "none" if n == 0 else "one or more"
    if i == "A-PDL1-MEL":
        return "none" if not (adm("mel_co", "PDL1") or adm("mel_ctl", "PDL1")) else "one or more"
    if i == "A-HLAE":
        return "none" if not (adm("mel_co", "HLAE") or adm("mel_ctl", "HLAE")) else "one or more"
    if kind == "A":
        return "yes" if s in adm(EV[s][1], EV[s][0]) else "no"
    if kind == "S":
        if s in ("MHC1-co", "MHC1-ctl"):
            st = R["statuses"].get("MHC1 mel_co" if s == "MHC1-co" else "MHC1 mel_ctl")
            return None if st is None else (st if st in BEFORE else ("the model's part closed (relative)" if st == CLOSED[1] else "another"))
        name = s.split("-")[1]
        e = ev(name)
        if e is None:
            return "not admitted" if name == "BRD4" else None
        st = R["targets"][f"PDL1 thp1"]["status_by_event"][name]
        return "open or undecided before the model" if st in BEFORE else ("the model's part closed" if st in CLOSED else "not closed, not tested or another")
    e = ev(s.split("-")[0])
    if e is None:
        return None
    if kind == "T":
        return e["as_trained"]["u"]
    if kind == "Z":
        return e["reference"]["u"]
    if kind == "A3":
        return "yes" if e["A3"]["holds"] else "no"
    if kind == "R":
        if not e["A3"]["holds"]:
            return None
        return e["arms"]["rows"]["u"] if s.endswith("-u") else ("yes" if e["closure"]["grade"] in ("relative", "absolute") else "no")
    if kind == "E":
        return "yes" if e["arms"]["earlier"]["beyond_reference"]["reading"] == "no added loss" else "no or unresolved"
    raise KeyError(i)


rows = []
for f in D:
    o = observed(f["id"])
    if o is None:
        rows.append({**f, "observed": None, "scored": False, "why": "its condition did not occur"})
    elif f["kind"] == "interval":
        rows.append({**f, "observed": o, "scored": True, "inside": bool(f["interval"][0] <= o <= f["interval"][1])})
    else:
        top = max(f["probabilities"], key=f["probabilities"].get)
        rows.append({**f, "observed": o, "scored": True, "top": top, "hit": bool(o == top), "probability_given_to_what_happened": f["probabilities"].get(o, 0.0)})
iv = [r for r in rows if r["scored"] and r["kind"] == "interval"]
oc = [r for r in rows if r["scored"] and r["kind"] == "outcome"]
res = {"intervals": {"scored": len(iv), "inside": sum(r["inside"] for r in iv)}, "outcomes": {"scored": len(oc), "hit": sum(r["hit"] for r in oc)},
       "not_scored": [r["id"] for r in rows if not r["scored"]], "rows": rows}
(RUN / "results" / "outcome.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
print(f"intervals inside: {res['intervals']['inside']} of {res['intervals']['scored']}; outcomes, the class given the highest probability happened: {res['outcomes']['hit']} of {res['outcomes']['scored']}")
print("| id | forecast | observed | |\n|---|---|---|---|")
for r in rows:
    fc = str(r.get("interval")) if r["kind"] == "interval" else ", ".join(f"{k} {v}" for k, v in r["probabilities"].items())
    ob = "not scored" if not r["scored"] else (f"{r['observed']:.3g}" if isinstance(r["observed"], float) else str(r["observed"]))
    print(f"| {r['id']} | {fc} | {ob} | {'' if not r['scored'] else ('inside' if r.get('inside') else 'outside') if r['kind'] == 'interval' else ('hit' if r['hit'] else 'miss')} |")
