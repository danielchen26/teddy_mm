#!/usr/bin/env python3
"""After the step: every target with its status, the open items derived again, what became of the 21 items the step
started from, the forecasts scored, the ledger rules, and package 2.

  python items.py            (reads results/*.json and the round's results; writes open_items.json, open_items.md,
                              outcome.json, package2.json; reads no cell)
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
RES = RUN / "results"
sys.path.insert(0, str(RUN.parents[1]))
from anm_diag import outer  # noqa: E402

TYPE = {"mel_ifn": {"cell_type": "melanoma", "stratum": "IFN-gamma"}, "mel_co": {"cell_type": "melanoma", "stratum": "co-culture"}, "thp1": {"cell_type": "THP-1", "stratum": "stimulated"}}
NAME = {"mel_ifn": "melanoma IFN-gamma", "mel_co": "melanoma co-culture", "thp1": "THP-1 stimulated"}
OTHER = {"mel_ifn": "mel_co", "mel_co": "mel_ifn", "thp1": "mel_ifn"}
pair = lambda g: " and ".join(sorted({NAME[g], NAME[OTHER[g]]}))
V = lambda x: {k: x[k] for k in ("u", "ci", "verdict")}
load = lambda f: json.loads((RES / f).read_text())


def from_bridge(key, e, readers, sel):
    readout, g = key.split()[:2]
    rows = e["arms"]["rows"]
    t = {"readout": readout, "group": g, "type": TYPE[g], "events_n": len(e["events"]), "admitted": True, "response": {"snr": None}, "attainable": e["closure"]["attainable"],
         "evidence": f"bridge.json:{key}", "reference": {"units": V(e["reference"]), "own": V(e["grid"]["own"]), "evidence": f"bridge.json:{key}.reference"},
         "model": {"own": V(rows), "added": rows["beyond_reference"], "reading": None, "rows_of_type": len(e["events"]), "evidence": f"bridge.json:{key}.arms.rows",
                   "events": {s: {**V(x), "added": x["beyond_reference"]} for s, x in rows["events"].items()}},
         "gates": {"evaluable": bool(e["gate3"]["passes"]), "why": [], "guard_kept": True}, "arms": {},
         "pool": {"more_units_of_the_type": 0 if (readout, g) == ("PDL1", "mel_ifn") else None}}
    if not e["gate3"]["passes"]:
        g3 = e["gate3"]["against_the_reference"]
        t["gates"]["why"].append({"question": "gate 3 as amended (3a) does not pass: on the cells the repair was fitted on, the interval of |u of the output| - |u of the reference| reaches above the step margin",
                                  "value": g3["added"], "evidence": f"bridge.json:{key}.gate3",
                                  "note": f"{g3['added']:.2f} [{g3['ci'][0]:.2f}, {g3['ci'][1]:.2f}]: on its own cells the output overshoots in some resamples (u below zero) and |u| counts an overshoot as a loss; "
                                          f"the interval rests on {len(e['events'])} genes. The audit half reads no added loss ({rows['beyond_reference']['added']:.2f}). An amendment for the next iteration"})
    for a in ("other", "rows_other"):
        if a in e["arms"]:
            t["arms"][a] = {**e["arms"][a]["beyond_reference"], "pair": pair(g), "evidence": f"bridge.json:{key}.arms.{a}"}
    r = readers.get(f"{readout} {g}")
    if r:
        iv, pc = r["input_view"], r["positive_control"]
        t["input_view"] = {"own": V(iv["own"]), "gap_to_later": iv["gap_to_later"], "classes": iv["classes"], "nonlinear": iv["nonlinear"], "evidence": f"readers.json:{readout} {g}",
                           "positive_control": pc["holds"], "readout_control": None,
                           "positive_control_note": f"of {pc['readings']} readings of the events and the pool (Z by ridge and boosting, 512 components by ridge and boosting, 4,000 genes by ridge) "
                                                    f"{pc['supported']} are supported, {pc['unresolved']} unresolved, {pc['rejected']} rejected",
                           "note": f"best reader of the input view: {iv['own']['view']}, {iv['own']['reader']}"}
    return t


def score(bridge, readers, sel, census, round2):
    """Every forecast of declared.json against what was read. An interval forecast is inside or outside; a probability forecast
    records the probability that was given to what happened."""
    rows = []

    def iv(h, what, lohi, got):
        rows.append({"hypothesis": h, "what": what, "forecast": lohi, "observed": float(got), "inside": bool(lohi[0] <= got <= lohi[1])})

    def pr(h, what, probs, got):
        p = probs.get(got, 0.0)
        rows.append({"hypothesis": h, "what": what, "forecast": probs, "observed": got, "probability_given": p, "most_likely_class_happened": bool(p == max(probs.values()))})

    yes = lambda p, b: ({"yes": p, "no": round(1 - p, 2)}, "yes" if b else "no")
    T = {"PDL1 mel_ifn": [-0.02, 0.04], "PDL1 mel_co": [-0.13, -0.07], "PDL1 thp1": [-0.11, -0.03], "MHC1 mel_co": [-0.07, -0.01]}
    for k, f in T.items():
        g = bridge[k]["gate3"]
        iv("P1", f"gate 3 against the reference, {k}", f, g["against_the_reference"]["added"])
        pr("P1", f"gate 3 passes, {k}", *yes(0.95, g["passes"]))
        pr("P1", f"control: the output as trained fails the gate, {k}", *yes(0.95, not g["control_the_output_as_trained"]["passes"]))
    P2 = {"PDL1 mel_ifn": ({"unresolved": 0.85, "refuted": 0.10, "confirmed": 0.05}, [0.05, 0.10], 0.05), "PDL1 mel_co": ({"unresolved": 0.65, "confirmed": 0.35}, [0.13, 0.19], 0.85),
          "PDL1 thp1": ({"unresolved": 0.85, "confirmed": 0.10, "refuted": 0.05}, [0.30, 0.42], 0.97), "MHC1 mel_co": ({"unresolved": 0.70, "confirmed": 0.30}, [0.24, 0.31], 0.99)}
    for k, (probs, f, p) in P2.items():
        e = bridge[k]
        pr("P2", f"outcome of the test with two controls, {k}", probs, e["test"]["outcome"])
        iv("P2", f"gain of rows minus gain of controls-only, {k}", f, e["rows_minus_controls"]["gain_difference"])
        pr("P2", f"that difference lies above the step margin, {k}", *yes(p, e["rows_minus_controls"]["reading"] == "a step"))
    R = sel["events_of_the_round"]
    ok = lambda k, genes=None: all(v["A1"] for s, v in R[k].items() if genes is None or s in genes)
    pr("P4", "every event of PDL1 mel_ifn passes on a third", *yes(0.80, ok("PDL1 mel_ifn")))
    pr("P4", "every event of PDL1 mel_co passes on a third", *yes(0.97, ok("PDL1 mel_co")))
    pr("P4", "IRF1 in THP-1 passes on a third", *yes(0.40, R["PDL1 thp1"]["IRF1"]["A1"]))
    pr("P4", "the other four of THP-1 pass on a third", *yes(0.97, ok("PDL1 thp1", ["IFNGR1", "IFNGR2", "JAK2", "STAT1"])))
    iv("P4", "library knockouts that pass for PD-L1 on a third, beyond the readout's own gene and the five", [0, 1], len(sel["PDL1"]["admitted"]) - 5)
    r = bridge["PDL1 mel_ifn"]["arms"]["rows"]
    iv("S1", "PDL1 mel_ifn, repaired output, type from the record", [-0.01, 0.05], r["u"])
    pr("S1", "PDL1 mel_ifn verdict", {"support": 0.60, "unresolved": 0.40}, r["verdict"])
    iv("S1", "PDL1 mel_co, repaired output", [0.10, 0.13], bridge["PDL1 mel_co"]["arms"]["rows"]["u"])
    iv("S1", "PDL1 thp1 (five events pooled), repaired output", [0.14, 0.17], bridge["PDL1 thp1"]["arms"]["rows"]["u"])
    pr("S1", "PDL1 mel_ifn closed (absolute) on development cells", *yes(0.55, bridge["PDL1 mel_ifn"]["closure"]["grade"] == "absolute"))
    A, Bt = bridge["PDL1 thp1 type A"], bridge["PDL1 thp1 type B"]
    iv("S2", "type A reference", [0.02, 0.12], A["reference"]["u"])
    pr("S2", "type A reference verdict", {"unresolved": 0.65, "support": 0.20, "reject": 0.15}, A["reference"]["verdict"])
    iv("S2", "type A repaired output", [0.06, 0.15], A["arms"]["rows"]["u"])
    pr("S2", "type A repaired output verdict", {"unresolved": 0.70, "support": 0.10, "reject": 0.20}, A["arms"]["rows"]["verdict"])
    pr("S2", "type A adds no loss beyond its reference", *yes(0.80, A["arms"]["rows"]["beyond_reference"]["reading"] == "no added loss"))
    iv("S2", "type B (IRF1, its own row) repaired output", [0.35, 0.55], Bt["arms"]["rows"]["u"])
    pr("S2", "type B verdict reject", *yes(0.90, Bt["arms"]["rows"]["verdict"] == "reject"))
    pr("S2", "type B adds no loss beyond a reader of Z with its own cells", *yes(0.60, Bt["arms"]["rows"]["beyond_reference"]["reading"] == "no added loss"))
    g = Bt["input_view"]["from_the_input_view_to_Z_own_cells"]
    iv("S2", "type B from the RNA view to Z, own cells", [0.10, 0.35], g["added"])
    pr("S2", "type B reading upstream", *yes(0.45, g["reading"] == "loss"))
    v = readers.get("PDL1 mel_co")
    if v:
        e = v["events"]["pool"]["by_stage"]
        iv("V1", "PDL1 mel_co, 512 components, boosting, own cells", [0.15, 0.30], e["RNA512"]["own"]["hgb"]["u"])
        iv("V1", "PDL1 mel_co, 4,000 genes, ridge, own cells", [0.08, 0.20], e["RNA4K"]["own"]["ridge"]["u"])
        pr("V1", "the input view is no longer dominated", *yes(0.60, not v["input_view"]["dominated"]))
        pr("V1", "positive control reached", *yes(0.15, v["positive_control"]["holds"]))
        pr("V1", "permuted view reads nothing", *yes(0.95, v["permuted_control"]["holds"]))
    v = readers.get("MHC1 mel_co")
    if v:
        e = v["events"]["pool"]["by_stage"]
        iv("V2", "MHC1 mel_co, Z, boosting, own cells", [0.20, 0.40], e["Z"]["own"]["hgb"]["u"])
        iv("V2", "MHC1 mel_co, 512 components, boosting, own cells", [0.22, 0.42], e["RNA512"]["own"]["hgb"]["u"])
        iv("V2", "MHC1 mel_co, 4,000 genes, ridge, own cells", [0.22, 0.36], e["RNA4K"]["own"]["ridge"]["u"])
        allrej = all(x["verdict"] == "reject" for st in ("Z", "RNA512", "RNA4K") for x in e[st]["own"].values())
        pr("V2", "every reader at every stage rejects (pool, own cells)", *yes(0.80, allrej))
        pr("V2", "positive control reached", *yes(0.10, v["positive_control"]["holds"]))
        b = v["events"]["pool"]["boosting_reads_better_than_ridge_by"]
        pr("V2", "boosting reads at least 0.10 better than ridge at a stage", *yes(0.20, any(x["value"] >= 0.10 for x in b.values())))
    adm = set(sel["MHC1"]["admitted"])
    pr("E1", "B2M and TAPBP both admitted on a third", *yes(0.95, {"B2M", "TAPBP"} <= adm))
    pr("E1", "exactly those two", *yes(0.70, adm == {"B2M", "TAPBP"}))
    for gname in ("SSR2", "STAT3", "CD44"):
        pr("E1", f"{gname} admitted on a third", *yes(0.15, gname in adm))
    if round2:
        ev = round2["events"]
        iv("E2", "head as trained, B2M", [0.85, 1.0], ev["B2M"]["as_trained"]["u"])
        iv("E2", "head as trained, TAPBP", [0.85, 1.0], ev["TAPBP"]["as_trained"]["u"])
        iv("E2", "Z, own cells, B2M", [0.30, 0.90], ev["B2M"]["reference"]["own"]["u"])
        pr("E2", "Z, own cells, B2M rejected", *yes(0.85, ev["B2M"]["reference"]["own"]["verdict"] == "reject"))
        iv("E2", "Z, own cells, TAPBP", [0.40, 0.95], ev["TAPBP"]["reference"]["own"]["u"])
        pr("E2", "Z, own cells, TAPBP rejected", *yes(0.85, ev["TAPBP"]["reference"]["own"]["verdict"] == "reject"))
        pr("E2", "the group passes A3", *yes(0.45, round2["pool_of_the_two_events"]["A3"] or any(x["admission"]["A3"] for x in ev.values())))
        v = readers.get("MHC1 mel_ifn")
        if v:
            iv("E2", "input view, own cells, B2M (best reader)", [0.20, 0.85], v["events"]["B2M"]["input_own"]["u"])
            for gname, p in (("B2M", 0.30), ("TAPBP", 0.25)):
                x = v["events"][gname]
                pr("E2", f"reading upstream, {gname}", *yes(p, x["from_the_input_view_to_Z"]["reading"] == "loss"))
            pr("E2", "positive control reached for the new target", *yes(0.15, v["positive_control"]["holds"]))
    for name, st in census["strata"].items():
        if st["unread"]:
            iv("E3", f"share of unread units set aside as final, {name}", [0.30, 0.37], st["final"] / st["unread"])
    n_iv = [x for x in rows if "inside" in x]
    n_pr = [x for x in rows if "probability_given" in x]
    return {"rows": rows, "intervals": {"n": len(n_iv), "inside": sum(x["inside"] for x in n_iv)},
            "probabilities": {"n": len(n_pr), "most_likely_class_happened": sum(x["most_likely_class_happened"] for x in n_pr),
                              "mean_probability_given_to_what_happened": round(sum(x["probability_given"] for x in n_pr) / max(len(n_pr), 1), 3),
                              "given_zero": [x["what"] for x in n_pr if x["probability_given"] == 0]},
            "by_hypothesis": {h: {"intervals_inside": f"{sum(x.get('inside', False) for x in rows if x['hypothesis'] == h and 'inside' in x)} of {sum('inside' in x for x in rows if x['hypothesis'] == h)}",
                                  "most_likely_happened": f"{sum(x.get('most_likely_class_happened', False) for x in rows if x['hypothesis'] == h and 'probability_given' in x)} of "
                                                          f"{sum('probability_given' in x for x in rows if x['hypothesis'] == h)}"} for h in sorted({x["hypothesis"] for x in rows})}}


def main():
    bridge, readers, sel, census, declared = load("bridge.json")["targets"], load("readers.json"), load("selection.json"), load("census.json"), load("declared.json")
    round2 = load("round2.json") if (RES / "round2.json").exists() else None
    before = json.loads((RUN.parent / "outer_loop_design_20261120" / "results" / "open_items.json").read_text())
    targets = [from_bridge(k, bridge[k], readers, sel) for k in ("PDL1 mel_ifn", "PDL1 mel_co", "PDL1 thp1", "MHC1 mel_co")]
    # the CD274 probe of the round stays what it was
    targets[0]["self_events"] = [{"site": "CD274", "u": 0.91, "ci": [0.87, 0.95], "evidence": "runs/teg_loop_20261120 round_development.json:PDL1.probe", "note": "not admitted (A2); no reader of Z reads it"}]
    # IRF1 read on its own (type B of the provisional split): its reference and its input view
    B = bridge["PDL1 thp1 type B"]
    thp = next(t for t in targets if t["group"] == "thp1")
    thp["model"]["events"]["IRF1"]["before"] = {"reference": {"units": V(B["reference"]), "evidence": "bridge.json:PDL1 thp1 type B.reference"},
                                                "input_view": {"own": V(B["input_view"]["own"]), "gap_to_later": B["input_view"]["from_the_input_view_to_Z_own_cells"], "classes": ["ridge"],
                                                               "nonlinear": False, "positive_control": True, "evidence": "bridge.json:PDL1 thp1 type B.input_view"}}
    if round2:
        ev = round2["events"]
        order = ["loss", "unresolved", "no added loss"]
        worst = min(ev, key=lambda s: order.index(ev[s]["arms"]["rows"]["beyond_reference"]["reading"]))
        w = ev[worst]
        rd = readers.get("MHC1 mel_ifn") or {}
        t = {"readout": "MHC1", "group": "mel_ifn", "type": TYPE["mel_ifn"], "events_n": len(ev), "admitted": True, "response": {"snr": None},
             "attainable": "relative" if any(x["admission"]["attainable"] == "relative" for x in ev.values()) else "absolute", "evidence": "round2.json",
             "reference": None, "control_adds_a_loss": bool(any(x["admission"]["A3"] for x in ev.values())),
             "model": {"own": V(w["arms"]["rows"]), "added": w["arms"]["rows"]["beyond_reference"], "reading": None, "rows_of_type": len(ev), "evidence": f"round2.json:events.{worst}.arms.rows", "events": {}},
             "gates": {"evaluable": all(x["gate3"]["passes"] for x in ev.values()), "why": [], "guard_kept": True}, "arms": {}, "pool": {"more_units_of_the_type": 0}}
        if not t["gates"]["evaluable"]:
            t["gates"]["why"].append({"question": "gate 3 as amended (3a) does not pass for an event of the new target", "evidence": "round2.json",
                                      "note": "; ".join(f"{s}: {x['gate3']['against_the_reference']['added']:.2f} {[round(v, 2) for v in x['gate3']['against_the_reference']['ci']]}" for s, x in ev.items())})
        for s, x in ev.items():
            r = (rd.get("events") or {}).get(s)
            b = {"reference": {"units": x["reference"]["own"], "evidence": f"round2.json:events.{s}.reference.own"}}
            if r:
                pc = rd["positive_control"]
                b["input_view"] = {"own": V(r["input_own"]), "gap_to_later": {k: r["from_the_input_view_to_Z"][k] for k in ("added", "ci", "reading")}, "classes": ["ridge", "hgb"], "nonlinear": True,
                                   "positive_control": pc["holds"], "evidence": f"readers.json:MHC1 mel_ifn.events.{s}",
                                   "positive_control_note": f"of {pc['readings']} readings of the two events and their pool, {pc['supported']} are supported, {pc['unresolved']} unresolved, {pc['rejected']} rejected",
                                   "note": f"best reader of the input view: {r['input_own']['view']}, {r['input_own']['reader']}"}
            t["model"]["events"][s] = {**V(x["arms"]["rows"]), "added": x["arms"]["rows"]["beyond_reference"], "before": b}
        targets.append(t)
    # what the claim declared at the start and this step did not turn into a target
    for readout, g, snr, why in (("HLAE", "mel_ifn", 2.21, "pooled 2.2 (wanted 3); on a selection third of the library no knockout passes for HLA-E (the largest is 1.6)"),
                                 ("HLAE", "mel_co", 3.12, "pooled 3.1, no event passes on its own")):
        targets.append({"readout": readout, "group": g, "type": TYPE[g], "events_n": 5, "admitted": False, "response": {"snr": snr}, "note": why,
                        "evidence": "runs/teg_loop_20261120 round_replication.json", "pool": {"more_units_of_the_type": None, "admissible_now": [] if g == "mel_ifn" else None}})
    final = census["roles"]["final"]
    world = {"typed_right": None, "selection_on_audited_cells": None, "final_units": len(final), "undeclared_cases": []}
    items = outer.open_items(targets, world)
    for i in items:                                          # an estimate that the pool cannot decide is reported as it is (a terminal state; read on development cells)
        if i["next_step_needs"] == "at the margin":
            i["state"] = "at the margin"
    c = outer.counts(items)
    diff = declared["diff"]
    for d in diff:                                           # the diff names the items of the list the step started from
        d["items_before"] = d.pop("items")
        d["items"] = []
    findings = [x for x in outer.outer_check(targets, items, diff=[], hypotheses=declared["hypotheses"], roles={"final": final}, reads={u: 0 for u in final})]
    out = {"after": "the first outer step (package 2), and one round under it", "targets": [{"readout": t["readout"], "group": t["group"], "status": outer.target_status(t), "attainable": t.get("attainable"),
                                                                                             "own_output": (t.get("model") or {}).get("own"), "added_beyond_the_reference": (t.get("model") or {}).get("added"),
                                                                                             "gates": t.get("gates")} for t in targets],
           "items": items, "count_by_lever": c["by_lever"], "count_by_reading": c["by_reading"], "count_by_what_the_next_step_needs": c["by_what_the_next_step_needs"],
           "coverage": outer.coverage(targets, items), "ledger_check": findings, "final_units": len(final), "items_before": len(before["items"])}
    (RES / "open_items.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float) + "\n")
    f = lambda v: "" if v is None else (f"{v:.2f}".replace("-", "−") if isinstance(v, float) else str(v))
    ci = lambda x: "" if not x else f"[{f(x[0])}, {f(x[1])}]"
    o = ["## Targets after the step", "", "| readout | group | status | own output | added beyond the reference |", "|---|---|---|---|---|"]
    for t in targets:
        m = t.get("model")
        o.append(f"| {t['readout']} | {t['group']} | {outer.target_status(t)} | {(f(m['own']['u']) + ' ' + ci(m['own']['ci']) + ' ' + m['own']['verdict']) if m else ''} | "
                 f"{(f(m['added']['added']) + ' ' + ci(m['added']['ci']) + ' ' + m['added']['reading']) if m else ''} |")
    o += ["", "## Open items after the step", "", "| id | readout | group | event | reading | lever | question | value | the next step needs |", "|---|---|---|---|---|---|---|---|---|"]
    for i in items:
        val = "; ".join(f"{f(x['value'])} {ci(x['ci'])}".strip() + (f" ({x['group']})" if len(i["observations"]) > 1 else "") for x in i["observations"])
        o.append(f"| {i['id']} | {i['readout']} | {i['group'] or ''} | {i['event']} | {i['reading']} | {i['lever']} | {i['question']} | {val} | {i.get('forecast') or i['next_step_needs']} |")
    o += ["", "| lever | open items |", "|---|---|"] + [f"| {k} | {v} |" for k, v in c["by_lever"].items()]
    (RES / "open_items.md").write_text("\n".join(o) + "\n")
    print("\n".join(o))
    print(out["coverage"], findings)
    # what became of the 21 items the step started from
    rd = lambda k: readers.get(k) or {}
    r_ifn = bridge["PDL1 mel_ifn"]["arms"]["rows"]
    co_p, co_m = rd("PDL1 mel_co"), rd("MHC1 mel_co")
    dom = (co_p.get("input_view") or {}).get("dominated")
    Bt = bridge["PDL1 thp1 type B"]
    amend = lambda k: ("closed", f"amendment {k} of protocol version 3 is in force (LOOP_PROTOCOL.md)")
    disp = {
        "O01": ("open", "not acted on: the owner's decision (a limit for knockouts of the readout's own gene, or a system that is told the perturbation)"),
        "O02": ("at the margin", f"bridge reading with the type from the record: {r_ifn['u']:.3f} [{r_ifn['ci'][0]:.3f}, {r_ifn['ci'][1]:.3f}], {r_ifn['verdict']}; on a selection third no further "
                                 f"knockout of the library passes for PD-L1, so the pool cannot decide it. Read on development cells"),
        "O03": ("open", "not acted on: needs units the system has not been read on"),
        "O04": ("open", "not acted on: needs units in which the condition varies; for PD-L1 the only further event of co-culture (NPC1) fell to the final units"),
        "O05": ("decided on development units" if dom is False else "open",
                None if not co_p else f"richer readers: the best reader of the input view leaves {co_p['input_view']['own']['u']:.2f}; from the input view to Z "
                f"{co_p['input_view']['gap_to_later']['added']:.2f} {[round(v, 2) for v in co_p['input_view']['gap_to_later']['ci']]}; dominated: {dom}"),
        "O06": ("open", None if not co_p else f"richer readers: of {co_p['positive_control']['readings']} readings {co_p['positive_control']['supported']} are supported. "
                + ("The positive control holds now" if co_p["positive_control"]["holds"] else "No positive control: the owner's decision (a tolerance for this target)")),
        "O07": ("at the margin", "unchanged: Z reads 0.12 [0.01, 0.19]; no further event of the type is left for development (NPC1 is final)"),
        "O08": ("open", f"read again without IRF1: the reference of the other four is {bridge['PDL1 thp1 type A']['reference']['u']:.2f} "
                        f"{[round(v, 2) for v in bridge['PDL1 thp1 type A']['reference']['ci']]}, {bridge['PDL1 thp1 type A']['reference']['verdict']}; new THP-1 events are in batch 1"),
        "O09": ("decided on development units", f"the split as proposed is not supported: with its own row IRF1 leaves {Bt['arms']['rows']['u']:.2f}, with the pooled rows "
                f"{Bt['with_the_pooled_rows']['u']:.2f}. Its loss is before the head: the RNA view reads {Bt['input_view']['own']['u']:.2f}, Z {Bt['reference']['u']:.2f}; "
                f"from the RNA view to Z {Bt['input_view']['from_the_input_view_to_Z_own_cells']['added']:.2f} "
                f"{[round(v, 2) for v in Bt['input_view']['from_the_input_view_to_Z_own_cells']['ci']]}. Read `upstream` (a new item)"),
        "O10": ("open", "not acted on: needs units the system has not been read on"),
        "O11": ("open", "not acted on: needs new THP-1 events; BRD4 and CMTM6 are in batch 1"),
        "O12": ("decided on development units", "B2M and TAPBP are admitted on a selection third; the target was read (round2.json). What it showed is in the new items"),
        "O13": ("open", None if not co_m else f"richer readers: of {co_m['positive_control']['readings']} readings {co_m['positive_control']['supported']} are supported. "
                + ("The positive control holds now" if co_m["positive_control"]["holds"] else "No positive control: the owner's decision (a tolerance for this target)")),
        "O14": amend("3a"), "O15": amend("3b"),
        "O16": ("open", "on a selection third no knockout of the library passes for HLA-E (the largest is 1.6): the owner's decision (not asked, or the experiment list)"),
        "O17": ("open", "not acted on: units of the pool without states, or the owner's decision"),
        "O18": amend("3e"), "O19": amend("3d"), "O20": amend("3c"),
        "O21": ("closed", f"{len(final)} units are set aside as final (census.json); none has states yet"),
    }
    assert set(disp) == {i["id"] for i in before["items"]}
    out["what_became_of_the_items_before"] = [{"id": i["id"], "readout": i["readout"], "group": i["group"], "reading": i["reading"], "lever": i["lever"], "question": i["question"],
                                               "state": disp[i["id"]][0], "evidence": disp[i["id"]][1]} for i in before["items"]]
    out["states_of_the_items_before"] = {s: sum(v[0] == s for v in disp.values()) for s in outer.STATES}
    (RES / "open_items.json").write_text(json.dumps(out, indent=1, ensure_ascii=False, default=float) + "\n")
    o2 = ["", "## What became of the 21 items the step started from", "", "| id | reading | lever | state after the step | evidence |", "|---|---|---|---|---|"]
    o2 += [f"| {x['id']} | {x['reading']} | {x['lever']} | {x['state']} | {x['evidence']} |" for x in out["what_became_of_the_items_before"]]
    (RES / "open_items.md").write_text("\n".join(o + o2) + "\n")
    print("\n".join(o2))
    print(out["states_of_the_items_before"])
    package = {"package": 2, "sealed": "22 November 2026", "against": "package 1: runs/teg_loop_20261120",
               "E": {"events": {"PDL1": {g: bridge[f"PDL1 {g}"]["events"] for g in ("mel_ifn", "mel_co", "thp1")}, "MHC1": {"mel_co": bridge["MHC1 mel_co"]["events"], "mel_ifn": sel["MHC1"]["admitted"]}},
                     "admission": {"MHC1 mel_ifn": "on the selection third (selection.json)", "the others": "as admitted in the round of 20 November, on all their cells (kept; amendment 3d)"},
                     "selection_third": sel["rule"], "roles": {"final": f"{len(final)} units (census.json: roles.final)", "development": "every unit the system has been read on",
                                                               "unread_not_final": f"{len(census['roles']['unread_not_final'])} units without states", "batch_1": census["batch_1"]["units"]}},
               "V": {"stages": ["input view: RNA as log1p of counts per 10,000 (512 components; 4,000 genes)", "Z: the embedding the head reads", "the system's own output"],
                     "readers": {"Z": ["ridge", "boosting"], "512 components": ["ridge", "boosting"], "4,000 genes": ["ridge"]}, "controls": ["the permuted view", "the positive control"],
                     "readouts": ["PDL1", "MHC1", "HLAE"]},
               "S": {"type": "cell type x condition, taken from the record (amendment 3e)", "pool": "the admitted events of a type",
                     "provisional": "the THP-1 pool as two types was operated once and is not kept: with its own row IRF1 reads worse than with the pooled rows",
                     "candidate_fields": ["the place of the gene in the pathway", "the relation of the gene to the readout (own gene; part of the measured complex; upstream of it)"]},
               "C": {"margin": 0.10, "step_margin": 0.10, "guard": 0.02, "targets": [{"readout": t["readout"], "group": t["group"], "status": outer.target_status(t)} for t in targets]},
               "P": "LOOP_PROTOCOL.md, version 3", "M": "the head with the repaired outputs of the round of 20 November; not changed in this step",
               "diff": [{"lever": d["lever"], "change": d["change"], "items_before": d["items_before"], "hypothesis": d["hypothesis"]} for d in diff]}
    (RES / "package2.json").write_text(json.dumps(package, indent=1, ensure_ascii=False) + "\n")
    sc = score(bridge, readers, sel, census, round2)
    (RES / "outcome.json").write_text(json.dumps(sc, indent=1, ensure_ascii=False, default=float) + "\n")
    print("forecasts:", sc["intervals"], sc["probabilities"])
    for h, v in sc["by_hypothesis"].items():
        print("  ", h, v)
    for x in sc["rows"]:
        if x.get("inside") is False or x.get("most_likely_class_happened") is False:
            print("   missed:", x["hypothesis"], x["what"], "| forecast", x["forecast"], "| observed", x["observed"] if isinstance(x["observed"], str) else round(x["observed"], 3))


if __name__ == "__main__":
    main()
