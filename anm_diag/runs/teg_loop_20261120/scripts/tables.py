#!/usr/bin/env python3
"""Print the tables of README.md from the result files.   python tables.py"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / "results"
A = json.loads((RES / "round_development.json").read_text())["PDL1"]
S = json.loads((RES / "round_development_soft.json").read_text())["PDL1"]
B = json.loads((RES / "round_replication.json").read_text())
AF = json.loads((RES / "after_replication.json").read_text())
NAME = {"mel_ifn": "melanoma, IFN-γ", "mel_co": "melanoma, co-culture", "thp1": "THP-1"}
f = lambda v: f"{v:.2f}".replace("-0.00", "0.00")
iv = lambda q: f"{f(q['u'])} [{f(q['ci'][0])}, {f(q['ci'][1])}] {q['verdict']}"
sp = lambda b: f"{f(b['added'])} [{f(b['ci'][0])}, {f(b['ci'][1])}] {b['reading']}"
ARM = {"other": "the repair of the other group, as it is", "controls": "this group's control cells only", "rows_other": "this group's controls, the other group's rows",
       "rows": "**this group's controls and admitted rows**", "cells": "this group's controls and every knockout cell (reference)"}


def closing(groups):
    print("| group | `units` at Z | attainable | own output as trained | with the repair by rows | added beyond Z | closure | test of the hypothesis |\n|---|---|---|---|---|---|---|---|")
    for g, e in groups.items():
        if not e["arms"]:
            continue
        r = e["arms"]["rows"]
        print(f"| {NAME[g]} | {iv(e['grid']['scopes']['units'])} | {e['admission']['attainable']} | {f(e['as_trained']['u'])} | {iv(r)} | {sp(r['beyond_units_at_Z'])} | **{e['closure']['grade']}** | {e['test']['outcome']} |")


def arms(g, e):
    print(f"\n**{NAME[g]}**\n\n| output | own-output *u* | gain over the head as trained | added beyond Z | by gene |\n|---|---|---|---|---|")
    print(f"| the head as trained | {iv(e['as_trained'])} | | {f(e['as_trained']['beyond_units_at_Z']['added'])} | |")
    for a, v in e["arms"].items():
        print(f"| {ARM[a]} | {iv(v)} | {f(v['gain'])} [{f(v['gain_ci'][0])}, {f(v['gain_ci'][1])}] | {sp(v['beyond_units_at_Z'])} | " + ", ".join(f"{s} {f(q['u'])}" for s, q in v["events"].items()) + " |")


print("## Part A\n")
closing(A["groups"])
for g, e in A["groups"].items():
    arms(g, e)
print("\n**The grid on the head as trained, PD-L1** (pooled events; `frozen` is the head's own output)\n\n| group | response (control SD) | frozen | train | other | units | own | reading of the grid |\n|---|---|---|---|---|---|---|---|")
for g, e in A["groups"].items():
    G = e["grid"]
    print(f"| {NAME[g]} | {f(G['response']['value'])} | " + " | ".join(iv(G["scopes"][s]) for s in ("frozen", "train", "other", "units", "own")) + f" | {G['reading']} |")
print("\n**The input view** (the same events read from RNA, 512 components)\n\n| group | from RNA, other genes | from RNA, own cells | added from RNA to Z (own cells) | what Z leaves is |\n|---|---|---|---|---|")
for g, e in A["groups"].items():
    v = e["admission"]["input_view"]["pool"]
    print(f"| {NAME[g]} | {iv(v['units'])} | {iv(v['own'])} | {sp(v['added_from_RNA_to_Z_own_scope'])} | {e['admission']['rest'] or 'nothing: Z is supported'} |")
print("\n**Gates 3 and 6, and the gate of the system**\n\n| group | on the cells it was fitted on | guard, validation cells (fall) | guard, test cells: readout, mean, worst protein | own class given to control cells | to knockout cells |\n|---|---|---|---|---|---|")
for g, e in A["groups"].items():
    q = [x for x in A["gate"] if x["group"] == g]
    print(f"| {NAME[g]} | {f(e['gate3']['on_the_cells_it_was_fitted_on'])} | {e['gate3']['guard_on_validation_cells']['fall']:.4f} | {e['gate6']['guard_on_test_cells']['fall']:.4f}, {e['gate6']['guard_on_test_cells']['mean_fall']:.5f}, "
          f"{e['gate6']['guard_on_test_cells']['worst_item_fall']:.4f} | {min(x['own_on_control_cells'] for x in q):.3f} | {min(x['own_on_knockout_cells'] for x in q):.3f} |")
bm = [x for x in A["gate"] if x["group"].startswith("bone")]
print(f"\nBone-marrow test cells given to a group by the gate: {max(x['cells_given_to_a_group'] for x in bm)} of 5,000.")
print("\n**With the readouts mixed by the gate's probabilities** (`--soft-gate`, a sensitivity)\n\n| group | with the repair by rows | on the cells it was fitted on | closure |\n|---|---|---|---|")
for g, e in S["groups"].items():
    print(f"| {NAME[g]} | {iv(e['arms']['rows'])} | {f(e['gate3']['on_the_cells_it_was_fitted_on'])} | {e['closure']['grade']} |")
q = A["probe"]["grid"]
print(f"\nProbe {A['probe']['site']} (response {f(q['response']['value'])} control SD): " + ", ".join(f"{k} {iv(v)}" for k, v in q["scopes"].items()) + f". {A['probe']['admission']['status']}: {A['probe']['admission']['reason']}.")
print("\nLedger, part A:", A["ledger"]["events"], "events;", [f"{c['code']} {c['severity']}" for c in A["ledger"]["checks"]], "; forecasts inside:", A["ledger"]["forecast_record"]["inside"], "of", A["ledger"]["forecast_record"]["n"])
fr = A["ledger"]["forecast_record"]["rows"]
print("\n| hypothesis | arm | forecast | measured gain | inside |\n|---|---|---|---|---|")
for x in fr:
    print(f"| {x['hypothesis']} | {x['arm']} | [{f(x['forecast'][0])}, {f(x['forecast'][1])}] | {f(x['gain'])} | {'yes' if x['inside_forecast'] else 'no'} |")

print("\n## Part B\n\n**Admission**\n\n| readout | group | pooled response (control SD) | times its half-width | events that pass A1 on their own | `units` at Z | admission |\n|---|---|---|---|---|---|---|")
for r, R in B.items():
    for g, e in R["groups"].items():
        G, a = e["grid"], e["admission"]
        n = sum(x["admission"]["status"] == "admitted" for x in e["events"])
        print(f"| {r} | {NAME[g]} | {f(G['response']['value'])} | {G['response']['snr']:.1f} | {n} of {len(e['events'])} | {iv(G['scopes']['units'])} | " + (f"**target**, attainable {a['attainable']}" if a["diagnostic"] == "ready" else "not a target: " + "; ".join(a["reasons"])) + " |")
for r, R in B.items():
    for g, e in R["groups"].items():
        if e["arms"]:
            print(f"\n{r}:")
            arms(g, e)
            print(f"\nGate 3: on the cells it was fitted on {f(e['gate3']['on_the_cells_it_was_fitted_on'])} (objective reached: {e['gate3']['objective_reached']}); guard on validation cells kept: {e['gate3']['guard_on_validation_cells']['kept']}. "
                  f"Gate 6: fall on test cells {e['gate6']['guard_on_test_cells']['fall']:.4f}. Closure: {e['closure']['grade']}. Test: {e['test']['outcome']} ({e['test']['why']}). Ceiling of the grid: {f(e['grid']['ceiling']['coverage']['removes'])}.")
    print(f"\n{r}: loop:", json.dumps({k: v for k, v in R["loop"].items() if k != "open"}), "| ledger checks:", [f"{c['code']} {c['severity']}: {c['message']}" for c in R["ledger"]["checks"]])
v = AF["mel_co_MHC1"]["pool"]
print(f"\nAfter the run, MHC-I in co-culture from RNA: other genes {iv(v['from_RNA']['units'])}, own cells {iv(v['from_RNA']['own'])}; added from RNA to Z {sp(v['added_from_RNA_to_Z_own_scope'])}.")
