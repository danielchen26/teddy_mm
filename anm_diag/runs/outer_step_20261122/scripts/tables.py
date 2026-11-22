#!/usr/bin/env python3
"""The tables of README.md, from results/*.json.   python tables.py > ../results/tables.md"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / "results"
L = lambda f: json.loads((RES / f).read_text())
f2 = lambda v: f"{v:.2f}".replace("-", "−")
f3 = lambda v: f"{v:.3f}".replace("-", "−")
ci = lambda x, f=f2: f"{f(x['u'] if 'u' in x else x['added'])} [{f(x['ci'][0])}, {f(x['ci'][1])}]"
NAME = {"PDL1 mel_ifn": "PD-L1, melanoma, IFN-γ", "PDL1 mel_co": "PD-L1, melanoma, co-culture", "PDL1 thp1": "PD-L1, THP-1 (five events)", "MHC1 mel_co": "MHC-I, melanoma, co-culture",
        "PDL1 thp1 type A": "PD-L1, THP-1 without IRF1", "PDL1 thp1 type B": "PD-L1, THP-1, IRF1 alone"}
old = json.loads((RES.parents[1] / "teg_loop_20261120" / "results" / "round_development.json").read_text())["PDL1"]["groups"]
old.update({"MHC1 " + k: v for k, v in json.loads((RES.parents[1] / "teg_loop_20261120" / "results" / "round_replication.json").read_text())["MHC1"]["groups"].items()})
B, S, C = L("bridge.json")["targets"], L("selection.json"), L("census.json")
R = L("readers.json") if (RES / "readers.json").exists() else {}
o = ["### Bridge reading: the repaired output, the round against package 2", "",
     "| target | reference at Z | output in the round (reader of types) | output with the type from the record | beyond the reference | grade |", "|---|---|---|---|---|---|"]
for k in ("PDL1 mel_ifn", "PDL1 mel_co", "PDL1 thp1", "MHC1 mel_co"):
    e = B[k]
    was = old[k.split()[1] if k.startswith("PDL1") else k]["arms"]["rows"]
    o.append(f"| {NAME[k]} | {ci(e['reference'], f3)} {e['reference']['verdict']} | {ci(was, f3)} {was['verdict']} | {ci(e['arms']['rows'], f3)} {e['arms']['rows']['verdict']} | "
             f"{ci(e['arms']['rows']['beyond_reference'], f3)} {e['arms']['rows']['beyond_reference']['reading']} | {e['closure']['grade']} |")
o += ["", "### Gate 3 as amended, and the test with two controls", "",
      "| target | gate 3: output − reference on fitted cells | passes | the output as trained at the same gate | rows − controls-only (gain) | reading | test |", "|---|---|---|---|---|---|---|"]
for k in ("PDL1 mel_ifn", "PDL1 mel_co", "PDL1 thp1", "MHC1 mel_co", "PDL1 thp1 type A"):
    e = B[k]
    g, d = e["gate3"], e["rows_minus_controls"]
    o.append(f"| {NAME[k]} | {ci(g['against_the_reference'])} | {'yes' if g['passes'] else '**no**'} | {ci(g['control_the_output_as_trained'])}: {'passes' if g['control_the_output_as_trained']['passes'] else 'fails'} | "
             f"{f3(d['gain_difference'])} [{f3(d['ci'][0])}, {f3(d['ci'][1])}] | {d['reading']} | {e['test']['outcome']} |")
o += ["", "### What each arm removes (gain: fall of \\|*u*\\| of the output)", "", "| target | as trained | controls only | rows | every cell | rows of the other group | the other group's repair |", "|---|---|---|---|---|---|---|"]
for k in ("PDL1 mel_ifn", "PDL1 mel_co", "PDL1 thp1", "MHC1 mel_co"):
    e = B[k]
    a = e["arms"]
    o.append(f"| {NAME[k]} | {f2(e['as_trained']['u'])} | " + " | ".join((f"{f2(a[x]['gain'])} (leaves {f2(a[x]['u'])})" if x in a else "") for x in ("controls", "rows", "cells", "rows_other", "other")) + " |")
A, Bt = B["PDL1 thp1 type A"], B["PDL1 thp1 type B"]
o += ["", "### The THP-1 pool, three ways", "", "| read as | reference | repaired output | beyond the reference |", "|---|---|---|---|"]
for k in ("PDL1 thp1", "PDL1 thp1 type A"):
    e = B[k]
    o.append(f"| {NAME[k]} | {ci(e['reference'])} {e['reference']['verdict']} | {ci(e['arms']['rows'])} {e['arms']['rows']['verdict']} | {ci(e['arms']['rows']['beyond_reference'])} {e['arms']['rows']['beyond_reference']['reading']} |")
o.append(f"| IRF1 with its own row | {ci(Bt['reference'])} {Bt['reference']['verdict']} (a reader of Z with its own cells) | {ci(Bt['arms']['rows'])} {Bt['arms']['rows']['verdict']} | "
         f"{ci(Bt['arms']['rows']['beyond_reference'])} {Bt['arms']['rows']['beyond_reference']['reading']} |")
p = Bt["with_the_pooled_rows"]
o.append(f"| IRF1 with the pooled rows | the same | {ci(p)} {p['verdict']} | {ci(p['beyond_reference'])} {p['beyond_reference']['reading']} |")
o.append(f"| IRF1, every cell of the fitted half | the same | {ci(Bt['arms']['cells'])} {Bt['arms']['cells']['verdict']} | {ci(Bt['arms']['cells']['beyond_reference'])} {Bt['arms']['cells']['beyond_reference']['reading']} |")
iv = Bt["input_view"]
o += ["", f"IRF1 before the head: the RNA view (512 components, ridge, own cells) leaves {ci(iv['own'])} {iv['own']['verdict']}; Z leaves {ci(Bt['reference'])}; from the RNA view to Z "
      f"{ci(iv['from_the_input_view_to_Z_own_cells'])}: {iv['from_the_input_view_to_Z_own_cells']['reading']}."]
for k, title in (("PDL1 mel_co", "PD-L1 in co-culture"), ("MHC1 mel_co", "MHC-I in co-culture"), ("MHC1 mel_ifn", "MHC-I under IFN-γ: B2M and TAPBP")):
    if k not in R:
        continue
    r = R[k]
    o += ["", f"### Readers: {title} (own cells)", "", "| event | Z, ridge | Z, boosting | 512 components, ridge | 512 components, boosting | 4,000 genes, ridge | permuted | from the input view to Z |", "|---|---|---|---|---|---|---|---|"]
    for u, e in r["events"].items():
        s = e["by_stage"]
        c = lambda st, cl: f"{f2(s[st]['own'][cl]['u'])} {s[st]['own'][cl]['verdict'][:3]}" if cl in s[st]["own"] else ""
        g = e["from_the_input_view_to_Z"]
        o.append(f"| {u} | {c('Z', 'ridge')} | {c('Z', 'hgb')} | {c('RNA512', 'ridge')} | {c('RNA512', 'hgb')} | {c('RNA4K', 'ridge')} | {f2(e['permuted']['u'])} | {ci(g)}{' dominated' if g['dominated'] else ''} |")
    pc = r["positive_control"]
    o += ["", f"Positive control: of {pc['readings']} readings (two scopes, five readers, every event and the pool) {pc['supported']} supported, {pc['unresolved']} unresolved, {pc['rejected']} rejected: "
          f"{'holds' if pc['holds'] else 'does not hold'}. Permuted control: {'holds' if r['permuted_control']['holds'] else 'does not hold'}. "
          f"The input view {'can' if r['input_view']['can_be_read'] else 'cannot'} be read."]
o += ["", "### The selection third", "", "| readout | knockouts that pass A1 on a third | admitted | the next one |", "|---|---|---|---|"]
for r in ("PDL1", "MHC1", "HLAE"):
    x = S[r]
    o.append(f"| {r} | {', '.join(x['pass_A1']) or 'none'} | {', '.join(x['admitted']) or 'none'} | {x['largest_not_passing']['event']} ({x['largest_not_passing']['snr']:.1f}) |")
o += ["", "| events of the round, on a third | |", "|---|---|"] + [f"| {k} | " + ", ".join(f"{s} {v['snr']:.1f}" for s, v in x.items()) + " |" for k, x in S["events_of_the_round"].items()]
if (RES / "round2.json").exists():
    r2 = L("round2.json")
    o += ["", "### The round under package 2: MHC-I in melanoma under IFN-γ", "",
          "| event | cells | change (control SD) | head as trained | Z, other cells of the type | Z, own cells | A3 | rows | controls only | every cell | test |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for s, e in r2["events"].items():
        a = e["arms"]
        o.append(f"| {s} | {e['n_cells']} | {f2(e['response']['value'])} | {ci(e['as_trained'])} | {ci(e['reference']['units'])} | {ci(e['reference']['own'])} | {'passes' if e['admission']['A3'] else 'fails'} "
                 f"({f2(e['admission']['gap_head_to_units']['value'])}) | {ci(a['rows'])} | {ci(a['controls'])} | {ci(a['cells'])} | {e['test']['outcome']} |")
o += ["", "### The census and the final units", "", "| stratum | units (at least 100 cells) | read by the system | unread | final | control cells | with states | unread |", "|---|---|---|---|---|---|---|---|"]
for k, v in C["strata"].items():
    o.append(f"| {k.replace(' | ', ', ')} | {v['units']} | {v['read']} | {v['unread']} | {v['final']} | {v['control_cells']} | {v['control_cells_with_states']} | {v['control_cells_unread']} |")
o += ["", "| stratum | batch 1: chosen for an open item | calibration | screened as clear, but final |", "|---|---|---|---|"]
for k, v in C["strata"].items():
    if v["unread"]:
        o.append(f"| {k.replace(' | ', ', ')} | {', '.join(v['batch_chosen']) or 'none'} | {', '.join(v['batch_calibration']) or 'none'} | {', '.join(v['screened_clear_and_final']) or 'none'} |")
if (RES / "outcome.json").exists():
    sc = L("outcome.json")
    o += ["", "### The forecasts", "", "| hypothesis | intervals inside | outcomes: the class given the highest probability happened |", "|---|---|---|"]
    o += [f"| {h} | {v['intervals_inside']} | {v['most_likely_happened']} |" for h, v in sc["by_hypothesis"].items()]
    o += [f"| all | {sc['intervals']['inside']} of {sc['intervals']['n']} | {sc['probabilities']['most_likely_class_happened']} of {sc['probabilities']['n']} |", "",
          "| missed | forecast | observed |", "|---|---|---|"]
    for x in sc["rows"]:
        if x.get("inside") is False or x.get("most_likely_class_happened") is False:
            fc = x["forecast"] if isinstance(x["forecast"], list) else ", ".join(f"{k} {v}" for k, v in x["forecast"].items())
            o.append(f"| {x['hypothesis']}: {x['what']} | {fc} | {x['observed'] if isinstance(x['observed'], str) else f3(x['observed'])} |")
print("\n".join(o))
