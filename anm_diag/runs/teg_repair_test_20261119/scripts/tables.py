"""Print the tables of README.md from the result files.   python tables.py"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / "results"
r = json.loads((RES / "reading.json").read_text())
p = json.loads((RES / "posthoc.json").read_text())
GENES = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
SET = [("1. same events, new cells", r["coculture"]["same_events"]), ("2. a new gene", r["coculture"]["new_gene"]),
       ("3. a new condition", r["from_ifn"]["same_genes"]), ("4. a new gene and a new condition", r["from_ifn"]["new_gene"])]
f = lambda v: f"{v:.2f}".replace("-0.00", "0.00")
iv = lambda q: f"{f(q['u'])} [{f(q['ci'][0])}, {f(q['ci'][1])}] {q['reading']}"

print("**The four settings** (pooled over the five genes; knockout cells resampled by gene)\n")
print("| setting | the head as trained | `controls` | `rows` | `cells` | outcome of `rows` |\n|---|---|---|---|---|---|")
for name, t in SET:
    print(f"| {name} | {iv(t['as_trained'])} | {iv(t['controls'])} | **{iv(t['rows'])}** | {iv(t['cells'])} | **{t['outcome']}** |")
print("\n**`rows` against the two others** (difference of |u| on the same resamples)\n")
print("| setting | `rows` − `controls` | `rows` − `cells` | `rows`, knockout cells resampled one by one |\n|---|---|---|---|")
for name, t in SET:
    a, b, q = t["rows"]["abs_u_minus_controls"], t["rows"]["abs_u_minus_cells"], t["rows"]
    print(f"| {name} | {a['difference']:.3f} [{a['ci'][0]:.3f}, {a['ci'][1]:.3f}] | {b['difference']:.3f} [{b['ci'][0]:.3f}, {b['ci'][1]:.3f}] | [{f(q['ci_by_cell'][0])}, {f(q['ci_by_cell'][1])}] {q['reading_by_cell']} |")
print("\n**`rows` by event** (cells resampled one by one)\n")
print("| setting | " + " | ".join(GENES) + " |\n|---|" + "---|" * len(GENES))
for name, t in SET:
    print(f"| {name} | " + " | ".join(f"{f(t['rows']['events'][g]['u'])} [{f(t['rows']['events'][g]['ci'][0])}, {f(t['rows']['events'][g]['ci'][1])}] {t['rows']['events'][g]['reading']}" for g in GENES) + " |")
print("| measured change of PD-L1, co-culture | " + " | ".join(f(r['coculture']['same_events']['rows']['events'][g]['measured_change']) for g in GENES) + " |")
print("\n**Other splits into halves** (`rows`, pooled point values; the rule reads seed 0)\n")
print("| setting | seed 0 | 1 | 2 | 3 | 4 |\n|---|---|---|---|---|---|")
for key, name in (("same_events", "1. same events, new cells"), ("new_gene", "2. a new gene")):
    print(f"| {name} | {r['coculture'][key]['rows']['u']:.3f} | " + " | ".join(f"{r['coculture']['other_splits'][str(s)][key]['rows']['u']:.3f}" for s in range(1, 5)) + " |")
print("\n**Guard and gate**\n")
print("| setting | PD-L1 on bone-marrow test cells, as trained → with `rows` | fall of the mean over proteins | gate on bone-marrow test cells (mean, highest) |\n|---|---|---|---|")
for name, t in SET:
    b = t["bone_marrow"]; g = b["gate_on_bone_marrow_test_cells"]
    print(f"| {name} | {b['rows']['pdl1_as_trained']:.4f} → {b['rows']['pdl1_now']:.4f} | {b['rows']['mean_fall']:.6f} | {g['mean']:.3f}, {g['highest']:.2f} |")
g1, g2 = r["coculture"]["gate"], r["from_ifn"]["gate"]
print(f"\nGate on the co-culture cells of the reading: {min(q['knockout_read'] for q in g1):.3f} (fitted in the co-culture condition), {g2['knockout_read']:.3f} (fitted in the IFN-γ condition); lowest cell {min(q['knockout_read_lowest'] for q in g1):.2f} and {g2['knockout_read_lowest']:.2f}.")
print("\n**The probe: CD274** (287 knockout cells, IFN-γ condition; measured change " + f(r["probe"]["arms"]["rows"]["measured_change"]) + ")\n")
print("| the head as trained | `controls` | `rows` | `cells` |\n|---|---|---|---|")
print("| " + " | ".join(iv(r["probe"]["arms"][k]) for k in ("as_trained", "controls", "rows", "cells")) + " |")
o = r["other_genes"]
print(f"\n**The {o['n_genes']} other genes** ({o['n_cells']} knockout cells; measured change: median size {o['measured']['median_abs']:.3f}, largest {o['measured']['largest_abs']:.2f})\n")
print("| output | median size of the predicted change | largest | genes beyond 0.2 | correlation with the measured changes |\n|---|---|---|---|---|")
for k in ("as_trained", "controls", "rows", "cells"):
    a = o["arms"][k]
    print(f"| {k} | {a['predicted_median_abs']:.3f} | {a['predicted_largest_abs']:.2f} ({a['largest_gene']}) | {a['genes_predicted_beyond_0.2']} | {a['correlation_with_measured']:.2f} |")
print("\n**Post hoc 1: the cells of the fit against the cells of the reading** (setting 1, pooled point values)\n")
print("| | `controls` | `rows` | `cells` |\n|---|---|---|---|")
h = p["fit_half_against_read_half"]
print("| on the half that was fitted | " + " | ".join(f(h[k]["fit_half"]) for k in ("controls", "rows", "cells")) + " |")
print("| on the other half (the reading) | " + " | ".join(f(h[k]["read_half"]) for k in ("controls", "rows", "cells")) + " |")
print("\n**Post hoc 2: by ridge weight** (setting 1, pooled point values; the rule chose 0.32 and 0.1 in the two halves)\n")
print("| ridge weight | `controls` | `rows` | `cells` |\n|---|---|---|---|")
for lam, v in p["by_ridge_weight"].items():
    print(f"| {lam} | " + " | ".join(f(v[k]) for k in ("controls", "rows", "cells")) + " |")
v5 = p["five_folds"]
print("\nPost hoc 3, five folds: " + "; ".join(f"`{k}` {iv(v5[k])}" for k in ("controls", "rows", "cells")))
t2 = p["two_conditions"]
print("\n**Post hoc 4: the two conditions**\n")
print("| gene | measured change, IFN-γ | co-culture | cosine of the two displacements | length, IFN-γ | co-culture |\n|---|---|---|---|---|---|")
for g in GENES:
    q = t2["by_gene"][g]
    print(f"| {g} | {f(q['measured_change_ifn'])} | {f(q['measured_change_coculture'])} | {f(q['cosine_of_the_displacements'])} | {q['length_ifn']:.3f} | {q['length_coculture']:.3f} |")
print(f"\nDistance between the control means of the two conditions: {t2['distance_between_the_control_means']:.2f}; mean length of a displacement: {t2['mean_length_of_a_displacement']:.2f}. Cosine of the readouts fitted in the two conditions: "
      + ", ".join(f"`{k}` {f(v)}" for k, v in t2["cosine_of_the_readouts"].items()) + ".")
m = p["mixed_repair"]
print(f"Mixed repair ({m['what']}): {iv(m)}; by event " + ", ".join(f"{g} {f(v)}" for g, v in m["events"].items()) + ".")
