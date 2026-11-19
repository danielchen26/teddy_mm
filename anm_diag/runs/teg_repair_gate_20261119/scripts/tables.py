"""Print the tables of README.md from results/*.json (markdown). Nothing is computed here.  python tables.py"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / "results"
Gt, Ph, Sc, Se = (json.load(open(RES / n)) for n in ("gate.json", "posthoc.json", "scale.json", "sensitivity.json"))
DS = {"frangieh2021": "melanoma", "papalexi2021": "THP-1"}
f = lambda v, d=2: f"{v:.{d}f}".replace("-", "−")
ci = lambda r: f"{f(r['u'])} [{f(r['ci'][0])}, {f(r['ci'][1])}] {r['reading']}"

print("**The gate, place A** (the PD-L1 output as a linear readout of the embedding; pooled)\n")
print("| given | melanoma | THP-1 | melanoma: \\|*u*\\| minus `cells` | THP-1: \\|*u*\\| minus `cells` |\n|---|---|---|---|---|")
for k in ("controls", "rows", "rows_hard", "rows_guides", "cells"):
    a = [Gt["data_sets"][ds]["places"]["A"][k] for ds in DS]
    d = ["" if k == "cells" else f"{f(x['abs_u_minus_cells']['difference'])} [{f(x['abs_u_minus_cells']['ci'][0])}, {f(x['abs_u_minus_cells']['ci'][1])}]" for x in a]
    print(f"| `{k}` | {ci(a[0])} | {ci(a[1])} | {d[0]} | {d[1]} |")
print("\n**By event, place A** (*u*; one gene, the fit that never saw it, all audit controls)\n")
ev = [(ds, g) for ds in DS for g in Gt["data_sets"][ds]["places"]["A"]["rows"]["events"]]
print("| given | " + " | ".join(f"{'mel' if ds == 'frangieh2021' else 'THP-1'} {g}" for ds, g in ev) + " |\n|---|" + "---|" * len(ev))
for k in ("controls", "rows", "cells"):
    print(f"| `{k}` | " + " | ".join(f(Gt["data_sets"][ds]["places"]["A"][k]["events"][g]["u"]) for ds, g in ev) + " |")
print("\n**`rows` by the number of events of the fit** (place A; mean \\|*u*\\| on the held-out event, every subset)\n")
print("| events of the fit | 0 (`controls`) | 1 | 2 | 3 | 4 |\n|---|---|---|---|---|---|")
for ds, n in DS.items():
    print(f"| {n} | " + " | ".join(f(Gt["data_sets"][ds]["rows_by_number_of_events"][str(k)]["mean_abs_u"]) for k in range(5)) + " |")
print("\n**The held-out displacement against the span of the events of the fit** (cosine)\n")
for ds, n in DS.items():
    print(f"- {n}: " + ", ".join(f"{g} {f(v['cosine_with_the_span'])}" for g, v in Gt["data_sets"][ds]["geometry"].items()))
print("\n**The gate, place B** (the head as trained plus a linear correction, bone-marrow term at weight 1; read as the predicted mean, as designed)\n")
print("| given | melanoma | THP-1 | bone marrow, PD-L1 score (0.408 as trained) |\n|---|---|---|---|")
for k in ("controls", "rows", "cells"):
    a = [Gt["data_sets"][ds]["places"]["B"][k] for ds in DS]
    print(f"| `{k}` | {ci(a[0])} | {ci(a[1])} | {f(a[0]['bone_marrow']['pdl1_corrected'], 3)} and {f(a[1]['bone_marrow']['pdl1_corrected'], 3)} |")
print("\n**The heads of `runs/teg_training_20261119` on two scales** (post hoc; pooled *u* of the own PD-L1 output)\n")
print("| arm | melanoma, as the predicted mean (as recorded) | melanoma, on the log scale | THP-1, as the predicted mean | THP-1, on the log scale |\n|---|---|---|---|---|")
for a in Sc["data_sets"]["frangieh2021"]["arms"]:
    r = [Sc["data_sets"][ds]["arms"][a] for ds in DS]
    w = lambda q: ci(q) + ("" if q["valid"] else ", wrong sign")
    print(f"| `{a}` | {w(r[0]['mean'])} | {w(r[0]['log'])} | {w(r[1]['mean'])} | {w(r[1]['log'])} |")
print("\n**A correction scoped to the cell type of its rows** (post hoc; pooled *u* on the log scale; bone marrow: PD-L1 score, 0.408 as trained)\n")
print("| correction | given | melanoma | THP-1 | bone marrow |\n|---|---|---|---|---|")
lab = {"0.0": "added to the head's output, no bone-marrow term", "1.0": "added to the head's output, bone-marrow term at weight 1 (place B)"}
for w in ("0.0", "1.0"):
    for k in ("controls", "rows", "cells"):
        r = [Ph["data_sets"][ds]["B_by_guard_weight"][w][k] for ds in DS]
        print(f"| {lab[w]} | `{k}` | {ci(r[0]['log'])} | {ci(r[1]['log'])} | {f(r[0]['bone_marrow']['pdl1_now'], 3)} and {f(r[1]['bone_marrow']['pdl1_now'], 3)} |")
for name, text in (("gate_add", "added to the head's output, only in this cell type"), ("gate_replace", "replaces the head's output, only in this cell type")):
    for k in ("controls", "rows", "cells"):
        r = [Ph["data_sets"][ds]["scoped"][name][k] for ds in DS]
        print(f"| {text} | `{k}` | {ci(r[0]['log'])} | {ci(r[1]['log'])} | {f(r[0]['bone_marrow']['pdl1_now'], 3)} and {f(r[1]['bone_marrow']['pdl1_now'], 3)} |")
print("\nplace A applied to bone marrow, PD-L1 score:", {DS[ds]: {k: round(v["pdl1_if_replaced"], 3) for k, v in Ph["data_sets"][ds]["A_on_bone_marrow"].items()} for ds in DS})
print("the gate between the cell type and bone marrow:", {DS[ds]: {k: round(v, 3) for k, v in Ph["data_sets"][ds]["gate"].items()} for ds in DS})
print("\n**The ridge weight** (post hoc; place A, pooled point values; the gate used 0.1, chosen by cross-validation on training cells)\n")
print("| ridge weight | melanoma: `controls` | `rows` | `cells` | THP-1: `controls` | `rows` | `cells` |\n|---|---|---|---|---|---|---|")
for lam in Se["data_sets"]["frangieh2021"]:
    print(f"| {lam} | " + " | ".join(f(Se["data_sets"][ds][lam][k]) for ds in DS for k in ("controls", "rows", "cells")) + " |")
print("\ndecision:", Gt["decision"])
