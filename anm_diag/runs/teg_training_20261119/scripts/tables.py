"""Print the tables of README.md from results/reading.json and results/fits.json (markdown). Nothing is computed here.

  python tables.py
"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parent.parent / "results"
R = json.load(open(RES / "reading.json"))
FITS = json.load(open(RES / "fits.json"))
NAME = {"frangieh2021": "Melanoma, IFN-γ", "papalexi2021": "THP-1"}
ORDER = ["as_trained", "data", "class", "edge05", "edge1", "edge2", "graph", "judge", "rows", "data_rows"]
LABEL = {"as_trained": "head as trained", "data": "`data`", "class": "`class`", "edge05": "`edge05`", "edge1": "`edge1`", "edge2": "`edge2`", "graph": "`graph`", "judge": "`judge`", "rows": "**`rows`**", "data_rows": "`data_rows`"}
STANDS = {"as_trained": "", "data": "the repair already made", "class": "an ontology", "edge05": "a knowledge graph, 0.5 SD", "edge1": "a knowledge graph, 1 SD", "edge2": "a knowledge graph, 2 SD",
          "graph": "a KG embedding", "judge": "a decision engine", "rows": "the TEG", "data_rows": "the repair + the TEG"}
f = lambda v, d=2: f"{v:.{d}f}".replace("-", "−")
ci = lambda r, d=2: f"{f(r['u'], d)} [{f(r['ci'][0], d)}, {f(r['ci'][1], d)}]"

for ds, e in R["data_sets"].items():
    print(f"\n**{NAME[ds]}** ({e['n_knockout']:,} knockout cells of 5 genes, each through the head that never saw its gene; {e['n_audit_controls']:,} audit controls)\n")
    print("| arm | stands for | *u*, own PD-L1 output | reading | slope among audit controls | output shifts the measured way | \\|*u*\\| minus that of `rows` | bone marrow, mean fall | trains the head well |")
    print("|---|---|---|---|---|---|---|---|---|")
    for a in ORDER:
        r = e["arms"][a]; d = r["abs_u_minus_reference"]
        diff = "" if a == "rows" else f"{f(d['difference'])} [{f(d['ci'][0])}, {f(d['ci'][1])}]"
        print(f"| {LABEL[a]} | {STANDS[a]} | {ci(r)} | {r['reading']} | {f(r['slope_among_audit_controls'])} | {'yes' if r['own_shift_has_the_measured_sign'] else '**no**'} | {diff} | "
              f"{f(r['bone_marrow']['mean_fall'], 3)} | {'yes' if R['trained_well'][a] else 'no'} |")

print("\n**The ten held-out events** (*u* of the own output; one gene, the head that never saw it, all audit controls)\n")
genes = [(ds, g) for ds, e in R["data_sets"].items() for g in e["arms"]["rows"]["events"]]
print("| arm | " + " | ".join(f"{'mel' if ds == 'frangieh2021' else 'THP-1'} {g}" for ds, g in genes) + " | inside ±0.10 |")
print("|---|" + "---|" * (len(genes) + 1))
for a in ORDER:
    print(f"| {LABEL[a]} | " + " | ".join(f(R["data_sets"][ds]["arms"][a]["events"][g]["u"]) for ds, g in genes) + f" | {R['family'][a]['events_inside_the_margin']} |")

print("\n**The events of the fit against the held-out event** (mean over the five heads of a data set)\n")
print("| arm | melanoma: *u* on the events of the fit | on the held-out event | THP-1: on the events of the fit | on the held-out event |")
print("|---|---|---|---|---|")
for a in ORDER[1:]:
    row = []
    for ds, e in R["data_sets"].items():
        ev = e["arms"][a]["events"]
        fit = [v for x in ev.values() for v in x["u_on_the_events_of_the_fit"].values()]
        row += [f(sum(fit) / len(fit)), f(sum(x["u"] for x in ev.values()) / len(ev))]
    print(f"| {LABEL[a]} | " + " | ".join(row) + " |")

print("\n**Without the statistic** (mean over the five held-out events of a data set)\n")
print("| arm | melanoma: correlation with the measured protein | shift of the own output / measured shift | THP-1: correlation | shift / measured shift |")
print("|---|---|---|---|---|")
for a in ORDER:
    row = []
    for ds, e in R["data_sets"].items():
        ev = list(e["arms"][a]["events"].values())
        row += [f(sum(x["plain"]["correlation_with_measured"] for x in ev) / len(ev)), f(sum(x["plain"]["own_shift_sd"] for x in ev) / sum(x["plain"]["measured_shift_sd"] for x in ev))]
    print(f"| {LABEL[a]} | " + " | ".join(row) + " |")

print("\n**The fits** (five per cell; kept step of 400; bone-marrow outputs moved, root mean square)\n")
print("| arm | melanoma: kept step | bone-marrow outputs moved | THP-1: kept step | bone-marrow outputs moved |")
print("|---|---|---|---|---|")
for a in ORDER[1:]:
    row = []
    for ds in R["data_sets"]:
        L = FITS[ds][a]; ks = [v["kept_step"] for v in L.values()]; mv = [v["bone_marrow_outputs_moved"] for v in L.values()]
        row += [f"{min(ks)}" if min(ks) == max(ks) else f"{min(ks)} to {max(ks)}", f"{f(min(mv), 3)} to {f(max(mv), 3)}"]
    print(f"| {LABEL[a]} | " + " | ".join(row) + " |")
s2 = R["statements"]["S2_rows_lower_than_the_repair_already_made"]
print("\nS2, |u| of `data` minus |u| of `rows`:", {ds: f"{f(v['difference'], 3)} [{f(v['ci'][0], 3)}, {f(v['ci'][1], 3)}]" for ds, v in s2.items()})
print("refit of `data` against the heads of runs/loop_guard_20261112, largest difference of the PD-L1 output:", {ds: round(e["arm_data_against_the_heads_of_loop_guard_20261112"]["largest_difference_of_the_pdl1_output"], 4) for ds, e in R["data_sets"].items()})
print("bone marrow, head as trained:", round(R["bone_marrow_as_trained"]["mean_pearson"], 4))
