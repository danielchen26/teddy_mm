#!/usr/bin/env python3
"""The general cycle on the cells of the first system: the numbers of the first outer step, computed again by
anm_diag.cycle and compared with what the step recorded.

  python regress.py --observe OBSERVE1 --observe2 OBSERVE2 --work2 WORK2       (writes results/regression.json)

Folders and environment as runs/outer_step_20261122 (bridge.py, readers.py, round2.py). The cells, the halves, the
seeds, the grids and the admitted events are those of that step: nothing is selected and nothing new is read. What
differs is the code that reads them. Compared, number by number:

  bridge.json    six readings: the reference, the output as trained, every arm, every event, rows minus controls, gate 3, the grade
  round2.json    two events read one by one: the same quantities, and A3
  readers.json   three placings: the reader of every view, the gap from the input view to Z, the positive control and the
                 control view, with the step's rule for choosing a reader ("smallest"); the tool's rule ("union") is listed beside it

and then the open items: the targets built by cycle.target and cycle.target_by_event and routed by outer.open_items,
against the items the step recorded. Everything this script knows about the first system is in the part marked
"the system"; the rest is anm_diag.cycle.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
RUNS = RUN.parent
STEP = RUNS / "outer_step_20261122" / "results"
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
from anm_diag import cycle, outer  # noqa: E402
from _common import SEED, forward_np, unit, weights  # noqa: E402     (the system: its head and its embedding)
from build import GROUPS, OUTPUT  # noqa: E402

# --- the system ----------------------------------------------------------------------------------------------------------
TYPE = {"mel_ifn": {"cell_type": "melanoma", "stratum": "IFN-gamma"}, "mel_co": {"cell_type": "melanoma", "stratum": "co-culture"}, "thp1": {"cell_type": "THP-1", "stratum": "stimulated"}}
NAME = {"mel_ifn": "melanoma IFN-gamma", "mel_co": "melanoma co-culture", "thp1": "THP-1 stimulated"}
pair = lambda g: " and ".join(sorted({NAME[g], NAME[GROUPS[g]["other"]]}))
TYPE_A, PATHWAY = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"], ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
WORDS = {"self event": "the knockout of the readout's own gene", "upstream": "a change of the encoder"}


def cells_of(path, readout, head):
    d = np.load(path)
    c = dict(ids=d["ids"].astype(str), x=unit(d["Z"]).astype(np.float64), treated=d["A"].astype(int) == 1, unit=d["unit"].astype(str), y=d[f"y_{readout}"].astype(np.float64))
    c["o"] = head(c["x"])
    c["first"] = cycle.halves(np.where(c["treated"], c["unit"], "control"), SEED)          # the halves of the round
    return c


sub = lambda c, m: {k: v[m] for k, v in c.items()}
whole = lambda c, sites: cycle.LinearRows(c["x"], c["y"], c["o"], np.where(~c["treated"])[0], [(s, np.where(c["treated"] & (c["unit"] == s))[0]) for s in sites])
other_of = lambda c, sites: (c["x"], c["y"], c["treated"].astype(int), c["unit"], sites)


# --- the comparison ------------------------------------------------------------------------------------------------------
class Compare:
    def __init__(self):
        self.rows, self.n, self.worst, self.unequal = [], 0, 0.0, []

    def walk(self, new, old, path):
        if isinstance(new, dict):
            for k, v in new.items():
                if isinstance(old, dict) and k in old:
                    self.walk(v, old[k], path + [str(k)])
        elif isinstance(new, (list, tuple)):
            for k, (a, b) in enumerate(zip(new, old)):
                self.walk(a, b, path + [str(k)])
        elif isinstance(new, (bool, str)) or new is None:
            self.n += 1
            if new != old:
                self.unequal.append({"where": ".".join(path), "general": new, "recorded": old})
        else:
            self.n += 1
            d = abs(float(new) - float(old))
            self.worst = max(self.worst, d)
            if d > 1e-6:
                self.unequal.append({"where": ".".join(path), "general": float(new), "recorded": float(old), "difference": d})

    def add(self, name, new, old):
        n0, u0, w0 = self.n, len(self.unequal), self.worst
        self.worst = 0.0
        self.walk(new, old, [name])
        self.rows.append({"what": name, "values_compared": self.n - n0, "largest_difference": self.worst, "unequal": len(self.unequal) - u0})
        self.worst = max(self.worst, w0)


def as_bridge(e):
    """A reading of cycle.read in the layout of bridge.json."""
    arms = {a: {**x, "events": {s: v["arms"][a] for s, v in e["by_event"].items()}} for a, x in e["arms"].items()}
    return {"events": e["events"], "reference": e["reference"], "as_trained": e["as_trained"], "arms": arms, "rows_minus_controls": e["rows_minus_controls"],
            "gate3": {k: e["gate3"][k] for k in ("on_the_cells_it_was_fitted_on", "against_the_reference", "passes", "control_the_output_as_trained")}, "closure": e["closure"]}


def as_readers(p):
    ev = {}
    for u, e in p["events"].items():
        v = lambda x: {k: x[k] for k in ("u", "ci", "verdict")}
        ev[u] = {"by_stage": e["by_stage"], "Z_own": {"reader": e["model_stage"]["reader"], **v(e["model_stage"])}, "input_own": {"view": e["input"]["stage"], "reader": e["input"]["reader"], **v(e["input"])},
                 "from_the_input_view_to_Z": e["from_the_input_to_the_model_stage"], "permuted": e["control"]}
    return {"events": ev, "positive_control": p["positive_control"], "permuted_control": p["control_view"]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--observe2", type=Path, required=True)
    ap.add_argument("--work2", type=Path, required=True)
    a = ap.parse_args()
    work = Path(os.environ["ANM_LOOP_WORK"])
    bridge, readers, round2 = (json.loads((STEP / f).read_text()) for f in ("bridge.json", "readers.json", "round2.json"))
    bridge = bridge["targets"]
    recorded = json.loads((STEP / "open_items.json").read_text())
    census = json.loads((STEP / "census.json").read_text())
    selection = json.loads((STEP / "selection.json").read_text())
    P0, names, _j = weights()
    probe = Compare()                                       # the comparison itself is checked first: it must see a difference where there is one
    probe.add("probe", {"a": {"u": 1.0, "ci": [0.0, 1.0], "verdict": "support"}}, {"a": {"u": 1.00001, "ci": [0.0, 1.0], "verdict": "reject"}})
    assert probe.n == 4 and len(probe.unequal) == 2, "the comparison does not see a planted difference"
    K, E, grids = Compare(), {}, {}
    for readout, groups in (("PDL1", ["mel_ifn", "mel_co", "thp1"]), ("MHC1", ["mel_co"])):
        j = names.index(OUTPUT[readout])
        head = lambda x, j=j: forward_np(P0, x)[1][:, j]
        cells = {g: cells_of(work / f"{g}.npz", readout, head) for g in GROUPS if readout in GROUPS[g]["readouts"]}
        sites = {g: bridge[f"{readout} {g}"]["events"] if f"{readout} {g}" in bridge else [] for g in cells}
        W = {g: whole(c, sites[g]) for g, c in cells.items() if sites[g]}
        for g in groups:
            c, og = cells[g], GROUPS[g]["other"]
            by, draws = grids[readout, g] = cycle.load(a.observe2 / "mel_co_MHC1_z2" if (g, readout) == ("mel_co", "MHC1") else a.observe / f"{g}_{readout}")
            cl = lambda u, sc, by=by: by[u]["Z"][sc]["class"]
            other, fixed = (other_of(cells[og], sites[og]), {"other": W[og]}) if sites.get(og) else (None, None)
            out, inside, _f = cycle.cross_fit(c["x"], c["y"], c["o"], c["treated"], c["unit"], sites[g], c["first"], other=other, fixed=fixed)
            refs = {s: draws[f"{s}|shift|Z|units|{cl(s, 'units')}"] for s in sites[g]}
            e = E[f"{readout} {g}"] = cycle.read(c["y"], c["o"], c["treated"], c["unit"], sites[g], out, inside, draws[f"pool|shift|Z|units|{cl('pool', 'units')}"], refs, seed=SEED)
            assert abs(e["as_trained"]["u"] - draws["pool|shift|system|frozen|system"][0]) < 2e-4, "the resamples do not pair"
            K.add(f"bridge: {readout} {g}", as_bridge(e), bridge[f"{readout} {g}"])
            if (g, readout) == ("thp1", "PDL1"):
                byA, dA = cycle.load(a.observe2 / "thp1a_PDL1")
                cA = sub(c, np.isin(c["unit"], TYPE_A + ["control"]))
                oA, iA, _f = cycle.cross_fit(cA["x"], cA["y"], cA["o"], cA["treated"], cA["unit"], TYPE_A, cA["first"], other=other, fixed=fixed)
                eA = cycle.read(cA["y"], cA["o"], cA["treated"], cA["unit"], TYPE_A, oA, iA, dA[f"pool|shift|Z|units|{byA['pool']['Z']['units']['class']}"],
                                {s: dA[f"{s}|shift|Z|units|{byA[s]['Z']['units']['class']}"] for s in TYPE_A}, seed=SEED)
                K.add("bridge: PDL1 thp1 type A", as_bridge(eA), bridge["PDL1 thp1 type A"])
                mb = np.isin(c["unit"], ["IRF1", "control"])
                cB, refB = sub(c, mb), draws[f"IRF1|shift|Z|own|{cl('IRF1', 'own')}"]
                oB, iB, _f = cycle.cross_fit(cB["x"], cB["y"], cB["o"], cB["treated"], cB["unit"], ["IRF1"], cB["first"])
                eB = cycle.read(cB["y"], cB["o"], cB["treated"], cB["unit"], ["IRF1"], oB, iB, refB, seed=SEED)
                up, _s = cycle.output_u(c["y"][mb], out["rows"][mb], c["treated"][mb], None, SEED)       # the same event with the pooled rows
                ub, _s = cycle.output_u(cB["y"], oB["rows"], cB["treated"], None, SEED)
                K.add("bridge: PDL1 thp1 type B", {**as_bridge(eB), "with_the_pooled_rows": {**cycle.val(up), "beyond_reference": cycle.diff(up, refB), "against_its_own_row": cycle.diff(up, ub)}},
                      bridge["PDL1 thp1 type B"])
    # the round under package 2: two events, each read on its own
    j = names.index(OUTPUT["MHC1"])
    head = lambda x: forward_np(P0, x)[1][:, j]
    c = cells_of(a.work2 / "mel_ifn_mhc.npz", "MHC1", head)
    co = cells_of(work / "mel_co.npz", "MHC1", head)
    new = selection["MHC1"]["admitted"]
    by2, d2 = cycle.load(a.observe2 / "mel_ifn_mhc_MHC1")
    out, inside, _f = cycle.cross_fit(c["x"], c["y"], c["o"], c["treated"], c["unit"], new, c["first"], other=other_of(co, PATHWAY), fixed={"other": whole(co, PATHWAY)})
    out["rows_covered"] = cycle.cross_fit(c["x"], c["y"], c["o"], c["treated"], c["unit"], PATHWAY, c["first"])[0]["rows"]
    R2 = {}
    for s in new:
        e = R2[s] = cycle.read(c["y"], c["o"], c["treated"], c["unit"], [s], out, inside, d2[f"{s}|shift|Z|units|{by2[s]['Z']['units']['class']}"], seed=SEED)
        old = round2["events"][s]
        K.add(f"round2: {s}", {"reference": {"units": e["reference"]}, "as_trained": e["as_trained"], "arms": e["arms"], "rows_minus_controls": e["rows_minus_controls"],
                                "gate3": {k: e["gate3"][k] for k in ("on_the_cells_it_was_fitted_on", "against_the_reference", "passes", "control_the_output_as_trained")},
                                "admission": {"A3": e["A3"]["holds"]}}, old)
    # the placings
    views = {"MHC1 mel_ifn": ("mel_ifn_mhc_MHC1", a.observe2 / "mel_ifn_mhc_MHC1"), "PDL1 mel_co": ("mel_co_PDL1", a.observe / "mel_co_PDL1"), "MHC1 mel_co": ("mel_co_MHC1", a.observe2 / "mel_co_MHC1_z2")}
    placed, rules = {}, {}
    for k, (stem, z) in views.items():
        g = cycle.merge(cycle.load(z), *(cycle.load(a.observe2 / f"{stem}_{v}") for v in ("in512", "in4k", "perm")))
        placed[k] = cycle.place(*g, "Z", inputs=["RNA512", "RNA4K"], control="PERM", rule="smallest")
        K.add(f"readers: {k}", as_readers(placed[k]), readers[k])
        un = cycle.place(*g, "Z", inputs=["RNA512", "RNA4K"], control="PERM", rule="union")
        rules[k] = [{"event": u, "view": w, "smallest": [e[w].get("stage"), e[w]["reader"], e[w]["verdict"]], "union": [x[w].get("stage"), x[w]["reader"], x[w]["verdict"]]}
                    for u, e in placed[k]["events"].items() for x in [un["events"][u]] for w in ("model_stage", "input") if (e[w].get("stage"), e[w]["reader"], e[w]["verdict"]) != (x[w].get("stage"), x[w]["reader"], x[w]["verdict"])]
    thp = cycle.place(*cycle.merge(grids["PDL1", "thp1"], cycle.load(a.observe / "thp1_PDL1_rna")), "Z", inputs=["RNA512"], rule="smallest")
    # the open items, by the general path
    ev = lambda k: f"general path: {k}"
    T = [cycle.target("PDL1", "mel_ifn", TYPE["mel_ifn"], E["PDL1 mel_ifn"], pair=pair("mel_ifn"), evidence=ev("PDL1 mel_ifn"), pool={"more_units_of_the_type": 0},
                      self_events=[{"site": "CD274", "u": 0.91, "ci": [0.87, 0.95], "evidence": "runs/teg_loop_20261120 round_development.json:PDL1.probe"}]),
         cycle.target("PDL1", "mel_co", TYPE["mel_co"], E["PDL1 mel_co"], placed["PDL1 mel_co"], pair=pair("mel_co"), evidence=ev("PDL1 mel_co"), pool={"more_units_of_the_type": None}),
         cycle.target("PDL1", "thp1", TYPE["thp1"], E["PDL1 thp1"], thp, events_alone=["IRF1"], pair=pair("thp1"), evidence=ev("PDL1 thp1"), pool={"more_units_of_the_type": None}),
         cycle.target("MHC1", "mel_co", TYPE["mel_co"], E["MHC1 mel_co"], placed["MHC1 mel_co"], pair=pair("mel_co"), evidence=ev("MHC1 mel_co"), pool={"more_units_of_the_type": None}),
         cycle.target_by_event("MHC1", "mel_ifn", TYPE["mel_ifn"], R2, placed["MHC1 mel_ifn"], evidence=ev("MHC1 mel_ifn"), pool={"more_units_of_the_type": 0})]
    del T[2]["input_view"]                                  # the step read no input view for the pool of this type, only for one of its events
    for readout, g, snr in (("HLAE", "mel_ifn", 2.21), ("HLAE", "mel_co", 3.12)):       # not targets since the round; no cell of theirs is read here
        T.append({"readout": readout, "group": g, "type": TYPE[g], "events_n": 5, "admitted": False, "response": {"snr": snr}, "pool": {"more_units_of_the_type": None, "admissible_now": [] if g == "mel_ifn" else None}})
    items = outer.open_items(T, {"typed_right": None, "selection_on_audited_cells": None, "final_units": len(census["roles"]["final"]), "undeclared_cases": [], "words": WORDS})
    key = lambda i: (i["readout"], i["group"] or "", i["event"], i["reading"])
    got, want = sorted(map(key, items)), sorted(map(key, recorded["items"]))
    status = {f"{t['readout']} {t['group']}": outer.target_status(t) for t in T}
    was = {f"{t['readout']} {t['group']}": t["status"] for t in recorded["targets"]}
    res = {"note": "the numbers of runs/outer_step_20261122, computed again by anm_diag.cycle on the same cells, halves, seeds and grids",
           "comparisons": K.rows, "values_compared": K.n, "largest_difference": max(r["largest_difference"] for r in K.rows), "unequal": K.unequal,
           "reader_rule": {"note": "where the tool's rule for choosing the reader of a view (union) gives another reader or verdict than the step's (smallest |u|)", "differences": rules},
           "statuses": {"general": status, "recorded": was, "equal": status == was},
           "items": {"general": len(got), "recorded": len(want), "equal": got == want, "only_general": [list(k) for k in got if k not in want], "only_recorded": [list(k) for k in want if k not in got]}}
    (RUN / "results").mkdir(exist_ok=True)
    (RUN / "results" / "regression.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("| what | values compared | largest difference | unequal |\n|---|---|---|---|")
    for r in K.rows:
        print(f"| {r['what']} | {r['values_compared']} | {r['largest_difference']:.1e} | {r['unequal']} |")
    print("values", K.n, "| unequal", len(K.unequal), "| statuses equal", status == was, "| items", len(got), "general,", len(want), "recorded, equal:", got == want)
    for x in K.unequal[:40]:
        print("  ", x)
    print("reader rule:", json.dumps(rules, default=float)[:1500])
    print("items only in the general path:", res["items"]["only_general"], "\nitems only in the record:", res["items"]["only_recorded"])
    print(status)


if __name__ == "__main__":
    main()
