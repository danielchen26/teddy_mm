#!/usr/bin/env python3
"""The round on batch 1: every admitted event read on its own, by the general code (anm_diag.cycle, anm_diag.outer).

  python round.py --work WORK --observe OBSERVE        (writes results/round.json; refuses to run when it exists)

WORK: the folder build.py wrote. OBSERVE: the grids of `anm_diag scope` for the SPEC files of WORK, one folder per SPEC.
Environment as runs/teg_loop_20261120.

A new event is a type of its own until a reading says otherwise: its reference is the reader of Z fitted with its own
cells (scope `own` of the grid), and resampling is by cell. Every cell is read by a fit that did not see it (two halves).
Arms: controls, rows (the event's own record), cells, earlier (the records of the stratum's earlier events only),
rows_other and other (the records, and the whole repair, of melanoma under IFN-gamma). DESIGN.md, Sections 4 and 5.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN, RUNS = HERE.parent, HERE.parents[1]
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
from anm_diag import cycle, outer  # noqa: E402
from _common import SEED, forward_np, unit, weights  # noqa: E402     (the system: its head and its embedding)
from build import OUTPUT  # noqa: E402

ROLES = json.loads((RUN / "results" / "roles.json").read_text())
TYPE = {g: {k: v[k] for k in ("cell_type", "stratum")} for g, v in ROLES["strata"].items()}
OTHER = {"PDL1": ("mel_ifn.npz", ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]), "MHC1": ("mel_ifn_mhc.npz", ["B2M", "TAPBP"])}
WORDS = {"self event": "the knockout of the readout's own gene", "upstream": "a change of the encoder"}
sub = lambda c, m: {k: v[m] for k, v in c.items()}


def cells_of(path, readout, head):
    d = np.load(path)
    c = dict(ids=d["ids"].astype(str), x=unit(d["Z"]).astype(np.float64), treated=d["A"].astype(int) == 1, unit=d["unit"].astype(str), y=d[f"y_{readout}"].astype(np.float64))
    c["o"] = head(c["x"])
    c["first"] = cycle.halves(np.where(c["treated"], c["unit"], "control"), SEED)
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--observe", type=Path, required=True)
    a = ap.parse_args()
    out_file = RUN / "results" / "round.json"
    if out_file.exists():
        raise SystemExit("round.json exists: this reading was made once and is not made again")
    sel = json.loads((RUN / "results" / "selection.json").read_text())["strata"]
    P0, names, _j = weights()
    res = {"margin": cycle.MARGIN, "step_margin": cycle.STEP, "note": "every event on its own, resampled by cell; its reference is a reader of Z with its own cells", "targets": {}, "not_targets": []}
    T, statuses = [], {}
    for g, S in sel.items():
        for r, x in S["readouts"].items():
            if not x["admitted"]:
                res["not_targets"].append({"stratum": g, "readout": r, "why": "no unit of the batch passes A1 and A2 on the selection third",
                                           "largest": {k: x["table"][0][k] for k in ("event", "cells", "change_in_control_sd", "snr")} if x["table"] else None})
                T.append({"readout": r, "group": g, "type": TYPE[g], "events_n": len(x["table"]), "admitted": False, "response": {"snr": x["table"][0]["snr"] if x["table"] else None},
                          "pool": {"more_units_of_the_type": None, "admissible_now": None}})
                continue
            j = names.index(OUTPUT[r])
            head = lambda z, j=j: forward_np(P0, z)[1][:, j]
            c = cells_of(a.work / f"b1_{g}.npz", r, head)
            by, draws = cycle.merge(cycle.load(a.observe / f"b1_{g}_{r}"), *(cycle.load(a.observe / f"b1_{g}_{r}_{v}") for v in ("in512", "in4k", "perm")))
            placed = cycle.place(by, draws, "Z", inputs=["RNA512", "RNA4K"], control="PERM")
            other = fixed = None
            if r in OTHER:
                oc = cells_of(a.work / OTHER[r][0], r, head)
                sites_o = OTHER[r][1]
                other = (oc["x"], oc["y"], oc["treated"].astype(int), oc["unit"], sites_o)
                fixed = {"other": cycle.LinearRows(oc["x"], oc["y"], oc["o"], np.where(~oc["treated"])[0], [(s, np.where(oc["treated"] & (oc["unit"] == s))[0]) for s in sites_o])}
            early = ROLES["fit_units"].get(g) or []
            if early:
                ec = cells_of(a.work / {"thp1": "thp1.npz", "mel_co": "mel_co.npz"}[g], r, head)       # the cells of the earlier rounds: only their events are used
                ec = sub(ec, np.isin(ec["unit"], early))
            reads, per_event = {}, {}
            for s in x["admitted"]:
                m = ~c["treated"] | (c["unit"] == s)
                cs = sub(c, m)
                out, inside, fits = cycle.cross_fit(cs["x"], cs["y"], cs["o"], cs["treated"], cs["unit"], [s], cs["first"], other=other, fixed=fixed)
                if early:
                    n = len(cs["y"])
                    xx, yy, oo = np.vstack([cs["x"], ec["x"]]), np.r_[cs["y"], ec["y"]], np.r_[cs["o"], ec["o"]]
                    tt, uu, ff = np.r_[cs["treated"], ec["treated"]], np.r_[cs["unit"], ec["unit"]], np.r_[cs["first"], ec["first"]]
                    out["earlier"] = cycle.cross_fit(xx, yy, oo, tt, uu, early, ff)[0]["rows"][:n]
                ref = draws[f"{s}|shift|Z|own|{by[s]['Z']['own']['class']}"]
                e = cycle.read(cs["y"], cs["o"], cs["treated"], cs["unit"], [s], out, inside, ref, seed=SEED)
                k = f"{s}|shift|system|frozen|system"
                assert k in draws and abs(e["as_trained"]["u"] - draws[k][0]) < 2e-4, "the resamples of the round and of the tool do not pair"
                e["cells"] = {"event": int(cs["treated"].sum()), "control": int((~cs["treated"]).sum())}
                e["selection_third"] = next(q for q in x["table"] if q["event"] == s)
                e["grid"] = {st: {sc: {k2: v[k2] for k2 in ("u", "ci", "verdict", "class") if k2 in v} for sc, v in by[s][st].items() if isinstance(v, dict) and "u" in v} for st in by[s]}
                e["fits"] = fits
                reads[s] = e
                t1 = cycle.target_by_event(r, g, TYPE[g], {s: e}, placed, evidence=f"round.json:{r} {g}.events.{s}", pool={"more_units_of_the_type": None})
                per_event[s] = {"status": outer.target_status(t1), "target": t1}
            t = cycle.target_by_event(r, g, TYPE[g], reads, placed, evidence=f"round.json:{r} {g}", pool={"more_units_of_the_type": None})
            T.append(t)
            statuses[f"{r} {g}"] = outer.target_status(t)
            res["targets"][f"{r} {g}"] = {"readout": r, "stratum": g, "type": TYPE[g], "events": reads, "placed": placed, "status": statuses[f"{r} {g}"],
                                          "status_by_event": {s: v["status"] for s, v in per_event.items()}, "target": t}
    items = outer.open_items(T, {"typed_right": None, "selection_on_audited_cells": None, "final_units": 133, "undeclared_cases": [], "words": WORDS})
    res["open_items"] = items
    res["statuses"] = statuses
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda v: f"{v:6.3f}"
    for k, t in res["targets"].items():
        print(f"\n== {k}: status {t['status']}")
        for s, e in t["events"].items():
            print(f"  {s}: cells {e['cells']} | as trained {f(e['as_trained']['u'])} {e['as_trained']['verdict']} | reference (Z, own) {f(e['reference']['u'])} {e['reference']['verdict']} | A3 {e['A3']['holds']} "
                  f"({f(e['A3']['added'])}) | gate 3 {e['gate3']['passes']} | grade {e['closure']['grade']} | rows-controls {e['rows_minus_controls']['reading']} | by event: {t['status_by_event'][s]}")
            for arm, v in e["arms"].items():
                print(f"     {arm:11s} {f(v['u'])} [{f(v['ci'][0])}, {f(v['ci'][1])}] {v['verdict']:10s} beyond reference {f(v['beyond_reference']['added'])} {v['beyond_reference']['reading']}")
            p = t["placed"]["events"][s]
            print(f"     input view {p['input']['stage']} {p['input']['reader']} {f(p['input']['u'])} {p['input']['verdict']} | input to Z {f(p['from_the_input_to_the_model_stage']['added'])} "
                  f"{p['from_the_input_to_the_model_stage']['reading']} | permuted {f(p['control']['u'])}")
    print("\nnot targets:", [(x["stratum"], x["readout"]) for x in res["not_targets"]])
    print("open items:", [(i["readout"], i["group"], i["event"], i["reading"], i["lever"]) for i in items])


if __name__ == "__main__":
    main()
