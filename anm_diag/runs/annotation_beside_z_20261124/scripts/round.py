#!/usr/bin/env python3
"""The round with the head reading Z, and with the head reading the first annotation position beside Z, on the same
cells, halves and resamples. The round of runs/input_beside_z_20261123 (the head reads the input view beside Z) is
made again beside them, for a paired comparison.

  python round.py --work WORK --observe OBSERVE [--results RESULTS]      (writes round.json; refuses to run when it exists)

WORK: the folder build_zc.py wrote. OBSERVE: the grids of `anm_diag scope` for its SPEC files, one folder per SPEC.
The reading code is anm_diag.cycle, unchanged; the repair is cycle.LinearRows. DESIGN.md, Sections 3 and 4.
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
sys.path.insert(0, str(RUNS / "input_beside_z_20261123" / "scripts"))
from anm_diag import cycle  # noqa: E402
from _common import SEED, unit  # noqa: E402
from build_zv import TARGETS  # noqa: E402

RECORDED = {"T1": "PDL1 thp1 type B", "T2": "MHC1 mel_co", "C1": "PDL1 thp1 type A", "C2": "PDL1 mel_co", "C3": "PDL1 mel_ifn"}
FULL = {"thp1a": "thp1"}                                       # a stratum whose halves are those of a larger one
WIDE = ["ZC", "C", "ZV"]                                       # ZC: the declared stage; C: the annotation position alone; ZV: the run before, again
ARMS = ("controls", "rows", "cells")
RANK = {"absolute": 2, "relative": 1}
STEP = cycle.STEP
step = lambda iv: "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")


def gain(a, b):
    """|a| - |b| on paired resamples: what the second stage removes that the first leaves."""
    d = np.abs(a) - np.abs(b)
    iv = cycle.interval(d[1:])
    return {"gain": float(d[0]), "ci": iv, "reading": step(iv), "above_zero": bool(iv[0] > 0), "no_harm": bool(iv[0] > -STEP)}


def cells(work, g, readout):
    d = np.load(work / f"{FULL.get(g, g)}_zc.npz")
    c = dict(ids=d["ids"].astype(str), Z=unit(d["Z"]).astype(np.float64), ZV=d["ZV"].astype(np.float64), C=unit(d["C"]).astype(np.float64), ZC=d["ZC"].astype(np.float64),
             treated=d["A"].astype(int) == 1, unit=d["unit"].astype(str), y=d[f"y_{readout}"].astype(np.float64), o=d[f"frozen_{readout}"].astype(np.float64))
    c["first"] = cycle.halves(np.where(c["treated"], c["unit"], "control"), SEED)      # the halves of the earlier rounds, made on the whole stratum
    if g in FULL:
        keep = np.isin(c["ids"], np.load(work / f"{g}_zc.npz")["ids"].astype(str))
        c = {k: v[keep] for k, v in c.items()}
    return c


def one(c, events, pooled, by, draws, stage):
    """The round for one target with the head reading `stage`. Returns the reading and the draws of every output."""
    if pooled:
        m, scope, key = np.ones(len(c["y"]), bool), "units", "pool"
    else:                                                       # an event on its own: its cells and the control cells, its own record
        m, scope, key = ~c["treated"] | (c["unit"] == events[0]), "own", events[0]
    x, y, o, tr, un, first = c[stage][m], c["y"][m], c["o"][m], c["treated"][m], c["unit"][m], c["first"][m]
    out, inside, fits = cycle.cross_fit(x, y, o, tr, un, events, first)
    ref = draws[f"{key}|shift|{stage}|{scope}|{by[key][stage][scope]['class']}"]
    refs = {s: draws[f"{s}|shift|{stage}|units|{by[s][stage]['units']['class']}"] for s in events} if pooled else None
    e = cycle.read(y, o, tr, un, events, out, inside, ref, refs, seed=SEED)
    k = f"{key}|shift|system|frozen|system"
    assert abs(e["as_trained"]["u"] - draws[k][0]) < 2e-4, "the resamples of the round and of the tool do not pair"
    cl = un if len(events) > 1 else None
    U = {a: cycle.output_u(y, p, tr, cl, SEED)[0] for a, p in out.items()}
    e["fits"] = fits
    e["reader_of_the_stage"] = {sc: {q: by[key][stage][sc][q] for q in ("u", "ci", "verdict", "class")} for sc in ("units", "own") if sc in by[key][stage]}
    return e, U, np.asarray(ref, np.float64), (key, scope)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results", help="another folder for a trial on made-up values")
    a = ap.parse_args()
    out_file = a.results / "round.json"
    if out_file.exists():
        raise SystemExit("round.json exists: this reading was made once and is not made again")
    rec = json.loads((RUNS / "outer_step_20261122" / "results" / "bridge.json").read_text())["targets"]
    rec3 = json.loads((RUNS / "input_beside_z_20261123" / "results" / "round.json").read_text())["targets"]
    res = {"margin": cycle.MARGIN, "step_margin": STEP, "note": "exploratory; development cells; the same cells, halves and resamples with the head reading each stage", "targets": {}}
    keep, worst, worst3 = {}, 0.0, 0.0
    for t, (g, study, readout, events, pooled, stages) in TARGETS.items():      # first the two stages that have a record, and the checks
        c = cells(a.work, g, readout)
        by, draws = cycle.load(a.observe / f"{t}_{g}_{readout}")
        E, U, REF = {}, {}, {}
        for st in ("Z", "ZV"):
            E[st], U[st], REF[st], (key, scope) = one(c, events, pooled, by, draws, st)
        r = rec[RECORDED[t]]
        d = max(abs(E["Z"]["arms"][arm]["u"] - r["arms"][arm]["u"]) for arm in ARMS)
        d3 = max(abs(E["ZV"]["arms"][arm]["u"] - rec3[t]["ZV"]["arms"][arm]["u"]) for arm in ARMS)
        worst, worst3 = max(worst, d), max(worst3, d3)
        keep[t] = (c, by, draws, E, U, REF, key, scope, d, d3)
        print(f"{t}: rows at Z {E['Z']['arms']['rows']['u']:.4f} (recorded {r['arms']['rows']['u']:.4f}), at ZV {E['ZV']['arms']['rows']['u']:.4f} (recorded {rec3[t]['ZV']['arms']['rows']['u']:.4f}); largest difference of an arm {d:.5f}, {d3:.5f}")
    res["checks_against_the_records"] = {"largest_difference_of_an_arm_at_Z": worst, "largest_difference_of_an_arm_at_ZV": worst3, "by_target": {t: [v[8], v[9]] for t, v in keep.items()}, "limit": 0.005}
    if max(worst, worst3) > 0.005:
        raise SystemExit(f"the readings at Z or at ZV do not repeat the recorded ones ({worst:.4f}, {worst3:.4f}): stopped, nothing written (DESIGN.md, Section 7)")
    for t, (g, study, readout, events, pooled, stages) in TARGETS.items():
        c, by, draws, E, U, REF, key, scope, _d, _d3 = keep[t]
        for st in ("ZC", "C"):
            E[st], U[st], REF[st], _k = one(c, events, pooled, by, draws, st)
        x = {"stratum": g, "readout": readout, "events": events, "read_as": "a pool" if pooled else "on its own, with its own record", "reference_scope": scope, **{st: E[st] for st in ("Z", "ZC", "C", "ZV")},
             "gain": {w: {arm: gain(U["Z"][arm], U[w][arm]) for arm in ARMS} for w in WIDE}, "gain_of_the_reference_reader": {w: gain(REF["Z"], REF[w]) for w in WIDE},
             "annotation_position_against_input_view": {**{arm: cycle.diff(U["ZC"][arm], U["ZV"][arm]) for arm in ARMS}, "reference_reader": cycle.diff(REF["ZC"], REF["ZV"])},
             "grade": {st: E[st]["closure"]["grade"] for st in ("Z", "ZC", "C", "ZV")}}
        x["grade_kept"] = {w: bool(RANK.get(x["grade"][w], 0) >= RANK.get(x["grade"]["Z"], 0)) for w in WIDE}
        if "V" in stages:                                       # the 512 components alone, in the same run and scope
            v = draws[f"{key}|shift|V|{scope}|{by[key]['V'][scope]['class']}"]
            x["input_view_alone"] = {**{q: by[key]["V"][scope][q] for q in ("u", "ci", "verdict", "class")}, **{f"from_the_input_view_to_{st}": cycle.diff(REF[st], v) for st in ("Z", "ZC", "C", "ZV")}}
        res["targets"][t] = x

    def hyp(t, w):
        x = res["targets"][t]
        closed = x[w]["closure"]["grade"] in ("relative", "absolute")
        r = x["gain"][w]["rows"]["reading"]
        return "confirmed" if closed and r == "a step" else ("refuted" if r == "no step" else "unresolved")
    res["hypotheses"] = {"H1": hyp("T1", "ZC"), "H2": hyp("T2", "ZC"), **{t: "holds" if res["targets"][t]["gain"]["ZC"]["rows"]["no_harm"] else "fails or unresolved" for t in ("C1", "C2", "C3")},
                         "grades_kept": {t: res["targets"][t]["grade_kept"]["ZC"] for t in TARGETS}}
    res["side"] = {"annotation_position_alone": {"H1": hyp("T1", "C"), "H2": hyp("T2", "C"), **{t: "holds" if res["targets"][t]["gain"]["C"]["rows"]["no_harm"] else "fails or unresolved" for t in ("C1", "C2", "C3")},
                                                 "grades_kept": {t: res["targets"][t]["grade_kept"]["C"] for t in TARGETS}},
                   "input_view_beside_Z_again": {"H1": hyp("T1", "ZV"), "H2": hyp("T2", "ZV")}}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda v: f"{v['u']:.3f} [{v['ci'][0]:.3f}, {v['ci'][1]:.3f}] {v['verdict']}"
    for t, x in res["targets"].items():
        print(f"\n== {t} {x['readout']} {x['stratum']} {x['events']} ({x['read_as']})")
        for st in ("Z", "ZC", "C", "ZV"):
            e = x[st]
            print(f"  head reads {st:2s}: reference {f(e['reference'])} | as trained {f(e['as_trained'])} | A3 {e['A3']['holds']} | gate 3 {e['gate3']['passes']} | grade {e['closure']['grade']}")
            for arm, v in e["arms"].items():
                print(f"       {arm:9s} {f(v)}  beyond reference {v['beyond_reference']['added']:.3f} {v['beyond_reference']['reading']}")
        for w in WIDE:
            print(f"  gain Z -> {w:2s}: " + " | ".join(f"{arm} {v['gain']:+.3f} [{v['ci'][0]:.3f}, {v['ci'][1]:.3f}] {v['reading']}" for arm, v in x["gain"][w].items())
                  + f" | reference {x['gain_of_the_reference_reader'][w]['gain']:+.3f} {x['gain_of_the_reference_reader'][w]['reading']}")
        q = x["annotation_position_against_input_view"]
        print("  ZC against ZV (|u| at ZC minus |u| at ZV): " + " | ".join(f"{k} {v['added']:+.3f} [{v['ci'][0]:.3f}, {v['ci'][1]:.3f}] {v['reading']}" for k, v in q.items()))
        if "input_view_alone" in x:
            i = x["input_view_alone"]
            print(f"  the 512 components alone {f(i)} ({i['class']}) | " + " | ".join(f"to {st} {i['from_the_input_view_to_' + st]['added']:+.3f} {i['from_the_input_view_to_' + st]['reading']}" for st in ("Z", "ZC", "C", "ZV")))
    print("\nhypotheses:", res["hypotheses"])
    print("side:", res["side"])


if __name__ == "__main__":
    main()
