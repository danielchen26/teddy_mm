#!/usr/bin/env python3
"""The round with the head reading Z and with the head reading ZV, on the same cells, halves and resamples.

  python round.py --work WORK --observe OBSERVE        (writes results/round.json; refuses to run when it exists)

WORK: the folder build_zv.py wrote. OBSERVE: the grids of `anm_diag scope` for its SPEC files, one folder per SPEC.
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
from anm_diag import cycle  # noqa: E402
from _common import SEED, unit  # noqa: E402
from build_zv import TARGETS  # noqa: E402

RECORDED = {"T1": "PDL1 thp1 type B", "T2": "MHC1 mel_co", "C1": "PDL1 thp1 type A", "C2": "PDL1 mel_co", "C3": "PDL1 mel_ifn"}
FULL = {"thp1a": "thp1"}                                       # a stratum whose halves are those of a larger one
STEP = cycle.STEP
step = lambda iv: "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")


def gain(a, b):
    """|a| - |b| on paired resamples: what the wider stage removes."""
    d = np.abs(a) - np.abs(b)
    iv = cycle.interval(d[1:])
    return {"gain": float(d[0]), "ci": iv, "reading": step(iv), "above_zero": bool(iv[0] > 0), "no_harm": bool(iv[0] > -STEP)}


def cells(work, g, readout):
    d = np.load(work / f"{FULL.get(g, g)}_zv.npz")
    c = dict(ids=d["ids"].astype(str), Z=unit(d["Z"]).astype(np.float64), ZV=d["ZV"].astype(np.float64), treated=d["A"].astype(int) == 1, unit=d["unit"].astype(str),
             y=d[f"y_{readout}"].astype(np.float64), o=d[f"frozen_{readout}"].astype(np.float64))
    c["first"] = cycle.halves(np.where(c["treated"], c["unit"], "control"), SEED)      # the halves of the earlier rounds, made on the whole stratum
    if g in FULL:
        keep = np.isin(c["ids"], np.load(work / f"{g}_zv.npz")["ids"].astype(str))
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
    a = ap.parse_args()
    out_file = RUN / "results" / "round.json"
    if out_file.exists():
        raise SystemExit("round.json exists: this reading was made once and is not made again")
    rec = json.loads((RUNS / "outer_step_20261122" / "results" / "bridge.json").read_text())["targets"]
    res = {"margin": cycle.MARGIN, "step_margin": STEP, "note": "exploratory; development cells; the same cells, halves and resamples with the head reading Z and reading ZV", "targets": {}}
    keep, worst = {}, 0.0
    for t, (g, study, readout, events, pooled, stages) in TARGETS.items():      # first the readings at Z, and the check against the record
        c = cells(a.work, g, readout)
        by, draws = cycle.load(a.observe / f"{t}_{g}_{readout}")
        eZ, UZ, refZ, (key, scope) = one(c, events, pooled, by, draws, "Z")
        r = rec[RECORDED[t]]
        d = max(abs(eZ["arms"][arm]["u"] - r["arms"][arm]["u"]) for arm in ("controls", "rows", "cells"))
        worst = max(worst, d)
        keep[t] = (c, by, draws, eZ, UZ, refZ, key, scope, d)
        print(f"{t} at Z: rows {eZ['arms']['rows']['u']:.4f} (recorded {r['arms']['rows']['u']:.4f}); reference {eZ['reference']['u']:.4f} (recorded {r['reference']['u']:.4f}); largest difference of an arm {d:.5f}")
    res["check_against_the_record"] = {"largest_difference_of_an_arm_at_Z": worst, "by_target": {t: v[8] for t, v in keep.items()}, "limit": 0.005}
    if worst > 0.005:
        raise SystemExit(f"the readings at Z do not repeat the recorded ones ({worst:.4f}): stopped, nothing written (DESIGN.md, Section 7)")
    for t, (g, study, readout, events, pooled, stages) in TARGETS.items():
        c, by, draws, eZ, UZ, refZ, key, scope, _d = keep[t]
        eV, UV, refV, _k = one(c, events, pooled, by, draws, "ZV")
        x = {"stratum": g, "readout": readout, "events": events, "read_as": "a pool" if pooled else "on its own, with its own record", "reference_scope": scope,
             "Z": eZ, "ZV": eV, "gain": {arm: gain(UZ[arm], UV[arm]) for arm in ("controls", "rows", "cells")}, "gain_of_the_reference_reader": gain(refZ, refV)}
        if "V" in stages:                                       # the 512 components alone, in the same run and scope
            v = draws[f"{key}|shift|V|{scope}|{by[key]['V'][scope]['class']}"]
            x["input_view_alone"] = {**{q: by[key]["V"][scope][q] for q in ("u", "ci", "verdict", "class")}, "from_the_input_view_to_Z": cycle.diff(refZ, v),
                                     "from_the_input_view_to_ZV": cycle.diff(refV, v)}
        res["targets"][t] = x
    def hyp(t):
        x = res["targets"][t]
        closed = x["ZV"]["closure"]["grade"] in ("relative", "absolute")
        r = x["gain"]["rows"]["reading"]
        return "confirmed" if closed and r == "a step" else ("refuted" if r == "no step" else "unresolved")
    res["hypotheses"] = {"H1": hyp("T1"), "H2": hyp("T2"), **{t: "holds" if res["targets"][t]["gain"]["rows"]["no_harm"] else "fails or unresolved" for t in ("C1", "C2", "C3")}}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda v: f"{v['u']:.3f} [{v['ci'][0]:.3f}, {v['ci'][1]:.3f}] {v['verdict']}"
    for t, x in res["targets"].items():
        print(f"\n== {t} {x['readout']} {x['stratum']} {x['events']} ({x['read_as']})")
        for st in ("Z", "ZV"):
            e = x[st]
            print(f"  head reads {st:2s}: reference {f(e['reference'])} | as trained {f(e['as_trained'])} | A3 {e['A3']['holds']} | gate 3 {e['gate3']['passes']} | grade {e['closure']['grade']}")
            for arm, v in e["arms"].items():
                print(f"       {arm:9s} {f(v)}  beyond reference {v['beyond_reference']['added']:.3f} {v['beyond_reference']['reading']}")
        for arm, v in x["gain"].items():
            print(f"  gain of {arm:9s} {v['gain']:.3f} [{v['ci'][0]:.3f}, {v['ci'][1]:.3f}] {v['reading']}")
        v = x["gain_of_the_reference_reader"]
        print(f"  gain of the reference reader {v['gain']:.3f} [{v['ci'][0]:.3f}, {v['ci'][1]:.3f}] {v['reading']}")
        if "input_view_alone" in x:
            i = x["input_view_alone"]
            print(f"  the 512 components alone {f(i)} ({i['class']}) | to Z {i['from_the_input_view_to_Z']['added']:.3f} {i['from_the_input_view_to_Z']['reading']} | to ZV {i['from_the_input_view_to_ZV']['added']:.3f} {i['from_the_input_view_to_ZV']['reading']}")
    print("\nhypotheses:", res["hypotheses"])


if __name__ == "__main__":
    main()
