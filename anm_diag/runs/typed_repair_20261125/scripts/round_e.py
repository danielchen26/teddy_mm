#!/usr/bin/env python3
"""Part E: the last layers trained on 24 November, read again with the repair of this run. results/round_e.json. Made once.

  python round_e.py --work WORK --observe OBSERVE --trained TRAINED [--results RESULTS]

WORK and OBSERVE as for round_a.py (the THP-1 targets). TRAINED: the folder of train.py of
runs/last_layer_trial_20261124 (<arm>_<fold>.npz: the state of every cell under the layer trained on that fold's half).

For each arm the stage is the state under the layer that was trained on the fitted half, for the cells that are fitted
and for the cells that are read alike: a readout is fitted on states of the layer it will be read with, and the cells
that are read were never seen by that layer or by the readout. (The trial's own reading gave each cell the state of the
layer that had not seen it, and fitted its readouts across the two layers.)

Read: IRF1 on its own with its record, by the old rule and by `Full`; IRF1 as an event without a protein record (none,
type, derived); the four pathway genes with their records. Every gain is over the released state Z on paired resamples.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from repair import SEED, STEP, cross, cycle, gain, u_of, unit
from round_a import RUN, TARGETS, cells, r_control
from round_b import sub, unrecorded

ARMS = {"rows4": "the records of the four pathway genes", "rows5": "and the record of IRF1", "mixed": "the event labels shuffled", "cells5": "the protein of every cell", "prop": "a derived record for IRF1: the mean of its type"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--trained", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    out_file = a.results / "round_e.json"
    if out_file.exists():
        raise SystemExit("round_e.json exists: this reading was made once and is not made again")
    res = {"margin": cycle.MARGIN, "step_margin": STEP, "note": "THP-1, development cells; no new training: the five last layers trained on 24 November, read with the repair of this run", "arms": ARMS, "targets": {}}
    worst = 0.0
    for t in ("T1", "C1"):
        g, readout, events, pooled, typ, _x = TARGETS[t]
        c = cells(a.work, g, readout)
        for arm in ARMS:
            F = [np.load(a.trained / f"{arm}_{f}.npz") for f in (0, 1)]
            for f, d in enumerate(F):
                assert np.array_equal(d["ids"].astype(str), c["ids"]), "the trained states are of other cells, or in another order"
                assert np.array_equal(d["fit"].astype(bool), c["first"] if f == 0 else ~c["first"]), "the half a layer was trained on is not the half the round fits on"
                worst = max(worst, float(np.abs(unit(d["Z0"]).astype(np.float64) - c["Z"]).max()))
            c[arm] = tuple(unit(d["ZT"]).astype(np.float64) for d in F)
        by, draws = cycle.load(a.observe / f"{t}_{g}_{readout}")
        if pooled:
            m, scope, key = np.isin(c["unit"], events + ["control"]), "units", "pool"
        else:
            m, scope, key = ~c["treated"] | (c["unit"] == events[0]), "own", events[0]
        s = sub(c, m)
        y, o, tr, un, first = s["y"], s["o"], s["treated"], s["unit"], s["first"]
        cl = un if len(events) > 1 else None
        ref = draws[f"{key}|shift|Z|{scope}|{by[key]['Z'][scope]['class']}"]
        refs = {e: draws[f"{e}|shift|Z|units|{by[e]['Z']['units']['class']}"] for e in events} if pooled else None
        x, U, UU = {"events": events, "stages": {}}, {}, {}
        for st in ["Z"] + list(ARMS):
            out, inside, fits = cross(s[st], y, o, tr, un, first, events)
            e = cycle.read(y, o, tr, un, events, out, inside, ref, refs, seed=SEED)
            U[st] = {arm: u_of(y, out[arm], tr, cl)[0] for arm in out}
            q = {"arms": e["arms"], "gate3": e["gate3"], "closure": e["closure"], "r_control": {arm: r_control(y, out[arm], tr) for arm in out}}
            if st != "Z":
                q["gain_over_Z"] = {arm: gain(U["Z"][arm], U[st][arm]) for arm in U[st]}
            if not pooled:
                q["without_a_protein_record"], UU[st] = unrecorded(c, st, events[0], typ)
                if st != "Z":
                    q["gain_over_Z_without_a_protein_record"] = {arm: gain(UU["Z"][arm], UU[st][arm]) for arm in UU[st]}
            x["stages"][st] = q
            print(f"{t} {st:6s}: " + " | ".join(f"{arm} {e['arms'][arm]['u']:+.3f}" for arm in ("controls", "rows", "cells")) + f" | {e['closure']['grade']}"
                  + ("" if pooled else " || " + " ".join(f"{k} {v['u']:+.2f}" for k, v in q["without_a_protein_record"].items())), flush=True)
        res["targets"][t] = x
    res["check"] = {"the_released_state_in_the_files_of_the_training_against_the_state_of_the_round": worst, "limit": 2e-3, "hold": bool(worst <= 2e-3)}
    if worst > 2e-3:
        raise SystemExit(f"the released state in the trained files is not the state of the round ({worst}): stopped, nothing written")
    T1 = res["targets"]["T1"]["stages"]
    steps = {arm: T1[arm]["gain_over_Z"]["rows"]["reading"] for arm in ("rows4", "rows5")}
    res["hypotheses"] = {"HE_a_trained_last_layer_adds_a_step_under_the_new_repair": {
        "reading": "confirmed" if "a step" in steps.values() else ("refuted" if all(v == "no step" for v in steps.values()) else "unresolved"), "IRF1_with_its_record_gain_over_Z": {arm: T1[arm]["gain_over_Z"]["rows"] for arm in ARMS}}}
    print(json.dumps({k: v["reading"] for k, v in res["hypotheses"].items()}))
    a.results.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("written", out_file.name, res["check"])


if __name__ == "__main__":
    main()
