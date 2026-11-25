#!/usr/bin/env python3
"""A check added after the reading of part A. It was not declared, and nothing is judged by it.

  python seam_check.py --work WORK [--results RESULTS]

Why. In the trial a cell takes its state from the fold that did not hold it, so the two halves of the cells have
their states from two different trained layers. The round and the tool fit a reader on one half and read the other:
a reader meets a seam between two layers. The reading of part A had the reference reader of the trained state
unchanged (0.40) while the round's arm fitted on every cell read a step more (0.25). This asks whether the seam is
why.

What. For IRF1 on its own (its cells and the control cells) and for the four pathway genes as a pool, inside each
half separately, where every state comes from one layer and no cell was seen by it: the half is split in two again
(anm_diag.cycle.halves, the seed of the rounds), a reader is fitted on one part and read on the other, and u is
computed as the tool computes it for an output. Two readers of a stage:
  tool    the tool's reference reader for the class ridge: standardised state, RidgeCV with the weights of the SPEC
  round   anm_diag.cycle.LinearRows, the arms of the round (controls, rows, cells), at unit length
Resamples are paired between stages (the same seed, the same cells). Half the cells, so wider intervals.

What it also shows. In the trial's own reading, the readout that reads a half is fitted on the control cells of the
other half: the cells the layer was trained with, and so nearly the readout the layer was trained under. Here the
readout is fitted on other control cells. The difference between the two is how much of a gain belongs to that one
readout.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml

RUN = Path(__file__).resolve().parent.parent
RUNS = RUN.parent
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
from anm_diag import cycle  # noqa: E402
from _common import SEED, unit  # noqa: E402
from sklearn.linear_model import RidgeCV  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

STEP = cycle.STEP
POOL = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
READERS = ("tool", "controls", "rows", "cells")
step = lambda iv: "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")
flat = lambda v: np.asarray(v, np.float64).ravel()


def one_half(d, stages, alphas, m, A, un, y, o, events):
    """One target inside one half: every reader at every stage, and the gain of every trained stage over Z70."""
    tr, uu, yy, oo = A[m], un[m], y[m], o[m]
    cl = uu if len(events) > 1 else None
    sub = cycle.halves(np.where(tr, uu, "control"), SEED)
    D = {}
    R = {"cells": int(m.sum()), "control_cells": int((~tr).sum()), "knockout_cells": int(tr.sum()),
         "the_layer_these_states_come_from": "the fold trained on the other half", "stages": {}}
    for st in stages:
        raw, x = d[st][m].astype(np.float64), unit(d[st][m]).astype(np.float64)
        pred = np.zeros(len(yy))
        for fit, read in ((sub, ~sub), (~sub, sub)):
            pred[read] = make_pipeline(StandardScaler(), RidgeCV(alphas=np.asarray(alphas, float))).fit(raw[fit], yy[fit]).predict(raw[read])
        out, _inside, _fits = cycle.cross_fit(x, yy, oo, tr, uu, events, sub)
        D[st] = {"tool": flat(cycle.output_u(yy, pred, tr, cl, SEED)[0]), **{arm: flat(cycle.output_u(yy, p, tr, cl, SEED)[0]) for arm, p in out.items()}}
        R["stages"][st] = {k: {"u": float(v[0]), "ci": cycle.interval(v[1:])} for k, v in D[st].items()}
    for st in stages[1:]:
        R["stages"][st]["gain_from_Z"] = {}
        for k in D[st]:
            g = np.abs(D["Z70"][k]) - np.abs(D[st][k])
            iv = cycle.interval(g[1:])
            R["stages"][st]["gain_from_Z"][k] = {"gain": float(g[0]), "ci": iv, "reading": step(iv)}
    return R


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    a = ap.parse_args()
    d = np.load(a.work / "thp1_tt.npz")
    alphas = yaml.safe_load((a.work / "T_thp1_IRF1.yaml").read_text())["readout_classes"]["ridge"]["alphas"]
    A, un, y, o = d["A"].astype(int) == 1, d["unit"].astype(str), d["y_PDL1"].astype(np.float64), d["frozen_PDL1"].astype(np.float64)
    first = cycle.halves(np.where(A, un, "control"), SEED)       # the halves of the rounds, which are the folds of the trial
    stages = [s for s in ("Z70", "T4", "T5", "TM", "TC", "TP") if s in d.files]
    res = {"what": "not declared, added after the reading of part A; nothing is judged by it. Read inside each half of the cells, where every state comes from one trained layer that did not see the cell, with readers fitted on cells of that half",
           "step_margin": STEP, "ridge_weights_of_the_tool": alphas, "targets": {}}
    for target, events in (("IRF1 on its own", ["IRF1"]), ("the four pathway genes, a pool", POOL)):
        res["targets"][target] = {"halves": {}}
        print(f"==== {target}")
        for name, half in (("first", first), ("second", ~first)):
            R = one_half(d, stages, alphas, (~A | np.isin(un, events)) & half, A, un, y, o, events)
            res["targets"][target]["halves"][name] = R
            print(f"== the {name} half: {R['cells']} cells ({R['control_cells']} control, {R['knockout_cells']} knockout)")
            for st in stages:
                s = R["stages"][st]
                print(f"  {st:<4}" + " | ".join(f"{k} {s[k]['u']:.3f} [{s[k]['ci'][0]:.3f}, {s[k]['ci'][1]:.3f}]" for k in READERS))
                if st != "Z70":
                    print("        gain from Z: " + " | ".join(f"{k} {g['gain']:+.3f} [{g['ci'][0]:.3f}, {g['ci'][1]:.3f}] {g['reading']}" for k, g in s["gain_from_Z"].items()))
    (a.results / "seam_check.json").write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
