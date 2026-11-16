#!/usr/bin/env python3
"""Gate 0c: the forecasts still to be read, from existing arrays. No head is fitted and no fitted head is read.

For every fit set of the design (four PBMC donors), a ridge readout of TEDDY's state (the stage the head reads) is fitted
  pbmc     on the 20,000 PBMC cells of the fit donors, and
  pooled   on those cells together with the 47,405 bone-marrow cells of seed 0's `coverage` arm, without a study label,
and read on the four donors it did not see. u is that of the round (variation target; scope.arms), the base is the head
in use. The gain of `pbmc` is the ceiling the protocol speaks of; `pooled` is what one map on the two studies' scales
can deliver. M1 and M2: the batch of the scored donors was seen; B1 and B2: it was not.

  python forecasts.py --bm-arrays cite_arrays.npz --bm-z z_rna.npy --pbmc-arrays ext_cite_arrays.npz --pbmc-z ext_z_rna.npy \\
      --head-weights head_weights.npz

Writes ../results/forecasts.json.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import common
import plan
import train_arm
from anm_diag import scope

DECISION = {"bootstrap": 1000, "seed": 0, "alpha": 0.05}


def ridge(X, Y):
    return make_pipeline(StandardScaler(), RidgeCV(alphas=np.asarray(common.ALPHAS))).fit(X, Y)


def main():
    p = argparse.ArgumentParser()
    for k in ("bm-arrays", "bm-z", "pbmc-arrays", "pbmc-z", "head-weights"):
        p.add_argument(f"--{k}", type=Path, required=True)
    a = p.parse_args()
    bm, pb, shared = common.load_studies(a.bm_arrays, a.bm_z, a.pbmc_arrays, a.pbmc_z)
    Y = pb["adt"][:, shared].astype(np.float64)
    head = common.head_forward(dict(np.load(a.head_weights, allow_pickle=True)), pb["z"])[:, shared]
    n_train = int((bm["split"] == "train").sum())
    (_, _, _, _), rec = train_arm.build_train({"seed": 0, "n_bm": n_train - plan.N_ADDED, "pbmc": {}}, bm, None)
    Xb, Yb = bm["z"][rec["bm_index"]], bm["adt"][rec["bm_index"]][:, shared]

    pred = {k: np.full(Y.shape, np.nan) for k in ("pbmc_seen", "pooled_seen", "pbmc_not_seen", "pooled_not_seen", "other_seven")}
    for fs in plan.FIT_SETS:
        fit = np.isin(pb["donors"], plan.fit_donors(fs))
        sc = np.isin(pb["donors"], plan.scored_donors(fs))
        tag = "seen" if plan.BATCH_SEEN[fs] else "not_seen"
        pred[f"pbmc_{tag}"][sc] = ridge(pb["z"][fit], Y[fit]).predict(pb["z"][sc])
        pred[f"pooled_{tag}"][sc] = ridge(np.vstack([Xb, pb["z"][fit]]), np.vstack([Yb, Y[fit]])).predict(pb["z"][sc])
        print(fs, tag, "fitted", int(fit.sum()), "scored", int(sc.sum()), flush=True)
    for d in plan.DONORS:                                           # as runs/scope_20261116 read it: the other seven donors
        me = pb["donors"] == d
        pred["other_seven"][me] = ridge(pb["z"][~me], Y[~me]).predict(pb["z"][me])
    assert not any(np.isnan(v).any() for v in pred.values())

    per_donor, U = {}, {k: [] for k in ["head"] + list(pred)}
    for d in plan.DONORS:
        me = pb["donors"] == d
        r = scope.arms(Y[me], {"head": head[me], **{k: v[me] for k, v in pred.items()}}, base="head", decision=DECISION)
        per_donor[d] = {k: {"u": x["u"], "ci": x["ci"], "gain": x["gain"], "gain_ci": x["gain_ci"]} for k, x in r["arms"].items()}
        per_donor[d]["score_head"] = common.score(Y[me], head[me])
        for k in U:
            U[k].append([r["arms"][k]["u"]])
    fam = scope.arms_family({k: np.array(v) for k, v in U.items()}, base="head")
    out = {"what": "ceilings read by ridge readouts of the state the head reads; base: the head in use; no fitted head of the round is read",
           "quantity": "fall of u (variation target, frozen scope) from the head in use, mean over the 8 donors; interval from resampling donors",
           "proteins": int(len(shared)), "decision": DECISION, "approximate": fam["approximate"],
           "forecast": {k: {"gain": v["gain"], "ci": v["ci"], "u": v["u"]} for k, v in fam["arms"].items()},
           "for": {"D1 (coverage on M1, M2)": "pbmc_seen (ceiling), pooled_seen", "D2 (coverage on B1, B2)": "pbmc_not_seen (ceiling), pooled_not_seen",
                   "check": "other_seven repeats the units scope of runs/scope_20261116 without its cap of 30,000 fit cells"},
           "per_donor": per_donor}
    (common.RUN / "results").mkdir(exist_ok=True)
    (common.RUN / "results" / "forecasts.json").write_text(json.dumps(out, indent=1))
    for k, v in out["forecast"].items():
        print(f"{k:16s} u {v['u']:.4f}  gain {v['gain']:+.4f} [{v['ci'][0]:+.4f}, {v['ci'][1]:+.4f}]")


if __name__ == "__main__":
    main()
