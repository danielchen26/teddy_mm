#!/usr/bin/env python3
"""The trained states as stages, cross-fitted by half, and the SPEC files of the reading. No readout enters.

  python build_tt.py --enc WORK_ENC --states OUT --work WORK [--results RESULTS] [--arms rows4,rows5]

WORK_ENC: the work folder of runs/encoder_stages_20261124 (thp1_enc.npz and the two SPEC files). OUT: the folder of
train.py, two files per arm (fold 0 fitted on the first half of the cells, fold 1 on the second).

**Every cell is given the state of the layer that was trained without it**: a cell of the first half takes its state
from fold 1, a cell of the second half from fold 0. So a stage here is a state of cells the trained layer never saw,
neither their tokens nor their protein. The halves are those of the earlier rounds (anm_diag.cycle.halves, seed 0).

Stages: V (the input view), Z70 (the state as released), T4, T5, TM, TC, TP (the arms rows4, rows5, mixed, cells5,
prop; with --arms only those named, for the first part of the reading), PERM. Writes WORK/thp1_tt.npz, thp1a_tt.npz, the two SPEC files, and RESULTS/trial_manifest.json: the training
logs (they hold the rows of the fitted halves) and, with no other protein value, how far the states moved.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import yaml

TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
ARMS = {"rows4": "T4", "rows5": "T5", "mixed": "TM", "cells5": "TC", "prop": "TP"}


def halves(groups, seed=0, share=0.5):
    """anm_diag.cycle.halves, copied."""
    rng, first = np.random.default_rng(seed), np.zeros(len(groups), bool)
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        first[idx[: int(round(share * len(idx)))]] = True
    return first


def unit(z):
    z = z.astype(np.float64)
    return z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enc", type=Path, required=True)
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=Path(__file__).resolve().parent.parent / "results")
    ap.add_argument("--arms", default=",".join(ARMS))
    a = ap.parse_args()
    arms = {k: ARMS[k] for k in a.arms.split(",")}
    STAGES = ["V", "Z70"] + list(arms.values()) + ["PERM"]
    a.work.mkdir(parents=True, exist_ok=True)
    K = np.load(a.enc / "thp1_enc.npz")
    ids, A, un = K["ids"].astype(str), K["A"].astype(int), K["unit"].astype(str)
    first = halves(np.where(A == 1, un, "control"), 0)
    ctl = A == 0
    S = {"V": K["V"].astype(np.float32), "Z70": K["L12"].astype(np.float32)}
    info = {"cells": int(len(ids)), "first_half": int(first.sum()), "arms": {}}
    for arm, st in arms.items():
        F = [np.load(a.states / f"{arm}_{f}.npz") for f in (0, 1)]
        logs = [json.loads(str(f["log"])) for f in F]
        for f, X in enumerate(F):
            assert np.array_equal(X["ids"].astype(str), ids), f"{arm}, fold {f}: other cells"
            assert np.array_equal(X["fit"], first if f == 0 else ~first), f"{arm}, fold {f}: the halves are not those of the rounds"
            assert float(np.abs(X["Z0"] - K["L12"]).max()) <= 5e-3, f"{arm}, fold {f}: the released state is not the recorded one"
        S[st] = np.where(first[:, None], F[1]["ZT"], F[0]["ZT"]).astype(np.float32)   # the state from the fold that did not hold the cell
        x0, xt = unit(K["L12"]), unit(S[st])
        v0 = float(((x0[ctl] - x0[ctl].mean(0)) ** 2).sum(1).mean())
        moved = ((xt - x0) ** 2).sum(1) / v0                      # in units of the spread of the control cells' states
        shift = {}                                                # the shift of each fold's own readout on the cells it did not hold: no protein value enters
        for f, X in enumerate(F):
            other = ~X["fit"]
            r = (unit(X["ZT"]) - X["xc"]) @ X["w"] + float(X["b"])
            r0 = (unit(X["Z0"]) - X["xc"]) @ X["w"] + float(X["b"])
            for g in sorted(set(un[A == 1])):
                m_ = other & (A == 1) & (un == g)
                shift.setdefault(g, []).append({"fold": f, "released": float(r0[m_].mean() - r0[other & ctl].mean()), "trained": float(r[m_].mean() - r[other & ctl].mean())})
        cos = float(F[0]["w"] @ F[1]["w"] / (np.linalg.norm(F[0]["w"]) * np.linalg.norm(F[1]["w"])))
        info["arms"][arm] = {"stage": st, "folds": [{k: L[k] for k in ("fold", "steps", "stopped_by", "rows", "events", "seconds", "largest_move_of_a_state")} | {"at_the_stop": L["curve"][-1], "at_the_start": L["curve"][0]} for L in logs],
                             "moved_cells_not_held": {"control_cells_mean": float(moved[ctl].mean()), **{g: float(moved[(A == 1) & (un == g)].mean()) for g in sorted(set(un[A == 1]))}},
                             "shift_of_the_folds_readout_on_cells_not_held": shift, "cosine_of_the_two_folds_readouts": cos}
    S["PERM"] = K["PERM"].astype(np.float32)
    assert list(S) == STAGES, list(S)
    rest = {k: K[k] for k in ("ids", "A", "unit", "y_PDL1", "frozen_PDL1")}
    np.savez(a.work / "thp1_tt.npz", **S, **rest)
    msk = np.isin(un, TYPE_A + ["control"])
    np.savez(a.work / "thp1a_tt.npz", **{k: v[msk] for k, v in S.items()}, **{k: v[msk] for k, v in rest.items()})
    for src, arrays in (("T_thp1_IRF1", "thp1_tt.npz"), ("C_thp1a_pool", "thp1a_tt.npz")):
        s = yaml.safe_load((a.enc / f"{src}.yaml").read_text())  # the same events, readers, decision and keys as the runs before; only the stages change
        s["keys"]["stages"] = {k: k for k in STAGES}
        s["reads"] = "Z70"
        s["arrays"] = {"audit": arrays}
        (a.work / f"{src}.yaml").write_text(yaml.safe_dump(s, sort_keys=False))
    a.results.mkdir(parents=True, exist_ok=True)
    (a.results / "trial_manifest.json").write_text(json.dumps({**info, "stages": STAGES}, indent=1) + "\n")
    print(json.dumps({arm: {"steps": [f["steps"] for f in v["folds"]], "stopped_by": [f["stopped_by"] for f in v["folds"]], "control_cells_moved": round(v["moved_cells_not_held"]["control_cells_mean"], 4)} for arm, v in info["arms"].items()}, indent=1))


if __name__ == "__main__":
    main()
