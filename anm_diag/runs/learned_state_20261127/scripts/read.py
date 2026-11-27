#!/usr/bin/env python3
"""The reading of the learned pooling: the loop judges it, on cells the pooling never saw. results/reading.json,
results/outcome.json. Committed before any trained state exists; the reading is made once.

  python read.py --enc THP1_ENC.npz --trained OUT_DIR --states-thp THP_STATES... [--results DIR] [--rehearsal]

OUT_DIR: the files of `learned_state_design_20261127/scripts/pool.py`, <arm>_<fold>.npz for the arms none, pool4,
pool5, mixed and the folds 0, 1. A cell is read at the state the pooling gives it that was trained on the OTHER half
of the cells (`account.reading` takes a stage as the pair of its two versions).

For IRF1 on PD-L1, `account.account` at: V the input view; Z the mean, as released; T the state typed by the declared
gene sets (run 40); P0, P4, P5, PM the state under the pooling of the arms none, pool4, pool5, mixed. The other four
pathway events are the events whose rows may teach a readout. A step is 0.10, as everywhere.

H1  the trained state holds more: `carried` at P5 over Z gains at least a step and its interval lies above zero.
    Refuted: the interval lies below 0.10.
H2  the gain belongs to the rows: `carried` at P5 over PM, the same rule.
H3  an event with no record gains: `carried` at P4 over Z, the same rule.
Guards: (1) for each of the four pathway events, `carried` at P5 is not a step below Z (the interval of the gain does
not lie below -0.10); (2) among control cells of the half a pooling never saw, the cross-validated error of a ridge
readout at P5 is at most 1.05 times that at Z; (3) the arm none moved nothing.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
sys.path.insert(0, str(RUNS / "final_audit_20261127" / "scripts"))
sys.path.insert(0, str(RUNS / "last_layer_trial_20261124" / "scripts"))
import _fa as F  # noqa: E402
from train import ALL5, TYPE_A, ridge  # noqa: E402
from anm_diag import account as AC  # noqa: E402

STEP, K, ARMS, EXPECTED = 0.10, 4, ("none", "pool4", "pool5", "mixed"), {"H1": 0.35, "H2": 0.25, "H3": 0.2, "the_guards_all_hold": 0.7}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def flat(B):
    """The state the pooling was trained under: its blocks, each at unit length, side by side."""
    return (B / (np.linalg.norm(B, axis=2, keepdims=True) + 1e-6)).reshape(len(B), -1).astype(np.float64) / B.shape[1] ** 0.5


def verdict(g):
    return "holds" if g["gain"] >= STEP and g["ci"][0] > 0 else "refuted" if g["ci"][1] < STEP else "unresolved"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enc", type=Path, required=True)
    ap.add_argument("--trained", type=Path, required=True)
    ap.add_argument("--states-thp", type=Path, nargs="+", required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    ap.add_argument("--rehearsal", action="store_true")
    a = ap.parse_args()
    assert not a.rehearsal or a.results.resolve() != (RUN / "results").resolve(), "a rehearsal does not write into the results of the run"
    out_file = a.results / "reading.json"
    if out_file.exists():
        raise SystemExit("reading.json exists: the reading was made once and is not made again")
    E = np.load(a.enc)
    ids, A = E["ids"].astype(str), E["A"].astype(int)
    unit = np.where(A == 1, E["unit"].astype(str), "control")
    y = E["y_PDL1"].astype(np.float64)
    P, first, logs = {}, None, {}
    for arm in ARMS:
        two = []
        for fold in (0, 1):
            d = np.load(a.trained / f"{arm}_{fold}.npz")
            assert np.array_equal(d["ids"].astype(str), ids) and d["PT"].shape == (len(ids), K, E["L12"].shape[1])
            fit = d["fit"].astype(bool)
            first = fit if first is None and fold == 0 else first
            assert np.array_equal(fit, first if fold == 0 else ~first), "the halves of the files differ"
            assert a.rehearsal or float(np.abs(d["Z0"] - E["L12"]).max()) < 5e-3, "the mean of a file is not the recorded state"
            two.append(flat(d["PT"]))
            logs[f"{arm}_{fold}"] = {k: v for k, v in json.loads(str(d["log"])).items() if k != "curve"} | {"curve_first_and_last": [json.loads(str(d["log"]))["curve"][0], json.loads(str(d["log"]))["curve"][-1]], "file_sha256": sha(a.trained / f"{arm}_{fold}.npz")}
            if arm == "none":
                logs[f"{arm}_{fold}"]["largest_move_of_a_block_from_the_mean"] = float(np.abs(d["PT"] - d["Z0"][:, None]).max())
        P[arm] = tuple(two)
    st = F.States(a.states_thp)
    raw, blocks, sets = st.rows(ids), st.rows(ids, "blocks"), st.sets
    assert float(np.abs(raw - E["L12"]).max()) < 5e-3, "the states file holds other states than the recorded ones"
    nz = np.linalg.norm(raw, axis=1, keepdims=True) + 1e-6
    blk = lambda name, f: blocks[:, sets.index(f"{name}|{f}")]
    T = tuple(np.hstack([blk("A", f), raw - blk("A", f)]) / nz for f in (0, 1))
    Z = raw / nz
    stages = {"V": E["V"].astype(np.float64), "Z": Z, "T": T, "P0": P["none"], "P4": P["pool4"], "P5": P["pool5"], "PM": P["mixed"]}
    pairs = (("Z", "P5"), ("PM", "P5"), ("Z", "P4"), ("Z", "PM"), ("Z", "P0"), ("P0", "P5"), ("Z", "T"), ("T", "P5"), ("Z", "V"), ("P5", "V"))
    acc = AC.account(stages, y, unit, first, "IRF1", others=TYPE_A, inputs=("V",), state="Z", pairs=pairs, seed=0)
    res = {"rehearsal_on_made_up_states": bool(a.rehearsal), "cells": int(len(ids)), "first_half": int(first.sum()), "trainings": logs, "IRF1": acc, "guards": {"pathway_events": {}}}
    print("IRF1 carried:", {s: round(acc["stages"][s]["shares"]["carried"]["share"], 3) for s in stages}, flush=True)
    for g in TYPE_A:
        o = AC.account({"Z": Z, "P5": P["pool5"]}, y, unit, first, g, others=[h for h in ALL5 if h != g], state="Z", pairs=(("Z", "P5"),), seed=0, controls=False)
        res["guards"]["pathway_events"][g] = {"carried": {s: o["stages"][s]["shares"]["carried"]["share"] for s in ("Z", "P5")}, "P5_over_Z": o["gains"]["P5 over Z"]["carried"]}
        print(g, res["guards"]["pathway_events"][g], flush=True)
    cv = {}
    for f, fit in ((0, first), (1, ~first)):                         # control cells of the half the pooling of fold f never saw
        m = ~fit & (unit == "control")
        t = (y[m] - y[m].mean()) / y[m].std()
        cv[f] = {"Z": ridge(Z[m], t)[4], "P5": ridge(P["pool5"][f][m], t)[4]}
    res["guards"]["cv_error_of_a_control_readout"] = {"by_fold": cv, "ratio_P5_to_Z": float(np.mean([cv[f]["P5"] for f in cv]) / np.mean([cv[f]["Z"] for f in cv]))}
    a.results.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    # --- the hypotheses
    G = {k: v["carried"] for k, v in acc["gains"].items()}
    H = {"H1_the_trained_state_holds_more": dict(G["P5 over Z"], outcome=verdict(G["P5 over Z"])),
         "H2_the_gain_belongs_to_the_rows": dict(G["P5 over PM"], outcome=verdict(G["P5 over PM"])),
         "H3_an_event_with_no_record_gains": dict(G["P4 over Z"], outcome=verdict(G["P4 over Z"]))}
    g1 = all(not v["P5_over_Z"]["ci"][1] < -STEP for v in res["guards"]["pathway_events"].values())
    g2 = res["guards"]["cv_error_of_a_control_readout"]["ratio_P5_to_Z"] <= 1.05
    g3 = all(v.get("largest_move_of_a_block_from_the_mean", 0.0) < 1e-4 for k, v in logs.items() if k.startswith("none"))
    out = {"rehearsal_on_made_up_states": bool(a.rehearsal), "carried": {s: acc["stages"][s]["shares"]["carried"] for s in stages}, "gains_in_carried": G, "hypotheses": H,
           "guards": {"the_pathway_events_are_read_no_worse_by_a_step": g1, "the_control_readout_is_no_worse_by_5_percent": g2, "the_arm_none_moved_nothing": g3, "all_hold": bool(g1 and g2 and g3)},
           "what_the_account_says": acc["chain"]["says"], "expected_before_any_trained_state": EXPECTED}
    (a.results / "outcome.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for k, v in H.items():
        print(f"  {k}: gain {v['gain']:+.3f} [{v['ci'][0]:+.3f}, {v['ci'][1]:+.3f}] -> {v['outcome']}")
    print("  guards:", out["guards"])


if __name__ == "__main__":
    main()
