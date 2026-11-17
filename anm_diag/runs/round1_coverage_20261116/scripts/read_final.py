#!/usr/bin/env python3
"""Gate 5: the heads of the test, read once on the final study F. Nothing is fitted here.

  python read_final.py --heads HEADS_DIR --bm-arrays cite_arrays.npz --bm-z z_rna.npy --pbmc-arrays ext_cite_arrays.npz \\
      --pbmc-z ext_z_rna.npy --final-arrays F_cite_arrays.npz --final-z F_z_rna.npy --head-weights head_weights.npz

Refuses to run unless ../results/gate5_declaration.json is there and every file it names by hash is the file given:
the pack and the states of F, the weights of the fifteen heads, and the scripts that do the reading. What is read, and
how the outcome is stated, is fixed in ../GATE5.md (part 1) and in the declaration.

Every head is run from its weights in numpy, as scripts/read.py does. For a donor of F, u of every arm on the donor's
cells and the proteins of the declaration (`scope.arms`); over donors, gains with intervals (`scope.arms_family`); the
two hypotheses with `loop.test`, over the three seeds together and for each seed. Writes ../results/final.json and
appends its events to ../results/ledger_final.jsonl.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import common
import plan
import read
import train_arm
from anm_diag import loop, scope

DECISION = {"bootstrap": 1000, "seed": 0, "alpha": 0.05}
ARMS = ("as_is", "coverage", "labels", "amount", "function")
HEAD_OF = {"as_is": "as_is", "coverage": "coverage_all", "labels": "labels_all", "amount": "amount", "function": "function"}
RUN = "round1_coverage_20261116/final"
SCRIPTS = ("read_final.py", "read.py", "common.py", "plan.py", "train_arm.py")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def check_declaration(decl, files):
    """files: {key: path}. Every key of the declaration's hashes has to be given and to match."""
    bad = [k for k, h in decl["sha256"].items() if k not in files or sha256(files[k]) != h]
    if bad:
        raise SystemExit("not the files of the declaration: " + ", ".join(sorted(bad)))


def outcome_of_round(together, by_seed):
    """GATE5.md, part 1: the outcome over the three seeds together, stated only if each seed alone reads the same."""
    o = together["outcome"]
    return o if all(t["outcome"] == o for t in by_seed) else "unresolved"


def group_gain(Ub, Ua, rows):
    """The gain of an arm over the donors `rows`, as scope.arms_family computes it."""
    r = scope.arms_family({"base": Ub[rows], "arm": Ua[rows]}, base="base", bootstrap=DECISION["bootstrap"], seed=DECISION["seed"],
                          alpha=DECISION["alpha"])["arms"]["arm"]
    return {"gain": r["gain"], "ci": r["ci"], "donors": int(len(rows))}


def group_difference(Ub, Ua, first, second, bootstrap=1000, seed=0, alpha=0.05):
    """Gain among the donors `first` minus gain among the donors `second`. Donors are resampled inside each group and
    the seeds inside each arm, as scope.arms_family resamples them."""
    Ub, Ua = np.asarray(Ub, np.float64), np.asarray(Ua, np.float64)
    rng = np.random.default_rng(seed + 3)
    d = np.empty(bootstrap)
    for b in range(bootstrap):
        g = Ub[:, rng.integers(0, Ub.shape[1], Ub.shape[1])].mean(axis=1) - Ua[:, rng.integers(0, Ua.shape[1], Ua.shape[1])].mean(axis=1)
        d[b] = g[rng.choice(first, len(first))].mean() - g[rng.choice(second, len(second))].mean()
    g = Ub.mean(axis=1) - Ua.mean(axis=1)
    return {"gain": float(g[first].mean() - g[second].mean()), "ci": [float(np.quantile(d, alpha / 2)), float(np.quantile(d, 1 - alpha / 2))]}


def load_final(arrays, z, names):
    e = np.load(arrays, allow_pickle=False)
    if [str(x) for x in e["adt_names"]] != [str(x) for x in names]:
        raise SystemExit("the final study does not list the proteins of the head in the head's order")
    f = {"z": train_arm.unit_length(np.load(z)), "adt": e["adt"], "donors": e["donors"].astype(str), "status": e["status"].astype(str),
         "barcodes": e["barcodes"].astype(str), "present": e["adt_present"].astype(bool)}
    if len(f["z"]) != len(f["adt"]):
        raise SystemExit(f"final study: {len(f['z'])} states for {len(f['adt'])} cells")
    return f


def tests_of(hyp, fam, admitted, arms_of_rule):
    """`loop.test` over the seeds together and for each seed. admitted: {arm: [bool per seed]}."""
    together = loop.test(hyp, {k: {"gain": fam["arms"][k]["gain"], "ci": fam["arms"][k]["ci"], "runs": [RUN], "admitted": all(admitted[k])}
                               for k in arms_of_rule})
    by_seed = [loop.test(hyp, {k: {"gain": fam["arms"][k]["per_fit"][j]["gain"], "ci": fam["arms"][k]["per_fit"][j]["ci"], "runs": [f"{RUN}:seed{s}"],
                                   "admitted": admitted[k][j]} for k in arms_of_rule}) for j, s in enumerate(plan.SEEDS)]
    return {"together": together, "by_seed": by_seed, "outcome_of_the_round": outcome_of_round(together, by_seed)}


def main():
    p = argparse.ArgumentParser()
    for k in ("heads", "bm-arrays", "bm-z", "pbmc-arrays", "pbmc-z", "final-arrays", "final-z", "head-weights"):
        p.add_argument(f"--{k}", type=Path, required=True)
    a = p.parse_args()
    res, here = common.RUN / "results", Path(__file__).resolve().parent
    decl = json.loads((res / "gate5_declaration.json").read_text())
    names = {f"s{s}_{HEAD_OF[k]}": (k, s) for s in plan.SEEDS for k in ARMS}
    files = {"final_arrays": a.final_arrays, "final_z": a.final_z, **{f"head:{n}": a.heads / n / "head.npz" for n in names},
             **{f"script:{s}": here / s for s in SCRIPTS}}
    check_declaration(decl, files)

    bm, pb, shared = common.load_studies(a.bm_arrays, a.bm_z, a.pbmc_arrays, a.pbmc_z)
    F = load_final(a.final_arrays, a.final_z, bm["names"])
    cols = np.array([j for j in shared if F["present"][j]])
    if [str(bm["names"][j]) for j in cols] != decl["proteins"]:
        raise SystemExit("the proteins F measures among the 112 are not those of the declaration")
    Y = F["adt"][:, cols].astype(np.float64)
    exact = set(decl["proteins_by_exact_name"])                  # a second reading, not part of the outcome: without the proteins matched by an alias
    sub = np.array([i for i, n in enumerate(decl["proteins"]) if n in exact])
    if np.isnan(Y).any():
        raise SystemExit("a protein of the declaration is not measured in every cell of F")
    donors = [d["donor"] for d in decl["donors"]]
    if sorted(set(F["donors"].tolist())) != sorted(donors):
        raise SystemExit("the donors of the pack are not the donors of the declaration")
    patient = np.array([d["status"] == "patient" for d in decl["donors"]])
    for d in decl["donors"]:                                    # the pack agrees with the declaration, donor by donor
        st = set(F["status"][F["donors"] == d["donor"]].tolist())
        if st != {d["status"]}:
            raise SystemExit(f"{d['donor']}: status {st} in the pack, {d['status']} in the declaration")
    fit_files = {sha256(a.bm_arrays), sha256(a.pbmc_arrays), sha256(a.bm_z), sha256(a.pbmc_z)}
    if {decl["sha256"]["final_arrays"], decl["sha256"]["final_z"]} & fit_files:
        raise SystemExit("a file of F is among the files the fits read")

    dev = json.loads((res / "development.json").read_text())["heads"]
    todo = {c["name"]: c for c in plan.heads(n_train=int((bm["split"] == "train").sum())) + plan.heads(n_train=int((bm["split"] == "train").sum()), final=True)}
    evalc = (bm["split"] == "val") | (bm["split"] == "test")

    def run(W, on_fit_side=True):
        out = {"final": common.head_forward(W, F["z"])[:, cols]}
        if on_fit_side:
            out["pbmc"] = common.head_forward(W, pb["z"])[:, shared]
            out["bm"] = np.zeros(bm["adt"].shape, np.float32)
            out["bm"][evalc] = common.head_forward(W, bm["z"][evalc])
        return out

    official = run(dict(np.load(a.head_weights, allow_pickle=False)), on_fit_side=False)
    U = {k: np.full((len(donors), len(plan.SEEDS)), np.nan) for k in ARMS + ("head_in_use",)}
    S = {k: np.full((len(donors), len(plan.SEEDS)), np.nan) for k in ARMS + ("head_in_use",)}
    U2 = {k: np.full((len(donors), len(plan.SEEDS)), np.nan) for k in ARMS}
    table, events, admitted, per_donor = {}, [], {k: [] for k in ARMS}, {}
    for j, s in enumerate(plan.SEEDS):
        P = {}
        for k in ARMS:
            n = f"s{s}_{HEAD_OF[k]}"
            d = a.heads / n
            P[k] = {**run(dict(np.load(d / "head.npz"))), "metrics": json.loads((d / "metrics.json").read_text()), "cfg": todo[n],
                    "cells": dict(np.load(d / "fit_cells.npz")), "name": n}
            if P[k]["metrics"]["head"] != todo[n]:
                raise SystemExit(f"{n}: the head on disk was fitted under another plan")
        base = P["as_is"]
        for k in ARMS:
            h, n = P[k], P[k]["name"]
            new = n not in dev
            g3 = read.gate3(h["cfg"], h["metrics"], h["cells"], h, base, bm, pb, shared) if new else dev[n]["gate3"]
            g6 = read.guard6(bm, h, base, shared)
            table[n] = {"arm": k, "seed": s, "fitted": "gate 5" if new else "development", "checkpoint_epoch": h["metrics"]["checkpoint_epoch"],
                        "n_fit": h["metrics"]["n_fit"], "n_fit_second_study": h["metrics"]["n_fit_second_study"], "gate3": g3, "gate6": g6}
            admitted[k].append(bool(g3["admitted"]))
            if new:
                c = h["cfg"]
                events.append(loop.intervention_event(
                    system=n, parent="teddy_mm phase-1 head recipe (scripts/04_train.py)", component="the cells the head is fitted on",
                    change={"arm": c["arm"], "bone_marrow_cells": h["metrics"]["n_fit_bone_marrow"], "second_study_cells": c["pbmc"],
                            "proteins_shuffled": c["shuffle_labels"], "hidden": c["hidden"], "epochs": c["epochs"], "seed": s},
                    role="located" if k == "coverage" else "control", held_fixed=["model", "optimizer", "loss per protein", "checkpoint rule", "TEDDY state"],
                    fitted_on={"fit_set": "all", "scored_donors": donors}))
        for i, donor in enumerate(donors):
            me = F["donors"] == donor
            r = scope.arms(Y[me], {**{k: P[k]["final"][me] for k in ARMS}, "head_in_use": official["final"][me]}, base="as_is", decision=DECISION)
            per_donor.setdefault(donor, {})[f"seed{s}"] = {k: {x: v[x] for x in ("u", "ci", "gain", "gain_ci")} for k, v in r["arms"].items()}
            for k in U:
                U[k][i, j] = r["arms"][k]["u"]
                S[k][i, j] = common.score(Y[me], (official["final"] if k == "head_in_use" else P[k]["final"])[me])
            r2 = scope.arms(Y[me][:, sub], {k: P[k]["final"][me][:, sub] for k in ARMS}, base="as_is", decision=DECISION)
            for k in ARMS:
                U2[k][i, j] = r2["arms"][k]["u"]

    fam = scope.arms_family({k: U[k] for k in ARMS}, base="as_is", bootstrap=DECISION["bootstrap"], seed=DECISION["seed"], alpha=DECISION["alpha"])
    fam2 = scope.arms_family({k: U2[k] for k in ARMS}, base="as_is", bootstrap=DECISION["bootstrap"], seed=DECISION["seed"], alpha=DECISION["alpha"])
    hyp = {h["id"]: h for h in decl["hypotheses"]}
    T = {"donors_hao_coverage": tests_of(hyp["donors_hao_coverage"], fam, admitted, ("coverage", "function", "amount")),
         "coverage_not_levels": tests_of(hyp["coverage_not_levels"], fam, admitted, ("coverage", "labels"))}
    ip, ih = np.where(patient)[0], np.where(~patient)[0]
    groups = {k: {"healthy": group_gain(U["as_is"], U[k], ih), "patients": group_gain(U["as_is"], U[k], ip),
                  "patients_minus_healthy": group_difference(U["as_is"], U[k], ip, ih, **DECISION)} for k in ARMS if k != "as_is"}
    guard_ok = bool(all(t["gate6"]["pass"] for t in table.values()))
    tests = [t for v in T.values() for t in [v["together"]] + v["by_seed"]]
    checks = loop.ledger_check(list(hyp.values()) + tests)
    out = {
        "what": "gate 5: the heads of the test on the final study, read once", "run": RUN, "declaration": decl["digest"],
        "donors": len(donors), "healthy": int(len(ih)), "patients": int(len(ip)), "cells": int(len(Y)), "proteins": int(len(cols)), "decision": DECISION,
        "approximate": fam["approximate"],
        "u": {k: {"mean": float(np.nanmean(v)), "by_seed": [float(x) for x in np.nanmean(v, axis=0)]} for k, v in U.items()},
        "score_within_donor": {k: {"mean": float(np.nanmean(v)), "by_seed": [float(x) for x in np.nanmean(v, axis=0)]} for k, v in S.items()},
        "gains_over_as_is": {k: {x: fam["arms"][k][x] for x in ("gain", "ci", "per_fit", "per_unit")} for k in ARMS if k != "as_is"},
        "tests": T, "outcome": {k: v["outcome_of_the_round"] for k, v in T.items()},
        "guard_holds": guard_ok, "admission": {k: dict(zip([f"seed{s}" for s in plan.SEEDS], v)) for k, v in admitted.items()},
        "exact_name_proteins_only": {"proteins": int(len(sub)), "u": {k: float(np.nanmean(v)) for k, v in U2.items()},
                                     "gains_over_as_is": {k: {x: fam2["arms"][k][x] for x in ("gain", "ci", "per_fit")} for k in ARMS if k != "as_is"}},
        "healthy_and_patients": groups, "forecasts": loop.forecast_record([T[k]["together"] for k in T]),
        "ledger_checks": checks, "heads": table, "per_donor": per_donor, "status_of_donor": {d["donor"]: d["status"] for d in decl["donors"]},
    }
    (res / "final.json").write_text(json.dumps(out, indent=1))
    with open(res / "ledger_final.jsonl", "a") as f:
        for e in events + tests:
            f.write(json.dumps(e, ensure_ascii=False, default=float) + "\n")
    print(f"donors {len(donors)} ({len(ih)} healthy, {len(ip)} patients); proteins {len(cols)}; guard holds {guard_ok}")
    for k, v in T.items():
        print(f"  {k}: {v['outcome_of_the_round']} (together: {v['together']['outcome']}; by seed: {[t['outcome'] for t in v['by_seed']]})")
    for k in ARMS[1:]:
        g = fam["arms"][k]
        print(f"  {k:10s} gain {g['gain']:+.4f} [{g['ci'][0]:+.4f}, {g['ci'][1]:+.4f}]  by seed {[round(x['gain'], 4) for x in g['per_fit']]}")


if __name__ == "__main__":
    main()
