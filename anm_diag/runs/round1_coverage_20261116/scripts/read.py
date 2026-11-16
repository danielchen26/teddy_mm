#!/usr/bin/env python3
"""Gates 3, 4 and 6 of the round, from the fitted heads. Nothing is fitted here.

  python read.py --heads HEADS_DIR --bm-arrays cite_arrays.npz --bm-z z_rna.npy --pbmc-arrays ext_cite_arrays.npz \\
      --pbmc-z ext_z_rna.npy --head-weights head_weights.npz [--allow-incomplete]

Every head is run from its weights in numpy, as scripts/05_eval.py scores a head (MLP and decoder, constant size
factor). A PBMC donor is scored only by heads that were not fitted on it. Writes ../results/development.json and
../results/ledger.jsonl. Refuses to run before all 45 heads are finished, unless told otherwise (then it decides nothing).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import common
import plan
from anm_diag import loop, scope

GUARD, FALL, OWN_OBJECTIVE, BOUND = 0.02, 0.005, 0.035, 0.03419135502693498     # the last: the rule of donors_hao_coverage
DECISION = {"bootstrap": 1000, "seed": 0, "alpha": 0.05}
ARMS = ("as_is", "coverage_seen", "coverage_not_seen", "amount", "function", "labels", "four_donors", "one_donor_own_batch",
        "one_donor_other_batch")


# ------------------------------------------------------------------ which head scores which donor


def heads_for(seed, donor):
    """The heads of a seed that may score a donor (none was fitted on it), by what they stand for."""
    m = next(fs for fs in ("M1", "M2") if donor in plan.scored_donors(fs))
    b = next(fs for fs in ("B1", "B2") if donor in plan.scored_donors(fs))
    own = plan.BATCH_OF[donor]
    one = {("own" if plan.BATCH_OF[f"hao2021_P{i}"] == own else "other"): f"s{seed}_one_donor_{m}_P{i}" for i in plan.ONE_DONOR[m]}
    return {"as_is": f"s{seed}_as_is", "coverage_seen": f"s{seed}_coverage_{m}", "coverage_not_seen": f"s{seed}_coverage_{b}",
            "amount": f"s{seed}_amount", "function": f"s{seed}_function", "labels": f"s{seed}_labels_{m}",
            "four_donors": f"s{seed}_four_donors_{m}", "one_donor_own_batch": one["own"], "one_donor_other_batch": one["other"]}


# ------------------------------------------------------------------ gate 3


def ran_without_falling(history):
    a = np.array([h["val_mlp_pearson"] for h in history])
    last, before = float(a[-5:].mean()), float(a[-10:-5].mean())
    return {"pass": bool(last >= before - FALL), "last5": last, "before5": before}


def within_donor_score(Y, P, donors):
    return float(np.mean([common.score(Y[donors == d], P[donors == d]) for d in np.unique(donors)]))


def gate3(cfg, metrics, cells, pred, base, bm, pb, shared):
    """The four checks for one head. pred / base: {"pbmc", "bm"} outputs of the head and of `as_is` of its seed."""
    out = {"ran": ran_without_falling(metrics["history"])}
    if cfg["arm"] == "coverage":
        i = cells["pbmc_index"]
        Y, d = pb["adt"][i][:, shared], pb["donors"][i]
        s, s0 = within_donor_score(Y, pred["pbmc"][i], d), within_donor_score(Y, base["pbmc"][i], d)
        out["own_objective"] = {"pass": bool(s - s0 >= OWN_OBJECTIVE), "score_on_fit_cells": s, "as_is_on_the_same_cells": s0, "gain": s - s0}
    elif cfg["arm"] == "function":
        l, l0 = metrics["history"][-1]["train_loss"], base["metrics"]["history"][-1]["train_loss"]
        out["own_objective"] = {"pass": bool(l < l0), "train_loss": l, "as_is": l0}
    fit_donors = set(pb["donors"][cells["pbmc_index"]].tolist())
    scored = set(plan.scored_donors(cfg["fit_set"])) if cfg["fit_set"] in plan.FIT_SETS else set()
    fit_cells = set(zip(pb["donors"][cells["pbmc_index"]].tolist(), pb["barcodes"][cells["pbmc_index"]].tolist()))
    sc = np.isin(pb["donors"], sorted(scored))
    scored_cells = set(zip(pb["donors"][sc].tolist(), pb["barcodes"][sc].tolist()))
    out["no_scored_cell_in_fit"] = {"pass": bool(not (fit_donors & scored) and not (fit_cells & scored_cells)
                                                 and np.all(bm["split"][cells["bm_index"]] == "train")
                                                 and fit_donors == set(cfg["pbmc"])),
                                    "fit_donors": sorted(fit_donors), "scored_donors": sorted(scored)}
    v = bm["split"] == "val"
    g, g0 = common.score(bm["adt"][v], pred["bm"][v]), common.score(bm["adt"][v], base["bm"][v])
    out["guard_validation"] = {"pass": bool(g >= g0 - GUARD), "score": g, "as_is": g0, "change": g - g0}
    out["admitted"] = bool(all(x["pass"] for x in out.values()))
    return out


def guard6(bm, pred, base, shared):
    """Gate 6, bone marrow site 4: the mean score, the worst single protein, and the proteins the second study lacks."""
    t = bm["split"] == "test"
    r, r0 = common.per_protein(bm["adt"][t], pred["bm"][t]), common.per_protein(bm["adt"][t], base["bm"][t])
    lack = np.setdiff1d(np.arange(len(r)), shared)
    d = r - r0
    return {"score": float(np.nanmean(r)), "as_is": float(np.nanmean(r0)), "change": float(np.nanmean(r) - np.nanmean(r0)),
            "worst_protein": {"name": str(bm["names"][int(np.nanargmin(d))]), "change": float(np.nanmin(d))},
            "proteins_measured_in_both": float(np.nanmean(d[shared])), "proteins_the_second_study_lacks": float(np.nanmean(d[lack])),
            "pass": bool(np.nanmean(r) >= np.nanmean(r0) - GUARD)}


# ------------------------------------------------------------------ main


def main():
    p = argparse.ArgumentParser()
    for k in ("heads", "bm-arrays", "bm-z", "pbmc-arrays", "pbmc-z", "head-weights"):
        p.add_argument(f"--{k}", type=Path, required=True)
    p.add_argument("--allow-incomplete", action="store_true")
    a = p.parse_args()
    bm, pb, shared = common.load_studies(a.bm_arrays, a.bm_z, a.pbmc_arrays, a.pbmc_z)
    todo = plan.heads(n_train=int((bm["split"] == "train").sum()))
    done = [c for c in todo if (a.heads / c["name"] / "DONE").exists()]
    complete = len(done) == len(todo)
    if not complete and not a.allow_incomplete:
        raise SystemExit(f"{len(done)} of {len(todo)} heads are finished; nothing is read before all are")
    seeds = [s for s in plan.SEEDS if all((a.heads / c["name"] / "DONE").exists() for c in todo if c["seed"] == s)]
    Ypb = pb["adt"][:, shared].astype(np.float64)
    evalc = (bm["split"] == "val") | (bm["split"] == "test")

    def run(W):
        out = {"pbmc": common.head_forward(W, pb["z"])[:, shared], "bm": np.zeros(bm["adt"].shape, np.float32)}
        out["bm"][evalc] = common.head_forward(W, bm["z"][evalc])
        return out

    official = run(dict(np.load(a.head_weights, allow_pickle=True)))
    table, U, events, per_donor = {}, {k: np.full((len(plan.DONORS), len(seeds)), np.nan) for k in ARMS + ("head_in_use",)}, [], {}
    scores = {k: np.full((len(plan.DONORS), len(seeds)), np.nan) for k in ARMS + ("head_in_use",)}
    for j, s in enumerate(seeds):
        P = {}
        for c in [c for c in todo if c["seed"] == s]:
            d = a.heads / c["name"]
            P[c["name"]] = {**run(dict(np.load(d / "head.npz"))), "metrics": json.loads((d / "metrics.json").read_text()), "cfg": c,
                            "cells": dict(np.load(d / "fit_cells.npz"))}
        base = P[f"s{s}_as_is"]
        for name, h in P.items():
            c, m = h["cfg"], h["metrics"]
            if m["head"] != c:
                raise SystemExit(f"{name}: the head on disk was fitted under another plan")
            table[name] = {"arm": c["arm"], "seed": s, "fit_set": c["fit_set"], "n_fit": m["n_fit"], "n_fit_second_study": m["n_fit_second_study"],
                           "checkpoint_epoch": m["checkpoint_epoch"], "machine": m["machine"], "minutes": m["seconds"] / 60,
                           "gate3": gate3(c, m, h["cells"], h, base, bm, pb, shared), "gate6": guard6(bm, h, base, shared)}
            events.append(loop.intervention_event(
                system=name, parent="teddy_mm phase-1 head recipe (scripts/04_train.py)", component="the cells the head is fitted on",
                change={"arm": c["arm"], "bone_marrow_cells": m["n_fit_bone_marrow"], "second_study_cells": c["pbmc"], "proteins_shuffled": c["shuffle_labels"],
                        "hidden": c["hidden"], "epochs": c["epochs"], "seed": s},
                role="baseline" if c["arm"] == "as_is" else ("located" if c["arm"] in ("coverage", "one_donor", "four_donors") else "control"),
                held_fixed=["model", "optimizer", "loss per protein", "checkpoint rule", "TEDDY state"],
                fitted_on={"fit_set": c["fit_set"], "scored_donors": plan.scored_donors(c["fit_set"]) if c["fit_set"] in plan.FIT_SETS else []}))
        for i, donor in enumerate(plan.DONORS):
            me = pb["donors"] == donor
            H = heads_for(s, donor)
            for name in H.values():                                   # the check of gate 3, again, where it matters
                assert donor not in P[name]["cfg"]["pbmc"], (donor, name)
            r = scope.arms(Ypb[me], {**{k: P[n]["pbmc"][me] for k, n in H.items()}, "head_in_use": official["pbmc"][me]}, base="as_is",
                           decision=DECISION)
            per_donor.setdefault(donor, {})[f"seed{s}"] = {k: {x: v[x] for x in ("u", "ci", "gain", "gain_ci")} for k, v in r["arms"].items()}
            for k in U:
                U[k][i, j] = r["arms"][k]["u"]
                scores[k][i, j] = common.score(Ypb[me], (official["pbmc"] if k == "head_in_use" else P[H[k]]["pbmc"])[me])

    fam = scope.arms_family({k: v for k, v in U.items() if k != "head_in_use"}, base="as_is")
    batch = scope.arms_family({"coverage_not_seen": U["coverage_not_seen"], "coverage_seen": U["coverage_seen"]}, base="coverage_not_seen")
    levels = scope.arms_family({"labels": U["labels"], "coverage_seen": U["coverage_seen"]}, base="labels")
    fc = json.loads((common.RUN / "results" / "forecasts.json").read_text())["forecast"]
    gain = lambda k: {x: fam["arms"][k][x] for x in ("gain", "ci", "per_fit", "per_unit")}
    share = lambda k, f: {"delivered": fam["arms"][k]["gain"], "ceiling": fc[f]["gain"], "share": fam["arms"][k]["gain"] / fc[f]["gain"]}
    adm = lambda arm: {n: t["gate3"]["admitted"] for n, t in table.items() if t["arm"] == arm}
    g6 = lambda arm: {n: t["gate6"]["pass"] for n, t in table.items() if t["arm"] == arm}
    cov_b = fam["arms"]["coverage_not_seen"]
    go = {"coverage_admitted_in_every_seed": bool(all(adm("coverage").values())), "gain_with_batch_not_seen_above_0": bool(cov_b["ci"][0] > 0),
          "guard_holds": bool(all(g6("coverage").values()) and all(t["gate3"]["guard_validation"]["pass"] for t in table.values() if t["arm"] == "coverage"))}
    out = {
        "complete": complete, "seeds": seeds, "heads_read": len(table), "proteins": int(len(shared)), "decision": DECISION,
        "numbers": {"guard": GUARD, "own_objective": OWN_OBJECTIVE, "bound_of_the_rule": BOUND},
        "u": {k: {"mean": float(np.nanmean(v)), "by_seed": [float(x) for x in np.nanmean(v, axis=0)]} for k, v in U.items()},
        "score_within_donor": {k: {"mean": float(np.nanmean(v)), "by_seed": [float(x) for x in np.nanmean(v, axis=0)]} for k, v in scores.items()},
        "gains_over_as_is": {k: gain(k) for k in ARMS if k != "as_is"}, "approximate": fam["approximate"],
        "D1_batch_seen": {**gain("coverage_seen"), "of_the_ceiling": share("coverage_seen", "pbmc_seen"), "one_linear_map_forecast": fc["pooled_seen"]["gain"]},
        "D2_batch_not_seen": {**gain("coverage_not_seen"), "of_the_ceiling": share("coverage_not_seen", "pbmc_not_seen"),
                              "one_linear_map_forecast": fc["pooled_not_seen"]["gain"],
                              "seen_minus_not_seen": {x: batch["arms"]["coverage_seen"][x] for x in ("gain", "ci")}},
        "D3_donors_or_cells": {k: gain(k) for k in ("one_donor_own_batch", "one_donor_other_batch", "four_donors", "coverage_seen")},
        "D4_levels_or_relation": {"labels": gain("labels"), "coverage_minus_labels": {x: levels["arms"]["coverage_seen"][x] for x in ("gain", "ci")}},
        "D5_cost": {arm: {"validation_change": float(np.mean([t["gate3"]["guard_validation"]["change"] for t in table.values() if t["arm"] == arm])),
                          "site4_change": float(np.mean([t["gate6"]["change"] for t in table.values() if t["arm"] == arm])),
                          "site4_worst_protein": min((t["gate6"]["worst_protein"] for t in table.values() if t["arm"] == arm), key=lambda w: w["change"]),
                          "site4_proteins_the_second_study_lacks": float(np.mean([t["gate6"]["proteins_the_second_study_lacks"] for t in table.values()
                                                                                  if t["arm"] == arm])),
                          "every_head_within_guard": bool(all(g6(arm).values()))}
                    for arm in ("amount", "function", "coverage", "labels", "one_donor", "four_donors")},
        "against_the_rule_of_the_hypothesis": {
            "note": "development: the same donors the forecast was read on, so this is not the test",
            **{k: {"gain": fam["arms"][k]["gain"], "ci": fam["arms"][k]["ci"],
                   "reads": ("above" if fam["arms"][k]["ci"][0] >= BOUND else "below" if fam["arms"][k]["ci"][1] < BOUND else "crosses") + " the bound"}
               for k in ("coverage_seen", "coverage_not_seen", "amount", "function", "labels")}},
        "admission": {arm: adm(arm) for arm in ("as_is", "amount", "function", "coverage", "labels", "one_donor", "four_donors")},
        "go_to_gate_5": {**go, "go": bool(complete and all(go.values()))},
        "heads": table, "per_donor": per_donor,
    }
    res = common.RUN / "results"
    res.mkdir(exist_ok=True)
    (res / ("development.json" if complete else "development_incomplete.json")).write_text(json.dumps(out, indent=1))
    if complete:
        loop.write_ledger(events, res / "ledger.jsonl")
    print(f"heads read {len(table)}; admitted {sum(t['gate3']['admitted'] for t in table.values())}; go {out['go_to_gate_5']}")
    for k in ARMS[1:]:
        g = fam["arms"][k]
        print(f"  {k:24s} gain {g['gain']:+.4f} [{g['ci'][0]:+.4f}, {g['ci'][1]:+.4f}]  by seed {[round(f['gain'], 4) for f in g['per_fit']]}")


if __name__ == "__main__":
    main()
