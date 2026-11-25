#!/usr/bin/env python3
"""Parts B and D: the round with the head reading a typed state. results/round_b.json. Made once.

  python round_b.py --work WORK --observe OBSERVE --states STATES [--results RESULTS] [--only T1,C1]

WORK and OBSERVE as for round_a.py. STATES: the folder of states.py (<stratum>_typed.npz).

Stages, each made for the half it is fitted on (the genes were chosen on that half; states.py):
  Z    the mean at unit length (read again here, so that every gain is on paired resamples)
  T    two blocks by the TEG: the genes of the recorded events of the type, and the rest
  R    the control: genes at random, and the rest
  P    four programs of the control cells, found with no record, and the rest
A block is the sum over its genes' positions divided by the number of positions and by the length of the mean state,
so the blocks of a stage add up to Z. The repair is `Full` for every stage; the reference of the grade is the tool's
reference reader of Z (no grid of the tool exists for a stage that differs by half).

For a target: the arms of the round at every stage and what a stage removes that Z leaves. For an event read on its
own or as an event without a protein record, the genes of T are those of the other recorded events of its type
(A-<event>); the arms are those of round_a.py.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from repair import SEED, STEP, cross, cycle, derive, gain, u_of, unit, val
from round_a import FULL, MEL, RANK, RUN, TARGETS, cells, r_control

STAGES = ("Z", "T", "R", "P")
ARMS = ("controls", "rows", "exact", "kalman", "cells")


def sub(c, m):
    f = lambda v: v[m] if isinstance(v, np.ndarray) and len(v) == len(m) else v
    return {k: (tuple(f(x) for x in v) if isinstance(v, tuple) else f(v)) for k, v in c.items()}


def typed(states, g, c):
    """The stages of a stratum from the file of states.py, in the order of the cells of `c`."""
    d = np.load(states / f"{FULL.get(g, g)}_typed.npz")
    ids = d["ids"].astype(str)
    assert np.array_equal(ids, c["ids"]), "the typed states are of other cells, or in another order"
    assert np.array_equal(d["first"], c["first"]), "the halves of the typed states are not the halves of the round"
    zc = d["Zc"].astype(np.float64)
    worst = float(np.abs(unit(d["Zc"]).astype(np.float64) - c["Z"]).max())
    assert worst <= 2e-3, f"the mean of the typed states is not the state of the round ({worst})"
    nz = np.linalg.norm(zc, axis=1, keepdims=True) + 1e-6
    names = [str(s) for s in d["sets"]]
    B = d["blocks"]
    blk = lambda name, f: B[:, names.index(f"{name}|{f}")].astype(np.float64)
    two = lambda name: tuple(np.hstack([blk(name, f), zc - blk(name, f)]) / nz for f in (0, 1))
    out = {"T": two("A"), "R": two("R")}
    for name in {n.split("|")[0] for n in names if n.startswith("A-")}:
        out[f"T:{name}"] = two(name)
    ks = sorted({n.split("|")[0] for n in names if n.startswith("P")})
    out["P"] = tuple(np.hstack([blk(k, f) for k in ks] + [zc - sum(blk(k, f) for k in ks)]) / nz for f in (0, 1))
    return out, worst


def unrecorded(c, stage, event, recorded):
    """As round_a.unrecorded, for a stage that may differ by half."""
    m = ~c["treated"] | np.isin(c["unit"], recorded + [event])
    s = sub(c, m)
    x, y, o, tr, un, first = s[stage], s["y"], s["o"], s["treated"], s["unit"], s["first"]
    logs = []

    def drv(f):
        fit = first if f == 0 else ~first
        change, se2, lg = derive(s["V"], y, np.where(fit & ~tr)[0], np.where(fit & (un == event))[0], [(h, np.where(fit & (un == h))[0]) for h in recorded])
        logs.append(lg)
        return {event: (change, se2)}

    typ, _i, _f = cross(x, y, o, tr, un, first, recorded, cells=False)
    der, _i, _f = cross(x, y, o, tr, un, first, [event], derived=drv, cells=False)
    own, _i, _f = cross(x, y, o, tr, un, first, [event], cells=False)
    rd = ~tr | (un == event)
    P = {"none": typ["controls"], "type": typ["rows"], "derived": der["rows"], "own_record": own["rows"]}
    U = {a: u_of(y[rd], p[rd], tr[rd]) for a, p in P.items()}
    return {a: {**val(u, sl), "r_control": r_control(y[rd], P[a][rd], tr[rd])} for a, (u, sl) in U.items()}, {a: u for a, (u, _s) in U.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    out_file = a.results / "round_b.json"
    if out_file.exists():
        raise SystemExit("round_b.json exists: this reading was made once and is not made again")
    only = [t for t in a.only.split(",") if t] or list(TARGETS)
    was = json.loads((RUN / "results" / "round_a.json").read_text())["targets"] if not a.only else None
    res = {"margin": cycle.MARGIN, "step_margin": STEP, "note": "development cells; the typed states are the encoder as released, kept by type; the genes of a type were chosen on the fitted half", "targets": {}}
    if a.only:
        res["TRIAL"] = "a trial of the script; not the reading"
    worst_z, worst_a = 0.0, 0.0
    for t in only:
        g, readout, events, pooled, typ, explored = TARGETS[t]
        c = cells(a.work, g, readout)
        st_all, wz = typed(a.states, g, c)
        worst_z = max(worst_z, wz)
        c.update(st_all)
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
        x = {"stratum": g, "readout": readout, "events": events, "read_as": "a pool" if pooled else "on its own, with its own record", "explored_before_the_design": explored, "stages": {}}
        U = {}
        for st in STAGES:
            out, inside, fits = cross(s[st], y, o, tr, un, first, events, cells=(st != "P"))
            e = cycle.read(y, o, tr, un, events, out, inside, ref, refs, seed=SEED)
            U[st] = {arm: u_of(y, out[arm], tr, cl)[0] for arm in out}
            x["stages"][st] = {"width": int((s[st][0] if isinstance(s[st], tuple) else s[st]).shape[1]), "arms": e["arms"], "gate3": e["gate3"], "closure": e["closure"],
                               "widening": [f["widening"] for f in fits], "share_of_each_row_closed": [[r["share_of_the_residual_closed"] for r in f["rows"]] for f in fits],
                               "r_control": {arm: r_control(y, out[arm], tr) for arm in out}}
            if st != "Z":
                x["stages"][st]["gain_over_Z"] = {arm: gain(U["Z"][arm], U[st][arm]) for arm in U[st] if arm in U["Z"]}
            print(f"{t} {st}: " + " | ".join(f"{arm} {e['arms'][arm]['u']:+.3f}" for arm in ("controls", "rows") if arm in e["arms"]) + f" | gate 3 {e['gate3']['passes']} | {e['closure']['grade']}", flush=True)
        x["T_over_R"] = {arm: gain(U["R"][arm], U["T"][arm]) for arm in ARMS if arm in U["T"]}
        x["P_against_T"] = {arm: gain(U["T"][arm], U["P"][arm]) for arm in ("controls", "rows")}
        if was is not None:
            worst_a = max(worst_a, abs(x["stages"]["Z"]["arms"]["rows"]["u"] - was[t]["stages"]["Z"]["with_the_rows_applied"]["arms"]["rows"]["u"]))
        x["events_one_by_one"] = {}
        for e in events:
            recorded = [h for h in typ if h != e]
            c["Te"] = c[f"T:A-{e}"] if e in typ else c["T"]
            q, UU = {}, {}
            for st in ("Z", "Te", "P"):
                q[st], UU[st] = unrecorded(c, st, e, recorded)
            x["events_one_by_one"][e] = {"Z": q["Z"], "T": q["Te"], "P": q["P"], "T_over_Z": {arm: gain(UU["Z"][arm], UU["Te"][arm]) for arm in UU["Z"]}, "P_over_Z": {arm: gain(UU["Z"][arm], UU["P"][arm]) for arm in UU["Z"]}}
            print(f"   {t} {e}: " + " || ".join(f"{st} " + " ".join(f"{arm} {q[st][arm]['u']:+.2f}" for arm in ("none", "type", "derived", "own_record")) for st in ("Z", "Te", "P")), flush=True)
        res["targets"][t] = x
    res["checks"] = {"the_mean_of_the_typed_states_against_the_state_of_the_round": worst_z, "limit": 2e-3, "rows_at_Z_against_part_A": worst_a, "limit_A": 0.005, "hold": bool(worst_z <= 2e-3 and worst_a <= 0.005)}
    if not res["checks"]["hold"]:
        raise SystemExit(f"a check does not hold ({res['checks']}): stopped, nothing written")
    if not a.only:
        T = res["targets"]
        three = lambda r: {"a step": "confirmed", "no step": "refuted"}.get(r, "unresolved")
        own = {e: q["T_over_Z"]["own_record"] for e, q in T["T2"]["events_one_by_one"].items()}
        n_step = sum(q["reading"] == "a step" for q in own.values())
        drop = T["T2"]["stages"]["T"]["r_control"]["controls"] - T["T2"]["stages"]["T"]["r_control"]["rows"]
        res["hypotheses"] = {
            "HB1_the_typed_state_where_the_state_is_the_limit": {"reading": three(T["T2"]["stages"]["T"]["gain_over_Z"]["rows"]["reading"]), "T2_rows_T_over_Z": T["T2"]["stages"]["T"]["gain_over_Z"]["rows"]},
            "HB2_it_is_the_type": {"reading": three(T["T2"]["T_over_R"]["rows"]["reading"]), "T2_rows_T_over_R": T["T2"]["T_over_R"]["rows"]},
            "HB3_one_event_its_record": {"reading": "confirmed" if n_step >= 3 else ("refuted" if n_step < 2 else "unresolved"), "MHC1_events_that_gain_a_step_at_T": n_step, "by_event": own},
            "HB4_the_guard_at_T": {"reading": "confirmed" if drop < 0.03 else "refuted", "fall_of_the_correlation_among_control_cells": drop},
            "HB5_no_harm": {"reading": "confirmed" if all(T[t]["stages"]["T"]["gain_over_Z"]["rows"]["no_harm"] for t in ("C2", "C3")) else "refuted",
                            "by_target": {t: T[t]["stages"]["T"]["gain_over_Z"]["rows"] for t in ("C2", "C3")}},
            "HD1_types_from_the_input_alone": {"reading": three(T["T2"]["stages"]["P"]["gain_over_Z"]["rows"]["reading"]), "T2_rows_P_over_Z": T["T2"]["stages"]["P"]["gain_over_Z"]["rows"], "P_against_T": T["T2"]["P_against_T"]}}
        print(json.dumps({k: v["reading"] for k, v in res["hypotheses"].items()}))
    a.results.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("written", out_file.name, res["checks"])


if __name__ == "__main__":
    main()
