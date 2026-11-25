#!/usr/bin/env python3
"""Part A: the round with the rows applied, on the stages that exist already. results/round_a.json. Made once.

  python round_a.py --work WORK --observe OBSERVE [--results RESULTS] [--only T1,C1]

WORK and OBSERVE are those of runs/input_beside_z_20261123: the cells of the five targets with the state Z, the input
view V and the two side by side (ZV), and the tool's grids for them. Nothing is computed by a model here.

For every target, with the head reading Z and reading ZV:
  the round as it was (cycle.LinearRows), which must repeat the record of 23 November;
  the round with `Full` (repair.py): the arms controls, rows, exact, kalman, cells; gate 3 and the grade by the new rows;
  on paired resamples, what the new rows remove that the old rows leave;
  the correlation of every output with the protein among control cells that the fit did not see.
For every event of a target: read on its own with its own record (old rows, new rows); and read as an event that has
no protein record: `none` (the readout of control cells), `type` (the readout repaired by the rows of the other events
of its type) and `derived` (the readout of control cells with a row derived from the event's RNA; repair.derive).
Two controls that must fail: the same repair on a state whose rows are permuted, and an event made of control cells.
--only is for a trial of the script on the two THP-1 targets, which were explored before the design was written.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from repair import SEED, STEP, Full, cross, cycle, derive, gain, u_of, unit, val

RUN, RUNS = Path(__file__).resolve().parents[1], Path(__file__).resolve().parents[2]
TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
PATHWAY = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
# target -> (stratum file, readout, events, pooled, the events of its type that have records, explored before the design)
TARGETS = {"T1": ("thp1", "PDL1", ["IRF1"], False, TYPE_A, True), "T2": ("mel_co", "MHC1", PATHWAY, True, PATHWAY, False),
           "C1": ("thp1a", "PDL1", TYPE_A, True, TYPE_A, True), "C2": ("mel_co", "PDL1", PATHWAY, True, PATHWAY, False), "C3": ("mel_ifn", "PDL1", PATHWAY, True, PATHWAY, False)}
FULL = {"thp1a": "thp1"}
STAGES = ("Z", "ZV")
ARMS_OLD, ARMS_NEW = ("controls", "rows", "cells"), ("controls", "rows", "exact", "kalman", "cells")


def cells(work, g, readout):
    d = np.load(work / f"{FULL.get(g, g)}_zv.npz")
    z = unit(d["Z"]).astype(np.float64)
    c = dict(ids=d["ids"].astype(str), Z=z, ZV_old=d["ZV"].astype(np.float64), ZV=np.hstack([z, d["V"].astype(np.float64)]), V=d["V"].astype(np.float64),
             treated=d["A"].astype(int) == 1, unit=d["unit"].astype(str), y=d[f"y_{readout}"].astype(np.float64), o=d[f"frozen_{readout}"].astype(np.float64))
    c["first"] = cycle.halves(np.where(c["treated"], c["unit"], "control"), SEED)
    return c


def sub(c, m):
    return {k: (v[m] if isinstance(v, np.ndarray) and len(v) == len(m) else v) for k, v in c.items()}


def r_control(y, p, treated):
    C = ~treated & np.isfinite(p)
    return float(np.corrcoef(p[C], y[C])[0, 1])


def unrecorded(c, stage, event, recorded):
    """`event` read as an event without a protein record; the other events of its type have records. Every cell is read
    by a fit that did not hold it; no protein value of `event` enters `none`, `type` or `derived`."""
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
    old, _i, _f = cross(x, y, o, tr, un, first, [event], repair=cycle.LinearRows)
    rd = ~tr | (un == event)
    P = {"none": typ["controls"], "type": typ["rows"], "derived": der["rows"], "own_record": own["rows"], "own_record_old_rule": old["rows"]}
    U = {a: u_of(y[rd], p[rd], tr[rd]) for a, p in P.items()}
    measured = float(y[tr & (un == event)].mean() - y[~tr].mean())
    out = {"measured_change": measured, "the_rule": {"by_fold": logs, "derived_over_measured": [lg["change"] / measured for lg in logs]},
           "arms": {a: {**val(u, sl), "r_control": r_control(y[rd], P[a][rd], tr[rd])} for a, (u, sl) in U.items()},
           "gain_over_none": {a: gain(U["none"][0], U[a][0]) for a in ("type", "derived", "own_record")},
           "derived_against_own_record": gain(U["derived"][0], U["own_record"][0]), "new_rule_over_old_with_its_record": gain(U["own_record_old_rule"][0], U["own_record"][0])}
    return out


def must_fail(c, stage, events, pooled, rng):
    """Two controls. A state with its rows permuted holds nothing of a cell: the new rows must not close there. An event
    made of control cells has no change: the widening must stay at one and the readout as it is."""
    m = np.ones(len(c["y"]), bool) if pooled else ~c["treated"] | (c["unit"] == events[0])
    s = sub(c, m)
    xp = s[stage][rng.permutation(len(s["y"]))]
    out, _i, _f = cross(xp, s["y"], s["o"], s["treated"], s["unit"], s["first"], events, cells=False)
    cl = s["unit"] if len(events) > 1 else None
    perm = val(*u_of(s["y"], out["rows"], s["treated"], cl))
    k = sub(c, ~c["treated"])
    n = len(k["y"])
    size = int(np.median([(c["unit"] == e).sum() for e in events]))
    fake = np.zeros(n, bool)
    fake[rng.choice(n, min(size, n // 3), replace=False)] = True
    un = np.where(fake, "made-up", "control")
    first = cycle.halves(un, SEED)
    o2, _i, fits = cross(k[stage], k["y"], k["o"], fake, un, first, ["made-up"], cells=False)
    spread = float(np.nanstd(o2["controls"][~fake]))
    return {"permuted_state": perm, "made_up_event": {"cells": int(fake.sum()), "widening": [f["widening"] for f in fits],
                                                      "mean_change_of_the_output_among_its_cells_in_spreads": float(abs(np.nanmean((o2["rows"] - o2["controls"])[fake])) / spread),
                                                      "largest_change_of_the_output_of_a_cell_in_spreads": float(np.nanmax(np.abs(o2["rows"] - o2["controls"])) / spread)}}


RANK = {"not a target for a repair": 0, "not tested (gate 3)": 0, "not closed": 0, "loss no longer resolved": 1, "relative": 2, "absolute": 3, "nothing to close": 3}
MEL = ("T2", "C2", "C3")


def hypotheses(T):
    """DESIGN.md, Section 6. Melanoma only; the THP-1 targets were explored before the design and are not counted."""
    three = lambda ok, bad: "confirmed" if ok else ("refuted" if bad else "unresolved")
    z = lambda t: T[t]["stages"]["Z"]
    h1 = z("T2")["new_rows_over_old_rows"]
    grades = {t: {st: (x["as_it_was"]["closure"]["grade"], x["with_the_rows_applied"]["closure"]["grade"]) for st, x in T[t]["stages"].items()} for t in MEL}
    kept = all(RANK[grades[t]["Z"][1]] >= RANK[grades[t]["Z"][0]] for t in MEL)
    unharmed = all(z(t)["new_rows_over_old_rows"]["no_harm"] for t in ("C2", "C3"))
    ev = lambda t: T[t]["events_one_by_one"]
    own = {t: {e: q["Z"]["new_rule_over_old_with_its_record"] for e, q in ev(t).items()} for t in MEL}
    steps3 = sum(q["reading"] == "a step" for q in own["T2"].values())
    harmed3 = sum(not q["no_harm"] for t in MEL for q in own[t].values())
    harmed_mhc, harmed_pdl1 = sum(not q["no_harm"] for q in own["T2"].values()), sum(not q["no_harm"] for t in ("C2", "C3") for q in own[t].values())
    perm = {f"{t} {st}": x["must_fail"]["permuted_state"]["u"] for t in MEL for st, x in T[t]["stages"].items()}
    made = {f"{t} {st}": x["must_fail"]["made_up_event"]["mean_change_of_the_output_among_its_cells_in_spreads"] for t in MEL for st, x in T[t]["stages"].items()}
    drop = {f"{t} {st}": x["r_control"]["with_the_rows_applied"]["controls"] - x["r_control"]["with_the_rows_applied"]["rows"] for t in MEL for st, x in T[t]["stages"].items()}
    der = {t: {e: q["Z"]["gain_over_none"]["derived"] for e, q in ev(t).items()} for t in MEL}
    steps6 = sum(q["reading"] == "a step" for q in der["T2"].values())
    harmed6 = sum(not q["no_harm"] for t in MEL for q in der[t].values())
    ratio = {f"{t} {e}": float(np.mean(q["Z"]["the_rule"]["derived_over_measured"])) for t in MEL for e, q in ev(t).items()}
    inside = sum(0.67 <= r <= 1.5 for r in ratio.values())
    return {"H1_the_rows_are_applied": {"reading": {"a step": "confirmed", "no step": "refuted"}.get(h1["reading"], "unresolved"), "T2_at_Z": h1},
            "H2_no_harm": {"reading": "confirmed" if kept and unharmed else "refuted", "grades_as_it_was_and_with_the_rows_applied": grades, "C2_and_C3_not_harmed": unharmed},
            "H3_one_event_its_record": {"reading": three(steps3 >= 3 and harmed_mhc == 0 and harmed_pdl1 <= 1, steps3 < 2 or harmed3 >= 2), "MHC1_events_that_gain_a_step": steps3, "harmed_of_fifteen": harmed3, "by_event": own},
            "H4_not_a_trick": {"reading": "confirmed" if min(perm.values()) >= 0.90 and max(made.values()) < 0.10 else "refuted", "permuted_state": perm, "made_up_event": made},
            "H5_the_guard": {"reading": "confirmed" if max(drop.values()) < 0.03 else "refuted", "fall_of_the_correlation_among_control_cells": drop},
            "H6_without_a_protein_record": {"reading": three(steps6 >= 4 and harmed6 <= 2, steps6 < 2 or harmed6 > 4), "MHC1_events_where_derived_gains_a_step_over_none": steps6, "harmed_of_fifteen": harmed6, "by_event": der},
            "H7_the_rule_on_the_records": {"reading": three(inside >= 12, inside < 9), "events_with_the_derived_change_between_0.67_and_1.5_of_the_measured": inside, "derived_over_measured": ratio}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results")
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    out_file = a.results / "round_a.json"
    if out_file.exists():
        raise SystemExit("round_a.json exists: this reading was made once and is not made again")
    only = [t for t in a.only.split(",") if t] or list(TARGETS)
    rec = json.loads((RUNS / "input_beside_z_20261123" / "results" / "round.json").read_text())["targets"]
    res = {"margin": cycle.MARGIN, "step_margin": STEP, "note": "development cells; nothing is computed by a model; the stages are those of runs/input_beside_z_20261123", "targets": {}}
    if a.only:
        res["TRIAL"] = "a trial of the script on targets explored before the design; not the reading"
    worst = 0.0
    rng = np.random.default_rng(20261125)
    for t in only:
        g, readout, events, pooled, typ, explored = TARGETS[t]
        c = cells(a.work, g, readout)
        by, draws = cycle.load(a.observe / f"{t}_{g}_{readout}")
        if pooled:
            m, scope, key = np.isin(c["unit"], events + ["control"]), "units", "pool"
        else:
            m, scope, key = ~c["treated"] | (c["unit"] == events[0]), "own", events[0]
        s = sub(c, m)
        y, o, tr, un, first = s["y"], s["o"], s["treated"], s["unit"], s["first"]
        cl = un if len(events) > 1 else None
        x = {"stratum": g, "readout": readout, "events": events, "read_as": "a pool" if pooled else "on its own, with its own record", "explored_before_the_design": explored, "stages": {}}
        for st in STAGES:
            ref = draws[f"{key}|shift|{st}|{scope}|{by[key][st][scope]['class']}"]
            refs = {e: draws[f"{e}|shift|{st}|units|{by[e][st]['units']['class']}"] for e in events} if pooled else None
            xo = s["ZV_old"] if st == "ZV" else s[st]              # the round as it was read ZV with its two blocks at one total variance
            out0, in0, _f = cycle.cross_fit(xo, y, o, tr, un, events, first)
            e0 = cycle.read(y, o, tr, un, events, out0, in0, ref, refs, seed=SEED)
            d = max(abs(e0["arms"][arm]["u"] - rec[t][st]["arms"][arm]["u"]) for arm in ARMS_OLD)
            worst = max(worst, d)
            out1, in1, fits = cross(s[st], y, o, tr, un, first, events)
            e1 = cycle.read(y, o, tr, un, events, out1, in1, ref, refs, seed=SEED)
            U0 = {arm: u_of(y, out0[arm], tr, cl)[0] for arm in ARMS_OLD}
            U1 = {arm: u_of(y, out1[arm], tr, cl)[0] for arm in ARMS_NEW}
            x["stages"][st] = {"as_it_was": {k: e0[k] for k in ("reference", "as_trained", "arms", "A3", "gate3", "closure")}, "difference_from_the_record": d,
                               "with_the_rows_applied": {**{k: e1[k] for k in ("arms", "gate3", "closure", "by_event")}, "fits": fits},
                               "new_rows_over_old_rows": gain(U0["rows"], U1["rows"]), "new_rows_over_old_cells": gain(U0["cells"], U1["rows"]),
                               "exact_over_new_rows": gain(U1["rows"], U1["exact"]), "kalman_on_one_scale_over_old_rows": gain(U0["rows"], U1["kalman"]),
                               "r_control": {"as_it_was": {arm: r_control(y, out0[arm], tr) for arm in ARMS_OLD}, "with_the_rows_applied": {arm: r_control(y, out1[arm], tr) for arm in ARMS_NEW}},
                               "must_fail": must_fail(c, st, events, pooled, rng)}
            print(f"{t} {st}: old rows {e0['arms']['rows']['u']:.3f} (record {rec[t][st]['arms']['rows']['u']:.3f}) -> new rows {e1['arms']['rows']['u']:.3f}; gate 3 {e0['gate3']['passes']} -> {e1['gate3']['passes']}; "
                  f"grade {e0['closure']['grade']} -> {e1['closure']['grade']}", flush=True)
        x["stages_gain_with_the_rows_applied"] = None
        singles = events if pooled else events
        x["events_one_by_one"] = {}
        for e in singles:
            recorded = [h for h in typ if h != e]
            x["events_one_by_one"][e] = {st: unrecorded(c, st, e, recorded) for st in STAGES}
            q = x["events_one_by_one"][e]["Z"]["arms"]
            print(f"   {t} {e} at Z: none {q['none']['u']:+.3f} | type {q['type']['u']:+.3f} | derived {q['derived']['u']:+.3f} | its record {q['own_record']['u']:+.3f} (old rule {q['own_record_old_rule']['u']:+.3f})", flush=True)
        del x["stages_gain_with_the_rows_applied"]
        res["targets"][t] = x
    res["check_against_the_record"] = {"largest_difference_of_an_arm": worst, "limit": 0.005, "hold": bool(worst <= 0.005)}
    if worst > 0.005:
        raise SystemExit(f"the round as it was does not repeat the record of 23 November ({worst:.4f}): stopped, nothing written")
    if not a.only:
        res["hypotheses"] = hypotheses(res["targets"])
        print(json.dumps({k: v["reading"] for k, v in res["hypotheses"].items()}))
    a.results.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("written", out_file.name, "| check against the record", worst)


if __name__ == "__main__":
    main()
