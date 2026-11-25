#!/usr/bin/env python3
"""The trial that trains the last layer: the reading, made once.

  python read.py --work WORK --observe OBSERVE --head HEAD_WEIGHTS.npz [--part A|B] [--results RESULTS] [--made-up]

Part A (results/trial.json): the arms rows4 and rows5, by which the three hypotheses are judged. Part B
(results/trial_B.json): all five arms, read after the others are trained; it adds the arms mixed, cells5 and prop and
the comparisons between arms. Each part is read once.

WORK: the folder build_tt.py wrote. OBSERVE: the grids of `anm_diag scope` for its two SPEC files, one folder per
SPEC. HEAD_WEIGHTS.npz: the weights of the prediction head as trained (for the guard and for the head's own output at
a trained state). The reading code is anm_diag.cycle, unchanged; the repair is cycle.LinearRows; the halves, cells and
resamples are those of the earlier rounds. DESIGN.md, Sections 5 and 6.

--made-up: for a trial of the script on made-up values only. It skips the check against the recorded rounds and
marks the result as made up.
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
sys.path.insert(0, str(RUNS / "scope_20261116" / "scripts"))
from anm_diag import cycle  # noqa: E402
from _common import SEED, unit  # noqa: E402
from build import head  # noqa: E402

TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
TARGETS = {"T": ("thp1", "T_thp1_IRF1", ["IRF1"], False, "PDL1 thp1 type B"), "C": ("thp1a", "C_thp1a_pool", TYPE_A, True, "PDL1 thp1 type A")}
ARM_OF = {"T4": "rows4", "T5": "rows5", "TM": "mixed", "TC": "cells5", "TP": "prop"}
TRAINED, ARMS, RANK, STEP, GUARD = [], ("controls", "rows", "cells"), {"absolute": 2, "relative": 1}, cycle.STEP, 0.10
step = lambda iv: "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")
v = lambda x: {k: x[k] for k in ("u", "ci", "verdict")}


def gain(a, b):
    """|a| - |b| on paired resamples: what the second stage removes that the first leaves."""
    d = np.abs(a) - np.abs(b)
    iv = cycle.interval(d[1:])
    return {"gain": float(d[0]), "ci": iv, "reading": step(iv), "no_harm": bool(iv[0] > -STEP)}


def against(a, b):
    """An output against another on paired resamples, both ways, by the tool's step reading: |u| of a minus |u| of b."""
    cost, back = cycle.diff(a, b), cycle.diff(b, a)
    word = ("brings back a step" if back["reading"] == "loss" else "costs a step" if cost["reading"] == "loss"
            else "no change of a step" if back["reading"] == cost["reading"] == "no added loss" else "unresolved")
    return {"change": cost["added"], "ci": cost["ci"], "reading": word}


def cells(work, g):
    d = np.load(work / "thp1_tt.npz")
    c = dict(ids=d["ids"].astype(str), treated=d["A"].astype(int) == 1, unit=d["unit"].astype(str), y=d["y_PDL1"].astype(np.float64), o=d["frozen_PDL1"].astype(np.float64),
             raw={st: d[st] for st in ["Z70"] + TRAINED}, **{st: unit(d[st]).astype(np.float64) for st in ["Z70"] + TRAINED})
    c["first"] = cycle.halves(np.where(c["treated"], c["unit"], "control"), SEED)      # the halves of the earlier rounds, made on the whole stratum
    if g == "thp1a":
        keep = np.isin(c["ids"], np.load(work / "thp1a_tt.npz")["ids"].astype(str))
        c = {k: ({s: x[keep] for s, x in val.items()} if isinstance(val, dict) else val[keep]) for k, val in c.items()}
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
    assert abs(e["as_trained"]["u"] - draws[f"{key}|shift|system|frozen|system"][0]) < 2e-4, "the resamples of the round and of the tool do not pair"
    cl = un if len(events) > 1 else None
    U = {a: cycle.output_u(y, p, tr, cl, SEED)[0] for a, p in out.items()}
    e["reader_of_the_stage"] = {sc: {q: by[key][stage][sc][q] for q in ("u", "ci", "verdict", "class")} for sc in ("units", "own") if sc in by[key][stage]}
    return e, U, np.asarray(ref, np.float64), (key, scope, m, cl)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--head", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results", help="another folder for a trial on made-up values")
    ap.add_argument("--made-up", action="store_true")
    ap.add_argument("--part", choices=("A", "B"), default="A")
    a = ap.parse_args()
    out_file = a.results / ("trial.json" if a.part == "A" else "trial_B.json")
    if out_file.exists():
        raise SystemExit(f"{out_file.name} exists: this reading was made once and is not made again")
    there = set(np.load(a.work / "thp1_tt.npz").files)
    TRAINED.extend(st for st in ARM_OF if st in there)
    assert TRAINED == (["T4", "T5"] if a.part == "A" else list(ARM_OF)), f"part {a.part} does not read the stages {TRAINED}"
    rec = json.loads((RUNS / "outer_step_20261122" / "results" / "bridge.json").read_text())["targets"]
    res = {"part": a.part, "margin": cycle.MARGIN, "step_margin": STEP, "note": "exploratory; development cells; every trained state is the state of a cell the trained layer did not hold", "targets": {}}
    if a.made_up:
        res["MADE_UP"] = "a trial of the script on made-up values; not a result"
    W = dict(np.load(a.head))
    keep, worst = {}, 0.0
    for t, (g, spec, events, pooled, name) in TARGETS.items():      # first the stage that has a record, and the check
        c = cells(a.work, g)
        by, draws = cycle.load(a.observe / spec)
        E, U, REF = {}, {}, {}
        E["Z70"], U["Z70"], REF["Z70"], extra = one(c, events, pooled, by, draws, "Z70")
        d = 0.0 if a.made_up else max(abs(E["Z70"]["arms"][arm]["u"] - rec[name]["arms"][arm]["u"]) for arm in ARMS)
        worst = max(worst, d)
        keep[t] = (c, by, draws, E, U, REF, extra, d)
        print(f"{t}: rows at Z {E['Z70']['arms']['rows']['u']:.4f}" + ("" if a.made_up else f" (recorded {rec[name]['arms']['rows']['u']:.4f}); largest difference of an arm {d:.5f}"))
    res["checks_against_the_records"] = {"largest_difference_of_an_arm_at_Z": worst, "by_target": {t: x[7] for t, x in keep.items()}, "limit": 0.005, "hold": bool(worst <= 0.005)}
    if worst > 0.005:
        raise SystemExit(f"the readings at Z do not repeat the recorded ones ({worst:.4f}): stopped, nothing written (DESIGN.md, Section 8)")
    # --- the guard: the head as trained, at the released state and at every trained state, on the control cells
    c0 = keep["T"][0]
    ctl = ~c0["treated"]
    H0 = head(W, c0["raw"]["Z70"])                                # HEAD_LAT, the log mean, the mean
    O0, M0 = H0[1].astype(np.float64), H0[2].astype(np.float64)   # an output of the head is its log mean: that is what is recorded
    j = int(np.abs(O0 - c0["o"][:, None]).max(0).argmin())
    own = float(np.abs(O0[:, j] - c0["o"]).max())
    own_ok = bool(a.made_up or own < 1e-3 * max(1.0, float(np.abs(c0["o"]).max())))
    if not own_ok:
        raise SystemExit(f"the head as trained does not repeat its recorded output at the released state (largest difference {own:.4f}): stopped, nothing written")
    sd, sdm = O0[ctl].std(0), M0[ctl].std(0)
    rest = np.arange(O0.shape[1]) != j
    OT, res["guard"] = {}, {"limit": GUARD, "the_heads_own_output_is_reproduced": own_ok, "largest_difference_from_the_recorded_output": own, "column_of_the_target": j,
                            "scale": "the log mean, which is what the head's output is recorded as; the same on the scale of the mean is given beside it and is not judged",
                            "other_outputs": int(rest.sum()), "by_arm": {}}
    for st in TRAINED:
        HT = head(W, c0["raw"][st])
        OT[st] = HT[1].astype(np.float64)
        per_cell = float((np.abs(OT[st][ctl] - O0[ctl])[:, rest] / sd[rest]).mean())
        of_mean = float((np.abs(OT[st][ctl].mean(0) - O0[ctl].mean(0))[rest] / sd[rest]).mean())
        per_cell_m = float((np.abs(HT[2].astype(np.float64)[ctl] - M0[ctl])[:, rest] / sdm[rest]).mean())
        res["guard"]["by_arm"][ARM_OF[st]] = {"stage": st, "other_outputs_on_control_cells_mean_change_of_a_cell": per_cell, "mean_change_of_the_mean": of_mean,
                                               "the_same_on_the_scale_of_the_mean": per_cell_m, "holds": bool(per_cell <= GUARD)}
    for t, (g, spec, events, pooled, name) in TARGETS.items():
        c, by, draws, E, U, REF, (key, scope, m, cl), _d = keep[t]
        for st in TRAINED:
            E[st], U[st], REF[st], _x = one(c, events, pooled, by, draws, st)
        vin = draws[f"{key}|shift|V|{scope}|{by[key]['V'][scope]['class']}"]
        x = {"events": events, "read_as": "a pool" if pooled else "on its own, with its own record", "reference_scope": scope,
             "input_view": {q: by[key]["V"][scope][q] for q in ("u", "ci", "verdict", "class")}, "stages": {}}
        for st in ["Z70"] + TRAINED:
            e = E[st]
            s = {"arm_of_the_trial": ARM_OF.get(st), "reference": e["reference"], "reader_of_the_stage": e["reader_of_the_stage"], "from_the_input_view": cycle.diff(REF[st], vin),
                 "as_trained": e["as_trained"], "arms": e["arms"], "A3": e["A3"], "gate3": e["gate3"], "grade": e["closure"]["grade"]}
            if st != "Z70":
                s["gain_from_Z"] = {**{arm: gain(U["Z70"][arm], U[st][arm]) for arm in ARMS}, "reference_reader": gain(REF["Z70"], REF[st])}
                s["grade_kept"] = bool(RANK.get(s["grade"], 0) >= RANK.get(E["Z70"]["closure"]["grade"], 0))
                ih = np.isin(c0["ids"], c["ids"][m])              # the head as trained, reading the trained state
                assert np.array_equal(c0["ids"][ih], c["ids"][m])
                uo = cycle.output_u(c["y"][m], OT[st][ih][:, j], c["treated"][m], cl, SEED)
                s["the_head_as_trained_at_this_state"] = {"u": float(uo[0][0]), "ci": cycle.interval(uo[0][1:]), "gain_from_Z": gain(np.asarray(draws[f"{key}|shift|system|frozen|system"], np.float64), uo[0])}
            x["stages"][st] = s
        pairs = {"T5_against_TM": ("T5", "TM"), "T5_against_TC": ("T5", "TC"), "T4_against_TP": ("T4", "TP"), "T4_against_T5": ("T4", "T5")}
        x["between_arms"] = {k: {**{arm: against(U[p][arm], U[q][arm]) for arm in ARMS}, "reference_reader": against(REF[p], REF[q])} for k, (p, q) in pairs.items() if p in TRAINED and q in TRAINED}
        res["targets"][t] = x
    T, C = res["targets"]["T"]["stages"], res["targets"]["C"]["stages"]
    word = lambda r: {"a step": "confirmed", "no step": "refuted"}.get(r, "unresolved")
    h2, r3 = T["T5"]["from_the_input_view"]["reading"], T["T5"]["gain_from_Z"]["rows"]["reading"]
    res["hypotheses"] = {
        "H1_without_its_own_record": {"reading": word(T["T4"]["gain_from_Z"]["controls"]["reading"]), "gain": T["T4"]["gain_from_Z"]["controls"],
                                      "the_output_at_the_trained_state": v(T["T4"]["arms"]["controls"]), "closed": bool(T["T4"]["arms"]["controls"]["verdict"] == "support")},
        "H2_the_encoders_loss_with_the_record": {"reading": {"no added loss": "confirmed", "loss": "refuted"}.get(h2, "unresolved") if T["Z70"]["from_the_input_view"]["reading"] == "loss" else "not testable: the released state is no loss from the input view",
                                                 "from_the_input_view": T["T5"]["from_the_input_view"], "the_released_state_from_the_input_view": T["Z70"]["from_the_input_view"],
                                                 "gain_of_the_reference_reader": T["T5"]["gain_from_Z"]["reference_reader"]},
        "H3_the_repair_with_the_record": {"reading": "confirmed" if (T["T5"]["grade"] in RANK and r3 == "a step") else ("refuted" if r3 == "no step" else "unresolved"),
                                          "gain": T["T5"]["gain_from_Z"]["rows"], "grade": T["T5"]["grade"]},
        "the_control_target_is_not_harmed": {ARM_OF[st]: bool(C[st]["gain_from_Z"]["rows"]["no_harm"]) for st in TRAINED},
        "the_control_target_keeps_its_grade": {ARM_OF[st]: C[st]["grade_kept"] for st in TRAINED},
        "the_guard_on_other_outputs_holds": {arm: g_["holds"] for arm, g_ in res["guard"]["by_arm"].items()},
    }
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda q: f"{q['u']:.3f} [{q['ci'][0]:.3f}, {q['ci'][1]:.3f}] {q['verdict']}"
    for t, x in res["targets"].items():
        print(f"\n== {t} {x['events']} ({x['read_as']}); the input view {f(x['input_view'])}")
        for st, s in x["stages"].items():
            print(f"  {st:3s} {str(s['arm_of_the_trial']):6s} reference {f(s['reference'])} | from the input {s['from_the_input_view']['added']:+.3f} {s['from_the_input_view']['reading']} | grade {s['grade']} | gate 3 {s['gate3']['passes']}")
            print("        " + " | ".join(f"{arm} {f(q)}" for arm, q in s["arms"].items()))
            if st != "Z70":
                print("        gain from Z: " + " | ".join(f"{k} {q['gain']:+.3f} [{q['ci'][0]:.3f}, {q['ci'][1]:.3f}] {q['reading']}" for k, q in s["gain_from_Z"].items())
                      + f" | the head as trained {s['the_head_as_trained_at_this_state']['u']:.3f}")
        for k, q in x["between_arms"].items():
            print(f"  {k}: " + " | ".join(f"{arm} {d['change']:+.3f} [{d['ci'][0]:.3f}, {d['ci'][1]:.3f}] {d['reading']}" for arm, d in q.items()))
    print("\nguard:", json.dumps(res["guard"]["by_arm"]))
    print("hypotheses:", json.dumps({k: (q["reading"] if isinstance(q, dict) and "reading" in q else q) for k, q in res["hypotheses"].items()}))


if __name__ == "__main__":
    main()
