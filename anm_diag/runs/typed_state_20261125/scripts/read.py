#!/usr/bin/env python3
"""The reading of the typed states: results/typed.json. Made once.

  python read.py --work WORK --observe OBSERVE [--results RESULTS] [--made-up]

WORK: the folder of build_ty.py. OBSERVE: the grids of `python -m anm_diag scope` for its two SPEC files.

For IRF1 on its own and for the four pathway genes as a pool, at every stage:
  the tool's reference reader of the stage (the linear class), and what it adds to the input view;
  the round with the head reading the stage (anm_diag.cycle, the arms controls, rows, cells), and its grade;
  against the state a head reads today (Z70), on paired resamples: the gain of every output.
What the round reads: at Z70 the state at unit length, as in every round before. A typed state is given to it
divided by the length of that same mean state, so its blocks add up to what the round reads at Z70 (the blocks of EA
are of what enters the first layer and are divided by the same length); a column that only says whether a gene is
among the tokens is left as it is.

Checks (the script stops, and writes nothing, if one fails): the round at Z70 repeats the first outer step, and the
reference readers of Z70 and of the input view repeat the last reading (runs/last_layer_trial_20261124), to 0.005.
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
from anm_diag import cycle  # noqa: E402
from _common import SEED, unit  # noqa: E402

TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
TARGETS = {"T": ("thp1", "T_thp1_IRF1", ["IRF1"], False, "PDL1 thp1 type B"), "C": ("thp1a", "C_thp1a_pool", TYPE_A, True, "PDL1 thp1 type A")}
TYPED, FLAGS = ["TK8", "TA", "TR", "G", "ZG", "EA"], {"G": 1, "ZG": 1}
WHAT = {"TK8": "8 types of genes from the model's gene embedding", "TA": "the 128 genes the pathway's events move most, and the rest", "TR": "128 genes at random, and the rest: the control of TA",
        "G": "the position of the measured protein's gene", "ZG": "the mean beside the position of the measured protein's gene", "EA": "the typing of TA, of what enters the first layer"}
ARMS, RANK, STEP = ("controls", "rows", "cells"), {"absolute": 2, "relative": 1}, cycle.STEP
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
    d = np.load(work / "thp1_ty.npz")
    nz = np.linalg.norm(d["Z70"].astype(np.float64), axis=1, keepdims=True) + 1e-6
    c = dict(ids=d["ids"].astype(str), treated=d["A"].astype(int) == 1, unit=d["unit"].astype(str), y=d["y_PDL1"].astype(np.float64), o=d["frozen_PDL1"].astype(np.float64),
             Z70=unit(d["Z70"]).astype(np.float64))
    for st in TYPED:
        x, k = d[st].astype(np.float64), FLAGS.get(st, 0)
        x[:, :x.shape[1] - k] /= nz                              # the blocks of a typed state add up to the state the round reads at Z70
        c[st] = x
    c["first"] = cycle.halves(np.where(c["treated"], c["unit"], "control"), SEED)      # the halves of the earlier rounds, made on the whole stratum
    if g == "thp1a":
        keep = np.isin(c["ids"], np.load(work / "thp1a_ty.npz")["ids"].astype(str))
        c = {k: val[keep] for k, val in c.items()}
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
    return e, U, np.asarray(ref, np.float64), (key, scope)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results", help="another folder for a trial on made-up values")
    ap.add_argument("--made-up", action="store_true")
    a = ap.parse_args()
    out_file = a.results / "typed.json"
    if out_file.exists():
        raise SystemExit("typed.json exists: this reading was made once and is not made again")
    res = {"margin": cycle.MARGIN, "step_margin": STEP, "note": "exploratory; development cells; nothing is trained; the typed states are the encoder as released, kept by type", "stages": WHAT, "targets": {}}
    if a.made_up:
        res["MADE_UP"] = "a trial of the script on made-up values; not a result"
        rec = last = None
    else:
        rec = json.loads((RUNS / "outer_step_20261122" / "results" / "bridge.json").read_text())["targets"]
        last = json.loads((RUNS / "last_layer_trial_20261124" / "results" / "trial.json").read_text())["targets"]
    keep, worst = {}, {"an_arm_of_the_round_at_Z": 0.0, "the_reference_reader_at_Z": 0.0, "the_input_view": 0.0}
    for t, (g, spec, events, pooled, name) in TARGETS.items():      # first the stage that has a record, and the checks
        c = cells(a.work, g)
        by, draws = cycle.load(a.observe / spec)
        E, U, REF = {}, {}, {}
        E["Z70"], U["Z70"], REF["Z70"], (key, scope) = one(c, events, pooled, by, draws, "Z70")
        if not a.made_up:
            worst["an_arm_of_the_round_at_Z"] = max(worst["an_arm_of_the_round_at_Z"], max(abs(E["Z70"]["arms"][arm]["u"] - rec[name]["arms"][arm]["u"]) for arm in ARMS))
            worst["the_reference_reader_at_Z"] = max(worst["the_reference_reader_at_Z"], abs(E["Z70"]["reference"]["u"] - last[t]["stages"]["Z70"]["reference"]["u"]))
            worst["the_input_view"] = max(worst["the_input_view"], abs(by[key]["V"][scope]["u"] - last[t]["input_view"]["u"]))
        keep[t] = (c, by, draws, E, U, REF, key, scope)
        print(f"{t}: rows at Z {E['Z70']['arms']['rows']['u']:.4f}, reference {E['Z70']['reference']['u']:.4f}, the input view {by[key]['V'][scope]['u']:.4f}")
    res["checks_against_the_records"] = {**worst, "limit": 0.005, "hold": bool(max(worst.values()) <= 0.005)}
    if max(worst.values()) > 0.005:
        raise SystemExit(f"the readings at Z or of the input view do not repeat the recorded ones ({json.dumps(worst)}): stopped, nothing written")
    for t, (g, spec, events, pooled, name) in TARGETS.items():
        c, by, draws, E, U, REF, key, scope = keep[t]
        for st in TYPED:
            E[st], U[st], REF[st], _x = one(c, events, pooled, by, draws, st)
            print(f"  {t} {st} read")
        vin = draws[f"{key}|shift|V|{scope}|{by[key]['V'][scope]['class']}"]
        x = {"events": events, "read_as": "a pool" if pooled else "on its own, with its own record", "reference_scope": scope,
             "input_view": {q: by[key]["V"][scope][q] for q in ("u", "ci", "verdict", "class")},
             "permuted_view": {q: by[key]["PERM"][scope][q] for q in ("u", "ci", "verdict")}, "stages": {}}
        for st in ["Z70"] + TYPED:
            e = E[st]
            s = {"what": WHAT.get(st, "the mean over the gene positions: what a head reads today"), "width": int(c[st].shape[1]), "reference": e["reference"], "reader_of_the_stage": e["reader_of_the_stage"],
                 "from_the_input_view": cycle.diff(REF[st], vin), "as_trained": e["as_trained"], "arms": e["arms"], "A3": e["A3"], "gate3": e["gate3"], "grade": e["closure"]["grade"]}
            if st != "Z70":
                s["gain_from_Z"] = {**{arm: gain(U["Z70"][arm], U[st][arm]) for arm in ARMS}, "reference_reader": gain(REF["Z70"], REF[st])}
                s["grade_kept"] = bool(RANK.get(s["grade"], 0) >= RANK.get(E["Z70"]["closure"]["grade"], 0))
            x["stages"][st] = s
        pairs = {"TR_against_TA": ("TR", "TA"), "G_against_ZG": ("G", "ZG"), "ZG_against_TA": ("ZG", "TA"), "TK8_against_TA": ("TK8", "TA"), "EA_against_TA": ("EA", "TA")}
        x["between_stages"] = {k: {**{arm: against(U[p][arm], U[q][arm]) for arm in ARMS}, "reference_reader": against(REF[p], REF[q])} for k, (p, q) in pairs.items()}
        res["targets"][t] = x
    T, C = res["targets"]["T"]["stages"], res["targets"]["C"]["stages"]
    word = lambda r: {"a step": "confirmed", "no step": "refuted"}.get(r, "unresolved")
    closes = lambda r: {"no added loss": "confirmed", "loss": "refuted"}.get(r, "unresolved")

    def four(st):
        g_, s = T[st]["gain_from_Z"], T[st]
        return {"the_state_keeps_it": {"reading": word(g_["reference_reader"]["reading"]), "gain_of_the_reference_reader": g_["reference_reader"], "the_reference_reader": s["reference"]},
                "the_encoders_loss_is_closed": {"reading": closes(s["from_the_input_view"]["reading"]) if T["Z70"]["from_the_input_view"]["reading"] == "loss" else "not testable: the state a head reads today is no loss from the input view",
                                                "from_the_input_view": s["from_the_input_view"]},
                "the_repair_with_the_record": {"reading": "confirmed" if (s["grade"] in RANK and g_["rows"]["reading"] == "a step") else ("refuted" if g_["rows"]["reading"] == "no step" else "unresolved"),
                                               "gain": g_["rows"], "the_output": v(s["arms"]["rows"]), "grade": s["grade"]},
                "without_its_own_record": {"reading": word(g_["controls"]["reading"]), "gain": g_["controls"], "the_output": v(s["arms"]["controls"]), "closed": bool(s["arms"]["controls"]["verdict"] == "support")}}

    res["hypotheses"] = {
        "TA": four("TA"), "ZG": four("ZG"),
        "the_choice_of_genes_matters": {"reading": res["targets"]["T"]["between_stages"]["TR_against_TA"]["reference_reader"]["reading"], "TR_against_TA": res["targets"]["T"]["between_stages"]["TR_against_TA"]},
        "the_control_target_is_not_harmed": {st: bool(C[st]["gain_from_Z"]["rows"]["no_harm"]) for st in TYPED},
        "the_control_target_keeps_its_grade": {st: C[st]["grade_kept"] for st in TYPED},
    }
    a.results.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda q: f"{q['u']:.3f} [{q['ci'][0]:.3f}, {q['ci'][1]:.3f}] {q['verdict']}"
    for t, x in res["targets"].items():
        print(f"\n== {t} {x['events']} ({x['read_as']}); the input view {f(x['input_view'])}; the permuted view {x['permuted_view']['u']:.3f}")
        for st, s in x["stages"].items():
            print(f"  {st:4s} [{s['width']:5d}] reference {f(s['reference'])} | from the input {s['from_the_input_view']['added']:+.3f} {s['from_the_input_view']['reading']} | grade {s['grade']} | gate 3 {s['gate3']['passes']}")
            print("        " + " | ".join(f"{arm} {f(q)}" for arm, q in s["arms"].items()))
            if st != "Z70":
                print("        gain from Z: " + " | ".join(f"{k} {q['gain']:+.3f} [{q['ci'][0]:.3f}, {q['ci'][1]:.3f}] {q['reading']}" for k, q in s["gain_from_Z"].items()))
        for k, q in x["between_stages"].items():
            print(f"  {k}: " + " | ".join(f"{arm} {d['change']:+.3f} [{d['ci'][0]:.3f}, {d['ci'][1]:.3f}] {d['reading']}" for arm, d in q.items()))
    print("\nhypotheses:", json.dumps({k: ({kk: qq["reading"] for kk, qq in q.items()} if k in ("TA", "ZG") else (q["reading"] if isinstance(q, dict) and "reading" in q else q)) for k, q in res["hypotheses"].items()}))


if __name__ == "__main__":
    main()
