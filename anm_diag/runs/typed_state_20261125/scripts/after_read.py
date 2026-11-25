#!/usr/bin/env python3
"""The reading of the check added after the reading: results/after.json. Made once. Not declared before the reading.

  python after_read.py --work WORK_AFTER --observe OBSERVE_AFTER [--results RESULTS] [--made-up]

Two questions, both about what the reading of the typed states means.

 1. The genes. TA holds IRF1 and JAK2, two knocked-out genes. Does the gain of TA at the reference reader rest on
    them? `TAx` (the rule of TA, the five knocked-out genes not eligible) against `TA`, and its gain over Z70; `EAx`
    against `EA` likewise.
 2. The scale. The round's linear reader penalises every column alike and does not bring the columns to one scale;
    the tool's reference reader does. Is that why the round did not use TA? The stages at one scale (`Zs`, `TAs`,
    `TRs`, `TAxs`, `EAs`): the gain of every output of the round over Z70 as it was read, and over `Zs`, the control
    of the scale.

As in read.py: the raw stages are given to the round divided by the length of the mean state; a stage at one scale is
given to the round and to the tool as it is. The cells, the halves and the resamples are those of the reading.

Checks (the script stops, and writes nothing, if one fails): Z70 and TA read here repeat results/typed.json to 0.005
(the same states; the grids of the tool are made anew).
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import read as R  # noqa: E402
from read import ARMS, SEED, against, cycle, gain, unit  # noqa: E402

RAW, SCALED = ["TA", "EA", "TAx", "EAx"], ["Zs", "TAs", "TRs", "TAxs", "EAs"]
WHAT = {"TA": "as read: the 128 genes the pathway's events move most, and the rest", "EA": "as read: the typing of TA, of what enters the first layer",
        "TAx": "the rule of TA, the five knocked-out genes not eligible", "EAx": "the typing of TAx, of what enters the first layer",
        "Zs": "the mean, at one scale: the control of the scale", "TAs": "TA at one scale", "TRs": "TR (128 genes at random) at one scale", "TAxs": "TAx at one scale", "EAs": "EA at one scale"}
PAIRS = {"TAx_against_TA": ("TAx", "TA"), "EAx_against_EA": ("EAx", "EA"), "TAxs_against_TAs": ("TAxs", "TAs"), "TRs_against_TAs": ("TRs", "TAs"), "EAs_against_TAs": ("EAs", "TAs"), "Zs_against_TAs": ("Zs", "TAs")}


def cells(work, g):
    d = np.load(work / "thp1_af.npz")
    nz = np.linalg.norm(d["Z70"].astype(np.float64), axis=1, keepdims=True) + 1e-6
    c = dict(ids=d["ids"].astype(str), treated=d["A"].astype(int) == 1, unit=d["unit"].astype(str), y=d["y_PDL1"].astype(np.float64), o=d["frozen_PDL1"].astype(np.float64),
             Z70=unit(d["Z70"]).astype(np.float64))
    for st in RAW:
        c[st] = d[st].astype(np.float64) / nz
    for st in SCALED:
        c[st] = d[st].astype(np.float64)
    c["first"] = cycle.halves(np.where(c["treated"], c["unit"], "control"), SEED)
    if g == "thp1a":
        keep = np.isin(c["ids"], np.load(work / "thp1a_af.npz")["ids"].astype(str))
        c = {k: val[keep] for k, val in c.items()}
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=R.RUN / "results", help="another folder for a trial on made-up values")
    ap.add_argument("--made-up", action="store_true")
    a = ap.parse_args()
    out_file = a.results / "after.json"
    if out_file.exists():
        raise SystemExit("after.json exists: this reading was made once and is not made again")
    res = {"note": "a check added after the reading of the typed states; not declared before it; exploratory", "margin": cycle.MARGIN, "step_margin": R.STEP, "stages": WHAT, "targets": {}}
    was = None
    if a.made_up:
        res["MADE_UP"] = "a trial of the script on made-up values; not a result"
    else:
        was = json.loads((R.RUN / "results" / "typed.json").read_text())["targets"]
    keep, worst = {}, 0.0
    for t, (g, spec, events, pooled, name) in R.TARGETS.items():
        c = cells(a.work, g)
        by, draws = cycle.load(a.observe / spec)
        E, U, REF = {}, {}, {}
        for st in ["Z70", "TA"]:
            E[st], U[st], REF[st], (key, scope) = R.one(c, events, pooled, by, draws, st)
            if was is not None:
                w = was[t]["stages"][st]
                worst = max(worst, abs(E[st]["reference"]["u"] - w["reference"]["u"]), max(abs(E[st]["arms"][arm]["u"] - w["arms"][arm]["u"]) for arm in ARMS))
        keep[t] = (c, by, draws, E, U, REF, key, scope)
        print(f"{t}: at Z rows {E['Z70']['arms']['rows']['u']:.4f} reference {E['Z70']['reference']['u']:.4f}; at TA rows {E['TA']['arms']['rows']['u']:.4f} reference {E['TA']['reference']['u']:.4f}")
    res["check_against_the_reading"] = {"largest_difference_at_Z70_and_TA": worst, "limit": 0.005, "hold": bool(worst <= 0.005)}
    if worst > 0.005:
        raise SystemExit(f"Z70 and TA do not repeat the reading ({worst}): stopped, nothing written")
    for t, (g, spec, events, pooled, name) in R.TARGETS.items():
        c, by, draws, E, U, REF, key, scope = keep[t]
        for st in ["EA", "TAx", "EAx"] + SCALED:
            E[st], U[st], REF[st], _x = R.one(c, events, pooled, by, draws, st)
            print(f"  {t} {st} read")
        vin = draws[f"{key}|shift|V|{scope}|{by[key]['V'][scope]['class']}"]
        x = {"events": events, "read_as": "a pool" if pooled else "on its own, with its own record", "input_view": {q: by[key]["V"][scope][q] for q in ("u", "ci", "verdict")}, "stages": {}}
        for st in ["Z70"] + RAW + SCALED:
            e = E[st]
            s = {"what": WHAT.get(st, "the mean over the gene positions: what a head reads today"), "width": int(c[st].shape[1]), "reference": e["reference"], "from_the_input_view": cycle.diff(REF[st], vin),
                 "arms": e["arms"], "grade": e["closure"]["grade"]}
            if st != "Z70":
                s["gain_from_Z"] = {**{arm: gain(U["Z70"][arm], U[st][arm]) for arm in ARMS}, "reference_reader": gain(REF["Z70"], REF[st])}
            if st in SCALED and st != "Zs":
                s["gain_from_Zs"] = {**{arm: gain(U["Zs"][arm], U[st][arm]) for arm in ARMS}, "reference_reader": gain(REF["Zs"], REF[st])}
            x["stages"][st] = s
        x["between_stages"] = {k: {**{arm: against(U[p][arm], U[q][arm]) for arm in ARMS}, "reference_reader": against(REF[p], REF[q])} for k, (p, q) in PAIRS.items()}
        res["targets"][t] = x
    T, B = res["targets"]["T"]["stages"], res["targets"]["T"]["between_stages"]
    res["questions"] = {
        "1_the_genes": {"TAx_against_TA_at_the_reference_reader": B["TAx_against_TA"]["reference_reader"], "the_gain_of_TAx_over_Z_at_the_reference_reader": T["TAx"]["gain_from_Z"]["reference_reader"],
                        "the_gain_of_TA_over_Z_at_the_reference_reader": T["TA"]["gain_from_Z"]["reference_reader"], "EAx_against_EA_at_the_reference_reader": B["EAx_against_EA"]["reference_reader"]},
        "2_the_scale": {arm: {"TAs_over_Z_as_read": T["TAs"]["gain_from_Z"][arm], "Zs_over_Z_as_read": T["Zs"]["gain_from_Z"][arm], "TAs_over_Zs": T["TAs"]["gain_from_Zs"][arm],
                              "TRs_over_Zs": T["TRs"]["gain_from_Zs"][arm], "TAxs_over_Zs": T["TAxs"]["gain_from_Zs"][arm], "EAs_over_Zs": T["EAs"]["gain_from_Zs"][arm]} for arm in ARMS},
        "the_control_target_is_not_harmed_with_its_records": {st: bool(res["targets"]["C"]["stages"][st]["gain_from_Z"]["rows"]["no_harm"]) for st in RAW + SCALED},
    }
    a.results.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda q: f"{q['u']:.3f} [{q['ci'][0]:.3f}, {q['ci'][1]:.3f}] {q['verdict']}"
    gq = lambda q: f"{q['gain']:+.3f} [{q['ci'][0]:.3f}, {q['ci'][1]:.3f}] {q['reading']}"
    for t, x in res["targets"].items():
        print(f"\n== {t} {x['events']} ({x['read_as']}); the input view {f(x['input_view'])}")
        for st, s in x["stages"].items():
            print(f"  {st:5s} [{s['width']:5d}] reference {f(s['reference'])} | from the input {s['from_the_input_view']['added']:+.3f} {s['from_the_input_view']['reading']} | grade {s['grade']}")
            print("        " + " | ".join(f"{arm} {f(q)}" for arm, q in s["arms"].items()))
            if "gain_from_Z" in s:
                print("        gain from Z : " + " | ".join(f"{k} {gq(q)}" for k, q in s["gain_from_Z"].items()))
            if "gain_from_Zs" in s:
                print("        gain from Zs: " + " | ".join(f"{k} {gq(q)}" for k, q in s["gain_from_Zs"].items()))
        for k, q in x["between_stages"].items():
            print(f"  {k}: " + " | ".join(f"{arm} {d['change']:+.3f} [{d['ci'][0]:.3f}, {d['ci'][1]:.3f}] {d['reading']}" for arm, d in q.items()))
    print("\nnot harmed (the control target, with its records):", json.dumps(res["questions"]["the_control_target_is_not_harmed_with_its_records"]))


if __name__ == "__main__":
    main()
