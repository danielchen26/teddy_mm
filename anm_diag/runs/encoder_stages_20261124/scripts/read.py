#!/usr/bin/env python3
"""Where between the input and Z the response of IRF1 stops being readable: the reading, made once.

  python read.py --observe OBSERVE        (writes results/stages.json; refuses to run when it exists)

OBSERVE: the grids of `anm_diag scope` for the two SPEC files of build_stages.py, one folder per SPEC. The reading code is
anm_diag.cycle, unchanged (`place` with the stages inside the model given in order). DESIGN.md, Section 4.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))
from anm_diag import cycle  # noqa: E402
from build_stages import CHAIN  # noqa: E402

RECORDED = {"V": 0.18370169645910994, "L12": 0.3958111239405937}      # runs/input_beside_z_20261123/results/round.json (T1: the input view alone; the reference at Z)
v = lambda x: {k: x[k] for k in ("u", "ci", "verdict")}


def chain(by, draws, key, scope):
    """Every stage of one event (or of the pool) at one scope: its reader and reading, what is added from the input view,
    and what is added from the stage before it, all on paired resamples."""
    e = {"by_stage": {st: {sc: cycle._classes(by[key][st][sc]) for sc in ("units", "own") if sc in by[key][st]} for st in by[key] if st != "system"}}
    dr = lambda st, c: draws[f"{key}|shift|{st}|{scope}|{c}"]
    R = {st: cycle._reader(e, [st], scope) for st in e["by_stage"]}
    out = {}
    for st, (_s, c, x) in R.items():
        out[st] = {"reader": c, **v(x), "by_class": e["by_stage"][st][scope], "from_the_input": cycle.diff(dr(st, c), dr("V", R["V"][1])) if st != "V" else None}
    for a, b in zip(CHAIN, CHAIN[1:]):
        out[b]["from_the_stage_before"] = {"stage": a, **cycle.diff(dr(b, R[b][1]), dr(a, R[a][1]))}
    return out, R, dr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results", help="another folder for a trial on made-up values")
    a = ap.parse_args()
    out_file = a.results / "stages.json"
    if out_file.exists():
        raise SystemExit("stages.json exists: this reading was made once and is not made again")
    res = {"margin": cycle.MARGIN, "step_margin": cycle.STEP, "chain": CHAIN, "note": "exploratory; development cells; one event (IRF1); pooled summaries of a stage are not nested"}
    # --- IRF1, read on its own, with its own cells
    by, draws = cycle.load(a.observe / "T_thp1_IRF1")
    S, R, dr = chain(by, draws, "IRF1", "own")
    res["IRF1"] = {"scope": "own", "stages": S}
    res["checks"] = {"V_repeats_the_record": abs(S["V"]["u"] - RECORDED["V"]), "L12_repeats_the_record": abs(S["L12"]["u"] - RECORDED["L12"]), "L12_equals_Z": abs(S["L12"]["u"] - S["Z"]["u"]),
                     "limit": 0.005, "control_view_reads_nothing": bool(abs(S["PERM"]["u"]) >= 0.9), "control_view": S["PERM"]["u"]}
    res["checks"]["hold"] = bool(max(res["checks"][k] for k in ("V_repeats_the_record", "L12_repeats_the_record", "L12_equals_Z")) <= 0.005)
    # --- the rule: where the loss appears
    read = [(st, S[st]["from_the_input"]["reading"]) for st in CHAIN[1:]]
    first = next((st for st, r in read if r == "loss"), None)
    before = [st for st, r in read[:[s for s, _r in read].index(first)] if r == "no added loss"] if first else []
    last_clear = (before[-1] if before else "V") if first else None
    res["located"] = {"first_stage_that_is_a_loss": first, "last_stage_before_it_with_no_added_loss": last_clear, "readings_from_the_input": dict(read),
                      "single_steps_that_are_a_loss": [[S[st]["from_the_stage_before"]["stage"], st] for st in CHAIN[1:] if S[st]["from_the_stage_before"]["reading"] == "loss"]}
    where = None if first is None else ("tokenisation" if first in ("SELV", "RANK") else ("the entry: one mean of embeddings" if first == "E0" else "in the layers"))
    res["located"]["statement"] = where
    # the same by the general code (place, with the stages inside the model in order)
    P = cycle.place(by, draws, "L12", inputs=["V"], between=CHAIN[1:-1], others=["SET", "W0", "W12"], control="PERM")["events"]["IRF1"]
    res["general_code"] = {"along_the_model": P["along_the_model"], "lost_between": P["lost_between"], "from_the_input_to_the_model_stage": P["from_the_input_to_the_model_stage"]}
    # --- the mean over positions, at the entry and at the end
    res["pooling"] = {"at_the_end": {"the_single_mean_loses_against_eight_windows": cycle.diff(dr("L12", R["L12"][1]), dr("W12", R["W12"][1])), "W12_from_the_input": S["W12"]["from_the_input"]},
                      "at_the_entry": {"the_single_mean_loses_against_eight_windows": cycle.diff(dr("E0", R["E0"][1]), dr("W0", R["W0"][1])), "W0_from_the_input": S["W0"]["from_the_input"]}}
    res["pooling"]["implicated_at_the_end"] = res["pooling"]["at_the_end"]["the_single_mean_loses_against_eight_windows"]["reading"]
    # --- what the layers do with what enters them
    res["layers"] = {"L12_against_E0": cycle.diff(dr("E0", R["E0"][1]), dr("L12", R["L12"][1])),      # above the step margin: the layers bring back at least a step of what the entry misses
                     "best_layer": min(CHAIN[4:], key=lambda st: abs(S[st]["u"])), "order_matters_at_the_input": cycle.diff(dr("SET", R["SET"][1]), dr("RANK", R["RANK"][1]))}
    # --- the control: four pathway genes, a pool, linear readers
    byc, drc = cycle.load(a.observe / "C_thp1a_pool")
    C, _R, _d = chain(byc, drc, "pool", "units")
    res["control"] = {"events": ["IFNGR1", "IFNGR2", "JAK2", "STAT1"], "scope": "units", "readers": "ridge", "stages": C,
                      "stages_that_are_a_loss_from_the_input": [st for st in CHAIN[1:] if C[st]["from_the_input"]["reading"] == "loss"]}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda x: f"{x['u']:.3f} [{x['ci'][0]:.3f}, {x['ci'][1]:.3f}] {x['verdict']}"
    g = lambda d: "" if d is None else f"{d['added']:+.3f} [{d['ci'][0]:.3f}, {d['ci'][1]:.3f}] {d['reading']}"
    print("checks:", res["checks"])
    for st in CHAIN + ["SET", "W0", "W12", "Z", "PERM"]:
        x = S[st]
        print(f"  {st:5s} {f(x):38s} ({x['reader']:5s}) from the input {g(x['from_the_input']):42s} | step {g(x.get('from_the_stage_before'))} || control {C[st]['u']:.3f} {g(C[st]['from_the_input'])}")
    print("located:", {k: res["located"][k] for k in ("first_stage_that_is_a_loss", "last_stage_before_it_with_no_added_loss", "statement", "single_steps_that_are_a_loss")})
    print("pooling at the end:", g(res["pooling"]["at_the_end"]["the_single_mean_loses_against_eight_windows"]), "| at the entry:", g(res["pooling"]["at_the_entry"]["the_single_mean_loses_against_eight_windows"]))
    print("layers: L12 against E0", g(res["layers"]["L12_against_E0"]), "| best layer", res["layers"]["best_layer"], "| order at the input", g(res["layers"]["order_matters_at_the_input"]))
    print("control: stages that are a loss from the input:", res["control"]["stages_that_are_a_loss_from_the_input"])


if __name__ == "__main__":
    main()
