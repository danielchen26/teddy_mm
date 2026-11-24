#!/usr/bin/env python3
"""Does the response of IRF1 come back when TEDDY is run with its annotation positions: the reading, made once.

  python read.py --observe OBSERVE        (writes results/modes.json; refuses to run when it exists)

OBSERVE: the grids of `anm_diag scope` for the two SPEC files of build_modes.py, one folder per SPEC. The reading code is
anm_diag.cycle, unchanged. DESIGN.md, Section 4.
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
sys.path.insert(0, str(HERE.parents[2]))
from anm_diag import cycle  # noqa: E402

RECORDED = {"V": 0.18370169645910994, "ZN": 0.3958111239405937}      # runs/encoder_stages_20261124/results/stages.json (V and L12)
MAIN, SIDE, POSITIONS = ["ZP", "ZT"], ["ZD", "ZW"], ["AP", "AT", "CP", "ZAP"]
STAGES = ["V", "ZN"] + MAIN + SIDE + POSITIONS + ["PERM"]
v = lambda x: {k: x[k] for k in ("u", "ci", "verdict")}


def against(s, ref):
    """A stage against another on paired resamples, both ways, by the tool's step reading: |u| at s minus |u| at ref."""
    cost, back = cycle.diff(s, ref), cycle.diff(ref, s)
    word = ("brings back a step" if back["reading"] == "loss" else "costs a step" if cost["reading"] == "loss"
            else "no change of a step" if back["reading"] == cost["reading"] == "no added loss" else "unresolved")
    return {"change": cost["added"], "ci": cost["ci"], "reading": word}


def stages(by, draws, key, scope):
    e = {"by_stage": {st: {sc: cycle._classes(by[key][st][sc]) for sc in ("units", "own") if sc in by[key][st]} for st in by[key] if st != "system"}}
    dr = lambda st, c: draws[f"{key}|shift|{st}|{scope}|{c}"]
    R = {st: cycle._reader(e, [st], scope) for st in e["by_stage"]}
    d = lambda st: dr(st, R[st][1])
    out = {}
    for st, (_s, c, x) in R.items():
        out[st] = {"reader": c, **v(x), "by_class": e["by_stage"][st][scope]}
        if st != "V":
            out[st]["from_the_input"] = cycle.diff(d(st), d("V"))
        if st not in ("V", "ZN", "PERM"):
            out[st]["against_no_annotations"] = against(d(st), d("ZN"))
    return out, d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results", help="another folder for a trial on made-up values")
    a = ap.parse_args()
    out_file = a.results / "modes.json"
    if out_file.exists():
        raise SystemExit("modes.json exists: this reading was made once and is not made again")
    man = json.loads((a.results / "modes_manifest.json").read_text())
    res = {"margin": cycle.MARGIN, "step_margin": cycle.STEP, "stages": STAGES, "note": "exploratory; development cells; one event (IRF1); every cell of the stratum has the same four annotation tokens"}
    by, draws = cycle.load(a.observe / "T_thp1_IRF1")
    assert sorted(st for st in by["IRF1"] if st != "system") == sorted(STAGES), sorted(by["IRF1"])
    S, d = stages(by, draws, "IRF1", "own")
    res["IRF1"] = {"scope": "own", "stages": S}
    res["checks"] = {"V_repeats_the_record": abs(S["V"]["u"] - RECORDED["V"]), "ZN_repeats_the_record": abs(S["ZN"]["u"] - RECORDED["ZN"]), "limit": 0.005,
                     "state_without_annotations_against_the_recorded_state": man["N_against_the_recorded_state"], "control_view_reads_nothing": bool(abs(S["PERM"]["u"]) >= 0.9), "control_view": S["PERM"]["u"]}
    res["checks"]["hold"] = bool(max(res["checks"]["V_repeats_the_record"], res["checks"]["ZN_repeats_the_record"]) <= 0.005 and res["checks"]["state_without_annotations_against_the_recorded_state"] <= 1e-4
                                 and res["checks"]["control_view_reads_nothing"])
    # --- the rule: what the annotation positions do to the state the head reads
    r = {st: S[st]["against_no_annotations"]["reading"] for st in MAIN}
    back = [st for st in MAIN if r[st] == "brings back a step"]
    if back:
        closed = [st for st in back if S[st]["from_the_input"]["reading"] == "no added loss"]
        still = [st for st in back if S[st]["from_the_input"]["reading"] == "loss"]
        word = "brings the response back" if closed else ("brings part of it back" if still else "brings back a step; what is left from the input is unresolved")
    elif any(x == "costs a step" for x in r.values()):
        word = "costs a step"
    elif all(x == "no change of a step" for x in r.values()):
        word = "changes nothing"
    else:
        word = "unresolved"
    res["annotations_in_the_input"] = {"against_no_annotations": r, "from_the_input": {st: S[st]["from_the_input"]["reading"] for st in ["ZN"] + MAIN}, "statement": word}
    res["content"] = {"labels_of_the_line_against_wrong_labels": against(d("ZT"), d("ZW")), "labels_of_the_line_against_the_official_defaults": against(d("ZT"), d("ZD")),
                      "labels_of_the_line_against_placeholders": against(d("ZT"), d("ZP")), "side": {st: S[st]["against_no_annotations"]["reading"] for st in SIDE}}
    res["annotation_positions"] = {st: {"against_no_annotations": S[st]["against_no_annotations"]["reading"], "from_the_input": S[st]["from_the_input"]["reading"]} for st in POSITIONS}
    res["annotation_positions"]["statement"] = {"brings back a step": "the model's own summary at the annotation positions reads a step more than the mean",
                                                "costs a step": "it reads a step less than the mean", "no change of a step": "it reads as much as the mean, within a step",
                                                "unresolved": "unresolved"}[S["AP"]["against_no_annotations"]["reading"]]
    # the same by the general code: what is lost from the input to the state, for every way of running
    res["general_code"] = {}
    for st in ["ZN"] + MAIN:
        P = cycle.place(by, draws, st, inputs=["V"], control="PERM")["events"]["IRF1"]
        res["general_code"][st] = P["from_the_input_to_the_model_stage"]
    # --- the control: four pathway genes, a pool, linear readers
    byc, drc = cycle.load(a.observe / "C_thp1a_pool")
    C, _d = stages(byc, drc, "pool", "units")
    bad = [st for st in MAIN if C[st]["from_the_input"]["reading"] == "loss" or C[st]["against_no_annotations"]["reading"] == "costs a step"]
    res["control"] = {"events": ["IFNGR1", "IFNGR2", "JAK2", "STAT1"], "scope": "units", "readers": "ridge", "stages": C, "ways_of_running_that_lose_for_the_control": bad, "holds": not bad}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda x: f"{x['u']:.3f} [{x['ci'][0]:.3f}, {x['ci'][1]:.3f}] {x['verdict']}"
    g = lambda q: "" if q is None else f"{q.get('added', q.get('change')):+.3f} [{q['ci'][0]:.3f}, {q['ci'][1]:.3f}] {q['reading']}"
    print("checks:", res["checks"])
    for st in STAGES:
        x = S[st]
        print(f"  {st:5s} {f(x):38s} ({x['reader']:5s}) from the input {g(x.get('from_the_input')):40s} | against no annotations {g(x.get('against_no_annotations')):44s} || control {C[st]['u']:.3f} {g(C[st].get('from_the_input'))}")
    print("annotations in the input:", res["annotations_in_the_input"])
    print("content:", {k: (g(x) if "ci" in x else x) for k, x in res["content"].items()})
    print("annotation positions:", res["annotation_positions"])
    print("control holds:", res["control"]["holds"], bad)


if __name__ == "__main__":
    main()
