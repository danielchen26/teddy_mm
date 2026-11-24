#!/usr/bin/env python3
"""Do the larger checkpoints of TEDDY-G lose the response of IRF1 as the small one does: the reading, made once.

  python read.py --observe OBSERVE        (writes results/sizes.json; refuses to run when it exists)

OBSERVE: the grids of `anm_diag scope` for the two SPEC files of build_sizes.py, one folder per SPEC. The reading code is
anm_diag.cycle, unchanged. DESIGN.md, Section 4.
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUN, RUNS = HERE.parent, HERE.parents[1]
sys.path.insert(0, str(RUNS.parent))
from anm_diag import cycle  # noqa: E402

BIG = (160, 400)
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
    return out, d


def size(S, d, m):
    """What is read about one of the larger checkpoints. A gain or a closure counts only if the state reduced to the
    width of the small model shows it too."""
    z, zr = f"Z{m}", f"Z{m}r"
    x = {"from_the_input": S[z]["from_the_input"], "from_the_input_512_components": S[zr]["from_the_input"],
         "against_the_small_model": against(d(z), d("Z70")), "against_the_small_model_512_components": against(d(zr), d("Z70")),
         "entry_from_the_input": S[f"E{m}"]["from_the_input"], "last_layer": against(d(z), d(f"P{m}")),
         "annotation_position_against_the_state": against(d(f"C{m}"), d(z)), "annotation_position_from_the_input": S[f"C{m}"]["from_the_input"],
         "annotation_position_against_the_small_model_s": against(d(f"C{m}"), d("C70")), "annotation_position_512_components_against_the_small_model_s": against(d(f"C{m}r"), d("C70")),
         "placeholders_against_none": against(d(f"ZP{m}"), d(z))}
    a, ar = x["against_the_small_model"]["reading"], x["against_the_small_model_512_components"]["reading"]
    f, fr = x["from_the_input"]["reading"], x["from_the_input_512_components"]["reading"]
    if f == "no added loss" and fr == "no added loss":
        word = "does not lose it"
    elif a == "brings back a step" and ar == "brings back a step":
        word = "loses less, by a step"
    elif a == "costs a step":
        word = "loses more, by a step"
    elif f == "loss" and a == "no change of a step":
        word = "loses it as the small model does"
    else:
        word = "unresolved"
    x["statement"] = word
    x["width_only"] = bool((f == "no added loss" and fr != "no added loss") or (a == "brings back a step" and ar != "brings back a step"))
    return x


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--results", type=Path, default=RUN / "results", help="another folder for a trial on made-up values")
    a = ap.parse_args()
    out_file = a.results / "sizes.json"
    if out_file.exists():
        raise SystemExit("sizes.json exists: this reading was made once and is not made again")
    rec_s = json.loads((RUNS / "encoder_stages_20261124" / "results" / "stages.json").read_text())["IRF1"]["stages"]
    rec_m = json.loads((RUNS / "annotation_input_20261124" / "results" / "modes.json").read_text())["IRF1"]["stages"]
    recorded = {"V": rec_s["V"]["u"], "E70": rec_s["E0"]["u"], "P70": rec_s["L11"]["u"], "Z70": rec_s["L12"]["u"], "C70": rec_m["CP"]["u"]}
    res = {"margin": cycle.MARGIN, "step_margin": cycle.STEP, "note": "exploratory; development cells; one event (IRF1); the larger checkpoints have no trained head, so only what a reader of a stage can read is asked"}
    by, draws = cycle.load(a.observe / "T_thp1_IRF1")
    S, d = stages(by, draws, "IRF1", "own")
    res["IRF1"] = {"scope": "own", "stages": S}
    res["checks"] = {"repeat_the_records": {k: abs(S[k]["u"] - r) for k, r in recorded.items()}, "limit": 0.005, "control_view_reads_nothing": bool(abs(S["PERM"]["u"]) >= 0.9), "control_view": S["PERM"]["u"]}
    res["checks"]["hold"] = bool(max(res["checks"]["repeat_the_records"].values()) <= 0.005 and res["checks"]["control_view_reads_nothing"])
    res["small_model"] = {"from_the_input": S["Z70"]["from_the_input"], "entry_from_the_input": S["E70"]["from_the_input"], "last_layer": against(d("Z70"), d("P70")),
                          "annotation_position_against_the_state": against(d("C70"), d("Z70"))}
    res["sizes"] = {str(m): size(S, d, m) for m in BIG}
    res["largest_against_middle"] = against(d("Z400"), d("Z160"))
    # the same by the general code: what is lost from the input to the state, for every size
    res["general_code"] = {st: cycle.place(by, draws, st, inputs=["V"], control="PERM")["events"]["IRF1"]["from_the_input_to_the_model_stage"] for st in ("Z70", "Z160", "Z400")}
    # --- the control: four pathway genes, a pool, linear readers
    byc, drc = cycle.load(a.observe / "C_thp1a_pool")
    C, _d = stages(byc, drc, "pool", "units")
    bad = [st for st in ("Z70", "Z160", "Z400") if C[st]["from_the_input"]["reading"] == "loss"]
    res["control"] = {"events": ["IFNGR1", "IFNGR2", "JAK2", "STAT1"], "scope": "units", "readers": "ridge", "stages": C, "states_that_lose_for_the_control": bad, "holds": not bad}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda x: f"{x['u']:.3f} [{x['ci'][0]:.3f}, {x['ci'][1]:.3f}] {x['verdict']}"
    g = lambda q: "" if q is None else f"{q.get('added', q.get('change')):+.3f} [{q['ci'][0]:.3f}, {q['ci'][1]:.3f}] {q['reading']}"
    print("checks:", res["checks"])
    for st in S:
        x = S[st]
        print(f"  {st:6s} {f(x):38s} ({x['reader']:5s}) from the input {g(x.get('from_the_input')):40s} || control {C[st]['u']:.3f} {g(C[st].get('from_the_input'))}")
    print("small model:", {k: g(q) for k, q in res["small_model"].items()})
    for m in BIG:
        x = res["sizes"][str(m)]
        print(f"{m}M: {x['statement']} | width only: {x['width_only']}")
        for k, q in x.items():
            if isinstance(q, dict):
                print(f"     {k:62s} {g(q)}")
    print("400M against 160M:", g(res["largest_against_middle"]))
    print("control holds:", res["control"]["holds"], bad)


if __name__ == "__main__":
    main()
