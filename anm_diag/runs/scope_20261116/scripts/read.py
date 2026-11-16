#!/usr/bin/env python3
"""Read the six scope runs against the three criteria of DESIGN.md, write the hypotheses and their tests into one
ledger, and set the corrected donor-shift numbers beside the old ones. Makes no fit.

  python read.py            (reads ../results and the records of earlier runs; writes ../results/summary.json, ledger.jsonl)
"""
from __future__ import annotations

import inspect
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
RUNS = RUN.parent
sys.path.insert(0, str(RUNS.parent))
from anm_diag import events as ev  # noqa: E402
from anm_diag import loop, scope  # noqa: E402
from anm_diag.stats import ci  # noqa: E402

RES = RUN / "results"
SPECS = ["ko_melanoma_pooled", "ko_thp1_pooled", "ko_melanoma_each", "ko_thp1_each", "donors_bm", "donors_hao"]
CT = {"melanoma": "frangieh2021", "thp1": "papalexi2021"}
R = {n: json.loads((RES / n / "scope.json").read_text()) for n in SPECS}
D = {n: dict(np.load(RES / n / "scope_draws.npz")) for n in SPECS}
ALPHA, STEP = 0.05, 0.10


def cell(r, target, stage, scope_):
    return r["targets"][target]["stages"][stage][scope_]


def draws(name, unit, target, stage, scope_):
    r = next(x for x in R[name]["results"] if x["event"]["unit"] == unit)
    cls = cell(r, target, stage, scope_)["class"]
    return np.abs(D[name][f"{unit}|{target}|{stage}|{scope_}|{cls}"].astype(np.float64))


def iv(x):
    return [float(v) for v in ci(np.asarray(x)[1:], ALPHA)]


out = {}

# ------------------------------------------------------------------ G1: one code, one format
fns = [scope.Study, scope.scope_grid, scope.read, scope.forecast, scope.family, scope._u_shift, scope._u_variation, scope._Side]
out["G1"] = {
    "specs_run": {n: bool((RES / n / "SPEC.yaml").exists() and (RES / n / "scope.json").exists()) for n in SPECS},
    "kinds": {n: R[n]["kind"] for n in SPECS},
    "functions_that_compute_without_reading_the_kind": {f.__name__: "kind" not in inspect.getsource(f) for f in fns},
    "checks": {n: [f"{c['code']} {c['severity']}" for c in R[n]["checks"]] for n in SPECS},
}
out["G1"]["holds"] = bool(all(out["G1"]["specs_run"].values()) and all(out["G1"]["functions_that_compute_without_reading_the_kind"].values()))

# ------------------------------------------------------------------ G2: what is already known comes back
sealed = {ct: json.loads((RUNS / f"{ds}_headfix" / "anm_record.json").read_text())["readouts"]["PDL1"]["states"] for ct, ds in CT.items()}
g2 = {"knockouts": {}, "donors": {}}
for ct in CT:
    r = R[f"ko_{ct}_pooled"]["results"][0]
    S = r["targets"]["shift"]["stages"]
    span = r["targets"]["shift"]["spans"]["units"]["Z->HEAD_OUT"]
    k = {"units": {st: {"u": S[st]["units"]["u"], "ci": S[st]["units"]["ci"], "verdict": S[st]["units"]["verdict"], "class": S[st]["units"]["class"]}
                   for st in ("Z", "HEAD_LAT", "HEAD_OUT")},
         "version_1_sealed": {st: {"u": sealed[ct][st]["u"], "ci": sealed[ct][st]["ci"], "verdict": sealed[ct][st]["verdict"]} for st in ("Z", "HEAD_LAT", "HEAD_OUT")},
         "added_Z_to_HEAD_OUT": {"value": span["added"], "ci": span["ci"]}, "response": r["delta"], "resampled": r["event"]["resampled"]}
    k["head_out_rejected"] = S["HEAD_OUT"]["units"]["verdict"] == "reject"
    k["added_loss_above_margin"] = bool(span["ci"][0] > STEP)
    k["z_not_rejected"] = S["Z"]["units"]["verdict"] != "reject"
    k["holds"] = bool(k["head_out_rejected"] and k["added_loss_above_margin"] and (k["z_not_rejected"] or ct != "melanoma"))
    g2["knockouts"][ct] = k
for name in ("donors_hao", "donors_bm"):
    rows = []
    for r in R[name]["results"]:
        g = r["targets"]["variation"]["gaps"]
        ft, tu, uo = g["frozen->train"]["value"], g["train->units"]["value"], g["units->own"]["value"]
        rows.append({"unit": r["event"]["unit"], "frozen->train": ft, "train->units": tu, "units->own": uo,
                     "a_head_not_the_bottleneck": bool(ft < tu), "b_coverage_positive": bool(tu > 0), "c_own_smaller": bool(uo < tu)})
    g2["donors"][name] = {"rows": rows, "n": len(rows), **{k: int(sum(x[k] for x in rows)) for k in ("a_head_not_the_bottleneck", "b_coverage_positive", "c_own_smaller")}}
h = g2["donors"]["donors_hao"]
g2["donors"]["holds"] = bool(h["a_head_not_the_bottleneck"] == h["n"] and h["b_coverage_positive"] == h["n"] and h["c_own_smaller"] == h["n"])
g2["holds"] = bool(all(k["holds"] for k in g2["knockouts"].values()) and g2["donors"]["holds"])
out["G2"] = g2

# ------------------------------------------------------------------ G3: the grid against interventions already made
ledger, tests, g3 = [], [], {}
for n in SPECS:
    ledger += [json.loads(line) for line in (RES / n / "scope_events.jsonl").read_text().splitlines() if line.strip()]
A_RUNS = RUNS / "loop_20261112"
for ct, ds in CT.items():
    name = f"ko_{ct}_pooled"
    r = R[name]["results"][0]
    base = draws(name, "pool", "shift", "HEAD_OUT", "units")
    fc = {arm: iv(base - draws(name, "pool", "shift", "Z", sc_)) for arm, sc_ in (("coverage", "units"), ("other study", "other"), ("function", "train"))}
    pt = {arm: float((base - draws(name, "pool", "shift", "Z", sc_))[0]) for arm, sc_ in (("coverage", "units"), ("other study", "other"), ("function", "train"))}
    hyp = loop.hypothesis_event(
        f"ko_{ct}_coverage", claim="coverage", rivals=["function", "other study"], retrospective=True, basis=[f"scope_20261116/{name}"],
        about={"event": f"knockouts of the IFN-gamma pathway, {ct}, pooled", "readout": "PDL1", "target": "shift",
               "quantity": "fall of |u| at HEAD_OUT in the audit of the changed head (version 1 statistic)"},
        forecast=fc, rule={"own": "coverage", "at_least": 0.20, "controls": {"other study": 0.20}},
        note="written from the grid after the outcomes were known; the rule was fixed in DESIGN.md before the grid was computed")
    b, m, t = (loop.Run.load(A_RUNS / f"pooled_{ds}_as_trained"), loop.Run.load(A_RUNS / f"pooled_{ds}_mlp_only"), loop.Run.load(A_RUNS / f"transfer_{ds}_mlp_only"))
    c = loop.contrast(b, m, "PDL1", ("Z", "HEAD_OUT"))
    tr = loop.stage(t, "PDL1", "HEAD_OUT")
    s0 = abs(sealed[ct]["HEAD_OUT"]["u"])
    arms = {"coverage": {"gain": c["removed"], "ci": c["removed_ci"], "runs": [f"loop_20261112/pooled_{ds}_as_trained", f"loop_20261112/pooled_{ds}_mlp_only"],
                         "admitted": True},
            "other study": {"gain": s0 - abs(tr["u"]), "ci": [s0 - max(abs(tr["ci"][0]), abs(tr["ci"][1])), s0 - min(abs(tr["ci"][0]), abs(tr["ci"][1]))],
                            "runs": [f"{ds}_headfix", f"loop_20261112/transfer_{ds}_mlp_only"], "admitted": True}}
    te = loop.test(hyp, arms)
    te["note"] = "the 'other study' interval takes the sealed baseline at its point value (its resamples were not kept)"
    ledger += [hyp, te]
    tests.append(te)
    rd = r["targets"]["shift"]["reading"]
    g3[ct] = {"reading_of_the_grid": rd["missing"], "detail": rd["detail"], "resolved_gaps": rd["resolved_gaps"],
              "scopes_at_Z": {s: {"u": cell(r, "shift", "Z", s)["u"], "ci": cell(r, "shift", "Z", s)["ci"], "verdict": cell(r, "shift", "Z", s)["verdict"],
                                  "class": cell(r, "shift", "Z", s)["class"]} for s in ("train", "other", "units", "own")},
              "frozen": {k: cell(r, "shift", "system", "frozen")[k] for k in ("u", "ci", "verdict")},
              "gaps": {k: v for k, v in r["targets"]["shift"]["gaps"].items() if k != "from_frozen"},
              "forecast": {a: {"point": pt[a], "ci": fc[a]} for a in fc}, "arms": te["arms"], "outcome": te["outcome"], "why": te["why"]}
    g3[ct]["grid_says_coverage"] = bool(rd["missing"] == "coverage")
    g3[ct]["coverage_gap_resolved"] = bool("train->units" in rd["resolved_gaps"])
# ---- an added check, not in DESIGN.md: the same hypothesis on the head's own output (frozen), measured by frozen_arms.py
added = {}
FA = json.loads((RES / "frozen_arms.json").read_text()) if (RES / "frozen_arms.json").exists() else {}
for ct in FA:
    name = f"ko_{ct}_pooled"
    fz = np.abs(D[name]["pool|shift|system|frozen|system"].astype(np.float64))
    fc = {arm: iv(fz - draws(name, "pool", "shift", "Z", sc_)) for arm, sc_ in (("coverage", "units"), ("other study", "other"), ("function", "train"))}
    pt = {arm: float((fz - draws(name, "pool", "shift", "Z", sc_))[0]) for arm, sc_ in (("coverage", "units"), ("other study", "other"), ("function", "train"))}
    hyp = loop.hypothesis_event(
        f"ko_{ct}_coverage_own_output", claim="coverage", rivals=["function", "other study"], retrospective=True, basis=[f"scope_20261116/{name}"],
        about={"event": f"knockouts of the IFN-gamma pathway, {ct}, pooled", "readout": "PDL1", "target": "shift",
               "quantity": "fall of |u| of the head's own PD-L1 output (frozen scope), from the saved outputs of the repaired heads"},
        forecast=fc, rule={"own": "coverage", "at_least": 0.20, "controls": {"other study": 0.20}},
        note="added after the design: scope.forecast is about the system's own output, which the design's quantity was not")
    loopc, full = FA[ct]["cells of the loop (half of the controls)"], FA[ct]["all cells of the sealed audit"]
    arms = {"coverage": {"gain": loopc["arms"]["coverage"]["fall"], "ci": loopc["arms"]["coverage"]["fall_ci"], "runs": [loopc["arms"]["coverage"]["run"]], "admitted": True}}
    if "other study" in full["arms"]:
        arms["other study"] = {"gain": full["arms"]["other study"]["fall"], "ci": full["arms"]["other study"]["fall_ci"],
                               "runs": [full["arms"]["other study"]["run"]], "admitted": True}
    te = loop.test(hyp, arms)
    ledger += [hyp, te]
    tests.append(te)
    added[ct] = {"as_trained": loopc["as_trained"], "forecast": {a_: {"point": pt[a_], "ci": fc[a_]} for a_ in fc}, "arms": te["arms"], "outcome": te["outcome"],
                 "elsewhere": loopc["arms"].get("elsewhere"), "share_of_forecast_delivered": float(arms["coverage"]["gain"] / pt["coverage"])}
out["G3_added"] = added
out["G3"] = {**g3, "holds": bool(all(v["outcome"] == "confirmed" and v["grid_says_coverage"] for v in g3.values())),
             "holds_by_gap": bool(all(v["outcome"] == "confirmed" and (v["grid_says_coverage"] or v["coverage_gap_resolved"]) for v in g3.values()))}

# ------------------------------------------------------------------ the open hypothesis for the donors (not retrospective)
ext = R["donors_hao"]["results"]
rng = np.random.default_rng(0)
idx = rng.integers(0, len(ext), size=(4000, len(ext)))


def fam(key):
    v = np.array([r["targets"]["variation"]["gaps"]["from_frozen"][key]["value"] for r in ext])
    return float(v.mean()), [float(x) for x in ci(v[idx].mean(axis=1), ALPHA)], [float(v.min()), float(v.max())]


cov, fun = fam("units"), fam("train")
bound = 0.5 * cov[1][0]
dh = loop.hypothesis_event(
    "donors_hao_coverage", claim="coverage", rivals=["function", "amount"], retrospective=False, basis=["scope_20261116/donors_hao", "generalization_20261112"],
    about={"event": "a donor the head was not trained on, external PBMC", "readout": "112 proteins", "target": "variation",
           "quantity": "fall of u (variation) at the head's own output, mean over held-out donors"},
    forecast={"coverage": cov[1], "function": fun[1], "amount": [-STEP / 2, STEP / 2]},
    rule={"own": "coverage", "at_least": bound, "controls": {"function": bound, "amount": bound}},
    note="arms of a training experiment that has not been run: PBMC cells of other donors added to the head's training; the same number of "
         "bone-marrow cells instead; another map on the same bone-marrow cells. The amount forecast is a chosen band, not read from the grid.")
ledger.append(dh)
out["donor_hypothesis"] = {"id": dh["id"], "forecast": {"coverage": {"mean": cov[0], "ci": cov[1], "range": cov[2]}, "function": {"mean": fun[0], "ci": fun[1], "range": fun[2]}},
                           "rule": dh["rule"], "status": "open: no arm has been run"}

# ------------------------------------------------------------------ the map, by family
def table(name, target):
    f = R[name]["family"][target]
    return {"n_events": f["n_events"], "approximate": f["approximate"], "missing": f["missing"],
            "stages": {st: {s: {k: x[k] for k in ("mean", "ci", "min", "max", "verdicts")} for s, x in by.items()} for st, by in f["stages"].items()},
            "gaps": f["gaps"]}


out["map"] = {name: {tg: table(name, tg) for tg in R[name]["family"]} for name in ("ko_melanoma_each", "ko_thp1_each", "donors_bm", "donors_hao")}
unseen_bm = [r for r in R["donors_bm"]["results"] if r["event"]["unit"] in ("13272:site4", "19593:site4", "15078:site4", "18303:site1")]
out["map"]["donors_bm_unseen"] = {tg: scope.family(unseen_bm, tg) for tg in ("shift", "variation")}

# ------------------------------------------------------------------ the correction to runs/generalization_20261112
old = {g["id"]: g for g in json.loads((RUNS / "generalization_20261112" / "results" / "gate0.json").read_text())["groups"]}
corr = []
for name, pre in (("donors_hao", "ext:hao2021_"), ("donors_bm", "bm:")):
    for r in R[name]["results"]:
        u = r["event"]["unit"]
        o = old[pre + u]["shift"]["by"]
        S = r["targets"]["shift"]["stages"]
        corr.append({"study": name, "unit": u, "seen_by_head": bool(old[pre + u]["seen_by_head"]),
                     "old_squared_share": {"head": o["A"]["u"], "state_other_donors": o["Z_L"]["u"], "state_training_cells": o["Z_T"]["u"]},
                     "old_on_the_scale_of_u": {k: float(np.sqrt(v)) for k, v in (("head", o["A"]["u"]), ("state_other_donors", o["Z_L"]["u"]),
                                                                                 ("state_training_cells", o["Z_T"]["u"]))},
                     "new_u": {"head": {k: S["system"]["frozen"][k] for k in ("u", "ci", "verdict")},
                               "state_other_donors": {k: S["Z"]["units"][k] for k in ("u", "ci", "verdict")},
                               "state_training_cells": {k: S["Z"]["train"][k] for k in ("u", "ci", "verdict")},
                               "state_own_cells": {k: S["Z"]["own"][k] for k in ("u", "ci", "verdict")}},
                     "reading": r["targets"]["shift"]["reading"]["missing"]})
out["correction"] = corr

# ------------------------------------------------------------------ the ledger
out["ledger"] = {"events": len(ledger), "by_type": {}, "checks": loop.ledger_check(ledger),
                 "conflicts": len(ev.ledger_conflicts([e for e in ledger if "verdict" in e and "run" in e])),
                 "forecast_record": loop.forecast_record(tests)}
for e in ledger:
    k = e.get("type") or ("verdict" if "verdict" in e else "source")
    out["ledger"]["by_type"][k] = out["ledger"]["by_type"].get(k, 0) + 1
loop.write_ledger(ledger, RES / "ledger.jsonl")
(RES / "summary.json").write_text(json.dumps(out, indent=1, default=float))
print("G1", out["G1"]["holds"], "| G2", out["G2"]["holds"], "| G3", out["G3"]["holds"], "(by gap:", out["G3"]["holds_by_gap"], ")")
for ct, v in g3.items():
    print(f"  {ct}: grid reads '{v['reading_of_the_grid']}', test {v['outcome']}; forecast coverage {v['forecast']['coverage']['ci']} measured "
          f"{v['arms']['coverage']['gain']:.3f}; other study {v['forecast']['other study']['ci']} measured {v['arms']['other study']['gain']:.3f}")
print("ledger:", out["ledger"]["by_type"], "| checks:", [f"{c['code']} {c['severity']}" for c in out["ledger"]["checks"]])
