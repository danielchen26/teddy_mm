"""Reading finished runs as a loop: admission, contrast of interventions, closure (LOOP_PROTOCOL.md).

A loop changes the audited system and audits it again. This module makes no fits. It reads the records of finished
runs (anm_record.json) and their bootstrap draws (anm_bootstrap.npz), and answers three questions, each with an interval:

  admit     Before any change: can this audit say anything about a change? The response is resolved, the reference stage
            is not rejected, the span to be repaired adds a resolved loss, and the interval is narrow enough for a
            repaired system to be recognized.
  diagnose  Is the located step the one that matters? Two changes are compared against the same baseline, one made at
            the located step and one made elsewhere: how much of the loss over a span each removes, and the difference.
            This is a comparison of interventions. It does not need the changed system to be free of loss.
  closure   Is the changed system free of the loss? Graded: "absolute" (the end stage is supported, two-sided, and the
            control's is rejected), "relative" (the span adds no loss, by its interval, and the control's does),
            "loss no longer resolved" (the span is unresolved: the loss is not shown and its absence is not shown
            either), "not closed".

Quantities. For a readout and a stage, u is the unexplained share of the response (stats.py). Over a span a -> b the
loss added is |u_b| - |u_a|, computed on every bootstrap draw, so its interval respects that both stages were measured
on the same cells. Runs of different systems on the same cells, folds and seed share their resamples draw for draw
(checked on the draws of the response); contrasts between such runs are paired. Otherwise the draws are combined as
independent and the result says so.

Stage verdicts here are two-sided: support needs the whole interval of some declared readout class inside
(-margin, margin). The one-sided verdicts of a record are left as they are.

  python -m anm_diag loop admit RUN --span Z:HEAD_OUT [--readout PDL1]
  python -m anm_diag loop diagnose --baseline RUN --located RUN --control RUN --span Z:HEAD_OUT [--readout PDL1]
  python -m anm_diag loop closure --arm RUN --control RUN --span Z:HEAD_OUT [--readout PDL1]
"""
import argparse
import json
from pathlib import Path

import numpy as np

from .stats import ci

DEFAULTS = {"min_clusters": 10, "min_delta_snr": 3.0}     # chosen numbers, reported with every admission


class LoopError(Exception):
    pass


class Run:
    """One finished audit: its record and, when present, its bootstrap draws."""

    def __init__(self, record, boots=None, label="run"):
        self.record, self.boots, self.label = record, boots, label
        d = record["declaration"]["decision"]
        self.margin = float(d["margin_u"])
        self.step_margin = float(d.get("step_margin", self.margin))
        self.alpha = float(d["alpha"])

    @classmethod
    def load(cls, path, label=None):
        p = Path(path)
        bf = p / "anm_bootstrap.npz"
        return cls(json.loads((p / "anm_record.json").read_text()), dict(np.load(bf)) if bf.exists() else None, label or p.name)

    def readouts(self):
        return [rid for rid, r in self.record["readouts"].items() if r["states"]]


def _iv(x, alpha):
    return [float(v) for v in ci(np.asarray(x, float), alpha)]


def stage(run, readout, sid, margin=None):
    """Two-sided verdict of one stage over the union of the declared readout classes, and the class it rests on."""
    m = run.margin if margin is None else float(margin)
    by = run.record["readouts"][readout]["states"][sid]["by_class"]
    sup = {c: max(abs(x["ci"][0]), abs(x["ci"][1])) for c, x in by.items() if -m < x["ci"][0] and x["ci"][1] < m}
    if sup:
        v, c = "support", min(sup, key=sup.get)
    else:
        c = min(by, key=lambda k: abs(by[k]["u"]))
        v = "reject" if all(x["ci"][0] > m or x["ci"][1] < -m for x in by.values()) else "unresolved"
    return {"stage": sid, "verdict": v, "class": c, "u": float(by[c]["u"]), "ci": [float(by[c]["ci"][0]), float(by[c]["ci"][1])]}


def _draws(run, readout, sid, cls):
    if run.boots is None:
        raise LoopError(f"{run.label}: no anm_bootstrap.npz, so no interval can be computed for a span or a contrast")
    return np.abs(run.boots[f"{readout}|{sid}/{cls}"].astype(float))


def added(run, readout, span):
    """Loss added over a span: the change in |u| from span[0] to span[1], with its interval and a three-way reading
    against the step margin ("loss": lower bound above it; "no added loss": upper bound below it; else "unresolved")."""
    a, b = span
    sa, sb = stage(run, readout, a), stage(run, readout, b)
    d = _draws(run, readout, b, sb["class"]) - _draws(run, readout, a, sa["class"])
    iv = _iv(d, run.alpha)
    reading = "loss" if iv[0] > run.step_margin else ("no added loss" if iv[1] < run.step_margin else "unresolved")
    return {"span": [a, b], "from": sa, "to": sb, "added": abs(sb["u"]) - abs(sa["u"]), "ci": iv,
            "half_width": (iv[1] - iv[0]) / 2, "reading": reading, "step_margin": run.step_margin, "draws": d}


def paired(run_a, run_b, readout):
    """True when two runs resampled the same cells in the same way: the draws of the response are identical."""
    if run_a.boots is None or run_b.boots is None:
        return False
    x, y = run_a.boots.get(f"{readout}|delta"), run_b.boots.get(f"{readout}|delta")
    return x is not None and y is not None and x.shape == y.shape and bool(np.array_equal(x, y))


def admit(run, readout, span, min_clusters=DEFAULTS["min_clusters"], min_delta_snr=DEFAULTS["min_delta_snr"]):
    """Gate 0. What this audit can and cannot decide about a change over the span, before any change is made."""
    r = run.record["readouts"][readout]
    lo, hi = r["delta_ci"]
    half = (hi - lo) / 2
    snr = abs(r["delta"]) / half if half > 0 else float("inf")
    out = {"run": run.label, "readout": readout, "span": list(span), "thresholds": {"min_clusters": min_clusters, "min_delta_snr": min_delta_snr,
                                                                                    "margin": run.margin, "step_margin": run.step_margin},
           "response": {"delta": float(r["delta"]), "ci": [float(lo), float(hi)], "resolved": not (lo <= 0 <= hi), "snr": float(snr), "ok": bool(not (lo <= 0 <= hi) and snr >= min_delta_snr)},
           "clusters": {"n": int(run.record["n_treated_clusters"]), "ok": bool(run.record["n_treated_clusters"] >= min_clusters)},
           "sealed": bool(not run.record.get("exploratory", True))}
    reasons = []
    if not out["response"]["resolved"] or not r["states"]:
        out.update({"diagnostic": "not ready", "certification": "not attainable", "reasons": ["the response is not resolved: there is nothing to follow"]})
        return out
    ref = stage(run, readout, span[0])
    end = stage(run, readout, span[1])
    base = added(run, readout, span)
    w = base["half_width"]
    end_half = (end["ci"][1] - end["ci"][0]) / 2
    out["reference"] = {**ref, "ok": ref["verdict"] != "reject"}
    out["baseline"] = {k: v for k, v in base.items() if k != "draws"}
    out["precision"] = {"span_half_width": w, "end_stage_half_width": end_half, "headroom": run.step_margin - w}
    if not out["response"]["ok"]:
        reasons.append(f"the response is small against its own interval (|delta| / half-width = {snr:.1f}, wanted {min_delta_snr:g}): u divides by it")
    if ref["verdict"] == "reject":
        reasons.append(f"the reference stage {span[0]} is rejected: the response is not in the span's input for the declared readout classes")
    if base["reading"] != "loss":
        reasons.append(f"the span adds no resolved loss as it stands ({base['reading']}): nothing is located here to repair")
    out["diagnostic"] = "ready" if (out["response"]["ok"] and ref["verdict"] != "reject" and base["reading"] == "loss") else "not ready"
    if out["diagnostic"] != "ready":
        out["certification"] = "not attainable"
    elif ref["verdict"] == "support" and end_half < run.margin:
        out["certification"] = "absolute"                 # a repaired end stage could be supported
    elif w <= run.step_margin / 2:
        out["certification"] = "relative"                 # a repaired span could be shown to add no loss
    elif w < run.step_margin:
        out["certification"] = "marginal"                 # only a span that adds almost nothing would be recognized
    else:
        out["certification"] = "not attainable"
    if out["diagnostic"] == "ready" and out["certification"] != "absolute":
        if ref["verdict"] != "support":
            reasons.append(f"the reference stage {span[0]} is {ref['verdict']} (two-sided), so no later stage can be called supported")
        elif end_half >= run.margin:
            reasons.append(f"the interval of {span[1]} is wider than the margin, so it could not be called supported")
        if out["certification"] in ("marginal", "not attainable"):
            reasons.append(f"the interval of the span's change has half-width {w:.3f} against a step margin of {run.step_margin:g}")
    if not out["clusters"]["ok"]:
        reasons.append(f"{out['clusters']['n']} treated clusters (wanted {min_clusters}): intervals are approximate")
    out["grade"] = "confirmatory possible" if (out["sealed"] and out["clusters"]["ok"]) else "exploratory"
    out["reasons"] = reasons
    return out


def contrast(base, arm, readout, span):
    """How much of the baseline's added loss an intervention removes over the span: absolute and as a share.
    Reading: "reduces the loss" (the interval of the removed part is above zero), "... by less than the margin" (above
    zero and wholly below the step margin), "adds to the loss" (wholly below zero: the arm adds more than the baseline),
    "does not reduce the loss" (not above zero, and below the margin), else "unresolved"."""
    is_paired = paired(base, arm, readout)
    b, c = added(base, readout, span), added(arm, readout, span)
    d0 = b["draws"] if is_paired else np.random.default_rng(0).permutation(b["draws"])
    n = min(len(d0), len(c["draws"]))
    rem = d0[:n] - c["draws"][:n]
    iv = _iv(rem, base.alpha)
    resolved, small = iv[0] > 0, iv[1] < base.step_margin          # a resolved reduction; a reduction below the step margin
    reading = (("reduces the loss by less than the margin" if small else "reduces the loss") if resolved
               else ("adds to the loss" if iv[1] < 0 else ("does not reduce the loss" if small else "unresolved")))
    out = {"baseline": base.label, "arm": arm.label, "readout": readout, "span": list(span), "paired": is_paired, "below_margin": bool(small),
           "added_baseline": {k: v for k, v in b.items() if k != "draws"}, "added_arm": {k: v for k, v in c.items() if k != "draws"},
           "removed": b["added"] - c["added"], "removed_ci": iv, "reading": reading, "share_removed": None, "draws": rem}
    if b["ci"][0] > 0:
        out["share_removed"] = {"point": 1 - c["added"] / b["added"], "ci": _iv(1 - c["draws"][:n] / d0[:n], base.alpha)}
    return out


def diagnose(base, located, control, readout, span):
    """The comparison of interventions. "located step confirmed": the baseline adds a resolved loss over the span, the
    change at the located step removes a resolved part of it that is not confined below the step margin, and it removes
    more than the change made elsewhere. "selective": whatever the control removes is below the step margin."""
    b = added(base, readout, span)
    out = {"readout": readout, "span": list(span), "baseline": {k: v for k, v in b.items() if k != "draws"}}
    if b["reading"] != "loss":
        return {**out, "outcome": "not testable", "why": f"the baseline adds no resolved loss over the span ({b['reading']})"}
    cl, cc = contrast(base, located, readout, span), contrast(base, control, readout, span)
    all_paired = cl["paired"] and cc["paired"]
    dc = cc["draws"] if all_paired else np.random.default_rng(1).permutation(cc["draws"])
    n = min(len(cl["draws"]), len(dc))
    diff = _iv(cl["draws"][:n] - dc[:n], base.alpha)
    if cl["reading"] == "reduces the loss" and diff[0] > 0:
        outcome = "located step confirmed"
    elif cl["below_margin"] or diff[1] < 0:
        outcome = "located step not confirmed"
    else:
        outcome = "unresolved"
    pub = lambda c: {k: v for k, v in c.items() if k != "draws"}
    return {**out, "located": pub(cl), "control": pub(cc), "located_minus_control": {"point": cl["removed"] - cc["removed"], "ci": diff},
            "paired": all_paired, "selective": cc["below_margin"], "outcome": outcome}


def closure(arm, control, readout, span):
    """Whether the changed system is free of the loss the control still shows, graded (see the module text)."""
    a_end, c_end = stage(arm, readout, span[1]), stage(control, readout, span[1])
    a, c = added(arm, readout, span), added(control, readout, span)
    if c["reading"] != "loss" and c_end["verdict"] != "reject":
        grade = "not testable"
    elif a_end["verdict"] == "support" and c_end["verdict"] == "reject":
        grade = "absolute"
    elif a["reading"] == "loss":
        grade = "not closed"
    elif a["reading"] == "no added loss" and c["reading"] == "loss":
        grade = "relative"
    elif a["reading"] == "unresolved" and c["reading"] == "loss":
        grade = "loss no longer resolved"
    else:
        grade = "not testable"
    return {"readout": readout, "span": list(span), "arm": arm.label, "control": control.label, "grade": grade,
            "arm_end": a_end, "control_end": c_end, "arm_span": {k: v for k, v in a.items() if k != "draws"},
            "control_span": {k: v for k, v in c.items() if k != "draws"}}


def guard(baseline_score, score, allowed_fall):
    """A guard on the system's original task: kept if the score fell by at most the allowed amount."""
    fall = float(baseline_score) - float(score)
    return {"baseline": float(baseline_score), "score": float(score), "fall": fall, "allowed_fall": float(allowed_fall), "kept": bool(fall <= allowed_fall)}


# ------------------------------------------------------------------ typed events of a loop (one line each in a ledger)
def intervention_event(system, parent, component, change, role, held_fixed=(), fitted_on=None, note=None):
    """A change to the system. role: "baseline", "located" (made at the located step) or "control" (made elsewhere,
    or the same procedure without the change)."""
    assert role in ("baseline", "located", "control"), role
    return {"type": "intervention", "system": system, "parent": parent, "component": component, "change": change, "role": role,
            "held_fixed": list(held_fixed), "fitted_on": fitted_on, "note": note}


def audit_event(run, system):
    r = run.record
    return {"type": "audit", "run": run.label, "system": system, "declaration_sha256": r.get("declaration_sha256"), "n_cells": r["n_cells"],
            "n_treated_clusters": r["n_treated_clusters"], "sealed": not r.get("exploratory", True), "anm_diag_version": r.get("anm_diag_version")}


def decision_event(kind, result, runs):
    """kind: "admission", "diagnosis" or "closure"; result: what admit(), diagnose() or closure() returned."""
    assert kind in ("admission", "diagnosis", "closure"), kind
    key = {"admission": "diagnostic", "diagnosis": "outcome", "closure": "grade"}[kind]
    return {"type": "decision", "kind": kind, "readout": result["readout"], "span": result["span"], "runs": list(runs), "outcome": result[key], "result": result}


def write_ledger(events, path):
    with open(path, "w") as f:
        for e in events:
            f.write(json.dumps(e, ensure_ascii=False, default=float) + "\n")
    return len(events)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="anm_diag loop", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("admit", "diagnose", "closure"):
        p = sub.add_parser(name)
        p.add_argument("--span", required=True, help="FROM:TO, two stages of the run")
        p.add_argument("--readout", help="default: every readout with a resolved response")
        p.add_argument("--out", help="write the result as JSON")
        if name == "admit":
            p.add_argument("run")
        elif name == "diagnose":
            p.add_argument("--baseline", required=True); p.add_argument("--located", required=True); p.add_argument("--control", required=True)
        else:
            p.add_argument("--arm", required=True); p.add_argument("--control", required=True)
    a = ap.parse_args(argv)
    span = tuple(a.span.split(":"))
    if len(span) != 2:
        raise SystemExit("--span takes FROM:TO")
    if a.cmd == "admit":
        run = Run.load(a.run)
        res = {rid: admit(run, rid, span) for rid in ([a.readout] if a.readout else list(run.record["readouts"]))}
        line = lambda r: f"diagnostic {r['diagnostic']}; certification {r['certification']}" + ("; " + "; ".join(r["reasons"]) if r.get("reasons") else "")
    elif a.cmd == "diagnose":
        b, l, c = Run.load(a.baseline), Run.load(a.located), Run.load(a.control)
        res = {rid: diagnose(b, l, c, rid, span) for rid in ([a.readout] if a.readout else b.readouts())}
        line = lambda r: r["outcome"] + (f": share of the loss removed {r['located']['share_removed']['point']:.2f} at the located step, "
                                         f"{r['control']['share_removed']['point']:.2f} by the control" if r.get("located") and r["located"]["share_removed"] else "")
    else:
        x, c = Run.load(a.arm), Run.load(a.control)
        res = {rid: closure(x, c, rid, span) for rid in ([a.readout] if a.readout else x.readouts())}
        line = lambda r: r["grade"]
    for rid, r in res.items():
        print(f"{rid}: {line(r)}")
    if a.out:
        Path(a.out).write_text(json.dumps(res, indent=1, default=float))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
