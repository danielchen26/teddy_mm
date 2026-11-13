"""The loop readings on hand-made runs: stage verdicts (two-sided), spans, admission, contrast of interventions, closure.

Runs are built from draws, with no fits. Draws of two runs made with the same seed share their resampling noise, as runs
of two systems on the same cells do."""
import json

import numpy as np
import pytest

from anm_diag import loop
from anm_diag.stats import ci

CHAIN = ["Z", "LAT", "OUT"]


def fake_run(u, label="run", seed=0, noise=0.02, own=0.004, delta=-0.9, delta_sd=0.03, clusters=12, sealed=False, B=600, u_hgb=None, boots=True):
    """u: {stage: value}. Resampling noise is drawn from `seed` and shared by every run with that seed."""
    shared = np.random.default_rng(seed)
    e = {s: shared.normal(0, noise, B) for s in CHAIN}
    d = delta + shared.normal(0, delta_sd, B)
    mine = np.random.default_rng(abs(hash(label)) % 2**32)
    bs, states = {"R|delta": d.astype(np.float32)}, {}
    for s in CHAIN:
        by = {}
        for cls, val in (("ridge", u[s]), ("hgb", (u_hgb or {}).get(s, u[s] + 0.03))):
            x = val + e[s] + mine.normal(0, own, B)
            bs[f"R|{s}/{cls}"] = x.astype(np.float32)
            by[cls] = {"u": float(val), "ci": ci(x, 0.05)}
        states[s] = {"role": "chain", "by_class": by}
    rec = {"declaration": {"decision": {"margin_u": 0.10, "alpha": 0.05}, "chain": CHAIN}, "n_cells": 3000, "n_treated": 1000,
           "n_treated_clusters": clusters, "exploratory": not sealed, "declaration_sha256": "x", "anm_diag_version": "test",
           "readouts": {"R": {"delta": delta, "delta_ci": ci(d, 0.05), "states": states}}}
    return loop.Run(rec, bs if boots else None, label)


BASE = {"Z": 0.03, "LAT": 0.22, "OUT": 0.45}


def test_stage_verdict_is_two_sided():
    r = fake_run({"Z": 0.0, "LAT": 0.3, "OUT": -0.3})
    assert loop.stage(r, "R", "Z")["verdict"] == "support"
    assert loop.stage(r, "R", "LAT")["verdict"] == "reject"
    assert loop.stage(r, "R", "OUT")["verdict"] == "reject"                      # far below zero is not support
    wide = fake_run({"Z": -0.6, "LAT": 0.0, "OUT": 0.0}, noise=0.3)
    assert wide.record["readouts"]["R"]["states"]["Z"]["by_class"]["ridge"]["ci"][1] < 0.10   # one-sided support
    assert loop.stage(wide, "R", "Z")["verdict"] == "unresolved"                 # not two-sided


def test_span_readings():
    assert loop.added(fake_run(BASE), "R", ("Z", "OUT"))["reading"] == "loss"
    assert loop.added(fake_run({"Z": 0.03, "LAT": 0.04, "OUT": 0.05}), "R", ("Z", "OUT"))["reading"] == "no added loss"
    assert loop.added(fake_run({"Z": 0.03, "LAT": 0.08, "OUT": 0.13}, noise=0.04), "R", ("Z", "OUT"))["reading"] == "unresolved"
    a = loop.added(fake_run(BASE), "R", ("Z", "OUT"))
    assert abs(a["added"] - 0.42) < 1e-9 and a["ci"][0] < 0.42 < a["ci"][1]
    # a negative u at the reference: the change is taken in |u|, so moving towards zero is not a loss
    assert loop.added(fake_run({"Z": -0.3, "LAT": -0.05, "OUT": -0.05}), "R", ("Z", "OUT"))["reading"] == "no added loss"


def test_admission():
    a = loop.admit(fake_run(BASE, sealed=True), "R", ("Z", "OUT"))
    assert a["diagnostic"] == "ready" and a["certification"] == "absolute" and a["grade"] == "confirmatory possible" and not a["reasons"]
    few = loop.admit(fake_run(BASE, clusters=5), "R", ("Z", "OUT"))
    assert few["diagnostic"] == "ready" and few["grade"] == "exploratory" and any("clusters" in x for x in few["reasons"])
    weak = loop.admit(fake_run(BASE, delta=-0.05, delta_sd=0.03), "R", ("Z", "OUT"))
    assert weak["diagnostic"] == "not ready" and weak["certification"] == "not attainable"
    small = loop.admit(fake_run(BASE, delta=-0.15, delta_sd=0.03), "R", ("Z", "OUT"))                 # resolved, but small against its interval
    assert small["response"]["resolved"] and not small["response"]["ok"] and small["diagnostic"] == "not ready"
    nothing = loop.admit(fake_run({"Z": 0.03, "LAT": 0.04, "OUT": 0.05}), "R", ("Z", "OUT"))
    assert nothing["diagnostic"] == "not ready" and any("nothing is located" in x for x in nothing["reasons"])
    # reference unresolved: a change can still be compared, and a span shown free of loss, but no stage called supported
    rel = loop.admit(fake_run({"Z": 0.08, "LAT": 0.25, "OUT": 0.45}, noise=0.012), "R", ("Z", "OUT"))
    assert rel["diagnostic"] == "ready" and rel["certification"] == "relative"
    wide = loop.admit(fake_run({"Z": 0.08, "LAT": 0.25, "OUT": 0.6}, noise=0.045), "R", ("Z", "OUT"))
    assert wide["diagnostic"] == "ready" and wide["certification"] in ("marginal", "not attainable")


def test_contrast_is_paired_on_shared_resamples():
    base, fixed = fake_run(BASE, "base"), fake_run({"Z": 0.03, "LAT": 0.06, "OUT": 0.08}, "fixed")
    other = fake_run({"Z": 0.03, "LAT": 0.06, "OUT": 0.08}, "other", seed=5)
    p, q = loop.contrast(base, fixed, "R", ("Z", "OUT")), loop.contrast(base, other, "R", ("Z", "OUT"))
    assert p["paired"] and not q["paired"]
    assert (p["removed_ci"][1] - p["removed_ci"][0]) < 0.5 * (q["removed_ci"][1] - q["removed_ci"][0])   # shared noise cancels
    assert p["reading"] == "reduces the loss" and abs(p["removed"] - 0.37) < 1e-9
    assert 0.8 < p["share_removed"]["point"] < 0.95 and p["share_removed"]["ci"][0] < p["share_removed"]["point"] < p["share_removed"]["ci"][1]
    same = loop.contrast(base, fake_run(BASE, "again"), "R", ("Z", "OUT"))
    assert same["reading"] == "does not reduce the loss" and same["below_margin"]
    little = loop.contrast(base, fake_run({"Z": 0.03, "LAT": 0.22, "OUT": 0.41}, "little"), "R", ("Z", "OUT"))
    assert little["reading"] == "reduces the loss by less than the margin"
    worse = loop.contrast(fake_run({"Z": 0.03, "LAT": 0.06, "OUT": 0.08}, "good"), fake_run({"Z": 0.03, "LAT": 0.08, "OUT": 0.13}, "worse"), "R", ("Z", "OUT"))
    assert worse["reading"] == "adds to the loss" and worse["removed_ci"][1] < 0


def test_diagnosis():
    base = fake_run(BASE, "base")
    mlp = fake_run({"Z": 0.03, "LAT": 0.06, "OUT": 0.08}, "located")
    dec = fake_run({"Z": 0.03, "LAT": 0.22, "OUT": 0.43}, "control")
    d = loop.diagnose(base, mlp, dec, "R", ("Z", "OUT"))
    assert d["outcome"] == "located step confirmed" and d["selective"] and d["paired"] and d["located_minus_control"]["ci"][0] > 0.2
    # the other way round: the change at the "located" step does nothing and the control removes the loss
    assert loop.diagnose(base, dec, mlp, "R", ("Z", "OUT"))["outcome"] == "located step not confirmed"
    # both changes remove the same part: nothing singles out the located step
    half = {"Z": 0.03, "LAT": 0.15, "OUT": 0.25}
    assert loop.diagnose(base, fake_run(half, "a"), fake_run(half, "b"), "R", ("Z", "OUT"))["outcome"] == "unresolved"
    healthy = fake_run({"Z": 0.03, "LAT": 0.04, "OUT": 0.05}, "healthy")
    assert loop.diagnose(healthy, mlp, dec, "R", ("Z", "OUT"))["outcome"] == "not testable"


def test_closure_grades():
    ctl = fake_run(BASE, "control")
    span = ("Z", "OUT")
    assert loop.closure(fake_run({"Z": 0.03, "LAT": 0.03, "OUT": 0.03}, "a"), ctl, "R", span)["grade"] == "absolute"
    assert loop.closure(fake_run({"Z": 0.03, "LAT": 0.2, "OUT": 0.4}, "b"), ctl, "R", span)["grade"] == "not closed"
    # reference near the margin: the end stage cannot be supported, the span can still be shown free of loss
    hi_ref, hi_ctl = {"Z": 0.09, "LAT": 0.10, "OUT": 0.11}, fake_run({"Z": 0.09, "LAT": 0.3, "OUT": 0.5}, "control2")
    assert loop.closure(fake_run(hi_ref, "c"), hi_ctl, "R", span)["grade"] == "relative"
    assert loop.closure(fake_run({"Z": 0.09, "LAT": 0.13, "OUT": 0.17}, "d", noise=0.03), hi_ctl, "R", span)["grade"] == "loss no longer resolved"
    healthy = fake_run({"Z": 0.03, "LAT": 0.03, "OUT": 0.04}, "healthy")
    assert loop.closure(healthy, healthy, "R", span)["grade"] == "not testable"


def test_needs_draws_and_guard():
    with pytest.raises(loop.LoopError):
        loop.added(fake_run(BASE, boots=False), "R", ("Z", "OUT"))
    assert loop.guard(0.602, 0.596, 0.02)["kept"] and not loop.guard(0.602, 0.56, 0.02)["kept"]


def test_ledger_and_command_line(tmp_path, capsys):
    runs = {"base": fake_run(BASE, "base"), "located": fake_run({"Z": 0.03, "LAT": 0.06, "OUT": 0.08}, "located"),
            "control": fake_run({"Z": 0.03, "LAT": 0.22, "OUT": 0.43}, "control")}
    for name, r in runs.items():
        (tmp_path / name).mkdir()
        (tmp_path / name / "anm_record.json").write_text(json.dumps(r.record, default=float))
        np.savez_compressed(tmp_path / name / "anm_bootstrap.npz", **r.boots)
    out = tmp_path / "d.json"
    assert loop.main(["diagnose", "--baseline", str(tmp_path / "base"), "--located", str(tmp_path / "located"), "--control", str(tmp_path / "control"),
                      "--span", "Z:OUT", "--out", str(out)]) == 0
    assert "located step confirmed" in capsys.readouterr().out and json.loads(out.read_text())["R"]["outcome"] == "located step confirmed"
    assert loop.main(["admit", str(tmp_path / "base"), "--span", "Z:OUT"]) == 0
    assert loop.main(["closure", "--arm", str(tmp_path / "located"), "--control", str(tmp_path / "base"), "--span", "Z:OUT"]) == 0
    d = loop.diagnose(runs["base"], runs["located"], runs["control"], "R", ("Z", "OUT"))
    ev = [loop.intervention_event("head-0", None, "head", "as trained", "baseline"),
          loop.intervention_event("head-1", "head-0", "head.mlp", "fine-tuned on the readout", "located", held_fixed=["head.decoder"]),
          loop.audit_event(runs["located"], "head-1"), loop.decision_event("diagnosis", d, ["base", "located", "control"])]
    assert loop.write_ledger(ev, tmp_path / "ledger.jsonl") == 4
    lines = [json.loads(x) for x in (tmp_path / "ledger.jsonl").read_text().splitlines()]
    assert [x["type"] for x in lines] == ["intervention", "intervention", "audit", "decision"] and lines[-1]["outcome"] == "located step confirmed"
