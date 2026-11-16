"""Contrast events, scoped verdict events, hypotheses and tests, and one file format for two kinds of event."""
import json

import numpy as np
import yaml

from anm_diag import events as ev
from anm_diag import loop
from anm_diag import scope as sc
from anm_diag.__main__ import main as cli


# ------------------------------------------------------------------ contrasts
def test_contrast_events_have_the_fields_of_a_source_event():
    spec = dict(name="d", kind="donor", treated=["P1"], control="others", study="hao", label_basis="design", confounded_with=["run"])
    es = ev.contrast_events(spec)
    assert {e["role"] for e in es} == {"treated", "control"}
    assert all(f in es[0] for f in ev.VERIFIER_FIELDS) and es[0]["channel"] == "donor" and es[0]["pool"] == "d:A=1"
    assert ev.check_contrast(spec) == []


def test_contrast_rules():
    codes = lambda **k: {(c["code"], c["severity"]) for c in ev.check_contrast(dict(name="d", kind="donor", treated=["P1"], study="hao", **k))}
    assert ("C2", "warning") in codes() and ("C3", "warning") in codes()
    assert ("C2", "error") in codes(label_basis="response", confounded_with=[])
    assert ("C4", "error") in codes(label_basis="design", confounded_with=[], control_study="bm")
    assert ("C4", "error") not in {(c["code"], c["severity"]) for c in
                                   ev.check_contrast(dict(name="s", kind="study", treated=["hao"], study="hao", control_study="bm",
                                                          label_basis="design", confounded_with=["tissue"]))}
    assert ("C1", "error") in {(c["code"], c["severity"]) for c in ev.check_contrast(dict(name="x", kind="colour", treated=["a"]))}
    pooled = ev.check_contrast(dict(name="d", kind="donor", treated=["P1", "P2", "P3"], study="hao", label_basis="design", confounded_with=[]))
    assert {c["code"] for c in pooled} == {"C5", "C6"}


def test_verdicts_of_different_scopes_do_not_conflict():
    mk = lambda run, scope, sign: dict(run=run, readout="PDL1", site="Z", scope=scope, target="shift", sign=sign,
                                       verdict={1: "support", -1: "reject"}[sign], amplitude=0.0, source_site="g", context="s")
    assert ev.ledger_conflicts([mk("a", "units", 1), mk("b", "train", -1)]) == []
    assert len(ev.ledger_conflicts([mk("a", "units", 1), mk("b", "units", -1)])) == 1
    old = lambda run, sign: dict(run=run, readout="PDL1", site="Z", sign=sign, verdict={1: "support", -1: "reject"}[sign], amplitude=0.0)
    assert len(ev.ledger_conflicts([old("a", 1), mk("b", "units", -1)])) == 1        # a v1 verdict is a units verdict


# ------------------------------------------------------------------ hypotheses and tests
def hyp(**k):
    base = dict(hid="h1", about={"event": "donor:P1", "readout": "pool", "target": "variation", "quantity": "score"}, claim="coverage",
                forecast={"coverage": [0.07, 0.10], "function": [-0.02, 0.03], "amount": [-0.02, 0.02]},
                rule={"own": "coverage", "at_least": 0.04, "controls": {"function": 0.04, "amount": 0.04}}, basis=["gate0"], rivals=["function"])
    return loop.hypothesis_event(**{**base, **k})


def arm(gain, half=0.01, runs=("r1",), admitted=True):
    return {"gain": gain, "ci": [gain - half, gain + half], "runs": list(runs), "admitted": admitted}


def test_confirmed_refuted_unresolved_not_testable():
    h = hyp()
    assert loop.test(h, {"coverage": arm(0.08), "function": arm(0.01), "amount": arm(0.0)})["outcome"] == "confirmed"
    assert loop.test(h, {"coverage": arm(0.01), "function": arm(0.01)})["outcome"] == "refuted"            # the own arm did not gain
    assert loop.test(h, {"coverage": arm(0.08), "function": arm(0.09)})["outcome"] == "refuted"            # the rival gained as much
    assert loop.test(h, {"coverage": arm(0.045, half=0.02), "amount": arm(0.0)})["outcome"] == "unresolved"
    assert loop.test(h, {"coverage": arm(0.08)})["outcome"] == "not testable"                               # no control arm
    assert loop.test(h, {"coverage": arm(0.08, admitted=False), "amount": arm(0.0)})["outcome"] == "not testable"


def test_forecast_is_scored_on_every_arm():
    t = loop.test(hyp(), {"coverage": arm(0.12), "function": arm(0.01), "amount": arm(0.0)})
    assert t["outcome"] == "confirmed" and t["arms"]["coverage"]["inside_forecast"] is False
    assert abs(t["arms"]["coverage"]["forecast_error"] - 0.035) < 1e-9
    rec = loop.forecast_record([t])
    assert rec["n"] == 3 and rec["inside"] == 2


def test_ledger_rules():
    h = hyp()
    ok = loop.test(h, {"coverage": arm(0.08), "amount": arm(0.0)})
    assert [c for c in loop.ledger_check([h, ok]) if c["severity"] == "error"] == []
    codes = lambda evs: {(c["code"], c["severity"]) for c in loop.ledger_check(evs)}
    assert ("H1", "error") in codes([ok])                                                                    # test without its hypothesis
    assert ("H2", "error") in codes([h, loop.test(h, {"coverage": arm(0.08, runs=("gate0",)), "amount": arm(0.0)})])
    assert ("H2", "warning") in codes([hyp(retrospective=True), loop.test(hyp(retrospective=True), {"coverage": arm(0.08), "amount": arm(0.0)})])
    assert ("H3", "error") in codes([h, loop.test(h, {"coverage": arm(0.08)})])
    assert ("H4", "warning") in codes([h, loop.test(h, {"coverage": arm(0.08), "amount": arm(0.0), "extra": arm(0.0)})])
    assert ("H5", "note") in codes([h])


# ------------------------------------------------------------------ one file format, two kinds of event
def _arrays(tmp_path, seed=0):
    """A study in which the state carries a shift that the training cells never showed (the coverage case), written
    twice: as a perturbation audit (treated units against unassigned control cells) and as donors of one study."""
    rng = np.random.default_rng(seed)
    n_u, per = 6, 300
    unit = np.repeat([f"u{i}" for i in range(n_u)], per)
    z = rng.normal(size=(n_u * per, 6))
    shift = np.repeat(rng.normal(scale=1.5, size=n_u), per)
    z[:, 1] += shift
    y = z[:, 2] + z[:, 1] + 0.1 * rng.normal(size=len(z))
    n_c = 1500
    zc = rng.normal(size=(n_c, 6))
    yc = zc[:, 2] + zc[:, 1] + 0.1 * rng.normal(size=n_c)
    zt = rng.normal(size=(2000, 6))
    yt = zt[:, 2] + 0.1 * rng.normal(size=2000)                         # the training cells: no dependence on z1
    np.savez(tmp_path / "train.npz", Z=zt, y=yt)
    np.savez(tmp_path / "donors.npz", Z=z, y=y, unit=unit, frozen=z[:, 2])
    np.savez(tmp_path / "perturb.npz", Z=np.vstack([z, zc]), y=np.concatenate([y, yc]), unit=np.concatenate([unit, ["control"] * n_c]),
             A=np.concatenate([np.ones(len(z), int), np.zeros(n_c, int)]), frozen=np.concatenate([z[:, 2], zc[:, 2]]))
    common = dict(study="synthetic", arrays=dict(train="train.npz"), reads="Z", decision=dict(bootstrap=200),
                  keys=dict(readouts=["y"], stages={"Z": "Z"}, unit="unit", frozen={"array": "frozen"}))
    donors = {**common, "name": "donors", "kind": "donor", "label_basis": "design", "confounded_with": ["run"],
              "arrays": {**common["arrays"], "audit": "donors.npz"}}
    perturb = {**common, "name": "perturb", "kind": "perturbation", "units_exclude": ["control"],
               "arrays": {**common["arrays"], "audit": "perturb.npz"}, "keys": {**common["keys"], "control": {"key": "A", "equals": 0}}}
    for name, s in (("donors", donors), ("perturb", perturb)):
        (tmp_path / f"{name}.yaml").write_text(yaml.safe_dump(s))
    return donors, perturb


def test_one_format_reads_a_perturbation_and_a_donor_contrast(tmp_path):
    _arrays(tmp_path)
    for name in ("donors", "perturb"):
        assert cli(["scope", str(tmp_path / f"{name}.yaml"), "--out", str(tmp_path / name)]) == 0
        res = json.loads((tmp_path / name / "scope.json").read_text())
        strong = [r for r in res["results"] if r["delta"]["resolved"] and abs(r["delta"]["value"]) > 0.5]
        assert len(strong) >= 3
        assert all(r["targets"]["shift"]["reading"]["missing"] == "coverage" for r in strong), name
        lines = [json.loads(l) for l in (tmp_path / name / "scope_events.jsonl").read_text().splitlines()]
        assert any(e.get("scope") == "units" and e.get("target") == "shift" for e in lines)
    donor_events = [json.loads(l) for l in (tmp_path / "donors" / "scope_events.jsonl").read_text().splitlines()]
    assert any(e.get("channel") == "donor" and e.get("role") == "treated" for e in donor_events)        # the contrast is in the ledger


def test_loop_test_command(tmp_path, capsys):
    (tmp_path / "h.json").write_text(json.dumps(hyp()))
    (tmp_path / "arms.json").write_text(json.dumps({"coverage": arm(0.08), "amount": arm(0.0)}))
    assert cli(["loop", "test", "--hypothesis", str(tmp_path / "h.json"), "--arms", str(tmp_path / "arms.json"), "--out", str(tmp_path / "t.json")]) == 0
    assert "confirmed" in capsys.readouterr().out and json.loads((tmp_path / "t.json").read_text())["outcome"] == "confirmed"


def test_units_set_aside_as_final_are_not_read(tmp_path):
    import pytest
    donors, _ = _arrays(tmp_path)
    locked = {**donors, "roles": {"final": ["u5"]}}
    with pytest.raises(sc.ScopeError, match="set aside as final"):
        sc.run_spec(locked, tmp_path)
    opened = {**locked, "final_audit": True, "hypotheses": ["h1"], "events": ["u5"]}
    assert len(sc.run_spec(opened, tmp_path)["results"]) == 1
