"""Events without control cells (GENERAL.md, Section 6; runs/unit_events_20261123): a unit against the other units, and a
group of units against another. The planted systems are in anm_diag/toy.py."""
import numpy as np
import pytest

from anm_diag import cycle, outer, toy
from anm_diag import scope as SC

SIZE = dict(n_ordinary=40, n_cell=200)


@pytest.fixture(scope="module")
def units(tmp_path_factory):
    return toy.run_units(tmp_path_factory.mktemp("units"), **SIZE)


@pytest.fixture(scope="module")
def disease(tmp_path_factory):
    return toy.run_groups(tmp_path_factory.mktemp("disease"))


def items(r, event):
    return {i["reading"]: i for i in r["items"] if i["event"] == event}


def test_the_other_side_counts_by_unit_at_admission():
    d = toy.planted_units(**SIZE)
    sel = cycle.select(d["ids"], d["unit"], d["y"], seed=1, control=None)
    assert 0.30 < sel["share"].mean() < 0.37
    assert "weak" not in sel["admitted"] and not [u for u in sel["admitted"] if u.startswith("u")]       # a unit that barely differs, and the ordinary units, are refused
    assert {"early", "enc", "modb", "post"} <= set(sel["admitted"])
    m = d["unit"] != "h1"
    by_cell, by_unit = cycle.screen(d["y"][~m], d["y"][m]), cycle.screen(d["y"][~m], d["y"][m], units=d["unit"][m])
    assert by_unit["change_in_control_sd"] == by_cell["change_in_control_sd"] and by_unit["snr"] < 0.5 * by_cell["snr"]      # the same change, against a wider interval


def test_a_unit_event_is_read_against_the_other_units_and_pairs_with_the_tool(units):
    assert units["read_as"] == "unit by unit"
    for s, e in units["reading"].items():
        assert e["resampled_by"] == "cell" and set(e["arms"]) == {"controls", "rows", "cells", "unseen"}
        assert "reference" in e["arms"]["unseen"] and "reference" not in e["arms"]["rows"]            # `unseen` is read against a reader that never saw the unit


def test_every_loss_before_the_head_is_put_where_it_was_planted_for_a_unit(units):
    ev = units["placed"]["events"]
    for name, (a, b) in toy.LOST_UNITS.items():
        assert ev[name]["input"]["u"] < 0.2 and ev[name]["model_stage"]["verdict"] == "reject" and ev[name]["lost_between"] == [[a, b]]
        assert not units["reading"][name]["A3"]["holds"]                                              # the head has nothing to remove
    for name in ("modb", "post"):
        assert ev[name]["input"]["verdict"] == "reject" and ev[name]["model_stage"]["verdict"] == "reject"
    assert ev["modb"]["other_inputs"]["view_b"]["u"] < 0.4 and ev["post"]["other_inputs"]["view_b"]["u"] > 0.8


def test_a_loss_in_the_head_is_removed_with_or_without_the_units_own_row(units):
    for h in ("h1", "h2", "h3"):
        if h not in units["reading"]:
            continue
        e = units["reading"][h]
        assert e["as_trained"]["verdict"] == "reject" and e["A3"]["holds"]
        assert abs(e["arms"]["rows"]["u"]) < 0.15 and abs(e["arms"]["unseen"]["u"]) < 0.15 and e["rows_minus_unseen"]["reading"] == "no step"
        assert "pooled, and differs" not in items(units, h)


def test_a_unit_that_only_its_own_row_reaches_is_flagged(units):
    e = units["reading"]["solo"]
    assert e["A3"]["holds"] and e["arms"]["unseen"]["u"] > 0.8 and e["arms"]["rows"]["u"] < 0.6 and e["rows_minus_unseen"]["reading"] == "a step"
    it = items(units, "solo")["pooled, and differs"]
    assert it["lever"] == "S" and it["question"] == "the rows of the other units do not carry to this unit: solo"
    assert not [i for i in units["items"] if i["reading"] == "pooled, and differs" and i["event"] != "solo"]


def test_nothing_false_is_said_of_a_unit(units):
    definite = {(i["event"], i["reading"]) for i in units["items"] if i["reading"] in ("upstream", "unreadable", "pooled, and differs")}
    assert definite <= set(toy.TRUTH_UNITS.items())


def test_a_disease_is_a_pool_of_patient_donors_read_by_unit_on_both_sides(disease):
    e = disease["reading"]
    assert disease["read_as"] == "a pool" and e["resampled_by"] == "event" and len(e["events"]) == 8
    assert e["as_trained"]["verdict"] == "reject" and e["A3"]["holds"] and e["arms"]["rows"]["verdict"] == "support" and e["gate3"]["passes"]
    assert disease["status"] == "closed (absolute)"
    d = toy.planted_groups()
    tr = d["treated"]
    by_cell, _a = cycle.output_u(d["y"], d["o"], tr, d["unit"])
    by_unit, _a = cycle.output_u(d["y"], d["o"], tr, d["unit"], control_cluster=d["unit"])
    assert by_cell[0] == by_unit[0] and by_unit[1:].std() >= by_cell[1:].std()                         # the same value; the healthy side now counts by donor


def test_a_disease_with_nothing_wrong_leaves_nothing_open(tmp_path):
    r = toy.run_groups(tmp_path, sound=True)
    assert r["status"] == "closed (absolute)" and r["items"] == [] and r["target"]["repair"] == "none: the output as trained is supported"


def test_the_repair_of_units_is_a_plug_in_of_the_same_shape():
    d = toy.planted_units(n_ordinary=12, n_cell=120)
    units_ = sorted(set(d["unit"]))
    R = cycle.UnitRows(d["z"], d["y"], d["o"], np.array([], int), [(u, np.where(d["unit"] == u)[0]) for u in units_])
    assert {"controls", "rows", "cells"} <= set(R.th) and all(f"unseen|{u}" in R.th for u in units_)
    assert R.out("rows", d["z"]).shape == d["y"].shape and len(R.log["rows"]) == len(units_)
    k = units_.index("solo")
    assert R.log["rows"][k]["share_of_the_residual_closed"] > 0.5                                      # its own row moves the readout toward it
    m = d["unit"] == "solo"
    shift = lambda arm: R.out(arm, d["z"][m]).mean() - R.out(arm, d["z"][~m]).mean()
    assert abs(shift("unseen|solo")) < 0.3 * abs(shift("rows"))                                         # and a fit that never saw the unit predicts almost none of its shift


def test_the_grid_of_each_readout_is_the_run_of_that_readout_alone(tmp_path):
    """One run of the tool on a panel gives, for every readout, what a run on that readout alone gives: the same readers
    (the ridge weight per readout), the same resamples."""
    rng = np.random.default_rng(0)
    n_u, n_c = 6, 150
    unit = np.repeat([f"d{k}" for k in range(n_u)], n_c)
    f = rng.normal(size=(len(unit), 3)) + np.repeat(rng.normal(size=(n_u, 3)), n_c, axis=0)
    z = np.hstack([f, rng.normal(size=(len(unit), 5))]).astype(np.float32)
    Y = (f @ np.array([[1.0, 0.2, 0.0], [0.0, 1.0, 0.5], [0.3, 0.0, 1.0]]) + 0.3 * rng.normal(size=f.shape)).astype(np.float32)
    np.savez(tmp_path / "a.npz", unit=unit, z=z, Y=Y, o=(Y + rng.normal(size=Y.shape)).astype(np.float32))
    base = dict(name="t", kind="donor", study="s", events=["d0", "d3"], reads="z", label_basis="design", confounded_with=[], keep_draws=True, arrays={"audit": "a.npz"},
                decision={"margin_u": 0.10, "step_margin": 0.10, "alpha": 0.05, "bootstrap": 200, "seed": 0, "min_clusters": 10},
                readout_classes={"ridge": {"regressor": "RidgeCV", "alphas": [0.1, 1.0, 10.0, 100.0], "alpha_per_target": True}})
    keys = dict(readouts="Y", stages={"z": "z"}, unit="unit", frozen={"array": "o"})
    panel = SC.run_spec({**base, "each_readout": True, "keys": keys}, tmp_path)
    assert "each_readout" not in SC.run_spec({**base, "keys": keys}, tmp_path)["results"][0]["targets"]["shift"]          # off unless a SPEC asks
    for j in range(3):
        alone = SC.run_spec({**base, "keys": {**keys, "readout_columns": [j], "frozen": {"array": "o", "columns": [j]}}}, tmp_path)
        for rp, ra in zip(panel["results"], alone["results"]):
            got, want = rp["targets"]["shift"]["each_readout"][j]["stages"], ra["targets"]["shift"]["stages"]
            for st in want:
                for sc in want[st]:
                    assert abs(got[st][sc]["u"] - want[st][sc]["u"]) < 1e-4 and np.allclose(got[st][sc]["ci"], want[st][sc]["ci"], atol=1e-4) and got[st][sc]["verdict"] == want[st][sc]["verdict"]
                    k = f"|{st}|{sc}|{want[st][sc]['class']}"
                    assert np.allclose(rp["_draws"][f"shift@{j}" + k], ra["_draws"]["shift" + k], atol=1e-4)
    by, draws = cycle.load({"json": {"results": [{k: v for k, v in r.items() if k != "_draws"} for r in panel["results"]]},
                            "draws": type("Z", (), {"files": [f"{r['event']['unit']}|{k}" for r in panel["results"] for k in r["_draws"]],
                                                    "__getitem__": lambda self, key: {f"{r['event']['unit']}|{k}": v for r in panel["results"] for k, v in r["_draws"].items()}[key]})()}, readout=1)
    assert set(by) == {"d0", "d3"} and "d0|shift|z|own|ridge" in draws and not [k for k in draws if "@" in k]


def test_an_item_when_the_rows_of_the_other_units_do_not_carry():
    V = lambda u, lo, hi, verdict: {"u": u, "ci": [lo, hi], "verdict": verdict}
    A = lambda a, lo, hi, reading: {"added": a, "ci": [lo, hi], "reading": reading}
    ev = lambda own: {**V(0.05, -0.02, 0.09, "support"), "added": A(0.0, -0.03, 0.03, "no added loss"), "in_pool": False, "reference": V(0.04, -0.02, 0.09, "support"),
                      "unseen": {**V(0.9, 0.8, 1.0, "reject"), "added": A(0.0, -0.05, 0.05, "no added loss"), "own_row_adds": own}}
    t = lambda own: {"readout": "r", "group": "g", "type": {}, "events_n": 1, "admitted": True, "reference": None, "attainable": "absolute", "gates": {"evaluable": True, "why": []}, "arms": {},
                     "model": {"own": V(0.05, -0.02, 0.09, "support"), "added": A(0.0, -0.03, 0.03, "no added loss"), "events": {"d1": ev(own)}}}
    step, none = {"gain_difference": 0.85, "ci": [0.7, 0.95], "reading": "a step"}, {"gain_difference": 0.01, "ci": [-0.02, 0.04], "reading": "no step"}
    got = outer.open_items([t(step)], {"final_units": 1})
    assert [(i["reading"], i["lever"], i["event"]) for i in got] == [("pooled, and differs", "S", "d1")]
    assert outer.open_items([t(none)], {"final_units": 1}) == []
