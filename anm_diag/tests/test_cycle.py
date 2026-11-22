"""The general cycle on a planted system with two modalities, two steps inside the model and two types of cell
(anm_diag/toy.py): every loss is put where it was planted, a system with nothing wrong leaves nothing open, and no name
of the first system is needed."""
import re
from pathlib import Path

import pytest

from anm_diag import cycle, toy


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    return toy.run(tmp_path_factory.mktemp("planted"))


@pytest.fixture(scope="module")
def sound(tmp_path_factory):
    return toy.run(tmp_path_factory.mktemp("sound"), sound=True)


def items(r, group, event):
    return {i["reading"]: i for i in r["items"] if i["group"] == group and i["event"] == event}


def test_an_event_that_barely_moves_the_readout_is_refused_on_the_selection_share(run):
    assert sorted(run["admitted"]["one"]) == ["early", "enc", "modb", "p1", "p2", "p3", "post"]          # `weak` is not among them
    assert sorted(run["admitted"]["two"]) == ["odd", "q1", "q2", "q3"]


def test_a_loss_in_the_head_is_closed_by_a_repair_of_the_head(run):
    e = run["types"]["one"]["reading"]
    assert e["as_trained"]["verdict"] == "reject" and e["as_trained"]["u"] > 0.9 and e["A3"]["holds"]
    assert e["reference"]["verdict"] == "support" and e["arms"]["rows"]["verdict"] == "support" and e["gate3"]["passes"]
    assert not e["gate3"]["control_the_output_as_trained"]["passes"]
    assert e["closure"]["grade"] == "absolute" and run["types"]["one"]["status"] == "closed (absolute)"


def test_every_loss_before_the_head_is_put_where_it_was_planted(run):
    got = {(i["group"], i["event"]): i for i in run["items"] if i["reading"] in ("upstream", "unreadable")}
    assert {k: v["reading"] for k, v in got.items()} == toy.TRUTH
    assert got["one", "enc"]["lever"] == "M" and got["two", "odd"]["lever"] == "M" and got["one", "modb"]["lever"] == "V" and got["one", "post"]["lever"] == "V"
    for k, (a, b) in toy.LOST.items():                                    # and, along the stages inside the model, in which step
        assert got[k]["note"] == f"a step is lost between {a} and {b}"
    assert "view_b" in (got["one", "modb"].get("forecast") or "")         # another measured input reads it
    assert "view_b" not in (got["one", "post"].get("forecast") or "")     # no input reads an event that acts on the readout itself


def test_the_stages_read_what_was_planted(run):
    ev = run["types"]["one"]["placed"]["events"]
    assert ev["enc"]["input"]["verdict"] == "support" and ev["enc"]["model_stage"]["u"] > 0.9
    assert ev["modb"]["input"]["u"] > 0.9 and ev["modb"]["other_inputs"]["view_b"]["u"] < 0.2
    assert all(ev["post"][k]["u"] > 0.9 for k in ("input", "model_stage")) and ev["post"]["other_inputs"]["view_b"]["u"] > 0.9
    for T in run["types"].values():
        assert T["placed"]["positive_control"]["holds"] and T["placed"]["control_view"]["holds"]


def test_the_repair_of_one_type_does_not_carry_to_the_other(run):
    for T in run["types"].values():
        assert T["reading"]["arms"]["other"]["u"] > 0.9
    across = [i for i in run["items"] if i["question"].startswith("the whole repair across")]
    assert len(across) == 1 and across[0]["reading"] == "pooled, and differs" and across[0]["lever"] == "S"
    assert not [i for i in run["items"] if i["reading"] == "pooled, may be the same"]


def test_an_event_that_differs_inside_its_pool_is_flagged_and_the_pool_waits_for_the_split(run):
    odd = items(run, "two", "odd")
    assert set(odd) == {"upstream", "pooled, and differs"}
    e = run["types"]["two"]["reading"]
    assert not e["A3"]["holds"] and e["closure"]["grade"] == "not a target for a repair"          # with the odd event inside, the pool does not decide that the head adds a loss
    assert run["types"]["two"]["status"] == "undecided in the model" and run["types"]["two"]["target"]["repair"].startswith("none")
    waiting = [i for i in run["items"] if i["group"] == "two" and i["event"] == "pool" and i["reading"] == "unresolved"]
    assert waiting and all(odd["pooled, and differs"]["id"] in i["after"] for i in waiting)


def test_an_event_read_on_its_own_is_not_called_pooled(run):
    assert not [i for i in run["items"] if i["reading"] == "pooled, and differs" and i["event"] in toy.TYPES["one"]["alone"]]


def test_nothing_else_is_said(run):
    definite = {(i["group"], i["event"], i["reading"]) for i in run["items"] if i["reading"] != "unresolved"}
    assert definite == {("one", "early", "upstream"), ("one", "enc", "upstream"), ("one", "modb", "unreadable"), ("one", "post", "unreadable"), ("one", "pool", "pooled, and differs"),
                        ("two", "odd", "upstream"), ("two", "odd", "pooled, and differs")}


def test_a_system_with_nothing_wrong_leaves_nothing_open_and_no_repair_is_read(sound):
    assert sound["items"] == []
    for T in sound["types"].values():
        e = T["reading"]
        assert e["as_trained"]["verdict"] == "support" and not e["A3"]["holds"] and e["closure"]["grade"] == "nothing to close"
        assert T["status"] == "closed (absolute)" and T["target"]["repair"] == "none: the output as trained is supported" and T["target"]["arms"] == {}


def test_a_status_before_the_model_needs_a_rejected_reference_stage():
    """Where no repair is read, the target is described by the output as trained; "open before the model" is said only where
    a stage before the model is rejected; where it is merely undecided the status says so (found on planted systems of
    medium size, where a system with nothing wrong was called open)."""
    from anm_diag import outer
    V = lambda u, lo, hi, verdict: {"u": u, "ci": [lo, hi], "verdict": verdict}
    A = lambda a, lo, hi, reading: {"added": a, "ci": [lo, hi], "reading": reading}

    def status(as_trained, reference, a3):
        e = {"events": ["a", "b", "c"], "reference": reference, "as_trained": {**as_trained, "beyond_reference": a3}, "A3": {**a3, "holds": a3["reading"] == "loss"},
             "arms": {"rows": {**V(0.01, -0.03, 0.05, "support"), "beyond_reference": A(0.0, -0.02, 0.02, "no added loss")}}, "by_event": {}, "gate3": {"passes": True},
             "closure": {"attainable": "relative"}}
        t = cycle.target("r", "g", {}, e)
        return outer.target_status(t), t["repair"]

    none, unr, loss = A(0.0, -0.03, 0.03, "no added loss"), A(0.2, 0.02, 0.4, "unresolved"), A(0.8, 0.7, 0.9, "loss")
    assert status(V(0.9, 0.8, 1.0, "reject"), V(0.9, 0.8, 1.0, "reject"), none)[0] == "open before the model"
    assert status(V(0.08, 0.02, 0.14, "unresolved"), V(0.07, 0.01, 0.13, "unresolved"), none)[0] == "undecided before the model"
    assert status(V(0.9, 0.8, 1.0, "reject"), V(0.4, 0.0, 0.9, "unresolved"), unr)[0] == "undecided in the model"
    assert status(V(0.02, -0.03, 0.07, "support"), V(0.01, -0.04, 0.06, "support"), none) == ("closed (absolute)", "none: the output as trained is supported")
    assert status(V(0.9, 0.8, 1.0, "reject"), V(0.01, -0.04, 0.06, "support"), loss) == ("closed (absolute)", "read")


def test_the_selection_share_is_a_third_and_enters_no_reading():
    d = toy.planted()["one"]
    sel = cycle.select(d["ids"], d["unit"], d["y"], seed=1)
    assert 0.30 < sel["share"].mean() < 0.37
    assert cycle.select(d["ids"], d["unit"], d["y"], seed=1, own=("post",))["admitted"] == [e for e in sel["admitted"] if e != "post"]


def test_a_repair_is_a_plug_in(run):
    """Any object with `th` (a name for every arm, `controls` and `rows` among them) and `out(arm, x)` can take the place of
    the linear readout."""
    import numpy as np

    class LeastSquares:
        def __init__(self, x, y, o, C, events, other=None):
            fit = lambda idx: np.linalg.lstsq(np.c_[x[idx], np.ones(len(idx))], y[idx], rcond=None)[0]
            self.th = {"controls": fit(C), "rows": fit(np.r_[C, np.concatenate([T for _n, T in events])])}

        def out(self, arm, x):
            return np.c_[x, np.ones(len(x))] @ self.th[arm]

    d = toy.planted(n_control=2000, n_event=400)["one"]
    sites = list(toy.TYPES["one"]["pool"])
    out, inside, _f = cycle.cross_fit(d["z"], d["y"], d["o"], d["treated"], d["unit"], sites, cycle.halves(d["unit"]), repair=LeastSquares)
    u, _slope = cycle.output_u(d["y"][np.isin(d["unit"], sites + ["control"])], out["rows"][np.isin(d["unit"], sites + ["control"])], d["treated"][np.isin(d["unit"], sites + ["control"])])
    assert set(out) == {"controls", "rows"} and np.isfinite(inside).all() and abs(u[0]) < 0.1


def test_the_general_modules_name_nothing_of_the_first_system():
    here = Path(cycle.__file__).parent
    names = re.compile(r"teddy|\brna\b|protein|\bgenes?\b|knockout|pd-?l1|mhc|melanoma|thp-?1|cite|guide|ifn", re.I)
    for f in ("cycle.py", "outer.py", "toy.py"):
        hits = [(k + 1, line.strip()) for k, line in enumerate((here / f).read_text().splitlines()) if names.search(line)]
        assert not hits, (f, hits)
