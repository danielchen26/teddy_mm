"""The outer step on planted targets: every reading of the routing table is produced by the case that should produce it,
nothing else is, and the ledger rules O1 to O6 catch what they are written for. No cell is read and nothing is fitted."""
import json
import sys
from pathlib import Path

import pytest

from anm_diag import outer

V = lambda u, lo, hi, verdict: {"u": u, "ci": [lo, hi], "verdict": verdict}
A = lambda a, lo, hi, reading: {"added": a, "ci": [lo, hi], "reading": reading}
SUP, REJ, UNR = V(0.02, -0.04, 0.08, "support"), V(0.40, 0.30, 0.50, "reject"), V(0.12, 0.01, 0.20, "unresolved")
NONE, LOSS = A(0.01, -0.03, 0.04, "no added loss"), A(0.30, 0.20, 0.40, "loss")


def target(**kw):
    """A target that is closed: admitted, the model adds nothing, the reference stage and the output are supported."""
    t = {"readout": "R", "group": "g", "type": {"cell_type": "c", "stratum": "s"}, "events_n": 12, "admitted": True, "response": {"snr": 9.0}, "attainable": "absolute",
         "reference": {"units": SUP, "own": SUP}, "model": {"own": SUP, "added": NONE, "reading": None, "rows_of_type": 5, "events": {}},
         "gates": {"evaluable": True, "why": [], "guard_kept": True}, "arms": {}, "pool": {"more_units_of_the_type": None}}
    t.update(kw)
    return t


def readings(targets, world=None):
    return [i["reading"] for i in outer.open_items(targets if isinstance(targets, list) else [targets], world)]


def test_a_closed_target_leaves_nothing():
    t = target()
    assert readings(t) == [] and outer.target_status(t) == "closed (absolute)"
    assert outer.coverage([t], [])["closed"]


CASES = {
    "invalid": target(reference={"units": REJ, "own": REJ}, input_view={"own": REJ, "gap_to_later": A(-0.09, -0.13, -0.04, "no added loss"), "classes": ["ridge", "hgb"], "nonlinear": True},
                      model={"own": REJ, "added": NONE, "events": {}}, attainable="relative"),
    "unresolved": target(reference={"units": UNR, "own": UNR}, model={"own": UNR, "added": NONE, "events": {}}, attainable="relative"),
    "small response": target(admitted=False, response={"snr": 2.0}),
    "unreadable": target(reference={"units": REJ, "own": REJ}, input_view={"own": REJ, "gap_to_later": A(0.0, -0.03, 0.03, "no added loss"), "classes": ["ridge", "hgb"], "nonlinear": True},
                         model={"own": REJ, "added": NONE, "events": {}}, attainable="relative"),
    "upstream": target(reference={"units": REJ, "own": REJ}, input_view={"own": SUP, "gap_to_later": A(0.38, 0.28, 0.48, "loss"), "classes": ["ridge", "hgb"], "nonlinear": True},
                       model={"own": REJ, "added": NONE, "events": {}}, attainable="relative"),
    "no rows": target(model={"own": REJ, "added": LOSS, "reading": "coverage", "rows_of_type": 0, "events": {}}),
    "pooled, and differs": target(model={"own": SUP, "added": NONE, "events": {"X": {**V(0.50, 0.40, 0.60, "reject"), "added": NONE}}}),
    "pooled, may be the same": target(arms={"rows_other": {**NONE, "pair": "a and b"}}),
    "mistyped": (target(), {"typed_right": 0.95, "typed_group": "g"}),
    "self event": target(self_events=[{"site": "OWN", "u": 0.9, "ci": [0.85, 0.95]}]),
    "rule": target(gates={"evaluable": False, "why": [{"question": "a gate that cannot be evaluated"}], "guard_kept": True}),
    "selection": (target(), {"selection_on_audited_cells": "screened on audited cells"}),
    "no final data": (target(), {"final_units": 0}),
    "confounded": target(confounded=[("study", "disease")]),
    "guard conflict": target(model={"own": REJ, "added": LOSS, "reading": "function", "rows_of_type": 5, "events": {}}, gates={"evaluable": True, "why": [], "guard_kept": False, "guard_holds_under_some_change": False}),
}


def test_the_cases_cover_the_table():
    assert set(CASES) == set(outer.ROUTES)


@pytest.mark.parametrize("reading", list(outer.ROUTES))
def test_each_reading_comes_from_its_case_and_only_from_it(reading):
    case = CASES[reading]
    t, world = case if isinstance(case, tuple) else (case, None)
    got = readings(t, world)
    assert got == [reading], got
    it = outer.open_items([t], world)[0]
    assert it["lever"] == outer.ROUTES[reading][0] and it["actions"] and it["state"] == "open" and it["id"] == "O01"
    assert not outer.outer_check([t], [it])


def test_a_loss_the_model_still_adds_is_the_inner_loops():
    t = target(model={"own": REJ, "added": LOSS, "reading": "function", "rows_of_type": 5, "events": {}})
    assert outer.target_status(t) == "open in the model" and readings(t) == []


def test_statuses():
    assert outer.target_status(target(admitted=False, response={"snr": 1.0})) == "not a target"
    assert outer.target_status(CASES["rule"]) == "not graded"
    assert outer.target_status(CASES["unresolved"]) == "the model's part closed (relative)"
    assert outer.target_status(target(model={"own": UNR, "added": A(0.08, 0.0, 0.16, "unresolved"), "events": {}})) == "undecided in the model"


def test_an_unreadable_verdict_on_one_linear_class_says_so():
    t = dict(CASES["unreadable"], input_view={"own": REJ, "gap_to_later": A(0.0, -0.03, 0.03, "no added loss"), "classes": ["ridge"], "nonlinear": False})
    it = outer.open_items([t])[0]
    assert it["reading"] == "unreadable" and "richer reader" in it["note"]


def test_a_dominated_view_gives_no_unreadable_or_upstream():
    assert readings(CASES["invalid"]) == ["invalid"]


def test_without_a_positive_control_no_unreadable_reading_is_drawn():
    t = dict(CASES["unreadable"], input_view=dict(CASES["unreadable"]["input_view"], positive_control=False))
    items = outer.open_items([t])
    assert [i["reading"] for i in items] == ["invalid"] and "positive control" in items[0]["question"]
    assert items[0]["next_step_needs"].startswith("no new data first")
    assert readings(dict(t, input_view=dict(t["input_view"], positive_control=True))) == ["unreadable"]


def test_a_loss_before_the_model_without_an_input_view_cannot_be_placed():
    t = target(reference={"units": REJ, "own": REJ}, model={"own": REJ, "added": NONE, "events": {}}, attainable="relative")
    assert readings(t) == ["invalid"]


def test_units_needed():
    assert outer.needed_units(0.043, [-0.021, 0.105], 5) == 6               # just over the line: one more unit, if the estimate holds
    assert outer.needed_units(0.117, [0.007, 0.193], 5) > 100               # next to the bound: not within reach
    assert outer.needed_units(0.02, [-0.04, 0.08], 5) == 5                  # already decided
    assert outer.needed_units(0.10, [0.0, 0.2], 5) is None                  # on the bound
    assert outer.needed_units(0.40, [0.30, 0.50], 5) == 5


def test_at_the_margin_and_what_the_pool_holds():
    near = target(reference={"units": V(0.117, 0.007, 0.193, "unresolved"), "own": UNR}, model={"own": UNR, "added": NONE, "events": {}}, attainable="relative", events_n=5)
    assert outer.open_items([near])[0]["next_step_needs"] == "at the margin"
    just = target(model={"own": V(0.043, -0.021, 0.105, "unresolved"), "added": NONE, "events": {}}, events_n=5, pool={"more_units_of_the_type": 0})
    it = outer.open_items([just])[0]
    assert it["units_needed"] == 6 and it["next_step_needs"] == "at the margin" and "no further unit" in it["forecast"]
    some = dict(just, pool={"more_units_of_the_type": 3})
    assert outer.open_items([some])[0]["next_step_needs"] == "units of the pool"


def test_a_split_inside_the_pool_comes_before_more_units():
    t = target(reference={"units": UNR, "own": UNR}, attainable="relative", model={"own": V(0.16, 0.08, 0.29, "unresolved"), "added": NONE, "events": {"X": {**V(0.54, 0.43, 0.65, "reject"), "added": NONE}}})
    items = outer.open_items([t])
    unresolved = [i for i in items if i["reading"] == "unresolved"][0]
    split = [i for i in items if i["reading"] == "pooled, and differs"][0]
    assert unresolved["after"] == [split["id"]]


def test_an_estimate_read_on_guessed_types_waits_for_the_bridge_reading():
    t = target(model={"own": V(0.043, -0.021, 0.105, "unresolved"), "added": NONE, "events": {}}, events_n=5, pool={"more_units_of_the_type": 0})
    items = outer.open_items([t], {"typed_right": 0.956, "typed_group": "g"})
    q = [i for i in items if i["reading"] == "unresolved"][0]
    m = [i for i in items if i["reading"] == "mistyped"][0]
    assert q["after"] == [m["id"]] and q["next_step_needs"] == "no new data" and q["forecast"].startswith("first the bridge reading")


def test_the_pool_at_hand_can_supply_events_for_a_small_response():
    has = target(admitted=False, response={"snr": 3.1}, pool={"more_units_of_the_type": 2, "admissible_now": ["B", "T"]})
    none = target(admitted=False, response={"snr": 2.2}, pool={"more_units_of_the_type": 0, "admissible_now": []})
    assert outer.open_items([has])[0]["next_step_needs"] == "units at hand"
    assert outer.open_items([none])[0]["next_step_needs"] == "the owner's decision"


def test_observations_on_one_question_are_one_item():
    a = target(group="a", arms={"other": {**LOSS, "pair": "a and b"}})
    b = target(group="b", arms={"other": {**A(0.25, 0.2, 0.3, "loss"), "pair": "a and b"}})
    items = outer.open_items([a, b])
    assert len(items) == 1 and len(items[0]["observations"]) == 2 and items[0]["standing"]


# --- the ledger rules ---------------------------------------------------------------------------------------------------
def _step():
    t = CASES["unresolved"]
    items = outer.open_items([t])
    h = outer.hypothesis("H1", items[0]["id"], "resolution", {"u": [0.05, 0.15]}, "the interval leaves the margin", "a split of the pool", "the new units read inside the old interval",
                         proposal_units=["e1", "e2"], test_units=["e3", "e4"])
    diff = [{"lever": "E", "change": "two more events of the type", "items": [items[0]["id"]], "hypothesis": "H1"}]
    return t, items, h, diff


def rules(findings, level="error"):
    return sorted({f["rule"] for f in findings if f["level"] == level})


def test_a_clean_step_passes():
    t, items, h, diff = _step()
    assert outer.outer_check([t], items, diff, [h], roles={"final": ["f1"]}, reads={"e1": 2, "f1": 0}) == []


def test_o1_one_status_per_target():
    t, items, h, diff = _step()
    assert rules(outer.outer_check([t, dict(t)], items)) == ["O1"]
    assert rules(outer.outer_check([{"readout": "R", "group": "x"}], [])) == ["O1"]


def test_o2_items_are_routed_and_limits_need_a_usable_reading():
    t, items, h, diff = _step()
    bad = [dict(items[0], state="limit", evidence="x")]
    assert rules(outer.outer_check([t], bad)) == ["O2"]                      # an unresolved item cannot become a limit
    assert rules(outer.outer_check([t], [dict(items[0], reading="hunch")])) == ["O2"]
    assert rules(outer.outer_check([t], [dict(items[0], state="at the margin", evidence=None)])) == ["O2"]


def test_o3_every_change_names_an_open_item():
    t, items, h, diff = _step()
    assert rules(outer.outer_check([t], items, [dict(diff[0], items=[])], [h])) == ["O3"]
    assert rules(outer.outer_check([t], items, [dict(diff[0], items=["O99"])], [h])) == ["O3"]


def test_o4_a_hypothesis_first_and_tested_on_other_units():
    t, items, h, diff = _step()
    assert rules(outer.outer_check([t], items, [dict(diff[0], hypothesis=None)], [h])) == ["O4"]
    assert rules(outer.outer_check([t], items, diff, [dict(h, control=None)])) == ["O4"]
    assert rules(outer.outer_check([t], items, diff, [dict(h, test_units=["e1", "e9"])])) == ["O4"]
    with pytest.raises(ValueError):
        outer.hypothesis("H", "O01", "hunch", 1, 1, 1, 1)


def test_o5_final_units_are_not_touched():
    t, items, h, diff = _step()
    assert rules(outer.outer_check([t], items, diff, [h], roles={"final": ["f1"]}, reads={"f1": 1})) == ["O5"]
    assert rules(outer.outer_check([t], items, diff, [dict(h, test_units=["f1"])], roles={"final": ["f1"]})) == ["O5"]
    assert rules(outer.outer_check([t], items, diff, [h], roles={"final": []}), "warning") == ["O5"]


def test_o6_a_change_of_the_claim_is_a_limit_and_recounts():
    t, items, h, diff = _step()
    ok = {"change": "drop the readout", "limit": {"item": "O01", "reading_passed_checks": True, "reopens": "a richer reader"}}
    assert outer.outer_check([t], items, claim_changes=[ok], cover=outer.coverage([t], items)) == []
    assert rules(outer.outer_check([t], items, claim_changes=[ok])) == ["O6"]
    assert rules(outer.outer_check([t], items, claim_changes=[{"change": "drop", "limit": None}], cover={})) == ["O6"]
    assert rules(outer.outer_check([t], items, claim_changes=[dict(ok, limit=dict(ok["limit"], reading_passed_checks=False))], cover={})) == ["O6"]


# --- the finished round of 20 November --------------------------------------------------------------------------------
ROUND = Path(__file__).resolve().parents[1] / "runs" / "teg_loop_20261120" / "results"
DESIGN = Path(__file__).resolve().parents[1] / "runs" / "outer_loop_design_20261120"


@pytest.mark.skipif(not (ROUND / "round_development.json").exists(), reason="the round's results are not in this checkout")
def test_the_round_of_20_november():
    sys.path.insert(0, str(DESIGN / "scripts"))
    import dryrun
    screen = json.loads((DESIGN / "results" / "pool_screen.json").read_text())
    targets, world = dryrun.targets_of_round(ROUND, screen)
    items = outer.open_items(targets, world)
    c = outer.counts(items)
    assert len(targets) == 7 and len(items) == 21
    assert c["by_lever"] == {"M": 0, "E": 8, "V": 3, "S": 5, "C": 1, "P": 4}
    assert c["by_reading"]["unreadable"] == 0 and c["by_reading"]["invalid"] == 3     # no positive control in co-culture: no reading of where the loss is
    assert [outer.target_status(t) for t in targets].count("the model's part closed (relative)") == 3
    assert rules(outer.outer_check(targets, items, roles={"final": []})) == []
    assert not outer.coverage(targets, items)["closed"]
    saved = json.loads((DESIGN / "results" / "open_items.json").read_text())
    assert [(i["id"], i["reading"], i["lever"]) for i in saved["items"]] == [(i["id"], i["reading"], i["lever"]) for i in items]
