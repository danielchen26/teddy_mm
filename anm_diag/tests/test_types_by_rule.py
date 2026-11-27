"""A type is proposed and tested as an edge is (LOOP_PROTOCOL.md, the candidate for version 7): flagged by a law that was
told of no type, grouped by the direction of the shift, and admitted or refused by rows that did not propose it."""
import numpy as np

from anm_diag import account as AC
from anm_diag import events as EV


def _condition(seed, plan, p=16, K=8, n=200, nc=2500, plain=20, misread=1.0):
    """Cells with a stage of p dimensions and K readouts the stage partly tells. plan: {event: (direction index, amplitude)}.
    An event of the plan moves the stage along its direction and changes two readouts in proportion, by more than a
    readout of control cells says (`misread`); the direction and what it does are the same in every condition."""
    rng, fixed = np.random.default_rng(seed), np.random.default_rng(7)
    Bm = fixed.normal(0, 0.25, (p, K))
    dirs = np.linalg.qr(fixed.normal(0, 1, (p, 6)))[0].T
    beta = fixed.normal(0, 1, (6, K)) * misread
    names = ["control"] + sorted(plan) + [f"n{i:02d}" for i in range(plain)]
    xs, Ys, un = [], [], []
    for u in names:
        m = nc if u == "control" else n
        x = rng.normal(0, 1, (m, p))
        y = x @ Bm + rng.normal(0, 0.5, (m, K))
        if u in plan:
            k, a = plan[u]
            x = x + a * dirs[k]
            y = y + a * (dirs[k] @ Bm + beta[k])
        elif u != "control":
            y = y + rng.normal(0, 0.02, K)
        xs.append(x); Ys.append(y); un += [u] * m
    x, y, un = np.vstack(xs), np.vstack(Ys), np.array(un)
    return x, {f"r{j}": y[:, j] for j in range(K)}, un


PLAN = {"a1": (0, 1.2), "a2": (0, 0.9), "a3": (0, 1.5), "a4": (0, 1.1), "b1": (1, 1.3), "b2": (1, 1.0), "b3": (1, 1.4), "m1": (2, 1.2), "m2": (3, 1.2)}


def test_a_law_told_of_no_type_proposes_the_types_that_are_there():
    x, Y, un = _condition(0, PLAN)
    law = AC.Law(x, Y, un, types={}).fit()
    flags = law.flags()
    assert {f["event"] for f in flags if f["flag"] == "a direction the law does not cover"} >= set(PLAN)
    got = EV.propose_types(flags, law.shifts())
    assert sorted(t["members"] for t in got) == [["a1", "a2", "a3", "a4"], ["b1", "b2", "b3"]]       # m1 and m2 each move the stage their own way: no type
    assert all(t["status"] == "proposed" and t["evidence"]["cosine_of_the_shifts"]["smallest"] > 0.8 for t in got)
    assert EV.propose_types(flags, law.shifts(), known=["a1", "a2", "a3", "a4"])[0]["members"] == ["b1", "b2", "b3"]   # an event that has a type is not proposed again
    assert EV.propose_types(flags, law.shifts(), minimum=5) == []
    # one flagged event seeds a type; the events that move the stage the same way join it, flagged or not
    seed = [f for f in flags if f["event"] in ("a1", "m1")]
    assert EV.propose_types(seed, law.shifts()) == []
    grown = EV.propose_types(seed, law.shifts(), moving=law.moving())
    assert [t["members"] for t in grown] == [["a1", "a2", "a3", "a4"]] and set(grown[0]["evidence"]["seeded_by"]) == {"a1"}
    assert grown[0]["evidence"]["joined_by_the_direction_of_their_shift"] == ["a2", "a3", "a4"] and set(law.moving()) >= set(PLAN)


def test_rows_that_did_not_propose_a_type_admit_it_or_refuse_it():
    x, Y, un = _condition(0, PLAN)
    law = AC.Law(x, Y, un, types={}).fit()
    proposed = EV.propose_types(law.flags(), law.shifts())
    made_up = {"type_id": "m1+m2", "kind": "type", "members": ["m1", "m2"], "basis": "flag", "status": "proposed", "stratum": None, "evidence": None, "provenance": "made up for the test"}
    x2, Y2, un2 = _condition(1, PLAN)                                # another condition: other cells, the same events
    res = {t["type_id"]: t for t in AC.test_types(proposed + [made_up], x2, Y2, un2)}
    for k in ("a1+a2+a3+a4", "b1+b2+b3"):
        t = res[k]["test"]
        on = t["on_the_pairs_it_teaches"]
        assert res[k]["status"] == "admitted" and on["error_with_the_type"] < 0.5 * on["error_without_the_type"] and t["events_with_a_pair_it_teaches"] == res[k]["members"]
        assert t["pairs_it_teaches"] >= 0.5 * t["pairs"] and on["closer"] >= 0.9 * t["pairs_it_teaches"]
    assert res["m1+m2"]["status"] == "refused"                       # two events that move the stage in different directions teach each other nothing
    assert res["m1+m2"]["test"]["on_the_pairs_it_teaches"]["share_of_the_squared_error_removed"] < 0.1
    x3, Y3, un3 = _condition(2, {k: (d, 0.0) for k, (d, _a) in PLAN.items()})     # a condition in which the events do nothing
    assert all(t["status"] == "undecided" for t in AC.test_types(proposed, x3, Y3, un3))
    again = AC.test_types(list(res.values()), x2, Y2, un2)
    assert [t["status"] for t in again] == [t["status"] for t in res.values()]    # a type that was tested is not tested again
