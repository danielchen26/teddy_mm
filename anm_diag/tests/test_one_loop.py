"""One loop for every shape of event (LOOP_PROTOCOL.md, version 4, rules 4i to 4l): the account where no cell is a control
cell, edges as records of the event graph, a law whose records teach by type and whose edges have a size, and the step
after a rejection."""
import numpy as np

from anm_diag import account as AC
from anm_diag import cycle, events as EV, outer, toy
from anm_diag import toy_account as T


def _units(seed=3):
    d = toy.planted_units(seed, n_ordinary=24, n_cell=200)
    un = d["unit"]
    return d, un, cycle.halves(un, 0)


def _acc(d, un, first, e, **kw):
    return AC.account({"view_a": d["view_a"], "z": d["z"]}, d["y"], un, first, e, others=[u for u in sorted(set(un.tolist())) if u != e], inputs=["view_a"], state="z", o=d["o"],
                      control=None, cluster=un, draws=80, **kw)


def test_a_unit_against_the_other_units_is_placed_where_it_was_planted():
    d, un, first = _units()
    part = lambda a, k: a["chain"]["parts"][k]["share"]
    h1 = _acc(d, un, first, "h1")
    assert h1["stages"]["z"]["read_against"] == "the other units" and h1["stages"]["z"]["cells"]["units_of_the_reference"] == 32
    assert part(h1, "read") > 0.85 and not h1["chain"]["open"]                        # a factor that varies inside units: a readout with no record reads it
    assert h1["stages"]["z"]["shares"]["as_trained"]["share"] < 0.2                  # and the head as trained does not
    for e in ("early", "enc"):                                                       # dropped by the model: in the input, not in the state
        a = _acc(d, un, first, e, controls=False)
        assert [p["part"] for p in a["chain"]["open"]] == ["lost_before_the_state"] and part(a, "lost_before_the_state") > 0.8
        assert a["plan"][0]["lever"].startswith("M: the model before the state")
    post = _acc(d, un, first, "post")                                                # the readout itself: no stage holds it, no unit forecasts it
    assert [p["part"] for p in post["chain"]["open"]] == ["not_in_the_input"] and post["chain"]["across_events"]["reading"] == "the other events do not forecast it"
    assert post["plan"][0]["lever"] == "G: the graph" and post["plan"][0]["reading"] == "not_in_the_input, not forecast"
    assert post["must_fail"]["made_up_event"]["no_change_to_account_for"] and abs(post["must_fail"]["permuted_stage"]["carried"]["share"]) < 0.2
    solo = _acc(d, un, first, "solo", controls=False)                                # told apart completely by a marker: held as its identity
    assert solo["stages"]["z"]["carried_by"] == "identity" and solo["plan"][0]["reading"] == "held_as_identity"
    assert abs(solo["stages"]["z"]["shares"]["others"]["share"]) < 0.3               # and the rows of units that do not differ teach nothing false
    weak = _acc(d, un, first, "weak", controls=False)
    assert weak["chain"]["says"].startswith("no change to account for") and weak["plan"] == []


def test_a_group_of_units_against_another_is_resampled_by_unit():
    d = toy.planted_groups(1, n_healthy=10, n_patients=6, n_cell=200)
    un = np.where(d["treated"], "disease", "control")
    first = cycle.halves(d["unit"], 0)
    kw = dict(inputs=["view_a"], state="z", o=d["o"], draws=80, controls=False)
    a = AC.account({"view_a": d["view_a"], "z": d["z"]}, d["y"], un, first, "disease", cluster=d["unit"], **kw)
    z = a["stages"]["z"]
    assert z["cells"]["units_of_the_event"] == 6 and z["cells"]["units_of_the_reference"] == 10
    assert z["shares"]["read"]["share"] > 0.9 and not a["chain"]["open"] and z["shares"]["as_trained"]["share"] < 0.1
    b = AC.account({"view_a": d["view_a"], "z": d["z"]}, d["y"], un, first, "disease", **kw)        # the same cells counted one by one
    assert b["stages"]["z"]["snr"] > 1.2 * z["snr"]                                  # donors, not cells, are the units: the change is less certain than cells say
    wa, wb = (np.diff(x["stages"]["z"]["shares"]["read"]["ci"])[0] for x in (a, b))
    assert wa > wb


def test_the_account_with_control_cells_is_what_it_was():
    d = T.planted(2)
    first = cycle.halves(d["unit"], 0)
    r = AC.reading(d["Z"], d["y"], d["unit"], first, "r1", others=["r2", "c1", "c2"], draws=60)
    assert r["read_against"] == "control cells" and "units_of_the_event" not in r["cells"]
    assert r["shares"]["others"]["share"] > 0.8 and r["shares"]["read"]["share"] < 0.4 and len(r["widening_by_the_rows_of_the_others"]) == 2


def test_edges_are_records_declared_proposed_and_tested():
    decl = {"readouts": [{"id": "PDL1", "features": ["CD274"]}, {"id": "MHC1", "features": ["HLA-A", "HLA-B", "HLA-C"]}]}
    own = EV.own_edges(decl, ["CD274", "HLA-B", "B2M", "JAK1"], partners={"MHC1": ["B2M"]}, scale="signal above background")
    assert {(e["site"], e["readout"]) for e in own} == {("CD274", "PDL1"), ("HLA-B", "MHC1"), ("B2M", "MHC1")}
    assert all(e["status"] == "declared" and e["class"] == "own" and e["basis"] == "declaration" for e in own) and not EV.check_edges(own)
    flags = [{"event": "TAPBP", "readout": "MHC1", "flag": "a missing edge", "change": -1.2, "forecast": -0.1, "z": -9.0},
             {"event": "STAT3", "readout": "MHC1", "flag": "a direction the law does not cover", "change": -0.4, "forecast": 0.0, "z": -4.0},
             {"event": "CD274", "readout": "PDL1", "flag": "a missing edge", "change": -2.0, "forecast": 0.0, "z": -20.0}]
    prop = EV.propose_edges(flags, known=own)
    assert [(e["site"], e["status"], e["basis"]) for e in prop] == [("TAPBP", "proposed", "flag")]          # a direction is no edge; a known pair is not proposed again
    assert EV.edge_table(own + prop) == {("CD274", "PDL1"): "own", ("HLA-B", "MHC1"): "own", ("B2M", "MHC1"): "own"}      # a proposed edge enters no forecast
    rows = [{"event": "TAPBP", "readout": "MHC1", "change": -1.9, "forecast": 0.05, "sd": 0.3}]
    tested = EV.test_edges(prop, rows)
    assert tested[0]["status"] == "admitted" and tested[0]["test"]["z"] < -1.96 and ("TAPBP", "MHC1") in EV.edge_table(own + tested)
    assert EV.test_edges(prop, [dict(rows[0], change=0.1)])[0]["status"] == "refused"
    assert EV.test_edges(prop, [dict(rows[0], change=1.9)])[0]["status"] == "undecided"                     # outside, on the other side
    assert EV.test_edges(prop, [])[0]["status"] == "undecided"
    own[0]["share"] = {"value": -0.56, "measured_in": "another condition"}                                  # what an edge carries from one condition to another
    assert EV.edge_shares(own + prop, 0.1) == {(own[0]["site"], own[0]["readout"]): (-0.56, 0.1)}
    bad = EV.check_edges([dict(prop[0], status="admitted"), EV.edge("X", "PDL1", None), own[0], own[0]])
    assert {b["code"] for b in bad} >= {"G1", "G2", "G5"}


def test_records_teach_within_a_type_and_no_other():
    d = T.planted(5)
    types = {"r1": "R", "r2": "R", "c1": "C", "c2": "C"}
    one = AC.Law(d["Z"], {"y": d["y"]}, d["unit"], edges=d["edges"]).fit()
    typed = AC.Law(d["Z"], {"y": d["y"]}, d["unit"], edges=d["edges"], types=types).fit()
    rows = {r["event"]: r for r in typed.forecast()}
    for e in ("r1", "r2"):                                                           # each is forecast from the other of its type
        assert abs(rows[e]["change"] - rows[e]["forecast"]) < 0.1 and abs(rows[e]["change"] - rows[e]["read"]) > 0.4 and rows[e]["type"] == "R"
    fill = [r for r in rows.values() if r["event"].startswith("q")]
    assert all(abs(r["from_the_other_events"]) < 1e-12 for r in fill)                # an event with no type is taught by no row
    rm = lambda rs: float(np.sqrt(np.mean([(r["change"] - r["forecast"]) ** 2 for r in rs])))
    assert rm(fill) <= rm([r for r in one.forecast() if r["event"].startswith("q")]) + 1e-9
    keep = d["unit"] != "r1"
    part = AC.Law(d["Z"][keep], {"y": d["y"][keep]}, d["unit"][keep], edges=d["edges"], types=types).fit()
    new = part.new(d["Z"][~keep], d["unit"][~keep], types={"r1": "R"})[0]
    none = part.new(d["Z"][~keep], d["unit"][~keep])[0]
    assert new["taught_by"] == 1 and abs(new["forecast"] - rows["r1"]["change"]) < 0.15
    assert none["taught_by"] == 0 and abs(none["from_the_recorded_events"]) < 1e-12 and abs(none["forecast"] - none["read"]) < 1e-12
    # the widening is a type's own: an event of no type keeps the control prior, so what it shows beyond it is flagged and not hidden
    assert set(typed.r["y"]["tau2"]) == {"R", "C"} and none["widening"] == 1.0 and new["widening"] == part.r["y"]["tau2"]["R"] > 1
    flags = typed.flags()
    assert {f["event"] for f in flags if f["flag"] == "a missing edge"} >= {"x"} and not any(f["event"].startswith("q") for f in flags)
    prop = EV.propose_edges(flags, known=[EV.edge(e, r, c) for (e, r), c in d["edges"].items()])          # the law's flags are what the graph reads
    assert "x" in {e["site"] for e in prop} and all(e["status"] == "proposed" and e["evidence"]["z"] is not None for e in prop)


def test_an_edge_has_a_size_its_class_share_of_a_scale():
    rng = np.random.default_rng(0)
    n_c, n_e, p = 3000, 300, 12
    th = rng.normal(size=p) * 0.3
    names = [f"q{k:02d}" for k in range(14)] + [f"k{k}" for k in range(8)]
    unit = np.array(["control"] * n_c + [e for e in names for _ in range(n_e)])
    x = rng.normal(size=(len(unit), p))
    y = x @ th + 0.5 * rng.normal(size=len(unit))
    sd = y[:n_c].std()
    share, scale, edges = -0.6, {}, {}
    for k, e in enumerate(names):
        m = unit == e
        x[m] += 0.15 * rng.normal(size=p)
        y[m] += (x[m].mean(0) - x[:n_c].mean(0)) @ th * 0.0                          # its shift is read by the control readout already
        if e.startswith("k"):                                                       # an edge: the readout loses a share of what it has to lose in that condition
            s = [0.4, 0.8, 1.2, 1.6, 2.0, 2.6, 3.2, 4.0][int(e[1:])]
            scale[(e, "y")], edges[(e, "y")] = s, "own"
            y[m] += share * s * sd
    flat = {r["event"]: r for r in AC.Law(x, {"y": y}, unit, edges=edges).fit().forecast()}
    law = AC.Law(x, {"y": y}, unit, edges=edges, scale=scale).fit()
    sized = {r["event"]: r for r in law.forecast()}
    ks = [e for e in names if e.startswith("k")]
    err = lambda rows: float(np.sqrt(np.mean([(rows[e]["change"] - rows[e]["forecast"]) ** 2 for e in ks])))
    assert abs(law.cl["own"]["mu"] - share) < 0.08 and err(sized) < 0.25 * err(flat)
    assert all(abs(v - share) < 0.15 for _e, _r, v in law.class_values()["own"])
    keep = unit != "k7"
    part = AC.Law(x[keep], {"y": y[keep]}, unit[keep], edges=edges, scale=scale).fit()
    f = part.new(x[~keep], unit[~keep], edges=edges, scale=scale)[0]
    change, _se = part.measure(y[~keep], unit[~keep], "y")["k7"]
    assert f["scale"] == 4.0 and abs(f["forecast"] - change) < 0.4 and f["interval"][0] < change < f["interval"][1]
    assert abs(part.new(x[~keep], unit[~keep], edges=edges)[0]["forecast"] - change) > 1.2          # with one size per class it is forecast at the class's mean


def test_a_pair_measured_in_one_condition_is_forecast_in_another_at_its_own_share():
    """Pairs of one class whose shares differ (one gene of several behind a protein removes less of it than the only
    gene): the class's mean share misses each; the pair's own share, measured where the scale is one value, forecasts it
    where the scale is another."""
    def condition(seed, scale_of):
        rng = np.random.default_rng(seed)
        n_c, n_e, p = 3000, 300, 10
        th = rng.normal(size=p) * 0.3
        names = [f"q{k:02d}" for k in range(12)] + [f"k{k}" for k in range(6)]
        unit = np.array(["control"] * n_c + [e for e in names for _ in range(n_e)])
        x = rng.normal(size=(len(unit), p))
        for e in names:
            x[unit == e] += 0.15 * rng.normal(size=p)                                # every event moves the stage, and the readout with it
        y = x @ th + 0.5 * rng.normal(size=len(unit))
        sd = y[:n_c].std()
        for k, sh in enumerate(SHARES):
            y[unit == f"k{k}"] += sh * scale_of * sd
        return x, y, unit
    SHARES = [-0.9, -0.6, -0.1, -0.15, -0.3, -0.75]
    edges = {(f"k{k}", "y"): "own" for k in range(6)}
    xa, ya, ua = condition(1, 2.0)                                                   # where the pairs were measured
    a = AC.Law(xa, {"y": ya}, ua, edges=edges, scale={k: 2.0 for k in edges}).fit()
    seen = {(e, r): (v, 0.05) for e, r, v in a.class_values()["own"]}
    assert all(abs(seen[(f"k{k}", "y")][0] - sh) < 0.08 for k, sh in enumerate(SHARES))
    xb, yb, ub = condition(2, 5.0)                                                   # another condition: the readout has more to lose
    new = np.isin(ub, [f"k{k}" for k in range(6)])
    b = AC.Law(xb[~new], {"y": yb[~new]}, ub[~new]).fit()
    kw = dict(edges=edges, scale={k: 5.0 for k in edges})
    own = b.new(xb[new], ub[new], shares=seen, **kw)
    cls = b.new(xb[new], ub[new], classes={"own": [v for v, _s in seen.values()]}, **kw)
    ch = {e: v[0] for e, v in b.measure(yb[new], ub[new], "y").items()}
    err = lambda rows: float(np.sqrt(np.mean([(ch[r["event"]] - r["forecast"]) ** 2 for r in rows])))
    assert err(own) < 0.25 * err(cls) and all(r["size_of_the_edge_from"].startswith("the pair's own share") for r in own)
    assert sum(r["interval"][0] < ch[r["event"]] < r["interval"][1] for r in own) >= 5 and cls[0]["size_of_the_edge_from"] == "its class"


def test_the_step_after_a_rejection_sends_each_open_part_to_its_lever():
    d, un, first = _units()
    accs = [{"readout": "y", "group": "units", "type": {"type": "units"}, "event": e, "account": _acc(d, un, first, e, controls=False)} for e in ("h1", "early", "post", "solo")]
    items = outer.open_parts(accs)
    by = {i["event"]: i for i in items}
    assert "h1" not in by                                                            # nothing open, no item
    assert by["early"]["lever"] == "M" and by["early"]["reading"] == "lost before the state" and by["early"]["state"] == "open"
    assert by["post"]["lever"] == "G" and by["post"]["next_step_needs"].startswith("no new data to propose")
    assert by["solo"]["lever"] == "E" and by["solo"]["reading"] == "held as the identity"
    assert outer.counts(items)["by_lever"]["G"] == 1 and [i["id"] for i in items] == ["P01", "P02", "P03"]
    with_edge = outer.open_parts(accs, edges=[EV.edge("post", "y", "own", provenance="declared for the test")])
    assert {i["event"]: i["state"] for i in with_edge}["post"] == "on the graph"
    assert set(outer.PARTS4.values()) == set(outer.ROUTES4) == set(outer.NEEDS4) and not set(outer.ROUTES4) & set(outer.ROUTES)      # the parts of an account have their own table
    assert {v[0] for v in outer.ROUTES4.values()} == {"M", "E", "G"} and all(k in outer.LEVERS for k in "MEG")


def test_a_forecast_that_overshoots_is_no_forecast():
    """`others` is a share of the change: 14 times the change misses it as surely as none of it (found by the planted
    study of units, where a readout of other units overshot a unit it had no row like)."""
    r = lambda lo, hi: AC.across_events({"shares": {"others": {"share": (lo + hi) / 2, "ci": [lo, hi]}}})["reading"]
    assert r(0.80, 1.10) == "the other events forecast it"
    assert r(-0.10, 0.40) == "the other events do not forecast it"
    assert r(3.2, 27.0) == "the other events do not forecast it"
    assert r(1.10, 1.42) == "unresolved" and r(0.40, 0.90) == "unresolved"


def test_a_unit_with_no_control_cells_is_forecast_from_the_other_units():
    """The law with no control cells: a unit left out is forecast from its state and the rows of the other units. A unit
    whose difference is a factor that varies inside units is forecast; an ordinary unit stays inside its interval; a
    unit whose difference is not in the state stands far outside, and does not bend the forecasts of the others."""
    d = toy.planted_units(0, n_ordinary=40, n_cell=300)
    un = d["unit"]
    out = {}
    for u in ["h1", "h2", "post", "modb"] + [f"u{k:02d}" for k in range(6)]:
        keep = un != u
        law = AC.Law(d["z"][keep], {"y": d["y"][keep]}, un[keep], control=None).fit()
        assert law.units_only and "control" not in law.events and len(law.events) == len(set(un[keep].tolist()))
        r = law.new(d["z"][~keep], un[~keep])[0]
        ch, _se = law.measure(d["y"][~keep], un[~keep], "y")[u]
        out[u] = (ch, r["forecast"], (ch - r["forecast"]) / r["sd"])
    assert all(abs(out[u][0]) > 1.5 and abs(out[u][0] - out[u][1]) < 0.25 for u in ("h1", "h2"))          # read from the state
    assert all(abs(out[u][2]) > 10 and abs(out[u][1]) < 0.3 for u in ("post", "modb"))                     # not in the state: no forecast of it, and it shows
    assert sum(abs(out[f"u{k:02d}"][2]) < 1.96 for k in range(6)) >= 5 and max(abs(out[f"u{k:02d}"][1]) for k in range(6)) < 0.2


def test_a_scale_that_measures_nothing_is_said_beforehand():
    """An edge sized by a scale keeps its class where the readout stands above its background by a spread or more; where
    it does not, the edge is of another class and has no scale. How well a share carries is asked of pairs measured twice."""
    table = {("A", "p"): "own", ("B", "q"): "own", ("C", "r"): "own"}
    out, scales = EV.by_scale(table, {"p": 3.1, "q": 0.46})
    assert out == {("A", "p"): "own", ("B", "q"): "own, at the background", ("C", "r"): "own"} and scales == {("A", "p"): 3.1}
    assert EV.carried_spread([(-0.56, -0.63), (-0.14, -0.13)]) is None
    assert abs(EV.carried_spread([(-0.5, -0.6), (-0.1, -0.2), (0.3, 0.2)]) - 0.1) < 1e-12
