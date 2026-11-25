"""The account of a change on the planted system of toy_account.py: every part where it was planted."""
import numpy as np

from anm_diag import account as AC
from anm_diag import cycle
from anm_diag import toy_account as T


def test_planted_system_as_planted():
    checks, rep = T.run(0, draws=80)
    wrong = [k for k, v in checks.items() if not v]
    assert len(wrong) <= 1, wrong                               # 73 checks, each with a 95% interval of 80 resamples
    assert rep["l"]["says"].startswith("lost before the state")
    assert rep["r1"]["says"].startswith("held not read")
    assert rep["x"]["says"].startswith("not in the input")
    assert rep["c1"]["says"].startswith(("read", "unresolved"))  # no part is open


def test_a_readout_made_to_follow_the_mean_is_not_a_measure():
    d = T.planted(0)
    un = d["unit"]
    m = np.isin(un, ["control", "m"])
    r = AC.reading(d["Z"][m], d["y"][m], un[m], cycle.halves(un, 0)[m], "m", draws=60)
    assert r["shares"]["found_again"]["share"] > 0.7            # the event's shift of the state is found again
    assert r["shares"]["carried"]["share"] < 0.2                # and the state holds almost none of its change
    assert r["carried_by"] == "cells"


def test_no_overlap_is_said_and_not_read_as_direct():
    d = T.planted(0, edited=1.0)
    un = d["unit"]
    m = np.isin(un, ["control", "r1"])
    r = AC.reading(d["Z"][m], d["y"][m], un[m], cycle.halves(un, 0)[m], "r1", draws=40)
    assert r["overlap"] < AC.MIN_OVERLAP and r["carried_by"] == "identity" and r["shares"]["carried"]["share"] == 1.0


def test_the_law_forecasts_by_an_edge_and_flags_a_missing_one():
    d = T.planted(1)
    law = AC.Law(d["Z"], {"y": d["y"]}, d["unit"], edges=d["edges"]).fit()
    rows = law.forecast()
    by = {r["event"]: r for r in rows}
    assert sum(abs(by[e]["z"]) < 1.96 for e in T.EDGES) >= 3
    flagged = {f["event"]: f["candidate"] for f in law.flags(rows)}
    assert flagged.get("x", "").startswith("a missing edge")
    assert not any(e.startswith("q") for e in flagged)
    cal = law.calibration(rows)
    assert cal["inside_95"] >= 0.85 and cal["by_class"]["own"]["error_forecast"] < cal["by_class"]["own"]["error_without_edges"]


def test_parts_sum_to_one_and_gains_are_paired():
    d = T.planted(2)
    un = d["unit"]
    m = np.isin(un, ["control", "l"])
    first = cycle.halves(un, 2)[m]
    R = {st: AC.reading(d[st][m], d["y"][m], un[m], first, "l", draws=60, seed=2) for st in ("V", "Z")}
    ch = AC.chain(R, ["V"], "Z")
    assert abs(sum(p["share"] for p in ch["parts"].values()) - 1) < 1e-9
    g = AC.gain(R["Z"], R["V"])
    assert g["reading"] == "a step" and abs(g["gain"] - ch["parts"]["lost_before_the_state"]["share"]) < 1e-9


def test_an_event_not_yet_measured_is_forecast_as_when_left_out():
    """Law.new on the cells of an event equals the forecast of that event left out, when the rest is held the same."""
    d = T.planted(5)
    full = AC.Law(d["Z"], {"y": d["y"]}, d["unit"], edges=d["edges"]).fit()
    rows = {r["event"]: r for r in full.forecast()}
    for e in ("c1", "l", "n2", "q07"):
        keep = d["unit"] != e
        part = AC.Law(d["Z"][keep], {"y": d["y"][keep]}, d["unit"][keep], edges=d["edges"]).fit()
        i = [k for k, h in enumerate(full.events) if h != e]
        part.sigma0, part.cl = full.sigma0, full.cl                          # the same spreads, weights and widening
        part.r["y"]["tau2"], part.r["y"]["w"] = full.r["y"]["tau2"], full.r["y"]["w"][i]
        vals = {c: [v for h, _r, v in rows_ if h != e] for c, rows_ in full.class_values().items()}
        f = part.new(d["Z"][~keep], d["unit"][~keep], edges=d["edges"], classes=vals)[0]
        r = rows[e]
        assert abs(f["read"] - r["read"]) < 1e-9 and abs(f["from_the_recorded_events"] - r["from_the_other_events"]) < 1e-6
        if r["no_forecast"] is None:
            assert abs(f["by_its_edge"] - r["by_its_edge"]) < 1e-9 and abs(f["forecast"] - r["forecast"]) < 1e-6
            assert abs(f["sd"] - r["sd"]) < 0.05 * r["sd"] + 0.02                # the row's error is taken at the control spread
        change, se = part.measure(d["y"][~keep], d["unit"][~keep], "y")[e]
        assert abs(change - r["change"]) < 1e-9 and abs(se - r["se"]) < 0.5 * r["se"]


def test_the_command_reads_a_spec_and_credits_a_repair_up_to_what_the_stage_carries(tmp_path):
    """`anm_diag account` on a SPEC file; a repair that follows the mean of an event no stage holds gets no credit."""
    d = T.planted(2)
    np.savez(tmp_path / "cells.npz", Z=d["Z"], V=d["V"], unit=d["unit"], A=(d["unit"] != "control").astype(int), y=d["y"])
    spec = {"name": "planted", "events": ["c1", "m"], "arrays": {"audit": "cells.npz"}, "account": {"unit_length": [], "others": ["c1", "c2", "r1", "r2"], "pairs": [["Z", "V"]]},
            "keys": {"readouts": ["y"], "stages": {"V": "V", "Z": "Z"}, "unit": "unit", "control": {"key": "A", "equals": 0}}}
    res = AC.run_spec(spec, tmp_path, draws=60)
    assert res["state"] == "Z" and res["inputs"] == ["V"] and set(res["events"]) == {"c1", "m"}
    m, c1 = res["events"]["m"], res["events"]["c1"]
    z = m["stages"]["Z"]["shares"]
    assert z["found_again"]["share"] > 0.4 and z["carried"]["share"] < 0.25
    rep = m["a_repair_that_follows_the_mean"]
    assert rep["credited"]["share"] < 0.25 and rep["beyond_what_the_stage_carries"]["reading"] == "a step"
    assert c1["stages"]["Z"]["shares"]["read"]["share"] > 0.8 and c1["a_repair_that_follows_the_mean"]["beyond_what_the_stage_carries"]["reading"] != "a step"
    assert "V over Z" in c1["gains"] and "_draws" not in c1["stages"]["Z"]
