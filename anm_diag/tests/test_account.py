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
    assert rep["x"]["says"].startswith("direct")
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
