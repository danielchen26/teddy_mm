"""Tests for bridge_anm/lib/v3_e1.py (E1 helpers) against v3_amend / v3_key and the ANM engine.

Synthetic tests always run; the engine test skips when ANM_ROOT is unavailable.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bridge_anm"))

from lib import v3_amend as va  # noqa: E402
from lib import v3_e1 as e1  # noqa: E402
from lib import v3_key as vk  # noqa: E402

ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))


@pytest.fixture(scope="module")
def reg():
    return va.load_registration_amended()


@pytest.fixture(scope="module")
def amend(reg):
    return reg["amendment_A1"]


# ----------------------------------------------------------------------------- matched coverage
def test_select_ordered_equals_coverage_select_on_materialised_replicates(amend):
    rng = np.random.default_rng(0)
    n = 300
    ids = np.sort(rng.choice(90261, size=n, replace=False))
    conf = np.round(rng.uniform(0, 1, n), 1)          # many ties
    conf[:40] = 1.0
    order = e1.ranked_order(conf, ids, amend)
    for rep in range(25):
        w = np.bincount(rng.integers(0, n, n), minlength=n) if rep else np.ones(n, dtype=np.int64)
        pos = np.repeat(np.arange(n), w)
        for c in (0.95, 0.7, 0.5, 0.333):
            mask = va.coverage_select(conf[pos], ids[pos], c, amend)
            want = np.bincount(pos[mask], minlength=n)
            s_o = e1.select_ordered(w[order], [e1.coverage_k(c, w.sum())])[0]
            got = np.zeros(n)
            got[order] = s_o
            np.testing.assert_array_equal(got, want)


def test_ranked_order_breaks_ties_by_registered_permutation(amend):
    ids = np.arange(1000, 1100)
    order = e1.ranked_order(np.zeros(100), ids, amend)
    rank = va.tie_break_rank(ids, amend)
    assert np.all(np.diff(rank[order]) > 0)
    order2 = e1.ranked_order(np.r_[np.ones(50), np.zeros(50)], ids, amend)
    assert set(order2[:50]) == set(range(50))


def test_select_lowest(amend):
    ids = np.arange(10)
    conf = np.array([0.5, 0.1, 0.1, 0.9, 0.2, 0.0, 0.3, 0.3, 0.4, 0.8])
    elig = np.array([1, 1, 1, 1, 0, 1, 1, 1, 0, 1], bool)
    o = e1.select_lowest(conf, elig, ids, amend)
    assert o[0] == 5 and set(o[1:3]) == {1, 2} and 4 not in o and 8 not in o
    assert np.all(np.diff(conf[o]) >= 0)


def test_aurc():
    grid = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5]
    assert e1.aurc([0.8] * 10, grid) == pytest.approx(0.8)
    lin = np.linspace(0.9, 1.0, 10)
    assert e1.aurc(lin, grid) == pytest.approx(0.95)


# ----------------------------------------------------------------------------- bootstrap and verdicts
def test_replicate_counts_two_and_one_stage():
    groups = [np.arange(0, 7), np.arange(7, 10)]
    reps = list(e1.replicate_counts(groups, 10, 200, seed=1))
    sizes = {int(w.sum()) for w in reps}
    assert sizes <= {6, 10, 14} and len(sizes) > 1          # donors drawn with replacement
    again = list(e1.replicate_counts(groups, 10, 200, seed=1))
    assert all(np.array_equal(x, y) for x, y in zip(reps, again))
    one = list(e1.replicate_counts([groups[1]], 10, 50, seed=1, two_stage=False))
    assert all(w.sum() == 3 and w[:7].sum() == 0 for w in one)


def test_verdict_rules():
    rng = np.random.default_rng(0)
    b = rng.normal(0.05, 0.005, 2000)
    assert e1.verdict(0.05, b, {"a": 0.04, "b": 0.06}, 0.01)["verdict"] == "win"
    assert e1.verdict(0.05, b, {"a": 0.004, "b": 0.06}, 0.01)["verdict"] == "inconclusive"   # a donor below margin/2
    assert e1.verdict(-0.05, -b, {"a": -0.04, "b": -0.06}, 0.01)["verdict"] == "loss"
    tight = rng.normal(0.0, 0.001, 2000)
    assert e1.verdict(0.0, tight, {"a": 0.0, "b": 0.0}, 0.01)["verdict"] == "equivalent"
    wide = rng.normal(0.0, 0.05, 2000)
    assert e1.verdict(0.0, wide, {"a": 0.0, "b": 0.0}, 0.01)["verdict"] == "inconclusive"
    assert e1.verdict(float("nan"), b, {"a": 0.0}, 0.01)["verdict"] == "inconclusive"


# ----------------------------------------------------------------------------- scores and classifiers
def test_entropy_confidence():
    S = np.array([[0.5, 0.5, 0.5, 0.5], [1.0, 0, 0, 0], [0, 0, 0, 0]])
    c = e1.entropy_confidence(S)
    assert c[0] == pytest.approx(0.0, abs=1e-12) and c[1] == pytest.approx(1.0) and c[2] == pytest.approx(0.0, abs=1e-12)


def test_tied_top_and_argmax_first():
    S = np.array([[1.0, 1.0, 0.2, 0.0], [0.1, 0.3, 0.2, 0.0]])
    assert e1.tied_top_classes(S).tolist() == [True, False]
    assert e1.argmax_calls(S).tolist() == ["B", "T"]


def test_fit_logreg_single_class_and_proba5():
    X = np.random.default_rng(0).normal(size=(20, 3))
    m = e1.fit_logreg(X, np.array(["T"] * 20), 1.0)
    P = e1.proba5(m, X)
    assert np.all(P[:, vk.CLASSES.index("T")] == 1.0)
    conf, call = e1.lineage_conf_call(P)
    assert np.all(call == "T") and np.all(conf == 1.0)
    y = np.array(["B", "OUT"] * 10)
    m2 = e1.fit_logreg(X, y, 1.0)
    P2 = e1.proba5(m2, X)
    assert np.allclose(P2.sum(axis=1), 1.0) and np.all(P2[:, vk.CLASSES.index("NK")] == 0)
    conf2, call2 = e1.lineage_conf_call(P2)
    assert set(call2) <= {"B", "T", "NK", "myeloid"}          # OUT is never called


def test_nested_draw():
    pool = np.arange(100, 200)
    a = e1.nested_draw(pool, 25, 3)
    b = e1.nested_draw(pool, 50, 3)
    assert a.size == 25 and np.array_equal(b[:25], a)
    assert np.array_equal(e1.nested_draw(pool[::-1], 25, 3), a)   # pool sorted by id first
    assert np.array_equal(e1.nested_draw(pool, "all", 3), pool)
    assert not np.array_equal(e1.nested_draw(pool, 25, 4), a)


# ----------------------------------------------------------------------------- E1.3
def test_top_sets_and_shares():
    v = np.array([[1.0, 1.0, 0.2], [0.3, 0.9, 0.1], [0.5, 0.5, 0.5]])
    m = e1.closed_form_top_sets(v)
    assert m.tolist() == [[True, True, False], [False, True, False], [True, True, True]]
    dec = v / 3.0
    assert np.array_equal(e1.engine_top_sets(dec), m)
    sh = e1.marker_shares(m, np.array([0, 0, 1]))
    np.testing.assert_allclose(sh["all"][0], [0.5 / 3, 1.5 / 3, 0.0])
    np.testing.assert_allclose(sh["within_class"][0], [0.25, 0.75, 0.0])
    np.testing.assert_allclose(sh["within_class"][1], [1 / 3] * 3)


def test_marker_share_null_mean_is_one_third():
    rng = np.random.default_rng(0)
    mask = np.zeros((600, 3), bool)
    mask[np.arange(600), rng.integers(0, 3, 600)] = True
    mask[:100] = True
    cls = rng.integers(0, 4, 600)
    null = e1.marker_share_null(mask, cls, 200, seed=5)
    np.testing.assert_allclose(np.nanmean(null["within_class"], axis=0), 1 / 3, atol=0.01)
    np.testing.assert_allclose(null["all"].sum(axis=(1, 2)), 1.0)
    obs = e1.marker_shares(np.tile([True, False, False], (600, 1)), cls)["within_class"][:, 0]
    assert all(e1.perm_p_two_sided(o, null["within_class"][:, k, 0]) < 0.01 for k, o in enumerate(obs))


# ----------------------------------------------------------------------------- ANM engine
def _bridge(reg):
    try:
        return e1.AnmBridge(reg, REPO / "bridge_anm/schemas/cite_lineage_finite_field_v0.json", ANM_ROOT)
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"ANM engine not importable ({exc})")


def test_engine_equals_recoded_rule_closure_and_leave_one_out(reg):
    br = _bridge(reg)
    names = reg["evidence"]["teddy_head"]["adt_names"]
    j = {n: i for i, n in enumerate(names)}
    p1, p3 = vk.question_panel(reg, "Q1"), vk.question_panel(reg, "Q3")
    rng = np.random.default_rng(11)
    v = rng.uniform(0, 1, size=(120, len(names)))
    v[rng.random(v.shape) < 0.3] = 1.0           # clipped ties, as in the real evidence
    v[:20] *= 0.35                                # some no-calls
    S1, S3 = vk.class_scores(v, names, reg, "Q1"), vk.class_scores(v, names, reg, "Q3")
    clr = vk.anm_action_scores(v, names, reg, "Q1", readout="closure")
    fld = reg["anm"]["field_representation"]
    for i in range(v.shape[0]):
        vals1 = {k: [v[i, j[p]] for p in p1[k]] for k in vk.LINEAGES}
        o = br.run(p1, vals1, readouts=("Q1", "Q2"), closure_readouts=("closure_Q1", "closure_Q2"))
        assert o["n_rejected"] == 0
        np.testing.assert_allclose(o["Q1"][0], vk.field_gain(3, fld) * 3 * S1[i], rtol=1e-12)
        np.testing.assert_allclose(o["closure_Q1"][0], clr[i], rtol=1e-12)
        for q in ("Q1", "Q2"):
            rc = vk.rule_calls(S1[i:i + 1], reg["questions"][q]["bar"])[0]
            assert (vk.LINEAGES[o[q][1]] if o[q][1] >= 0 else "") == rc
            cc = vk.rule_calls(clr[i:i + 1], reg["anm"][f"closure_bar_{q}"])[0]
            assert (vk.LINEAGES[o[f"closure_{q}"][1]] if o[f"closure_{q}"][1] >= 0 else "") == cc
        o3 = br.run(p3, {k: [v[i, j[p]] for p in p3[k]] for k in vk.LINEAGES}, readouts=("Q3",))
        rc3 = vk.rule_calls(S3[i:i + 1], reg["questions"]["Q3"]["bar"])[0]
        assert (vk.LINEAGES[o3["Q3"][1]] if o3["Q3"][1] >= 0 else "") == rc3
        c = o["Q1"][1]
        if c < 0:
            continue
        k = vk.LINEAGES[c]
        vc = np.asarray(vals1[k])
        dec, newc = [], []
        for p in p1[k]:
            od = br.run(p1, vals1, readouts=("Q1",), drop=(k, p))
            dec.append(o["Q1"][0][c] - od["Q1"][0][c])
            newc.append(od["Q1"][1])
        eng = e1.engine_top_sets(np.asarray(dec)[None])[0]
        assert np.array_equal(eng, e1.closed_form_top_sets(vc[None])[0])
        flip_engine = newc[int(np.argmax(eng))] != c
        flip_cf = bool(va.loo_flip_exact(S1[i:i + 1], vc[None], np.array([c]), reg["questions"]["Q1"]["bar"], fld)[0])
        assert flip_engine == flip_cf
        assert len({newc[t] for t in np.where(eng)[0]}) == 1     # tied top events give the same call
