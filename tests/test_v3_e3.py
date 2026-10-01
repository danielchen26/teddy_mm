"""Tests for E3 (bridge_anm/lib/v3_e3.py, v3_e3_readouts.py, scripts/v3_e3_nkt_repair.py).

Synthetic tests always run; tests that need the processed data or the registration skip without them.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bridge_anm"))
sys.path.insert(0, str(REPO))

from lib import v3_amend as va  # noqa: E402
from lib import v3_e3 as e3  # noqa: E402
from lib.knn_cosine import knn_cosine  # noqa: E402

DATA = REPO / "data" / "processed"
REG = REPO / "registration" / "registration_v3.json"


def _builder_pairs(z, k):
    """The registration builder's construction (v3_build_registration.py, E3 flag)."""
    nn_idx, nn_sim = knn_cosine(z, k)
    pairs = {}
    for i in range(nn_idx.shape[0]):
        for jj, s in zip(nn_idx[i], nn_sim[i]):
            a_, b_ = (i, int(jj)) if i < jj else (int(jj), i)
            pairs[(a_, b_)] = float(s)
    return pairs


# ----------------------------------------------------------------------------- pairs
def test_knn_edges_equal_builder_dict():
    rng = np.random.default_rng(0)
    z = rng.normal(size=(300, 16)).astype(np.float32)
    E, s = e3.knn_edges(z, 10)
    ref = _builder_pairs(z, 10)
    got = {(int(a), int(b)): float(c) for (a, b), c in zip(E, s)}
    assert got.keys() == ref.keys()
    assert all(abs(got[k] - ref[k]) < 1e-7 for k in ref)
    assert np.all(E[:, 0] < E[:, 1])


def test_nkt_pairs_within_donors_never_cross_donor_and_match_per_donor_build():
    rng = np.random.default_rng(1)
    n = 400
    z = rng.normal(size=(n, 8)).astype(np.float32)
    donors = np.where(np.arange(n) < 250, "d1", "d2")
    key = rng.choice(["NK", "T", "B", "unscored"], size=n, p=[0.3, 0.4, 0.2, 0.1])
    cells = np.arange(20, n)  # cells 0-19 are not in the evaluated split
    P = e3.nkt_pairs_within_donors(z, cells, donors, {"all": key}, 10, 0.5)
    pv = P["variants"]["all"]
    assert np.all(donors[pv["nk"]] == donors[pv["t"]])
    assert np.all(pv["donor"] == donors[pv["nk"]])
    assert np.all(key[pv["nk"]] == "NK") and np.all(key[pv["t"]] == "T")
    assert np.all(pv["nk"] >= 20) and np.all(pv["t"] >= 20)
    assert np.array_equal(pv["flag"], pv["cos"] >= 0.5)
    # per-donor reference
    want = set()
    for d in ("d1", "d2"):
        rows = np.sort(cells[donors[cells] == d])
        for (a_, b_), s in _builder_pairs(z[rows], 10).items():
            ga, gb = rows[a_], rows[b_]
            if {key[ga], key[gb]} == {"NK", "T"}:
                want.add((int(ga if key[ga] == "NK" else gb), int(gb if key[ga] == "NK" else ga)))
    assert set(zip(pv["nk"].tolist(), pv["t"].tolist())) == want


def test_pairs_sha256_is_order_free():
    nk, t = np.array([5, 1, 9]), np.array([2, 7, 3])
    o = np.array([2, 0, 1])
    assert e3.pairs_sha256(nk, t) == e3.pairs_sha256(nk[o], t[o])
    assert e3.pairs_sha256(nk, t) == e3.pairs_sha256(t, nk)  # stored as (lo, hi) codes, as E5


# ----------------------------------------------------------------------------- gap ratio and D
def test_gap_ratio_known_values():
    e_nk = np.array([[1.0, 0.5], [0.8, 0.5], [0.6, 0.5]])
    e_t = np.array([[0.0, 0.5], [0.4, 0.3], [0.5, 0.1]])
    m_nk = np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]])
    m_t = np.zeros((3, 2))
    gr, per = e3.gap_ratio(e_nk, e_t, m_nk, m_t)
    assert np.allclose(per, [0.4, 0.2])  # medians of |diff|: (1, .4, .1) -> .4 ; (0, .2, .4) -> .2
    assert math.isclose(gr, 0.3)
    gr0, _ = e3.gap_ratio(e_nk[:0], e_t[:0], m_nk[:0], m_t[:0])
    assert math.isnan(gr0)
    gz, perz = e3.gap_ratio(e_nk, e_t, m_nk * 0, m_t)
    assert math.isnan(gz) and np.all(np.isnan(perz))


def test_d_statistic_algebra_head_cancels():
    rng = np.random.default_rng(2)
    gr = {r: {"flagged": float(rng.uniform()), "unflagged": float(rng.uniform())} for r in ("head", "R1", "null")}
    d = e3.d_statistic(gr, "R1")
    alt = (gr["R1"]["flagged"] - gr["null"]["flagged"]) - (gr["R1"]["unflagged"] - gr["null"]["unflagged"])
    assert math.isclose(d, alt, abs_tol=1e-12)
    gr2 = {**gr, "head": {"flagged": 9.0, "unflagged": -3.0}}
    assert math.isclose(e3.d_statistic(gr2, "R1"), d, abs_tol=1e-12)


# ----------------------------------------------------------------------------- matched coverage
def test_select_top_equals_coverage_select():
    amend = {"matched_coverage": {"tie_break": {"seed": 29, "n_cells": 2000}}}
    rng = np.random.default_rng(3)
    ids = rng.choice(2000, 500, replace=False)
    conf = np.round(rng.uniform(size=500), 1)  # many ties
    conf[:40] = np.nan
    rank = va.tie_break_rank(ids, amend)
    for c in (0.95, 0.9, 0.8, 0.5):
        assert np.array_equal(e3.select_top(conf, rank, c), va.coverage_select(conf, ids, c, amend))
    # repeated ids (a bootstrap replicate) keep their rank
    rep = np.concatenate([ids[:100], ids[:100]])
    assert np.array_equal(e3.select_top(conf[np.r_[0:100, 0:100]], va.tie_break_rank(rep, amend), 0.7),
                          va.coverage_select(conf[np.r_[0:100, 0:100]], rep, 0.7, amend))


def test_selective_accuracy_at():
    calls = np.array(["NK", "T", "T", "B", "NK"])
    truth = np.array(["NK", "T", "NK", "T", "NK"])
    conf = np.array([0.9, 0.8, 0.7, 0.6, 0.1])
    rank = np.arange(5)
    assert e3.selective_accuracy_at(calls, truth, conf, rank, 0.4) == 1.0  # top 2
    assert math.isclose(e3.selective_accuracy_at(calls, truth, conf, rank, 0.8), 0.5)  # top 4: NK T (T!=NK) (B!=T)


# ----------------------------------------------------------------------------- scores
def test_rule_scores_and_confidences():
    names = ["a1", "a2", "b1", "b2", "c1", "c2", "d1", "d2"]
    panel = {"B": ["a1", "a2"], "T": ["b1", "b2"], "NK": ["c1", "c2"], "myeloid": ["d1", "d2"]}
    ev = np.array([[1, 0, .5, .5, 0, 0, 0, 0], [.25, .25, .25, .25, .25, .25, .25, .25]], dtype=float)
    S = e3.class_scores(ev, names, panel)
    assert np.allclose(S, [[.5, .5, 0, 0], [.25, .25, .25, .25]])
    assert list(e3.argmax_calls(S)) == ["B", "B"]  # ties -> first lineage (numpy argmax)
    assert np.allclose(e3.top_score(S), [.5, .25])
    assert np.allclose(e3.margin_score(S), [0, 0])
    ent = e3.entropy_score(np.array([[1, 0, 0, 0], [1, 1, 1, 1], [0, 0, 0, 0]], dtype=float))
    assert np.allclose(ent, [1.0, 0.0, 0.0])


def test_knn_label_disagreement():
    zr = np.array([[1, 0], [1, 0.01], [0, 1], [0.01, 1]], dtype=np.float32)
    lab = np.array(["NK", "NK", "T", "OUT"])
    ze = np.array([[1, 0.001], [0.001, 1]], dtype=np.float32)
    d = e3.knn_label_disagreement(ze, zr, lab, np.array(["NK", "T"]), k=2)
    assert np.allclose(d, [0.0, 0.5])


# ----------------------------------------------------------------------------- bootstrap and verdicts
def test_two_stage_indices_resample_within_donors():
    groups = np.array(["a"] * 5 + ["b"] * 3)
    reps = list(e3.two_stage_indices(groups, ["a", "b"], 50, 1))
    reps2 = list(e3.two_stage_indices(groups, ["a", "b"], 50, 1))
    assert all(np.array_equal(x, y) for x, y in zip(reps, reps2))
    for ix in reps:
        # each drawn donor contributes its own size, from its own units only
        assert ix.size in (6, 8, 10)
        na = np.sum(groups[ix] == "a")
        assert na in (0, 5, 10)


def test_verdict_rules():
    v = e3.verdict
    assert v(0.06, 0.01, 0.1, {"x": 0.04, "y": 0.03}, 0.05) == "win"
    assert v(0.06, 0.01, 0.1, {"x": 0.04, "y": 0.02}, 0.05) == "inconclusive"  # a donor below margin/2
    assert v(0.06, -0.01, 0.1, {"x": 0.06, "y": 0.06}, 0.05) == "inconclusive"  # lower bound <= 0
    assert v(-0.06, -0.1, -0.01, {"x": -0.03, "y": -0.04}, 0.05) == "loss"
    assert v(0.0, -0.04, 0.04, {"x": 0.0, "y": 0.1}, 0.05) == "equivalent"
    assert v(float("nan"), None, None, {}, 0.05) == "inconclusive"


def test_boot_stat_paired_draws():
    groups = np.array(["a"] * 30 + ["b"] * 20)
    x = np.arange(50.0)
    out = e3.boot_stat(50, groups, ["a", "b"], 20, 1, lambda ix: {"m": x[ix].mean(), "m2": 2 * x[ix].mean()})
    assert np.allclose(out["m2"], 2 * out["m"])


# ----------------------------------------------------------------------------- readouts
torch = pytest.importorskip("torch")


def test_r2_zero_query_is_gene_mean_and_padding_invariant():
    from lib.v3_e3_readouts import R2AttnPool, pad_states
    rng = np.random.default_rng(4)
    chunks = [rng.normal(size=(n, 6)).astype(np.float16) for n in (3, 7, 1)]
    H, M = pad_states(chunks)
    m = R2AttnPool(6, 2)
    mu, sd = rng.normal(size=6), rng.uniform(0.5, 2, size=6)
    m.set_standardisation(mu, sd)
    with torch.no_grad():
        y, a = m(H, M)
        for i, c in enumerate(chunks):
            gm = ((c.astype(np.float32) - mu) / sd).mean(0)
            ref = m.out(torch.from_numpy(gm.astype(np.float32)))
            assert torch.allclose(y[i], ref, atol=1e-5)
            assert torch.allclose(a[i, :c.shape[0]], torch.full((c.shape[0],), 1.0 / c.shape[0]), atol=1e-6)
            assert float(a[i, c.shape[0]:].abs().sum()) == 0.0
        # extra padding does not change the output
        H2 = torch.cat([H, torch.full((3, 4, 6), 7.0)], dim=1)
        M2 = torch.cat([M, torch.zeros((3, 4), dtype=torch.bool)], dim=1)
        m.q.data = torch.from_numpy(rng.normal(size=6).astype(np.float32))
        y1, _ = m(H, M)
        y2, _ = m(H2, M2)
        assert torch.allclose(y1, y2, atol=1e-5)


def test_r1_standardisation_and_early_stopping_restores_best():
    from lib.v3_e3_readouts import R1MLP, fit_with_early_stopping
    torch.manual_seed(0)
    X = torch.randn(64, 5)
    Y = X[:, :2] * 0.5
    m = R1MLP(5, 2, hidden=8, dropout=0.0)
    m.set_standardisation(np.zeros(5), np.ones(5))
    vals = iter([0.5, 0.2, 0.3, 0.4, 0.6])
    snaps = []

    def val_mse(mod):
        snaps.append({k: v.clone() for k, v in mod.state_dict().items()})
        return next(vals)

    def batches(ep):
        yield X, Y

    fit = fit_with_early_stopping(m, batches, lambda mod, b: torch.mean((mod(b[0]) - b[1]) ** 2), val_mse,
                                  lr=1e-2, max_epochs=5, patience=2, log=lambda s: None)
    assert fit["best_epoch"] == 2 and fit["epochs_run"] == 4
    for k, v in m.state_dict().items():
        assert torch.equal(v, snaps[1][k])


# ----------------------------------------------------------------------------- script / data
def _load_script():
    spec = importlib.util.spec_from_file_location("v3_e3_nkt_repair", REPO / "scripts" / "v3_e3_nkt_repair.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_declared_targets_cover_panel_and_gap_proteins():
    s = _load_script()
    if not REG.exists():
        pytest.skip("registration not built")
    reg = va.load_registration_amended()
    pan = s.panel_in_targets(reg)
    assert sum(len(v) for v in pan.values()) == 12
    assert all(p in s.DECLARED["targets"] for p in e3.GAP_PROTEINS)
    assert len(s.DECLARED["targets"]) == 13


def test_addendum_matches_script_when_present():
    f = REPO / "registration" / "addenda" / "E3.json"
    if not f.exists():
        pytest.skip("E3 addendum not written")
    s = _load_script()
    add = json.loads(f.read_text())
    assert add["declared"] == s.DECLARED
    assert add["computed_train_val"]["val_pairs"]["reproduces_registration"] is True
    assert add["leakage_check"]["site4_rows_poisoned_numbers_identical"] is True
    assert add["leakage_check"]["val_rows_poisoned_numbers_change"] is True
    hashes = (REPO / "registration" / "addenda" / "HASHES.txt").read_text().splitlines()
    assert f"{va.vk.sha256_file(f)}  E3.json" in hashes


def test_val_pairs_reproduce_registration():
    if not (DATA / "cite" / "cite_arrays.npz").exists() or not REG.exists():
        pytest.skip("data or registration absent")
    from lib import v3_key as vk
    reg = va.load_registration_amended()
    npz = np.load(DATA / "cite" / "cite_arrays.npz", allow_pickle=False)
    split = npz["split"].astype(str)
    vr = np.where(split == "val")[0]
    keys = vk.build_keys(npz["cell_types"].astype(str)[vr], npz["adt"][vr], [str(x) for x in npz["adt_names"]], reg)
    full = np.full(split.size, "unscored", dtype=object)
    full[vr] = keys["primary"]
    z = np.load(DATA / "cite_official" / "z_rna.npy", mmap_mode="r")
    P = e3.nkt_pairs_within_donors(z, vr, npz["donors"].astype(str), {"all": full.astype(str)}, 10, reg["e3"]["flag_cosine"])
    pv = P["variants"]["all"]
    assert pv["nk"].size == reg["e3"]["n_val_nkt_pairs"] == 171
    assert round(float(np.median(pv["cos"])), 6) == reg["e3"]["flag_cosine"]


def test_registered_d_rewards_uniform_shrinkage_but_dlog_does_not():
    """Documents why the addendum adds Dlog (secondary_H3a): shrinking every difference by s < 1 gives D > 0."""
    head = {"flagged": 0.4, "unflagged": 0.8}
    s = 0.5
    gr = {"head": head, "null": dict(head), "R1": {k: s * v for k, v in head.items()}}
    assert math.isclose(e3.d_statistic(gr, "R1"), (1 - s) * (0.8 - 0.4))
    assert abs(e3.d_log_statistic(gr, "R1")) < 1e-12
    gr["R1"]["flagged"] = 0.6  # more of the gap on flagged pairs only
    assert e3.d_log_statistic(gr, "R1") > 0
    gr["R1"]["flagged"] = 0.0
    assert math.isnan(e3.d_log_statistic(gr, "R1"))
