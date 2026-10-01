"""Tests for E4 (two-channel fusion): bridge_anm/lib/v3_e4.py and the runner's helpers.

Synthetic tests always run; tests that need the ANM engine skip without ANM_ROOT, and tests that need
the registration skip when it is absent.
"""
from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bridge_anm"))
sys.path.insert(0, str(REPO))

from lib import v3_amend as va  # noqa: E402
from lib import v3_e4 as e4  # noqa: E402
from lib import v3_key as vk  # noqa: E402

REG_PATH = REPO / "registration" / "registration_v3.json"
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))
FIELD = {"kind": "finite_graph_scalar", "steps": 4, "retention": 0.82, "diffusion": 0.16, "source_scale": 1.0}
AMEND = {"matched_coverage": {"tie_break": {"seed": 29, "n_cells": 90261}}}


def _reg_stub():
    """Minimal registration: the primary panel and the registered ANM field."""
    pan = {"B": ["CD20", "CD22", "CD268"], "T": ["CD3", "CD2", "CD5"], "NK": ["CD122", "CD94", "CD56"],
           "myeloid": ["CD172a", "CD11c", "CD62P"]}
    return {"panels": {"primary": {"proteins": pan}}, "anm": {"field_representation": dict(FIELD)}}


@pytest.fixture(scope="module")
def reg():
    if REG_PATH.exists():
        return va.load_registration_amended()
    return _reg_stub()


def _anm():
    try:
        return e4.import_anm(ANM_ROOT)
    except Exception:
        pytest.skip("ANM engine not importable (set ANM_ROOT)")


def _schema(reg):
    base = json.loads((REPO / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
    return e4.anm_schema(base, reg, 0.3)


# ----------------------------------------------------------------------------- closed form vs engine
def test_star_unit_response_is_the_registered_field_gain():
    for n in (1, 3, 6, 8):
        assert e4.star_unit_response(n, FIELD["steps"] + 1, FIELD) == pytest.approx(vk.field_gain(n, FIELD), rel=1e-12)
    assert e4.star_unit_response(6, 5, FIELD, propagate=False) == 0.0


@pytest.mark.parametrize("t2,contra,prop", [(0, True, True), (2, True, True), (0, False, True), (0, True, False), (2, False, True)])
def test_closed_form_equals_anm_engine(reg, t2, contra, prop):
    ffr = _anm()
    schema = _schema(reg)
    rng = np.random.default_rng(1)
    n = 60
    e1 = np.clip(rng.uniform(-0.3, 1.3, (n, 12)), 0, 1)
    e2 = np.clip(rng.uniform(-0.3, 1.3, (n, 12)), 0, 1)
    e1[:5] = 1.0  # ties at the clip: top class taken first in LINEAGES order
    e2[5:8] = 0.0
    C = e4.anm_fusion_closed(e1, e2, 0.83, 1.0, reg, t2=t2, contradiction=contra, propagate=prop)
    for i in range(n):
        ev = e4.fusion_events(e1[i], e2[i], 0.83, 1.0, reg, t2=t2, contradiction=contra)
        A, rec, nrej = e4.anm_fusion_engine_cell(ffr, schema, ev, propagate=prop)
        assert nrej == 0
        np.testing.assert_allclose(A, C[i], rtol=0, atol=1e-12)
        best = int(np.argmax(A))
        assert rec == (vk.LINEAGES[best] if A[best] >= 0.3 else None)


def test_events_count_values_and_contradiction_target(reg):
    e1 = np.linspace(0, 1, 12)
    e2 = np.zeros(12)
    e2[:3] = 0.9  # channel 2 says B
    ev = e4.fusion_events(e1, e2, 0.5, 1.0, reg, t2=2, contradiction=True)
    assert len(ev) == 24 + 6
    assert all(0.0 <= x["value"] <= 1.0 for x in ev)
    con = [x for x in ev if x["polarity"] == "contradiction"]
    ch2 = [x for x in con if x["modality"] == e4.CHANNELS[1]]
    assert sorted(x["action"] for x in ch2) == ["NK", "T", "myeloid"]
    assert all(x["time"] == 2 for x in ch2)
    assert ch2[0]["value"] == pytest.approx(1.0 * 0.5 * 0.9)
    ch1 = [x for x in con if x["modality"] == e4.CHANNELS[0]]
    assert sorted(x["action"] for x in ch1) == ["B", "NK", "T"]  # channel 1's top class is myeloid
    assert ch1[0]["value"] == pytest.approx(0.5 * 0.5 * e1[9:12].mean())


def test_no_contradiction_is_proportional_to_trust_weighted_average(reg):
    rng = np.random.default_rng(2)
    e1, e2 = rng.uniform(0, 1, (200, 12)), rng.uniform(0, 1, (200, 12))
    t1, t2 = 0.71, 0.86
    s = e4.source_scales({e4.CHANNELS[0]: t1, e4.CHANNELS[1]: t2})
    A = e4.anm_fusion_closed(e1, e2, s[e4.CHANNELS[0]], s[e4.CHANNELS[1]], reg, contradiction=False)
    F2 = e4.fuse_trust_weighted(e4.class_scores12(e1, reg), e4.class_scores12(e2, reg), t1, t2)
    ratio = A / F2
    np.testing.assert_allclose(ratio, ratio.flat[0], rtol=1e-12)  # exact source-scale ratio: proportional to float rounding
    G6 = vk.field_gain(6, FIELD)
    assert ratio.flat[0] == pytest.approx(G6 * 3 * (t1 + t2) / max(t1, t2), rel=1e-12)


def test_c3_and_f2_select_and_call_the_same_cells(reg):
    """C3 (F5 without contradiction events) is a constant multiple of F2, so with the A1.1 tie-break both select the
    same cells at every coverage and make the same calls (needs the exact source-scale ratio)."""
    rng = np.random.default_rng(12)
    n = 3000
    e1 = np.clip(rng.uniform(-0.2, 1.3, (n, 12)), 0, 1)
    e2 = np.clip(rng.uniform(-0.2, 1.3, (n, 12)), 0, 1)
    e1[:400, 3:6] = 1.0  # many cells tied at the clip
    ids = rng.choice(90261, size=n, replace=False)
    for t1, t2 in ((0.828938, 0.813822), (0.754129, 0.452151), (0.611111, 0.777777)):
        s = e4.source_scales({e4.CHANNELS[0]: t1, e4.CHANNELS[1]: t2})
        C3 = e4.anm_fusion_closed(e1, e2, s[e4.CHANNELS[0]], s[e4.CHANNELS[1]], reg, contradiction=False)
        F2 = e4.fuse_trust_weighted(e4.class_scores12(e1, reg), e4.class_scores12(e2, reg), t1, t2)
        assert np.array_equal(e4.calls_from_scores(C3), e4.calls_from_scores(F2))
        for cov in (0.5, 0.8, 0.9):
            assert np.array_equal(e4.select_at_coverage(C3.max(1), ids, cov, AMEND),
                                  e4.select_at_coverage(F2.max(1), ids, cov, AMEND))


# ----------------------------------------------------------------------------- matched coverage and bootstrap
def test_weighted_selective_accuracy_matches_explicit_coverage_select():
    rng = np.random.default_rng(3)
    n = 120
    ids = rng.choice(90261, size=n, replace=False)
    conf = np.round(rng.uniform(0, 1, n), 1)  # many ties
    correct = rng.uniform(size=n) < 0.7
    scored = rng.uniform(size=n) < 0.9
    order = e4.ranking_order(conf, ids, AMEND)
    W = np.stack([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(30)])
    for cov in (0.5, 0.8, 0.9):
        got = e4.weighted_selective_accuracy(order, W, correct & scored, scored, cov)
        for b in range(W.shape[0]):
            rep = np.repeat(np.arange(n), W[b])
            sel = va.coverage_select(conf[rep], ids[rep], cov, AMEND)
            ns = int(np.sum(sel & scored[rep]))
            want = np.sum(sel & correct[rep] & scored[rep]) / ns
            assert got[b] == pytest.approx(want, abs=1e-12)
    ones = np.ones((1, n), dtype=np.int64)
    sel = va.coverage_select(conf, ids, 0.8, AMEND)
    assert e4.weighted_selective_accuracy(order, ones, correct & scored, scored, 0.8)[0] == pytest.approx(
        np.sum(sel & correct & scored) / np.sum(sel & scored))


def test_select_at_coverage_rounds_float_noise_into_ties():
    ids = np.array([5, 9, 11, 40])
    conf = np.array([0.5, 0.5 + 1e-15, 0.2, 0.9])
    a = e4.select_at_coverage(conf, ids, 0.5, AMEND)
    b = e4.select_at_coverage(np.array([0.5, 0.5, 0.2, 0.9]), ids, 0.5, AMEND)
    assert np.array_equal(a, b)


def test_two_stage_weights():
    donor = np.array(["a"] * 7 + ["b"] * 3)
    W1 = e4.two_stage_weights(donor, 50, 1)
    W2 = e4.two_stage_weights(donor, 50, 1)
    assert np.array_equal(W1, W2)
    tot = W1.sum(axis=1)
    assert set(tot.tolist()) <= {6, 10, 14}  # (b, b), (a, b) or (a, a)
    for b in range(50):
        na, nb = W1[b, :7].sum(), W1[b, 7:].sum()
        assert na in (0, 7, 14) and nb in (0, 3, 6)


def test_decide():
    m = 0.01
    assert e4.decide(0.02, 0.001, 0.04, [0.02, 0.006], m) == "win"
    assert e4.decide(0.02, 0.001, 0.04, [0.02, 0.004], m) == "inconclusive"  # one donor below margin/2
    assert e4.decide(-0.02, -0.04, -0.001, [-0.03, -0.006], m) == "loss"
    assert e4.decide(0.001, -0.005, 0.008, [0.0, 0.002], m) == "equivalent"
    assert e4.decide(0.005, -0.012, 0.02, [0.0, 0.01], m) == "inconclusive"


# ----------------------------------------------------------------------------- noise, key, channel 2, trust
def test_dropout_is_per_cell_and_order_free():
    rng = np.random.default_rng(4)
    w = rng.uniform(0.1, 2.0, (40, 122))
    ids = rng.choice(90261, size=40, replace=False)
    d = e4.dropout_inputs(w, ids, 0.5, 17)
    assert np.all((d == 0).sum(axis=1) == 61)
    perm = rng.permutation(40)
    d2 = e4.dropout_inputs(w[perm], ids[perm], 0.5, 17)
    assert np.array_equal(d2, d[perm])
    assert not np.array_equal(e4.dropout_inputs(w, ids, 0.5, 18), d)


def test_count_recovery_and_thinning():
    counts = np.array([1, 3, 7, 2, 40])
    stored = counts / 3.7  # stored rows are counts over one per-cell factor
    rc, dev, ok = e4.recover_counts(stored.astype(np.float32))
    assert ok and dev < 1e-3 and np.array_equal(rc, counts)
    t1 = e4.thin_counts(counts, 0.2, 17, 123)
    assert np.array_equal(t1, e4.thin_counts(counts, 0.2, 17, 123))
    assert np.all(t1 <= counts) and np.all(t1 >= 0)
    big = np.full(20000, 5)
    assert abs(e4.thin_counts(big, 0.2, 17, 1).mean() / 5 - 0.2) < 0.01


def test_e4_key_and_tau(reg):
    m = np.zeros((6, 12))
    m[0, 0:3] = 0.9          # B
    m[1, 3:6] = 0.3          # T, weak
    m[2, 6:9] = 1.0          # NK
    m[2, 9:12] = 1.0         # tie NK / myeloid -> NK (first in order)
    m[3, 9:12] = 0.6         # myeloid
    key = e4.e4_key(m, reg, 0.5)
    assert key.tolist() == ["B", "OUT", "NK", "myeloid", "OUT", "OUT"]
    annot = np.array(["B", "T", "NK", "myeloid", "OUT", "OUT"])
    tau, rows = e4.choose_tau(m, annot, reg, grid=(0.2, 0.25, 0.5))
    assert tau == 0.2  # 0.2 and 0.25 tie at kappa 1 -> the smaller
    assert [r["tau"] for r in rows] == [0.2, 0.25, 0.5]


def test_measured_evidence_and_normalised_evidence():
    adt = np.array([[0.0, 2.0, 4.0], [1.0, 1.0, 9.0]])
    v = e4.measured_evidence(adt, ["a", "b", "c"], ["c", "a"], [4.0, 0.5])
    np.testing.assert_allclose(v, [[1.0, 0.0], [1.0, 1.0]])
    np.testing.assert_allclose(e4.normalised_evidence(np.array([[-1.0, 0.5]]), [1.0, 1.0]), [[0.0, 0.5]])


def test_trust_and_source_scales():
    rng = np.random.default_rng(5)
    meas = rng.uniform(size=(500, 3))
    t, r = e4.channel_trust(np.c_[meas[:, 0], -meas[:, 1], np.full(500, 0.3)], meas)
    assert r[0] == pytest.approx(1.0) and r[1] == pytest.approx(-1.0) and r[2] == 0.0
    assert t == pytest.approx(1.0 / 3)
    s = e4.source_scales({e4.CHANNELS[0]: 0.6, e4.CHANNELS[1]: 0.8})
    assert s == {e4.CHANNELS[0]: 0.6 / 0.8, e4.CHANNELS[1]: 1.0}  # the exact ratio, not rounded
    assert s[e4.CHANNELS[0]] == pytest.approx(0.75)
    assert e4.source_scales({e4.CHANNELS[0]: 0.0, e4.CHANNELS[1]: 0.0}) == {e4.CHANNELS[0]: 0.0, e4.CHANNELS[1]: 0.0}


def test_ridge_alpha_selection_prefers_the_smaller_on_ties():
    rng = np.random.default_rng(6)
    X = rng.normal(size=(300, 5))
    Y = X @ rng.normal(size=(5, 2)) + 0.01 * rng.normal(size=(300, 2))
    alpha, rows, model = e4.ridge_fit_select(X[:200], Y[:200], X[200:], Y[200:], alphas=(0.1, 1.0, 1000.0))
    assert alpha in (0.1, 1.0) and len(rows) == 3
    np.testing.assert_allclose(e4.linear_predict(X[200:], model.coef_, model.intercept_), model.predict(X[200:]))


def _clr_float32(x):
    """Stored ADT form: per-cell CLR log1p(x / g_c), g_c from every protein (panel included), stored as float32."""
    g = np.exp(np.log1p(x).mean(axis=1, keepdims=True))
    return np.log1p(x / g).astype(np.float32).astype(np.float64)


def _toy_channel2(reg, n=400, seed=0):
    rng = np.random.default_rng(seed)
    flat = e4.panel_flat(reg)
    names = flat[:6] + [f"np{i}" for i in range(40)] + flat[6:]
    x = rng.poisson(rng.gamma(2.0, 6.0, size=(n, len(names)))).astype(np.float64)
    x[:, 7] += 1
    adt = _clr_float32(x)
    w, _ = va.channel2_inputs(adt, names, flat)
    m = e4.measured_evidence(adt, names, flat, np.percentile(adt[:, [names.index(p) for p in flat]], 95, axis=0))
    _, _, ridge = e4.ridge_fit_select(w[: n // 2], m[: n // 2], w[n // 2:], m[n // 2:])
    q95 = np.percentile(e4.linear_predict(w, ridge.coef_, ridge.intercept_), 95, axis=0)
    return names, flat, x, adt, ridge.coef_, ridge.intercept_, q95, rng.choice(90261, size=n, replace=False)


def test_channel2_evidence_is_invariant_to_panel_counts_end_to_end(reg):
    """A1.2 end to end: inputs -> (L2 dropout) -> ridge -> normaliser. Changing only the panel counts changes the
    stored non-panel values (the CLR factor includes the panel) but not the channel-2 evidence (float32 residue)."""
    names, flat, x, adt, coef, icpt, q95, ids = _toy_channel2(reg)
    pj = [names.index(p) for p in flat]
    x2 = x.copy()
    x2[:, pj] = np.random.default_rng(1).poisson(80.0, size=(x.shape[0], len(pj)))
    adt2 = _clr_float32(x2)
    nj = [i for i in range(len(names)) if i not in pj]
    assert np.max(np.abs(adt2[:, nj] - adt[:, nj])) > 0.1  # the stored non-panel values do move
    for kw in ({}, {"cell_ids": ids, "dropout_fraction": 0.5, "seed": 17}):
        a = e4.channel2_evidence(adt, names, flat, coef, icpt, q95, **kw)
        b = e4.channel2_evidence(adt2, names, flat, coef, icpt, q95, **kw)
        assert np.max(np.abs(a - b)) < 1e-5
        assert np.array_equal(e4.top_class(e4.class_scores12(a, reg)), e4.top_class(e4.class_scores12(b, reg)))


def test_channel2_panel_invariance_check_passes_and_detects_a_leak(reg):
    names, flat, x, adt, coef, icpt, q95, ids = _toy_channel2(reg)
    inv = e4.channel2_panel_invariance(adt, names, flat, coef, icpt, q95, reg, ids, 0.5, 17)
    assert inv["passed"] and inv["L0"]["panel_overwritten_bit_identical"] and inv["L2"]["panel_overwritten_bit_identical"]
    assert inv["stored_nonpanel_max_abs_change"] > 0.1
    # positive control: a path that lets one panel protein through as an input must fail the check
    leaky_panel = flat[1:]  # flat[0] would then be read as a channel-2 input
    w, cols = va.channel2_inputs(adt, names, leaky_panel)
    assert flat[0] in cols
    from sklearn.linear_model import Ridge

    m = e4.measured_evidence(adt, names, flat, np.percentile(adt[:, [names.index(p) for p in flat]], 95, axis=0))
    rl = Ridge(alpha=0.1).fit(w, m)
    ql = np.percentile(e4.linear_predict(w, rl.coef_, rl.intercept_), 95, axis=0)
    bad = e4.channel2_panel_invariance(adt, names, flat, rl.coef_, rl.intercept_, ql, reg, ids, 0.5, 17,
                                       input_exclude=leaky_panel)
    assert not bad["passed"] and not bad["L0"]["panel_overwritten_bit_identical"]


def test_prepare_refuses_on_engine_mismatch_or_a1_2_failure():
    run = _runner()
    ok_cell = {"max_abs_engine_minus_closed_lt_1e-9": True, "n_cells_argmax_differs": 0, "n_rejected_events": 0}
    comp = {"anm_engine_check_val": {L: {v: dict(ok_cell) for v in run.ANM_VARIANTS} for L in e4.LEVELS},
            "channel2": {"a1_2_val_panel_invariance": {"passed": True}}}
    assert run.prepare_refusals(comp) == []
    comp["anm_engine_check_val"]["L2"]["F5"]["n_cells_argmax_differs"] = 1
    assert len(run.prepare_refusals(comp)) == 1
    comp["anm_engine_check_val"]["L2"]["F5"]["n_cells_argmax_differs"] = 0
    comp["channel2"]["a1_2_val_panel_invariance"]["passed"] = False
    assert len(run.prepare_refusals(comp)) == 1
    del comp["channel2"]["a1_2_val_panel_invariance"]
    assert len(run.prepare_refusals(comp)) == 1


def test_names_say_declared_rule_and_closed_form():
    run = _runner()
    assert "declared rule" in run.ARM_NAMES["F1"] and "declared rule" in run.ARM_NAMES["F2"]
    assert "ANM" not in run.ARM_NAMES["F1"] and "ANM" not in run.ARM_NAMES["F2"]
    assert "closed form" in run.ARM_NAMES["F5"] and "contradiction" in run.ARM_NAMES["F5"]
    assert "by construction" not in run.ARM_NAMES["C3"]
    txt = run.ADDENDUM_TEXT["open_choices_fixed"]
    assert "not rounded" in txt["trust"] and "F5_closed_form" in txt and "interpretation" in txt


# ----------------------------------------------------------------------------- registration-level checks
def test_registered_channel2_inputs_are_the_122_non_panel_proteins(reg):
    if not REG_PATH.exists():
        pytest.skip("registration not built")
    names = list(reg["evidence"]["teddy_head"]["adt_names"])
    x = np.zeros((2, len(names)))
    _, cols = va.channel2_inputs(x, names, e4.panel_flat(reg))
    assert len(cols) == 122
    assert set(cols) == set(reg["experiments"]["E4"]["channels"].get("channel2_input_proteins", cols))
    assert not set(cols) & set(e4.panel_flat(reg))


def _runner():
    spec = importlib.util.spec_from_file_location("v3_e4_fusion", REPO / "bridge_anm" / "v3_e4_fusion.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_update_hash_lines_keeps_other_builders_lines(tmp_path):
    run = _runner()
    hf = tmp_path / "HASHES.txt"
    hf.write_text("aaa  E5.json\nold  E4.json\n")
    run.update_hash_lines(tmp_path, {"E4.json": "new", "E4_models.json": "mm"})
    assert hf.read_text().splitlines() == ["aaa  E5.json", "new  E4.json", "mm  E4_models.json"]


def test_lineage_confidence_ignores_out_probability():
    run = _runner()
    P = np.array([[0.1, 0.05, 0.7, 0.1, 0.05], [0.2, 0.5, 0.1, 0.1, 0.1]])
    conf, call = run.lineage_conf_calls(P, ["B", "NK", "OUT", "T", "myeloid"])
    np.testing.assert_allclose(conf, [0.1, 0.5])
    assert call.tolist() == ["B", "NK"]
    p2 = run.softmax_proba(np.array([[1.0, 2.0]]), np.array([[1.0, 0.0], [0.0, 1.0]]), np.array([0.0, 0.0]))
    assert p2[0, 1] == pytest.approx(math.exp(2) / (math.exp(1) + math.exp(2)))


def test_commit_addendum_commits_only_e4_files(tmp_path):
    import subprocess
    import types

    run = _runner()
    repo = tmp_path / "repo"
    add_dir = repo / "registration" / "addenda"
    add_dir.mkdir(parents=True)

    def g(*args):
        return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, check=True, text=True).stdout

    g("init", "-q", "-b", "exp")
    g("config", "user.email", "t@example.com")
    g("config", "user.name", "t")
    (add_dir / "HASHES.txt").write_text("e5hash  E5.json\n")
    g("add", ".")
    g("commit", "-q", "-m", "base")
    # another builder: an uncommitted HASHES line and a staged file
    (add_dir / "HASHES.txt").write_text("e5hash  E5.json\ne2hash  E2.json\n")
    (repo / "other.txt").write_text("other builder\n")
    g("add", "other.txt")
    models = b'{"m": 1}\n'
    core = "c" * 64
    (add_dir / "E4_models.json").write_bytes(models)
    (add_dir / "E4.json").write_text(json.dumps({"core_sha256": core, "models_sha256": run.sha_bytes(models),
                                                 "computed": {"key": {"tau_K": 0.5}, "channel2": {"alpha": 1.0},
                                                              "best_comparator": "F2", "trust_untestable": False},
                                                 "smoke": None}))
    (add_dir / "E4_leakage_check.json").write_text(json.dumps({"passed": True, "addendum_core_sha256": core}))
    old_root = run.ROOT
    try:
        run.ROOT = repo
        a = types.SimpleNamespace(registration_dir=repo / "registration", smoke=None, out_dir=tmp_path / "out")
        run.stage_commit_addendum(a)
    finally:
        run.ROOT = old_root
    changed = sorted(g("show", "--name-only", "--format=", "HEAD").split())
    assert changed == ["registration/addenda/E4.json", "registration/addenda/E4_leakage_check.json",
                       "registration/addenda/E4_models.json", "registration/addenda/HASHES.txt"]
    committed_hashes = g("show", "HEAD:registration/addenda/HASHES.txt").splitlines()
    assert committed_hashes[0] == "e5hash  E5.json"
    assert not any("E2.json" in ln for ln in committed_hashes)
    assert {ln.split("  ")[1] for ln in committed_hashes[1:]} == {"E4.json", "E4_models.json", "E4_leakage_check.json"}
    assert "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" in g("log", "-1", "--format=%B")
    status = g("status", "--porcelain")
    assert "A  other.txt" in status  # the other builder's staged file is still staged, not committed
    assert "e2hash  E2.json" in (add_dir / "HASHES.txt").read_text()  # its uncommitted line is untouched
    assert "E4.json" not in status
