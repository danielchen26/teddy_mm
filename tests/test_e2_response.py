"""Tests for E2 (C3 response decomposition, C5 linearity): bridge_anm/lib/e2_response.py and the
pure parts of scripts/e2_response_decomposition.py.

Synthetic tests always run. The ANM-engine test skips without ANM_ROOT; the row-discipline test of the
register stage skips without the data pack (data/processed/cite)."""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from bridge_anm.lib import e2_response as er  # noqa: E402
from bridge_anm.lib import v3_amend as va  # noqa: E402
from bridge_anm.lib import v3_key as vk  # noqa: E402

ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))
DATA = REPO / "data/processed/cite/cite_arrays.npz"


def _script():
    spec = importlib.util.spec_from_file_location("e2_script", REPO / "scripts/e2_response_decomposition.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def reg():
    return va.load_registration_amended()


# ----------------------------------------------------------------------------- counts and perturbations
def test_recover_counts_from_scaled_row():
    rng = np.random.default_rng(0)
    counts = rng.poisson(0.8, size=500)
    counts[3] = 1  # make sure a count of 1 exists
    for factor in (0.37, 1.0, 2.913):
        x = (counts * factor).astype(np.float32)
        c, unit, dev = er.recover_counts(x)
        assert np.array_equal(c, counts) and dev <= 1.0 and unit == pytest.approx(factor, rel=1e-6)
    big = np.array([1, 2, 7, 150000, 999999])  # float32 rounding of large counts stays inside the tolerance
    for factor in (0.1937, 0.7, 1.33):
        c, unit, dev = er.recover_counts((big * factor).astype(np.float32))
        assert np.array_equal(c, big) and dev <= 1.0
    c, unit, dev = er.recover_counts(np.zeros(10, dtype=np.float32))
    assert c.sum() == 0 and unit == 0.0 and dev == 0.0
    _, _, dev = er.recover_counts(np.array([1.0, 1.5, 0.0]))
    assert dev > 1.0  # not a scaled count vector


def test_thinning_is_binomial_and_seeded():
    c = np.array([0, 1, 5, 100, 1000])
    assert np.array_equal(er.thin_counts(c, 1.0, er.thinning_rng(11, 7, 0)), c)
    a = er.thin_counts(c, 0.5, er.thinning_rng(11, 7, 3))
    b = er.thin_counts(c, 0.5, er.thinning_rng(11, 7, 3))
    d = er.thin_counts(c, 0.5, er.thinning_rng(12, 7, 3))
    assert np.array_equal(a, b) and np.all(a <= c) and a[0] == 0
    assert not np.array_equal(a, d) or c.sum() < 5
    big = np.full(20000, 10)
    t = er.thin_counts(big, 0.2, er.thinning_rng(11, 1, 4))
    assert t.mean() == pytest.approx(2.0, abs=0.05)


def test_gene_scaling_moves_or_removes_the_token():
    from teddy_mm.teddy_encoder import rank_encode_official

    vals = np.array([[5.0, 4.0, 3.0, 2.0, 1.0, 0.0]], dtype=np.float32)
    tok = np.arange(10, 16)
    base, _ = rank_encode_official(vals, tok, max_len=6, pad_id=0, pad_to=6)
    masked, _ = rank_encode_official(er.scale_gene(vals, 1, 0.0), tok, max_len=6, pad_id=0, pad_to=6)
    half, _ = rank_encode_official(er.scale_gene(vals, 1, 0.6), tok, max_len=6, pad_id=0, pad_to=6)
    assert list(base[0][:5]) == [10, 11, 12, 13, 14]
    assert list(masked[0][:4]) == [10, 12, 13, 14] and masked[0][4] == 0  # removed, others keep their order
    assert list(half[0][:5]) == [10, 12, 11, 13, 14]  # 4.0 * 0.6 = 2.4: moved below 3.0
    assert vals[0, 1] == 4.0  # input untouched


# ----------------------------------------------------------------------------- scores and linearity
def test_margins_and_extrapolation():
    S = np.array([[0.9, 0.5, 0.1, 0.0], [0.2, 0.6, 0.6, 0.1]])
    np.testing.assert_allclose(er.top_margin(S), [0.4, 0.0])
    np.testing.assert_allclose(er.signed_margin(S, np.array([0, 0])), [0.4, -0.4])
    np.testing.assert_allclose(er.extrapolate(np.array([1.0]), np.array([0.9]), 0.05, 0.8), [1.0 - 1.6])


def test_linearity_classifies_linear_quadratic_and_floor():
    rng = np.random.default_rng(1)
    d = rng.normal(size=(50, 8))
    floor = np.full(50, 1e-3)
    lin = er.linearity(0.2 * d, 0.1 * d, floor)
    assert lin["linear"].all() and np.allclose(lin["rho"], 1.0) and np.allclose(lin["cos"], 1.0)
    quad = er.linearity(0.04 * d, 0.01 * d, floor)
    assert not quad["linear"].any() and np.allclose(quad["rho"], 2.0)
    tiny = er.linearity(2e-4 * d / np.linalg.norm(d, axis=1, keepdims=True),
                        1e-4 * d / np.linalg.norm(d, axis=1, keepdims=True), floor)
    assert tiny["below_floor"].all() and not tiny["linear"].any()
    zero = er.linearity(d, np.zeros_like(d), floor)
    assert np.isnan(zero["rho"]).all() and not zero["linear"].any()
    rot = er.linearity(0.2 * d, 0.1 * np.roll(d, 1, axis=1), floor)
    assert np.mean(rot["linear"]) < 0.2  # direction changes -> cosine fails


# ----------------------------------------------------------------------------- decomposition
def test_lost_and_labels():
    key = np.array(["T", "T", "NK", "B", "OUT", "T", "T"])
    base = np.array(["T", "T", "NK", "", "", "NK", "T"])
    pert = np.array(["NK", "", "NK", "B", "T", "T", ""])
    lost = er.lost_cases(base, pert, key)
    assert lost.tolist() == [True, True, False, False, False, False, True]
    probe = np.array(["NK", "T", "NK", "B", "OUT", "T", ""])
    argmax = np.array(["NK", "B", "NK", "B", "T", "T", "T"])
    lab = er.label_lost(probe, argmax, key)
    assert lab[0] == "representation" and lab[1] == "head" and lab[6] == "representation"
    lab2 = er.label_lost(np.array(["T"]), np.array(["T"]), np.array(["T"]))
    assert lab2[0] == "decision"


def test_unit_tables_and_stats():
    # 3 cells; cell 0 has 2 cases, cell 1 one, cell 2 one (unscored)
    cc = np.array([0, 0, 1, 2])
    key = np.array(["T", "T", "OUT", "unscored"])
    cb = np.array(["T", "T", "", "B"])
    cp = np.array(["T", "", "NK", "B"])
    lab = np.array(["decision", "decision", "head", "decision"])
    T = er.unit_cell_tables(cc, 3, key=key, call_base=cb, call_pert=cp, label=lab)
    st = er.unit_stats(T)
    assert st["decision_acc_base"] == pytest.approx(1.0)  # T, T, OUT declined
    assert st["decision_acc_pert"] == pytest.approx(1 / 3)  # T ok; no call on T wrong; NK on OUT wrong
    assert st["decision_loss"] == pytest.approx(2 / 3)
    assert st["n_lost_cases"] == 1 and st["n_lost_cells"] == 1 and st["share_decision"] == pytest.approx(1.0)
    assert st["coverage_pert"] == pytest.approx(3 / 4)
    W = np.array([[2, 1, 0], [0, 1, 1]])
    bs = er.unit_stats(T, W)
    assert bs["decision_acc_pert"][0] == pytest.approx(2 / 5)  # cell 0 x2 (1 of 2 ok each) + cell 1 (wrong)
    assert np.isnan(bs["share_decision"][1])  # no lost case in replicate 2


def test_two_stage_weights_respect_donors():
    donor = np.array(["a"] * 5 + ["b"] * 3)
    W = er.two_stage_weights(donor, 200, 1)
    assert W.shape == (200, 8)
    for row in W:
        wa, wb = row[:5].sum(), row[5:].sum()
        assert (wa, wb) in {(10, 0), (5, 3), (0, 6)}
    W2 = er.two_stage_weights(donor, 200, 1)
    assert np.array_equal(W, W2)


def test_auroc_matches_sklearn_with_ties_and_weights():
    from sklearn.metrics import roc_auc_score

    rng = np.random.default_rng(3)
    s = rng.integers(0, 6, size=300).astype(float)
    y = rng.random(300) < 0.3 + 0.08 * s
    w = rng.integers(0, 3, size=300).astype(float)
    assert er.auroc(s, y) == pytest.approx(roc_auc_score(y, s))
    assert er.auroc(s, y, w) == pytest.approx(roc_auc_score(y, s, sample_weight=w))
    g = er.auroc_groups(s)
    assert er.auroc(s, y, w, g) == pytest.approx(roc_auc_score(y, s, sample_weight=w))
    assert np.isnan(er.auroc(s, np.zeros(300, dtype=bool)))


def _unit(loss, n_lost, shares):
    d = {"decision_loss": loss, "selective_loss": loss, "n_lost_cells": n_lost}
    d.update({f"share_{k}": v for k, v in zip(er.LABELS, shares)})
    return d


def test_criterion_i_needs_matched_loss_counts_margin_and_interval():
    B = 500
    rng = np.random.default_rng(5)
    units = {"u": _unit(0.05, 60, (0.8, 0.1, 0.1)), "w": _unit(0.06, 50, (0.6, 0.3, 0.1)),
             "far": _unit(0.30, 200, (0.0, 0.0, 1.0)), "few": _unit(0.05, 10, (0.0, 1.0, 0.0))}
    boot = {u: {f"share_{k}": units[u][f"share_{k}"] + rng.normal(0, 0.05, B) for k in er.LABELS} for u in units}
    res = er.criterion_i(units, boot)
    assert res["n_pairs_compared"] == 1 and res["holds"]
    rows = {r["label"]: r for r in res["rows"]}
    assert rows["representation"]["passes"] and not rows["decision"]["passes"]  # decision 0.1 vs 0.1: no difference
    noisy = {u: {f"share_{k}": units[u][f"share_{k}"] + rng.normal(0, 0.5, B) for k in er.LABELS} for u in units}
    assert not er.criterion_i(units, noisy)["holds"]  # interval covers 0


def test_win_rule():
    assert er.win_rule(0.05, [0.01, 0.09], {"a": 0.04, "b": 0.03}, 0.02) == "win"
    assert er.win_rule(0.05, [0.01, 0.09], {"a": 0.04, "b": 0.0}, 0.02) == "inconclusive"
    assert er.win_rule(0.0, [-0.01, 0.01], {"a": 0.0, "b": 0.0}, 0.02) == "equivalent"
    assert er.win_rule(-0.05, [-0.09, -0.01], {"a": -0.04, "b": -0.03}, 0.02) == "loss"


def test_flatten_genes_matches_registration(reg):
    G = reg["experiments"]["E2"]["perturbations"]["genes"]
    syms = er.flatten_genes(G["panel_coding_genes"], G["nk_t_genes"])
    assert len(syms) == len(set(syms)) == 21
    assert syms[:4] == ["MS4A1", "CD22", "TNFRSF13C", "CD3E"] and syms[-1] == "CD28"


# ----------------------------------------------------------------------------- ANM engine
def test_anm_engine_equals_closed_form_and_rule(reg):
    sys.path.insert(0, str(ANM_ROOT))
    try:
        from active_neural_matter.field import finite_field_runner as ffr
    except Exception:
        pytest.skip("ANM engine not importable (set ANM_ROOT)")
    base = json.loads((REPO / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
    thr = vk.anm_readout_threshold(reg, "Q1")
    schema = er.anm_schema(base, reg["anm"]["field_representation"], thr)
    rng = np.random.default_rng(9)
    v12 = np.clip(rng.uniform(-0.2, 1.2, size=(300, 12)), 0, 1)
    ok = np.ones(300, dtype=bool)
    ok[5] = False
    A, calls = er.anm_engine(v12, ok, ffr, schema)
    S = v12.reshape(300, 4, 3).mean(axis=2)
    G3 = vk.field_gain(3, reg["anm"]["field_representation"])
    np.testing.assert_allclose(A[ok], 3 * G3 * S[ok], rtol=1e-9)
    rule = vk.rule_calls(S, reg["questions"]["Q1"]["bar"])
    near = np.abs(S.max(axis=1) - reg["questions"]["Q1"]["bar"]) < 1e-12
    assert np.all(calls[ok & ~near] == rule[ok & ~near]) and calls[5] == ""


# ----------------------------------------------------------------------------- script pieces
def test_variants_for_cell_counts_and_identities(reg):
    m = _script()
    rng = np.random.default_rng(2)
    G = 60
    counts = rng.poisson(2.0, size=G)
    counts[:3] = [40, 30, 20]
    row = (counts * 0.61).astype(np.float32)
    token_ids = np.arange(100, 100 + G)
    genes = [{"symbol": "g1", "column": 1, "vocab_id": 101}, {"symbol": "absent", "column": 59, "vocab_id": 159}]
    row[59] = 0.0
    V = m.variants_for_cell(5, row, np.ones(G), token_ids, 0, genes, reg, seq_len=2048, norm=1e4)
    n_thin = 6 * 2
    assert len(V["ids"]) == 1 + n_thin + 3  # floor, 12 thinning draws, 3 eps for the one present gene
    assert V["identical"][0] and V["desc"][0][0] == m.KIND_FLOOR
    gene_rows = [i for i, d in enumerate(V["desc"]) if d[0] == m.KIND_GENE]
    masked = gene_rows[0]  # eps 1.0 first
    assert V["desc"][masked][1] == 1.0 and 101 not in V["ids"][masked].tolist()
    assert V["count_dev"] <= 1.0 and V["ntok_base"] == int(np.sum(row > 0))
    V2 = m.variants_for_cell(5, row, np.ones(G), token_ids, 0, genes, reg, seq_len=2048, norm=1e4)
    assert all(np.array_equal(x, y) for x, y in zip(V["ids"], V2["ids"]))  # seeded


def test_site4_run_refused_without_committed_addendum(tmp_path):
    m = _script()
    rd = tmp_path / "registration"
    rd.mkdir()
    for f in ("registration_v3.json", "registration_v3.json.sha256", "amendment_A1.json", "amendment_A1.json.sha256"):
        shutil.copy(REPO / "registration" / f, rd / f)
    a = SimpleNamespace(registration_dir=rd, cell_pool="e2_subset", stage="embed")
    with pytest.raises(SystemExit, match="site4 run refused"):
        m.check_registration(a)
    a.cell_pool = "val"
    info = m.check_registration(a)["info"]
    assert info["addendum_exists"] is False


def test_spec_hash_is_stable():
    m = _script()
    assert m.spec_sha() == _script().spec_sha()


def test_register_core_reads_no_site4_rows():
    """Poisoning every site4 row of z, ADT, cell type and RNA leaves the addendum core unchanged; poisoning
    val rows changes it (the positive control). Uses the real metadata and ADT with an 8-d random z."""
    if not DATA.exists():
        pytest.skip("data pack not available")
    m = _script()
    reg_ = va.load_registration_amended()
    npz = np.load(DATA, allow_pickle=False)
    split, sites, donors = (npz[k].astype(str) for k in ("split", "sites", "donors"))
    ct = npz["cell_types"].astype(str)
    adt = npz["adt"].astype(np.float32)
    names = [str(x) for x in npz["adt_names"]]
    rng = np.random.default_rng(0)
    z = rng.normal(size=(split.size, 8)).astype(np.float32)

    def rna_row(i):
        return (np.random.default_rng(int(i)).poisson(1.0, size=50) * 0.5).astype(np.float32)

    core = m.register_core(z, adt, ct, rna_row, split, sites, donors, names, reg_)
    test = np.where(split == "test")[0]
    z2, a2, c2 = z.copy(), adt.copy(), ct.copy()
    z2[test] = rng.normal(size=(test.size, 8))
    a2[test] = rng.normal(size=(test.size, adt.shape[1]))
    c2[test] = rng.choice(sorted(set(ct.tolist())), size=test.size)
    ts = set(test.tolist())
    core_t = m.register_core(z2, a2, c2, lambda i: np.full(50, 0.3, np.float32) if i in ts else rna_row(i),
                             split, sites, donors, names, reg_)
    assert core_t == core
    vl = np.where(split == "val")[0]
    z3 = z.copy()
    z3[vl] = rng.normal(size=(vl.size, 8))
    assert m.register_core(z3, adt, ct, rna_row, split, sites, donors, names, reg_) != core
