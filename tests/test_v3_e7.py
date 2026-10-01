"""Tests for E7 (bridge_anm/lib/v3_e7.py, scripts/mode_a_inverse_e7.py).

Synthetic tests only (small random transformer layers, synthetic outputs); nothing here reads val or site4 data.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bridge_anm"))
sys.path.insert(0, str(REPO))

from lib import v3_e7 as v7  # noqa: E402


def _script():
    spec = importlib.util.spec_from_file_location("mode_a_inverse_e7", REPO / "scripts" / "mode_a_inverse_e7.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ----------------------------------------------------------------------------- patches
def test_pairs_for_fraction_and_codes():
    assert v7.pairs_for_fraction(1 / 32, 1000, 4) == 15
    assert v7.pairs_for_fraction(1 / 2, 1000, 4) == 249
    assert v7.pairs_for_fraction(1 / 32, 10, 1) == 1  # at least one pair
    assert [v7.frac_code(f) for f in (1 / 32, 1 / 8, 1 / 4, 1 / 2)] == [32, 128, 256, 512]


def test_draw_patch_sets_matched_disjoint_non_g_and_twins_add_only_g_pairs():
    g, ng = [3, 10, 50], [i for i in range(300) if i not in (3, 10, 50)]
    M, T = v7.draw_patch_sets(np.random.default_rng(0), g, ng, 20, 8)
    assert len(M) == len(T) == 8
    for pm, pt in zip(M, T):
        assert pm.shape == (20, 2) and np.unique(pm).size == 40 and not set(pm.ravel()) & set(g)
        assert np.array_equal(pt[:20], pm) and pt.shape == (23, 2)
        assert sorted(pt[20:, 0].tolist()) == g and not set(pt[20:, 1]) & set(pm.ravel()) and not set(pt[20:, 1]) & set(g)
        assert np.unique(pt).size == pt.size  # disjoint pairs
    M2, T2 = v7.draw_patch_sets(np.random.default_rng(0), g, ng, 20, 8)
    assert all(np.array_equal(x, y) for x, y in zip(M + T, M2 + T2))
    assert v7.draw_patch_sets(np.random.default_rng(0), [], ng, 5, 3)[1] == []
    with pytest.raises(ValueError):
        v7.draw_patch_sets(np.random.default_rng(0), g, ng[:40], 20, 1)


def test_g_only_patches_average_every_g_token_with_a_distinct_non_g_partner():
    g, ng = [1, 5], list(range(10, 60))
    P = v7.draw_g_only(np.random.default_rng(1), g, ng, 4)
    assert len(P) == 4 and all(sorted(p[:, 0].tolist()) == g and set(p[:, 1]) <= set(ng) and np.unique(p).size == 4 for p in P)
    assert v7.draw_g_only(np.random.default_rng(1), [], ng, 4) == []


def test_average_patch_preserves_the_gene_mean_and_one_sided_moves_it():
    rng = np.random.default_rng(2)
    H = rng.standard_normal((40, 6))
    pairs = np.array([[0, 7], [3, 9], [11, 20]])
    Ha = v7.apply_average(H, pairs)
    assert np.abs(Ha.mean(0) - H.mean(0)).max() < 1e-14
    assert np.allclose(Ha[0], Ha[7]) and np.allclose(Ha[0], 0.5 * (H[0] + H[7]))
    H1 = v7.apply_average(H, pairs, one_sided=True)
    want = (H[pairs[:, 1]] - H[pairs[:, 0]]).sum(0) / (2 * 40)
    assert np.allclose(H1.mean(0) - H.mean(0), want)
    torch = pytest.importorskip("torch")
    Ht = v7.apply_average(torch.from_numpy(H), pairs)
    assert np.allclose(Ht.numpy(), Ha)


def test_region_shift_moves_the_mean_by_the_shift_through_the_patch_tokens_only():
    rng = np.random.default_rng(7)
    H, s = rng.standard_normal((40, 5)), rng.standard_normal(5)
    pairs = np.array([[0, 7], [3, 9], [11, 20]])
    Hr = v7.apply_region_shift(H, pairs, s)
    assert np.allclose(Hr.mean(0) - H.mean(0), s)
    changed = np.where(np.any(Hr != H, axis=1))[0]
    assert changed.tolist() == sorted(pairs.ravel().tolist())
    torch = pytest.importorskip("torch")
    Ht = v7.apply_region_shift(torch.from_numpy(H), pairs, torch.from_numpy(s))
    assert np.allclose(Ht.numpy(), Hr)


def test_clamp_moves_the_mean_by_the_shift_and_keeps_the_deviations():
    rng = np.random.default_rng(3)
    H, s = rng.standard_normal((30, 5)), rng.standard_normal(5)
    Hc = v7.apply_clamp(H, s)
    assert np.allclose(Hc.mean(0) - H.mean(0), s) and np.allclose(Hc - Hc.mean(0), H - H.mean(0))


# ----------------------------------------------------------------------------- consumer architecture facts (synthetic)
def _layer(torch, d=16, heads=2, seed=0):
    torch.manual_seed(seed)
    return torch.nn.TransformerEncoderLayer(d, heads, 2 * d, dropout=0.0, batch_first=True, activation="gelu").eval().double()


def test_consumer_without_positional_input_is_permutation_invariant_and_blind_to_swaps():
    torch = pytest.importorskip("torch")
    m = _script()
    layer = _layer(torch)
    H = torch.randn(1, 25, 16, dtype=torch.float64)
    with torch.no_grad():
        z = m.mld.layer_fn(layer, H).mean(1)
        perm = torch.randperm(25)
        assert torch.allclose(m.mld.layer_fn(layer, H[:, perm]).mean(1), z, atol=1e-12)
        Hs = H.clone()
        Hs[0, [2, 9]] = H[0, [9, 2]]
        assert torch.allclose(m.mld.layer_fn(layer, Hs).mean(1), z, atol=1e-12)
        Ha = v7.apply_average(H[0], np.array([[2, 9]]))[None]
        assert float((m.mld.layer_fn(layer, Ha).mean(1) - z).abs().max()) > 1e-6  # an average patch is visible


def test_padding_position_patch_is_exactly_invisible_with_a_bool_mask():
    torch = pytest.importorskip("torch")
    layer = _layer(torch, seed=1).float()
    L, P = 20, 5
    H = torch.randn(1, L + P, 16)
    mask = torch.zeros(1, L + P, dtype=torch.bool)
    mask[0, L:] = True
    with torch.no_grad():
        y0 = layer(H, src_key_padding_mask=mask)[0, :L]
        H2 = H.clone()
        H2[0, L + 2] = 7.0 * H2[0, L + 2]
        y1 = layer(H2, src_key_padding_mask=mask)[0, :L]
        yu = layer(H[:, :L])[0]
    assert float((y1 - y0).abs().max()) == 0.0
    assert float((yu - y0).abs().max()) < 1e-5


def test_layer12_output_mean_preserving_patch_leaves_the_pooled_state_unchanged():
    rng = np.random.default_rng(4)
    Y = rng.standard_normal((50, 8))
    Ya = v7.apply_average(Y, np.array([[1, 30]]))
    assert np.abs(Ya.mean(0) - Y.mean(0)).max() < 1e-14


# ----------------------------------------------------------------------------- per-cell statistics
def _cell(O0, Om, Ot, Og, Oc, code=128, Or=None):
    Or = np.zeros((0, 4)) if Or is None else Or
    O = np.concatenate([Om, Ot, Og, Oc, Or])
    n_m, n_t, n_g, n_c, n_r = len(Om), len(Ot), len(Og), len(Oc), len(Or)
    kind = np.array([v7.KIND_MATCHED] * n_m + [v7.KIND_TWIN] * n_t + [v7.KIND_G_ONLY] * n_g + [v7.KIND_CLAMP] * n_c
                    + [v7.KIND_REGION] * n_r)
    codes = np.array([code] * (n_m + n_t) + [0] * (n_g + n_c + n_r))
    base = np.array([-1] * n_m + list(range(n_t)) + [-1] * (n_g + n_c) + [0] * n_r)
    return O, np.asarray(O0, float), kind, codes, base


def test_cell_summary_separates_gene_mean_from_gene_mean_plus_g():
    u = np.ones(4)

    def dfun(a, b):
        return v7.deltas(a, b, u)
    O0 = np.zeros(4)
    Om = np.full((3, 4), 0.01)  # non-G rearrangements: small
    Ot = Om + np.array([0.2, 0, 0, 0])  # G twins: only the G states differ from the matched patch -> large
    Og = np.full((2, 4), 0.03)
    Oc = np.array([[0.5, 0, 0, 0], [0, 0.4, 0, 0]])
    Or = np.array([[0, 0, 0.45, 0], [0, 0, 0, 0.35]])
    s = v7.cell_summary(*_cell(O0, Om, Ot, Og, Oc, Or=Or), [128], dfun)
    assert np.isclose(s["D_gene_mean_plus_G"], 0.01) and np.isclose(s["D_gene_mean"], 0.21)
    assert np.isclose(s["W_context"], 0.2) and np.isclose(s["W_pure"], 0.03) and np.isclose(s["W"], 0.2)
    assert np.isclose(s["P"], 0.5) and np.isclose(s["R"], 0.45) and np.isnan(s["D_all_token_states"])
    assert s["n_region"] == 2 and s["n_clamp"] == 2  # the shift controls never enter a matched family
    s2 = v7.cell_summary(*_cell(O0, Om, Ot, Og, Oc, Or=Or), [32], dfun)  # other magnitude: only code-0 patches count
    assert np.isnan(s2["D_gene_mean_plus_G"]) and np.isclose(s2["W"], 0.03) and np.isclose(s2["P"], 0.5)
    assert np.isclose(s2["R"], 0.45)


def test_score_and_deltas_units():
    O = np.array([[1.0, 2.0, 3.0, 1.0]])
    assert np.isclose(v7.score(O)[0], 1.0)
    assert np.isclose(v7.deltas(O, np.zeros(4), np.array([1.0, 4.0, 1.0, 2.0]))[0], 3.0)
    assert np.isclose(v7.score_deltas(O, np.zeros(4), 2.0)[0], 0.5)
    assert np.isclose(v7.observable_deltas(O, O, "panel", {"panel": np.ones(4)})[0], 0.0)


# ----------------------------------------------------------------------------- selection rule
def test_dev_selection_minimal_retained_candidate_or_none():
    rng = np.random.default_rng(5)
    small, big = rng.uniform(0, 0.04, 200), rng.uniform(0.06, 0.2, 200)
    nan = np.full(200, np.nan)
    r = v7.dev_selection({"gene_mean": small, "gene_mean_plus_G": small, "all_token_states": nan}, 0.05)
    assert r["selected"] == "gene_mean"
    r = v7.dev_selection({"gene_mean": big, "gene_mean_plus_G": small, "all_token_states": nan}, 0.05)
    assert r["selected"] == "gene_mean_plus_G" and not r["candidates"]["gene_mean"]["retained"]
    r = v7.dev_selection({"gene_mean": big, "gene_mean_plus_G": big, "all_token_states": nan}, 0.05)
    assert r["selected"] is None and r["candidates"]["all_token_states"]["reason"].startswith("no matched pair")
    x = small.copy()
    x[0] = 0.3  # one cell above tol: the q95 rule keeps the candidate, the literal reading would not
    r = v7.dev_selection({"gene_mean": x, "gene_mean_plus_G": x, "all_token_states": nan}, 0.05)
    assert r["selected"] == "gene_mean" and r["literal_rule_selected"] is None


# ----------------------------------------------------------------------------- bootstrap, statuses, verdict
def test_two_stage_draws_are_seeded_and_resample_donors_then_cells():
    donor = np.array(["a"] * 6 + ["b"] * 4)
    cls = np.array(["NK", "T"] * 5)
    d1, d2 = v7.two_stage_draws(donor, cls, 50, 1), v7.two_stage_draws(donor, cls, 50, 1)
    assert all(np.array_equal(x, y) for x, y in zip(d1, d2))
    assert any(set(donor[x]) == {"a"} for x in d1)  # a replicate can draw one donor twice


def test_bounds_support_and_reject():
    rng = np.random.default_rng(6)
    n = 400
    donor = np.array(["13272"] * 200 + ["19593"] * 200)
    cls = np.array(["NK", "T"] * 200)
    ix = np.arange(n)
    small = rng.uniform(0, 0.02, n)
    big = rng.uniform(0.1, 0.3, n)
    for x, want in ((small, "supported"), (big, "rejected")):
        bm = v7.bounds(lambda ii, x=x: v7.finite_median(x[ii]), ix, donor, cls, 300, 1)
        bq = v7.bounds(lambda ii, x=x: v7.q95(x[ii]), ix, donor, cls, 300, 1)
        assert v7.candidate_status(bm, bq, 0.05) == want
        assert set(bm["per_donor"]) == {"13272", "19593"}
    mid = rng.uniform(0.0, 0.1, n)
    bm = v7.bounds(lambda ii: v7.finite_median(mid[ii]), ix, donor, cls, 300, 1)
    bq = v7.bounds(lambda ii: v7.q95(mid[ii]), ix, donor, cls, 300, 1)
    assert v7.candidate_status(bm, bq, 0.05) == "unresolved"


def test_verdict_table():
    sup, rej, unr = "supported", "rejected", "unresolved"
    def V(sel, st, pc, sep, impl, f64, det=True):
        return v7.verdict(sel, st, pc, sep, impl, f64, det)
    assert V("gene_mean_plus_G", {"gene_mean": rej, "gene_mean_plus_G": sup}, True, True, True, True)["code"] == "LOOP_COMPLETE"
    assert V("gene_mean_plus_G", {"gene_mean": unr, "gene_mean_plus_G": sup}, True, True, True, True)["code"] == \
        "SUPPORTED_REVISION_NOT_SHOWN_NECESSARY"
    assert V("gene_mean_plus_G", {"gene_mean": rej, "gene_mean_plus_G": sup}, True, False, True, True)["code"] == \
        "SUPPORTED_REVISION_NOT_SHOWN_NECESSARY"
    assert V("gene_mean", {"gene_mean": sup, "gene_mean_plus_G": sup}, True, False, True, True)["code"] == "GENE_MEAN_SUPPORTED"
    assert V("gene_mean", {"gene_mean": sup, "gene_mean_plus_G": sup}, False, False, True, True)["code"] == "UNINFORMATIVE"
    # the patch-region detectability witness gates every support claim, as the clamp does
    assert V("gene_mean", {"gene_mean": sup, "gene_mean_plus_G": sup}, True, False, True, True, False)["code"] == "UNINFORMATIVE"
    assert V("gene_mean_plus_G", {"gene_mean": rej, "gene_mean_plus_G": sup}, True, True, True, True, False)["code"] == \
        "UNINFORMATIVE"
    assert V("gene_mean", {"gene_mean": rej, "gene_mean_plus_G": rej}, True, True, True, True)["code"] == "SELECTED_REJECTED"
    # a rejection needs no control: a failed control never masks it
    assert V("gene_mean", {"gene_mean": rej, "gene_mean_plus_G": rej}, False, False, True, True, False)["code"] == \
        "SELECTED_REJECTED"
    assert V("gene_mean", {"gene_mean": unr, "gene_mean_plus_G": unr}, False, False, True, True, False)["code"] == \
        "SELECTED_UNRESOLVED"
    assert V(None, {"gene_mean": rej, "gene_mean_plus_G": rej}, False, False, True, True, False)["code"] == \
        "REJECTED_NO_REVISION"
    assert V(None, {"gene_mean": rej, "gene_mean_plus_G": unr}, True, True, True, True)["code"] == \
        "REJECTED_GENE_MEAN_NO_SELECTION"
    assert V(None, {"gene_mean": unr, "gene_mean_plus_G": unr}, True, True, True, True)["code"] == "NO_SELECTION_UNRESOLVED"
    assert V("gene_mean", {"gene_mean": sup}, True, True, False, True)["code"] == "NOT_VALIDATED"
    assert V("gene_mean", {"gene_mean": sup}, True, True, True, False)["code"] == "UNRESOLVED_PRECISION"
    assert set(v7.VERDICTS) >= {"LOOP_COMPLETE", "NOT_VALIDATED", "UNRESOLVED_PRECISION", "UNINFORMATIVE"}


# ----------------------------------------------------------------------------- script: design block, guards, histories
def test_declared_block_is_json_stable_and_states_the_design():
    m = _script()
    d = json.loads(json.dumps(m.DECLARED))
    assert d == json.loads(json.dumps(m.DECLARED, sort_keys=True))
    assert d["candidates"]["order"] == ["gene_mean", "gene_mean_plus_G", "all_token_states"]
    assert set(d["candidates"]["G"]["genes"]) == {"NCAM1", "KLRD1", "NCR1", "CD3E", "CD3D", "CD3G"}
    assert d["tolerance"]["tol"] == 0.05 and "no_swaps" in d["histories"]
    assert m.LADDER == (1 / 32, 1 / 8, 1 / 2) and m.FRESH == (1 / 4,) and not set(m.FRESH) & set(m.LADDER)
    assert set(m.G_GENES) <= set(m.mld.TOKEN_PROBE_IDS)


def test_histories_are_seeded_per_phase_and_fresh_at_confirmation():
    m = _script()
    ids = np.arange(1000, 1400)
    gtok = {1003, 1100, 1250}
    h1, g, ng = m.build_histories(77, ids, gtok, m.LADDER, "dev")
    h2, _, _ = m.build_histories(77, ids, gtok, m.LADDER, "dev")
    h3, _, _ = m.build_histories(77, ids, gtok, m.LADDER, "site4")
    assert g == [3, 100, 250] and len(ng) == 397
    assert all(np.array_equal(a[2], b[2]) for a, b in zip(h1, h2))
    kinds = [h[0] for h in h1]
    assert kinds.count(v7.KIND_MATCHED) == 3 * m.N_MATCHED and kinds.count(v7.KIND_TWIN) == 3 * m.N_MATCHED
    assert kinds.count(v7.KIND_G_ONLY) == m.N_G_ONLY and kinds.count(v7.KIND_CLAMP) == 2
    reg = [h for h in h1 if h[0] == v7.KIND_REGION]
    assert len(reg) == m.N_REGION and sorted(h[4] for h in reg) == [-1, 1]
    top = h1[reg[0][3]]  # the region is matched_1 of the largest magnitude: non-G tokens only
    assert top[0] == v7.KIND_MATCHED and top[1] == v7.frac_code(max(m.LADDER)) and np.array_equal(reg[0][2], top[2])
    assert not set(v7.region_positions(reg[0][2])) & set(g)
    pm = [h for h in h1 if h[0] == v7.KIND_MATCHED]
    p3 = [h for h in h3 if h[0] == v7.KIND_MATCHED]
    assert not any(np.array_equal(a[2], b[2]) for a, b in zip(pm, p3))  # fresh pairs at confirmation
    for h in h1:
        if h[0] == v7.KIND_TWIN:
            assert np.array_equal(h[2][: len(h1[h[3]][2])], h1[h[3]][2])


def test_smoke_flags_and_registered_cells():
    m = _script()
    with pytest.raises(SystemExit):
        m.parse_args(["--stage", "develop", "--cell-pool", "smoke"])
    with pytest.raises(SystemExit):
        m.parse_args(["--stage", "develop", "--smoke", "t"])  # registered cells are never smoked
    with pytest.raises(SystemExit):
        m.parse_args(["--stage", "develop", "--smoke", "t", "--cell-pool", "smoke", "--n-smoke", "6"])
    a = m.parse_args(["--stage", "confirm", "--smoke", "t", "--cell-pool", "smoke"])
    assert a.n_smoke == 5 and a.family == "site4"


def _fake_registration(tmp_path, m, sel=None):
    reg = tmp_path / "registration"
    (reg / "addenda").mkdir(parents=True)
    add = {"declared": json.loads(json.dumps(m.DECLARED)), "leakage_check": {"passed": True},
           "inputs": {"code_sha256": m.code_hashes(), "head_sha256": "x"}, "amendments_sha256": {}}
    (reg / "addenda" / "E7.json").write_text(json.dumps(add))
    if sel is not None:
        (reg / "addenda" / "E7_selection.json").write_text(json.dumps(sel))
    return SimpleNamespace(registration_dir=reg, smoke=None, out_dir=tmp_path / "out", head_ckpt=tmp_path / "none.pt")


def test_development_and_site4_refused_without_committed_records(tmp_path):
    m = _script()
    a = _fake_registration(tmp_path, m)
    for need in ("develop", "select", "followup", "confirm", "report"):
        with pytest.raises(SystemExit, match="refused"):
            m.check_registration(a, need)
    a.smoke = "t"
    info = m.check_registration(a, "confirm")  # smoke: recorded, not refused
    assert info["declared_equals_script"] and info["code_matches_addendum"] and info["addendum_committed"] is False


def test_confirm_refused_without_a_final_committed_selection_record(tmp_path):
    m = _script()
    a = _fake_registration(tmp_path, m, sel={"final": False, "addendum_sha256": "y", "primary": {"selected": None}})
    with pytest.raises(SystemExit, match="final development selection record"):
        m.check_registration(a, "confirm")


def test_run_cell_gates_flag_a_no_op_patch_and_require_the_controls():
    m = _script()
    n = 6
    kind = np.array([v7.KIND_MATCHED, v7.KIND_TWIN, v7.KIND_G_ONLY, v7.KIND_CLAMP, v7.KIND_CLAMP, v7.KIND_REGION])
    arr = {"O0": np.zeros(4), "O_pad0": np.zeros((1, 4)), "O_pad1": np.zeros((1, 4)), "O_l12": np.zeros((1, 4)),
           "O_perm": np.zeros((1, 4)), "dmean": np.zeros(n), "rows_mismatch": np.zeros(n, np.int64), "kind": kind}
    mc = {"controls": {"module_vs_explicit_abs": 0.0, "padded_vs_unpadded_abs": 0.0, "official_layer11_gene_mean_abs": 0.0}}
    units = {"panel": np.ones(4), "score": 1.0}
    g = m.cell_gates(mc, {**arr, "kind": np.append(kind, v7.KIND_REGION), "dmean": np.zeros(n + 1),
                          "rows_mismatch": np.zeros(n + 1, np.int64)}, units)
    assert g["ok"] and g["controls_present"]
    assert not m.cell_gates(mc, arr, units)["ok"]  # one patch-region history missing
    full = {**arr, "kind": np.append(kind, v7.KIND_REGION), "dmean": np.zeros(n + 1)}
    bad = {**full, "rows_mismatch": np.array([0, 0, 2, 0, 0, 0, 0])}  # a pair patch that changed no row
    assert not m.cell_gates(mc, bad, units)["ok"]
    shift = {**full, "rows_mismatch": np.zeros(n + 1, np.int64), "dmean": np.array([0, 0, 0, 0, 1e-3, 0, 0.0])}
    assert not m.cell_gates(mc, shift, units)["ok"]  # the clamp did not move the gene-mean by its vector


def test_select_writes_its_record_once_and_a_stopped_record_ends_e7(tmp_path):
    m = _script()
    a = _fake_registration(tmp_path, m, sel={"final": True, "primary": {"selected": None, "stopped": None}})
    a.roster, a.cell_pool = "dev", "registered"
    with pytest.raises(SystemExit, match="written once"):
        m.stage_select(a, {})
    (a.out_dir / "site4" / "cells").mkdir(parents=True)
    (a.out_dir / "site4" / "cells" / "cell_1.npz").write_bytes(b"x")
    with pytest.raises(SystemExit, match="confirmation outputs already exist"):
        m.stage_select(a, {})
    with pytest.raises(SystemExit, match="E7 ends with no verdict"):
        m.refuse_if_stopped({"primary": {"stopped": "implementation check failed on a development cell"}}, "confirm")
    with pytest.raises(SystemExit, match="E7 ends with no verdict"):
        m.refuse_if_stopped({"primary": {"stopped": None}, "followup": {"stopped": "float64"}}, "report")
    m.refuse_if_stopped({"primary": {"stopped": None}, "followup": {"stopped": None}}, "confirm")


def test_committed_addendum_matches_the_code_on_disk():
    m = _script()
    f = REPO / "registration" / "addenda" / "E7.json"
    if not f.exists():
        pytest.skip("no E7 addendum")
    add = json.loads(f.read_text())
    assert add["declared"] == json.loads(json.dumps(m.DECLARED)), "DECLARED changed after the E7 addendum"
    assert add["inputs"]["code_sha256"] == m.code_hashes(), "reused E5 functions or v3_e7.py changed after the addendum"
    assert add["leakage_check"]["passed"] and not add["smoke"]
