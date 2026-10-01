"""Tests for E5-M (bridge_anm/lib/v3_e5m.py, scripts/mode_a_observer_sufficiency.py).

Synthetic tests always run; the ANM test skips without ANM_ROOT; nothing here reads site4 data.
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bridge_anm"))
sys.path.insert(0, str(REPO))

from lib import v3_e5m as e5m  # noqa: E402


def _script():
    spec = importlib.util.spec_from_file_location("mode_a_observer_sufficiency", REPO / "scripts" / "mode_a_observer_sufficiency.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ----------------------------------------------------------------------------- per-cell linear algebra
def test_residual_zero_for_exact_factorisation_and_one_for_zero_map():
    rng = np.random.default_rng(0)
    Kz, H = rng.standard_normal((20, 5)), rng.standard_normal((3, 20))
    assert e5m.factorization_residual(H @ Kz, H, Kz) < 1e-12
    assert abs(e5m.factorization_residual(H @ Kz, np.zeros_like(H), Kz) - 1.0) < 1e-12
    assert np.isnan(e5m.factorization_residual(np.zeros((3, 5)), H, Kz))


def test_same_cell_fit_is_automatic_when_kz_has_full_column_rank():
    rng = np.random.default_rng(1)
    Kz, KO = rng.standard_normal((50, 6)), rng.standard_normal((4, 6))  # K_O unrelated to K_z
    assert np.linalg.matrix_rank(Kz) == 6
    assert e5m.same_cell_fit_residual(KO, Kz) < 1e-10


def test_near_null_share_finds_the_weak_source_combination():
    rng = np.random.default_rng(2)
    U, _ = np.linalg.qr(rng.standard_normal((30, 3)))
    Kz = U @ np.diag([1.0, 0.5, 1e-4])  # column 3 barely moves z
    KO = np.zeros((2, 3))
    KO[:, 2] = [1.0, 2.0]  # the observer responds only to the near-null column
    nn = e5m.near_null_share(Kz, KO, 0.05)
    assert nn["near_null_dim"] == 1 and nn["rank"] == 2 and abs(nn["share"] - 1.0) < 1e-10
    KO2 = np.zeros((2, 3))
    KO2[:, 0] = 1.0
    assert e5m.near_null_share(Kz, KO2, 0.05)["share"] < 1e-12


def test_quotient_diagnostic_rank_and_participation_ratio():
    KO = np.diag([1.0, 1.0, 0.0])[:2]  # rank 2, equal singular values
    q = e5m.quotient_diagnostic(KO, 0.05)
    assert q["rank"] == 2 and abs(q["participation_ratio"] - 2.0) < 1e-12
    assert e5m.quotient_diagnostic(np.array([[1.0, 0.0], [0.0, 1e-3]]), 0.05)["rank"] == 1


def test_ridge_select_recovers_a_linear_map_and_reports_lambda():
    rng = np.random.default_rng(3)
    H = rng.standard_normal((4, 10))
    Xtr, Xva = rng.standard_normal((400, 10)), rng.standard_normal((100, 10))
    f = e5m.ridge_select(Xtr, Xtr @ H.T + 0.5, Xva, Xva @ H.T + 0.5, (1e-4, 1e-2, 1.0))
    assert f["lambda"] == 1e-4 and np.max(np.abs(f["H"] - H)) < 1e-3
    assert np.allclose(f["intercept"], 0.5, atol=1e-3) and min(f["val_r2"]) > 0.999


def test_richardson_removes_the_quadratic_truncation_term():
    f, x0, h = (lambda x: np.sin(3 * x)), 0.3, 0.02
    D = lambda hh: (f(x0 + hh) - f(x0 - hh)) / (2 * hh)  # noqa: E731
    exact = 3 * np.cos(3 * x0)
    assert abs(e5m.richardson(D(h), D(2 * h)) - exact) < abs(D(h) - exact) / 100


def test_fd_check():
    J = np.array([1.0, 2.0, 3.0])
    assert e5m.fd_check(J, J * 1.01, 0.05, 0.99)[2]
    assert not e5m.fd_check(J, J * 1.2, 0.05, 0.99)[2]
    assert not e5m.fd_check(J, -J, 0.05, 0.99)[2]


# ----------------------------------------------------------------------------- bootstrap and decisions
def test_two_stage_draws_are_seeded_and_resample_donors_then_cells():
    donor = np.array(["a"] * 6 + ["b"] * 4)
    role = np.array(["pair", "random"] * 5)
    d1 = e5m.two_stage_draws(donor, role, 50, 1)
    d2 = e5m.two_stage_draws(donor, role, 50, 1)
    assert all(np.array_equal(x, y) for x, y in zip(d1, d2))
    sizes = {x.size for x in d1}
    assert sizes <= {8, 10, 12}  # (b, b), (a, b) or (a, a)
    for x in d1:  # every draw keeps whole donor x role strata sizes
        for dn in ("a", "b"):
            k = np.sum(donor[x] == dn)
            assert k % (6 if dn == "a" else 4) == 0


def test_contrast_and_stratum_median_diff():
    x = np.array([0.5, 0.1, 0.6, 0.2, 0.7, 0.3])
    role = np.array(["pair", "random"] * 3)
    donor = np.array(["a", "a", "a", "b", "b", "b"])
    f = e5m.stratum_median_diff(x, role)
    assert abs(f(np.arange(6)) - 0.4) < 1e-12
    c = e5m.contrast({"d": f}, np.arange(6), donor, e5m.two_stage_draws(donor, role, 200, 1))
    assert abs(c["d"]["value"] - 0.4) < 1e-12 and set(c["d"]["per_donor"]) == {"a", "b"}
    assert c["d"]["ci_two_stage"][0] > 0


def test_common_rule():
    assert e5m.common_rule(0.1, [0.02, 0.2], {"a": 0.08, "b": 0.03}, 0.05) == "win"
    assert e5m.common_rule(0.1, [0.02, 0.2], {"a": 0.08, "b": 0.02}, 0.05) == "inconclusive"  # one donor < margin / 2
    assert e5m.common_rule(-0.1, [-0.2, -0.02], {"a": -0.08, "b": -0.03}, 0.05) == "loss"
    assert e5m.common_rule(0.0, [-0.03, 0.04], {"a": 0.0, "b": 0.0}, 0.05) == "equivalent"
    assert e5m.common_rule(0.0, None, {}, 0.05) == "inconclusive"


def _c(delta, ci, pdn, dd_lb, level_ub):
    return {"delta_tok": {"value": delta, "ci_two_stage": ci, "per_donor": pdn},
            "dd_tok_minus_rand": {"value": dd_lb, "ci_two_stage": [dd_lb, dd_lb + 0.1], "per_donor": {}},
            "level_tok_pair": {"value": level_ub, "ci_two_stage": [0.0, level_ub], "per_donor": {}}}


def test_h5m_verdict_table():
    win = (0.1, [0.02, 0.2], {"a": 0.08, "b": 0.06})
    eq = (0.0, [-0.03, 0.03], {"a": 0.0, "b": 0.0})
    assert e5m.h5m_verdict(_c(*win, 0.01, 0.3), 0.05, False, True)["code"] == "NOT_VALIDATED"
    assert e5m.h5m_verdict(_c(*win, 0.01, 0.3), 0.05, True, False)["code"] == "NOT_VALIDATED"
    assert e5m.h5m_verdict(_c(*win, 0.01, 0.3), 0.05, True, True)["code"] == "REJECT_SUFFICIENCY"
    assert e5m.h5m_verdict(_c(*win, -0.01, 0.3), 0.05, True, True)["code"] == "INCONCLUSIVE"  # fails the O_rand null
    assert e5m.h5m_verdict(_c(*eq, 0.0, 0.04), 0.05, True, True)["code"] == "PASS_SUFFICIENT"
    assert e5m.h5m_verdict(_c(*eq, 0.0, 0.2), 0.05, True, True)["code"] == "H5M_FALSIFIED"


# ----------------------------------------------------------------------------- ANM
@pytest.mark.skipif(not os.environ.get("ANM_ROOT"), reason="ANM_ROOT not set")
def test_anm_status_not_informative_for_injective_kz_and_residual_decides():
    f = Path(os.environ["ANM_ROOT"]) / "benchmarks" / "paper1_retained_response" / "analysis.py"
    if not f.exists():
        pytest.skip("ANM analysis.py not found")
    spec = importlib.util.spec_from_file_location("anm_analysis", f)
    anm = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(anm)
    rng = np.random.default_rng(4)
    Kz, H = rng.standard_normal((16, 4)), rng.standard_normal((3, 16))
    ok = e5m.anm_factorization_status(anm.analyze_factorization, Kz, H @ Kz, H, atol=0.0, rtol=1e-6)
    assert ok["source_kernel_status"] == "NOT_INFORMATIVE" and ok["factorization_status"] == "PASS"
    assert ok["status"] == "NOT_INFORMATIVE"
    bad = e5m.anm_factorization_status(anm.analyze_factorization, Kz, H @ Kz + 0.5, H, atol=0.0, rtol=1e-6)
    assert bad["factorization_status"] == "REJECTED" and bad["status"] == "REJECTED"


# ----------------------------------------------------------------------------- clamp map, observers, layers (torch)
def test_clamp_map_residual_is_the_pooling_discarded_response():
    torch = pytest.importorskip("torch")
    from torch.func import jvp

    from lib.v3_e3_readouts import R2AttnPool
    m = _script()
    torch.manual_seed(0)
    L, d = 9, 8
    Hs = torch.randn(1, L, d)
    obs = R2AttnPool(d, 3)
    with torch.no_grad():
        obs.q.copy_(torch.randn(d) * 3)
    f = lambda X: obs(X, torch.ones(1, L, dtype=torch.bool))[0][0]  # noqa: E731
    Lam = m.clamp_map_tokens(f, Hs)  # [3, d]
    dH = torch.randn(1, L, d)
    KO = jvp(f, (Hs,), (dH,))[1]
    dz = dH[0].mean(0)
    dev = dH - dH.mean(dim=1, keepdim=True)
    direct = jvp(f, (Hs,), (dev,))[1]
    assert torch.allclose(KO - Lam @ dz, direct, atol=1e-5)
    # for R2's form the clamp leaves the attention unchanged: Lambda = W diag(1/sd) at every input
    assert torch.allclose(Lam, obs.out.weight / obs.sd[None, :], atol=1e-5)
    assert torch.allclose(m.clamp_map_tokens(f, Hs * 3 + 1), Lam, atol=1e-5)
    # an observer of the gene-mean only (query 0: uniform attention) has clamp residual 0
    obs0 = R2AttnPool(d, 3)
    f0 = lambda X: obs0(X, torch.ones(1, L, dtype=torch.bool))[0][0]  # noqa: E731
    Lam0 = m.clamp_map_tokens(f0, Hs)
    K0 = jvp(f0, (Hs,), (dH,))[1]
    assert float((K0 - Lam0 @ dz).detach().norm() / K0.detach().norm()) < 1e-5


def test_e5_layer_fn_jvp_matches_module_finite_differences():
    torch = pytest.importorskip("torch")
    from torch.func import jvp
    m = _script()
    torch.manual_seed(1)
    layer = torch.nn.TransformerEncoderLayer(16, 2, 32, dropout=0.0, batch_first=True, activation="gelu").eval().double()
    x, t = torch.randn(1, 7, 16, dtype=torch.float64), torch.randn(1, 7, 16, dtype=torch.float64)
    with torch.no_grad():
        assert torch.allclose(m.mld.layer_fn(layer, x), layer(x), atol=1e-10)
        J = jvp(lambda h: m.mld.layer_fn(layer, h), (x,), (t,))[1]
        e = 1e-5
        D = (layer(x + e * t) - layer(x - e * t)) / (2 * e)
    assert float((D - J).norm() / J.norm()) < 1e-6


def test_rand_params_are_seeded_and_row_norms_matched():
    torch = pytest.importorskip("torch")
    from lib.v3_e3_readouts import R2AttnPool
    m = _script()
    r2 = R2AttnPool(12, 5)
    a, b = m.rand_params(r2, [3, 1], [23, 5]), m.rand_params(r2, [3, 1], [23, 5])
    assert all(np.array_equal(a[k], b[k]) for k in a)
    assert abs(np.linalg.norm(a["v_rand"]) - 1) < 1e-12 and abs(np.linalg.norm(a["q_rand"]) - 1) < 1e-12
    rows = r2.out.weight.detach().numpy()[[3, 1]]
    assert np.allclose(np.linalg.norm(a["W_rand"], axis=1), np.linalg.norm(rows, axis=1))
    del torch


def test_rowview_poisons_only_the_chosen_rows():
    m = _script()
    A = np.arange(40, dtype=np.float32).reshape(10, 4)
    v = m.RowView(A, np.array([2, 7]), 101)
    X = v.rows(np.array([1, 2, 3]))
    assert np.array_equal(X[[0, 2]], A[[1, 3]]) and not np.array_equal(X[1], A[2])
    assert np.array_equal(m.RowView(A).rows(np.array([2, 7])), A[[2, 7]])


def test_declared_block_is_json_stable_and_names_the_framing():
    m = _script()
    d = json.loads(json.dumps(m.DECLARED))
    assert d == json.loads(json.dumps(m.DECLARED, sort_keys=True))
    assert "NOT_INFORMATIVE" in d["anm_framing"] and "diagnostic" in d["anm_framing"]
    assert set(d["hypothesis"]["verdicts"]) == {"REJECT_SUFFICIENCY", "PASS_SUFFICIENT", "H5M_FALSIFIED", "INCONCLUSIVE",
                                                "NOT_VALIDATED"}


def test_site4_refused_without_a_committed_addendum(tmp_path):
    m = _script()
    reg = tmp_path / "registration"
    (reg / "addenda").mkdir(parents=True)
    a = SimpleNamespace(registration_dir=reg)
    with pytest.raises(SystemExit, match="site4 run refused"):
        m.check_registration(a, site4=True)
    info = m.check_registration(a, site4=False)  # smoke / select: reported, not refused
    assert info["addendum_sha256"] is None


def test_smoke_flags_need_a_smoke_note():
    m = _script()
    with pytest.raises(SystemExit):
        m.parse_args(["--stage", "respond", "--out-dir", "/tmp/x", "--cell-pool", "val"])
    a = m.parse_args(["--stage", "respond", "--out-dir", "/tmp/x", "--cell-pool", "val", "--smoke", "t"])
    assert a.cell_pool == "val" and a.jvp_chunk == 8


# ----------------------------------------------------------------------------- guard and report labels (verifier, before site4)
def test_guard_records_code_and_amendments_against_the_addendum(tmp_path):
    m = _script()
    reg = tmp_path / "registration"
    (reg / "addenda").mkdir(parents=True)
    am = reg / "amendment_A1.json"
    am.write_text("{}")
    import hashlib
    sha_am = hashlib.sha256(am.read_bytes()).hexdigest()
    now = {"e5_function_source_sha256": m.code_hashes(),
           "lib_v3_e5m_sha256": m.sha256_file(REPO / "bridge_anm" / "lib" / "v3_e5m.py")}
    add = {"declared": {}, "calibration": {"H_file_sha256": "x"}, "leakage_check": {"passed": True},
           "amendments_sha256": {"amendment_A1.json": sha_am}, "inputs": now}
    (reg / "addenda" / "E5M.json").write_text(json.dumps(add))
    info = m.check_registration(SimpleNamespace(registration_dir=reg), site4=False)
    assert info["code_matches_addendum"] and info["addendum_amendments_match_disk"]
    add["inputs"] = {**now, "lib_v3_e5m_sha256": "0" * 64}
    add["amendments_sha256"] = {"amendment_A1.json": "0" * 64}
    (reg / "addenda" / "E5M.json").write_text(json.dumps(add))
    info = m.check_registration(SimpleNamespace(registration_dir=reg), site4=False)
    assert not info["code_matches_addendum"] and not info["addendum_amendments_match_disk"]
    with pytest.raises(SystemExit, match="versions the addendum records"):
        m.check_registration(SimpleNamespace(registration_dir=reg), site4=True)


def test_committed_addendum_matches_the_code_and_amendments_on_disk():
    m = _script()
    if not (REPO / "registration" / "addenda" / "E5M.json").exists():
        pytest.skip("no E5M addendum")
    info = m.check_registration(SimpleNamespace(registration_dir=REPO / "registration"), site4=False)
    assert info["code_matches_addendum"], "E5's reused functions or v3_e5m.py changed after the E5M addendum"
    assert info["addendum_amendments_match_disk"]
    assert info["declared_equals_script"]


def test_smoke_on_site4_respond_or_report_is_refused():
    m = _script()
    for st in ("respond", "report", "all"):
        with pytest.raises(SystemExit):
            m.parse_args(["--stage", st, "--out-dir", "/tmp/x", "--smoke", "t"])
    assert m.parse_args(["--stage", "select", "--out-dir", "/tmp/x", "--smoke", "t"]).cell_pool == "site4"


def test_registered_block_marks_interim_and_corrects_the_loss_wording():
    m = _script()
    base = {"code": "H5M_FALSIFIED", "verdict": "x", "delta_tok_rule": "loss", "level_condition": True}
    b = m.registered_block({}, base, "site4", 200, 200)
    assert b["final"] and b["code"] == "H5M_FALSIFIED" and "level condition holds" in b["verdict"]
    b = m.registered_block({}, {**base, "level_condition": False}, "site4", 120, 200)
    assert not b["final"] and b["verdict"].startswith("INTERIM (120 of 200") and b["verdict"].endswith("x")
    b = m.registered_block({}, {"code": "PASS_SUFFICIENT", "verdict": "y", "level_condition": True}, "val", 3, 12)
    assert not b["final"] and b["verdict"].startswith("(smoke on val cells")
