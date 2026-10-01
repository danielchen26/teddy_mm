"""Tests for bridge_anm/lib/v3_amend.py (amendments A1 and A2 to the v3 registration).

Synthetic tests always run; tests that need the ANM engine skip without ANM_ROOT, tests that need
the built amendments skip when registration/amendment_A1.json or amendment_A2.json is absent, and
tests that need the data pack skip without data/processed/cite.
"""
from __future__ import annotations

import copy
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
sys.path.insert(0, str(REPO / "bridge_anm"))
sys.path.insert(0, str(REPO))

from lib import v3_amend as va  # noqa: E402
from lib import v3_key as vk  # noqa: E402

REG_PATH = REPO / "registration" / "registration_v3.json"
AMEND_PATH = REPO / "registration" / "amendment_A1.json"
A2_PATH = REPO / "registration" / "amendment_A2.json"
CITE = REPO / "data" / "processed" / "cite" / "cite_arrays.npz"
REG_FILES = ("registration_v3.json", "registration_v3.json.sha256", "amendment_A1.json", "amendment_A1.json.sha256",
             "amendment_A2.json", "amendment_A2.json.sha256")
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))


@pytest.fixture(scope="module")
def reg():
    if not REG_PATH.exists():
        pytest.skip("registration not built")
    return vk.load_registration(REG_PATH)


@pytest.fixture(scope="module")
def amended():
    if not AMEND_PATH.exists() or not A2_PATH.exists():
        pytest.skip("amendment not built")
    return va.load_registration_amended()


@pytest.fixture(scope="module")
def a2():
    if not A2_PATH.exists():
        pytest.skip("amendment A2 not built")
    return va.load_amendment_A2()


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _copy_registration(dst: Path) -> None:
    for f in REG_FILES:
        shutil.copy(REPO / "registration" / f, dst / f)


def _clr_per_cell(x):
    """The per-cell CLR form of the stored ADT: log1p(x / g_c), g_c = exp(mean_p log1p(x_p))."""
    g = np.exp(np.log1p(x).mean(axis=1, keepdims=True))
    return np.log1p(x / g)


# ----------------------------------------------------------------------------- A1.2 channel 2
def test_channel2_inputs_do_not_depend_on_panel_counts():
    rng = np.random.default_rng(0)
    names = [f"P{i}" for i in range(30)]
    panel = names[:6]
    x = rng.poisson(rng.gamma(2.0, 5.0, size=(200, 30))).astype(np.float64)
    x[:, 10] += 1  # every cell has a non-panel count
    w1, cols = va.channel2_inputs(_clr_per_cell(x), names, panel)
    x2 = x.copy()
    x2[:, :6] = rng.poisson(50.0, size=(200, 6))  # change only the panel counts
    w2, _ = va.channel2_inputs(_clr_per_cell(x2), names, panel)
    assert cols == names[6:]
    np.testing.assert_allclose(w1, w2, rtol=1e-9, atol=1e-9)
    # equals log1p(1e4 * x_q / sum_r x_r) over the non-panel counts
    xn = x[:, 6:]
    np.testing.assert_allclose(w1, np.log1p(1e4 * xn / xn.sum(axis=1, keepdims=True)), rtol=1e-9, atol=1e-9)
    # while the stored non-panel values do change with the panel counts (the A1.2 finding)
    assert not np.allclose(_clr_per_cell(x)[:, 6:], _clr_per_cell(x2)[:, 6:])


# ----------------------------------------------------------------------------- A1.1 matched coverage
def test_coverage_select_ties_by_permutation_not_index():
    amend = {"matched_coverage": {"tie_break": {"seed": 29, "n_cells": 1000}}}
    ids = np.arange(100, 300)
    conf = np.ones(ids.size)          # all tied
    conf[:20] = 2.0                   # 20 clear winners
    m = va.coverage_select(conf, ids, 0.5, amend)
    assert m.sum() == 100 and m[:20].all()
    # the tied part is not "lowest index first"
    assert not m[20:100].all()
    # deterministic and invariant to the order in which cells are passed
    perm = np.random.default_rng(1).permutation(ids.size)
    m2 = va.coverage_select(conf[perm], ids[perm], 0.5, amend)
    assert set(ids[perm][m2]) == set(ids[m])
    # ceil rule
    assert va.coverage_select(conf, ids, 0.501, amend).sum() == int(np.ceil(0.501 * ids.size))


def test_tie_break_rank_is_a_permutation():
    amend = {"matched_coverage": {"tie_break": {"seed": 29, "n_cells": 50}}}
    r = va.tie_break_rank(np.arange(50), amend)
    assert sorted(r.tolist()) == list(range(50))


# ----------------------------------------------------------------------------- A1.3 leave-one-out
def test_top_marker_sets():
    v = np.array([[1.0, 1.0, 0.2], [0.3, 0.9, 0.1]])
    s = va.top_marker_sets(v)
    assert s[0].tolist() == [0, 1] and s[1].tolist() == [1]


def _anm():
    sys.path.insert(0, str(ANM_ROOT))
    try:
        from active_neural_matter.field import finite_field_runner as ffr
    except Exception:
        pytest.skip("ANM engine not importable (set ANM_ROOT)")
    return ffr


def _engine_scores(ffr, schema, vals, drop=None, zero=None):
    actions = [{"id": k, "label": k} for k in vk.LINEAGES]
    events = []
    for i, k in enumerate(vk.LINEAGES):
        for jj in range(vals.shape[1]):
            if drop == (i, jj):
                continue
            value = 0.0 if zero == (i, jj) else float(vals[i, jj])
            events.append({"event_id": f"{k}{jj}", "time": 0, "action": k, "modality": "evidence",
                           "polarity": "support", "value": value, "provenance": "test"})
    inst = {"instance_id": "x", "question": "q", "actions": actions, "proposed_source_events": events}
    fe = ffr.validate_source_events(schema, inst)["field_events"]
    st = ffr.evolve_field(schema, ffr.build_graph(inst, fe), fe)["state"]
    return np.array([st[f"action:{k}"] for k in vk.LINEAGES])


def test_loo_flip_exact_matches_engine_deletion(reg):
    ffr = _anm()
    schema = json.loads((REPO / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
    fld = reg["anm"]["field_representation"]
    schema["field_representation"].update({k: v for k, v in fld.items() if k != "kind"})
    schema["source_event_schema"]["allowed_modalities"] = ["evidence"]
    n = 3
    G = vk.field_gain(n, fld)
    rng = np.random.default_rng(4)
    n_flip = 0
    for _ in range(150):
        vals = rng.uniform(0, 1, size=(4, n))
        vals[rng.integers(4)] *= 0.4
        S = vals.mean(axis=1)
        bar = float(rng.uniform(0.25, 0.6))
        thr = G * n * bar
        A = _engine_scores(ffr, schema, vals)
        if A.max() < thr:
            continue
        c = int(A.argmax())
        p = int(vals[c].argmax())
        A_del = _engine_scores(ffr, schema, vals, drop=(c, p))
        call_del = int(A_del.argmax()) if A_del.max() >= thr else -1
        flip_engine = call_del != c
        flip_exact = bool(va.loo_flip_exact(S[None], vals[c][None], np.array([c]), bar, fld)[0])
        assert flip_engine == flip_exact
        n_flip += flip_engine
        A_zero = _engine_scores(ffr, schema, vals, zero=(c, p))
        call_zero = int(A_zero.argmax()) if A_zero.max() >= thr else -1
        flip_zero = bool(va.loo_flip_registered_form(S[None], vals[c][None], np.array([c]), bar)[0])
        assert (call_zero != c) == flip_zero  # the old form is exact for zeroing
    assert n_flip > 0


def test_loo_ratio(reg):
    fld = reg["anm"]["field_representation"]
    assert va.loo_ratio(3, fld) == pytest.approx(vk.field_gain(2, fld) / vk.field_gain(3, fld))
    assert va.loo_ratio(1, fld) == 1.0


# ----------------------------------------------------------------------------- loading and overrides
def test_amendment_hashes_and_overrides(amended, tmp_path):
    a = amended["amendment_A1"]
    assert a["amends"]["registration_sha256"] == vk.sha256_file(REG_PATH)
    e = amended["experiments"]
    assert "panel-free" in e["E4"]["channels"]["channel2"]
    assert "rho" in e["E1"]["exp3_why_this_call"]["rule"]
    assert e["common"]["matched_coverage"]["tie_break"]["seed"] == a["matched_coverage"]["tie_break"]["seed"]
    assert "within each donor" in e["E3"]["cells"]
    # tamper detection
    bad = tmp_path / "amendment_A1.json"
    shutil.copy(AMEND_PATH, bad)
    shutil.copy(AMEND_PATH.with_name("amendment_A1.json.sha256"), tmp_path / "amendment_A1.json.sha256")
    bad.write_text(bad.read_text().replace('"A1.1"', '"A1.1x"', 1))
    with pytest.raises(ValueError):
        va.load_amendment(bad)
    # every override path exists in the frozen experiments (no silent typo creates a new branch)
    reg = vk.load_registration(REG_PATH)
    for path in a["experiments_overrides"]:
        cur = reg["experiments"]
        for k in path.split("/")[:-1]:
            assert k in cur, path
            cur = cur[k]


# ============================================================================= amendment A2
# ----------------------------------------------------------------------------- loading
def test_loader_applies_A1_then_A2(amended, a2):
    assert [x["amendment_id"] for x in amended["amendments"]] == ["teddy_mm_v3_A1", "teddy_mm_v3_A2"]
    assert amended["amendment_A2"] == a2
    assert a2["amends"]["registration_sha256"] == vk.sha256_file(REG_PATH)
    assert a2["amends"]["amendment_A1_sha256"] == vk.sha256_file(AMEND_PATH)
    e = amended["experiments"]
    assert "A2.1" in e["E3"]["endpoints"]["E3.H3a"] and "scale artefact" in e["E3"]["endpoints"]["E3.H3a"]
    assert "A2 win" in e["E3"]["falsification"]
    assert "A2.2" in e["E3"]["gap_ratio"]
    assert e["E2"]["falsification"].startswith(amended["amendment_A1"]["experiments_overrides"]["E2/falsification"])
    assert "each primary donor" in e["E2"]["falsification"]
    assert "within each donor" in e["E3"]["cells"]          # A1 stays applied underneath A2
    assert e["E3"]["endpoints"]["E3.H3b"] == vk.load_registration(REG_PATH)["experiments"]["E3"]["endpoints"]["E3.H3b"]
    # A2 only replaces existing fields (no typo can create a new branch)
    reg = vk.load_registration(REG_PATH)
    for path in a2["experiments_overrides"]:
        cur = reg["experiments"]
        for k in path.split("/"):
            assert k in cur, path
            cur = cur[k]
    assert a2["leakage_check"]["passed"] is True
    assert a2["leakage_check"]["poison_test"]["identical_to_real"] is True
    assert a2["leakage_check"]["poison_train"]["identical_to_real"] is False


def test_A2_tamper_chain_and_missing_file(tmp_path, a2):
    _copy_registration(tmp_path)
    reg_p, a1_p, a2_p = tmp_path / "registration_v3.json", tmp_path / "amendment_A1.json", tmp_path / "amendment_A2.json"
    assert va.load_registration_amended(reg_p, a1_p)["amendment_A2"]["amendment_id"] == va.A2_ID
    txt = a2_p.read_text()
    a2_p.write_text(txt.replace('"A2.1"', '"A2.1x"', 1))           # altered A2, old hash file
    with pytest.raises(ValueError):
        va.load_registration_amended(reg_p, a1_p)
    d = json.loads(txt)
    d["amends"]["amendment_A1_sha256"] = "0" * 64                  # self-consistent A2 naming another A1
    a2_p.write_text(json.dumps(d))
    (tmp_path / "amendment_A2.json.sha256").write_text(f"{vk.sha256_file(a2_p)}  amendment_A2.json\n")
    with pytest.raises(ValueError, match="A1"):
        va.load_registration_amended(reg_p, a1_p)
    a2_p.unlink()
    with pytest.raises(FileNotFoundError):
        va.load_registration_amended(reg_p, a1_p)


def test_amendment_A2_status_needs_commit_and_match(tmp_path, a2):
    _copy_registration(tmp_path)
    assert va.amendment_A2_status(tmp_path, lambda p: True)["ok"]
    assert not va.amendment_A2_status(tmp_path, lambda p: p.name != "amendment_A2.json")["ok"]
    assert not va.amendment_A2_status(tmp_path, lambda p: p.name != "amendment_A2.json.sha256")["ok"]
    assert not va.amendment_A2_status(tmp_path, lambda p: None)["ok"]
    (tmp_path / "amendment_A1.json").write_text((tmp_path / "amendment_A1.json").read_text() + " ")
    assert not va.amendment_A2_status(tmp_path, lambda p: True)["ok"]          # A2 names another A1
    _copy_registration(tmp_path)
    (tmp_path / "amendment_A2.json.sha256").write_text("0" * 64 + "  amendment_A2.json\n")
    assert not va.amendment_A2_status(tmp_path, lambda p: True)["ok"]
    (tmp_path / "amendment_A2.json").unlink()
    st = va.amendment_A2_status(tmp_path, lambda p: True)
    assert not st["ok"] and not st["amendment_A2_exists"]


def test_e3_site4_guard_requires_A2(monkeypatch, a2):
    s = _load("v3_e3_nkt_repair_a2test", "scripts/v3_e3_nkt_repair.py")
    a = SimpleNamespace(registration_dir=REPO / "registration")
    monkeypatch.setattr(s, "_git_committed", lambda p: True)
    monkeypatch.setattr(s, "_git_head_text", lambda p: Path(p).read_text())
    info = s.check_registration(a, site4=True)                     # everything committed and matching: allowed
    assert info["amendment_A2"]["ok"] and info["amendment_A2_sha256"] == vk.sha256_file(A2_PATH)
    monkeypatch.setattr(s, "_git_committed", lambda p: Path(p).name != "amendment_A2.json")
    with pytest.raises(SystemExit, match="amendment_A2"):
        s.check_registration(a, site4=True)
    assert s.check_registration(a, site4=False)["amendment_A2"]["ok"] is False   # train/val stages still run


def test_e2_site4_guard_requires_A2(monkeypatch, a2):
    s = _load("e2_script_a2test", "scripts/e2_response_decomposition.py")
    a = SimpleNamespace(registration_dir=REPO / "registration", cell_pool="e2_subset", stage="report")
    monkeypatch.setattr(s, "_git_committed", lambda p: True)
    rc = s.check_registration(a)
    assert rc["info"]["amendment_A2"]["ok"] and "amendment_A2" in rc["reg"]
    monkeypatch.setattr(s, "_git_committed", lambda p: Path(p).name != "amendment_A2.json")
    for stage in ("embed", "report", "all"):
        with pytest.raises(SystemExit, match="amendment_A2"):
            s.check_registration(SimpleNamespace(registration_dir=a.registration_dir, cell_pool="e2_subset", stage=stage))
    s.check_registration(SimpleNamespace(registration_dir=a.registration_dir, cell_pool="val", stage="report"))  # smoke


def test_A2_leaves_addendum_parts_of_scripts_unchanged():
    """Stages already run were gated on the addenda; A2 must not change DECLARED (E3) or SPEC (E2)."""
    e3f, e2f = REPO / "registration/addenda/E3.json", REPO / "registration/addenda/E2.json"
    if not (e3f.exists() and e2f.exists()):
        pytest.skip("addenda not written")
    s3 = _load("v3_e3_nkt_repair_a2decl", "scripts/v3_e3_nkt_repair.py")
    s2 = _load("e2_script_a2spec", "scripts/e2_response_decomposition.py")
    assert json.loads(e3f.read_text())["declared"] == s3.DECLARED
    add2 = json.loads(e2f.read_text())
    assert add2["spec"] == s2.SPEC and add2["spec_sha256"] == s2.spec_sha()


# ----------------------------------------------------------------------------- A2.2 R2 normaliser sample
def test_a2_draw_is_label_free_seeded_and_train_only():
    split = np.array(["train", "val", "test"] * 400 + ["train"] * 300)
    s1 = va.a2_draw_r2_q95_sample(split, 31, 200)
    assert s1.size == 200 and np.all(split[s1] == "train") and np.all(np.diff(s1) > 0)
    assert np.array_equal(s1, va.a2_draw_r2_q95_sample(split.copy(), 31, 200))
    assert not np.array_equal(s1, va.a2_draw_r2_q95_sample(split, 32, 200))
    # only the positions of split == train matter: renaming the other splits changes nothing
    other = np.where(split == "train", "train", "x")
    assert np.array_equal(s1, va.a2_draw_r2_q95_sample(other, 31, 200))


def test_a2_registered_sample_reproduces(a2):
    if not CITE.exists():
        pytest.skip("data pack not available")
    split = np.load(CITE, allow_pickle=False)["split"].astype(str)
    s = va.a2_r2_q95_sample(split, a2)
    spec = a2["e3"]["r2_normaliser"]
    assert s.size == spec["n_cells"] == 10_000 and np.all(split[s] == "train")
    assert vk.index_hash(np.where(split == "train")[0]) == spec["population_sha256"]
    bad = copy.deepcopy(a2)
    bad["computed"]["e3_r2_normaliser_sample"]["cells"][0] += 1
    with pytest.raises(ValueError):
        va.a2_r2_q95_sample(split, bad)
    bad = copy.deepcopy(a2)
    bad["e3"]["r2_normaliser"]["seed"] = 32
    with pytest.raises(ValueError):
        va.a2_r2_q95_sample(split, bad)


def test_e3_predict_plan_uses_A2_sample(a2, amended):
    if not CITE.exists():
        pytest.skip("data pack not available")
    s = _load("v3_e3_nkt_repair_a2plan", "scripts/v3_e3_nkt_repair.py")
    split = np.load(CITE, allow_pickle=False)["split"].astype(str)
    D = SimpleNamespace(split=split, rows=lambda name: np.where(split == name)[0])
    plan = dict(s.predict_plan(SimpleNamespace(smoke=None), D, amended))
    assert set(plan) == {"val", "trainA2", "test"}
    assert np.array_equal(plan["trainA2"], va.a2_r2_q95_sample(split, a2))
    assert plan["val"].size == int(np.sum(split == "val")) and plan["test"].size == int(np.sum(split == "test"))


# ----------------------------------------------------------------------------- A2.1 E3 H3a
def test_a2_positive_condition_and_h3a_verdict():
    d2 = {"13272": 0.05, "19593": 0.2}
    assert va.a2_positive_condition(0.1, 0.02, d2)["holds"]
    assert not va.a2_positive_condition(0.1, -0.01, d2)["holds"]
    assert not va.a2_positive_condition(0.1, 0.0, d2)["holds"]
    assert not va.a2_positive_condition(-0.1, 0.02, d2)["holds"]
    assert not va.a2_positive_condition(0.1, 0.02, {"13272": -0.01, "19593": 0.2})["holds"]
    assert not va.a2_positive_condition(0.1, 0.02, {"13272": float("nan"), "19593": 0.2})["holds"]
    assert not va.a2_positive_condition(0.1, 0.02, {"13272": None, "19593": 0.2})["holds"]
    assert not va.a2_positive_condition(float("nan"), 0.02, d2)["holds"]
    assert not va.a2_positive_condition(0.1, None, d2)["holds"]
    assert not va.a2_positive_condition(0.1, 0.02, {})["holds"]
    assert va.a2_h3a_verdict("win", True, True) == "win"
    assert va.a2_h3a_verdict("win", False, True) == va.A2_SCALE_ARTEFACT == "not supported (scale artefact)"
    assert va.a2_h3a_verdict("win", True, False) == va.A2_SCALE_ARTEFACT
    for v in ("loss", "equivalent", "inconclusive"):
        assert va.a2_h3a_verdict(v, True, True) == v
    assert va.a2_e3_falsification("win", "inconclusive")["rejected"] is False
    assert va.a2_e3_falsification(va.A2_SCALE_ARTEFACT, "win")["rejected"] is False
    assert va.a2_e3_falsification(va.A2_SCALE_ARTEFACT, "inconclusive")["rejected"] is True


def _synthetic_pairs(keep: dict[str, tuple[float, float]], n_per_donor: int = 120, seed: int = 0):
    """Pairs of two donors, half flagged; measured NK-T difference ~1 on the 4 gap proteins; readout r keeps the
    fraction keep[r] = (flagged, unflagged) of the measured difference (plus small noise)."""
    rng = np.random.default_rng(seed)
    n = 2 * n_per_donor
    nk, t = np.arange(n), np.arange(n, 2 * n)
    donor = np.repeat(["13272", "19593"], n_per_donor)
    flag = np.tile(np.arange(n_per_donor) % 2 == 0, 2)
    meas = np.zeros((2 * n, 13))
    meas[nk] = 1.0 + 0.05 * rng.standard_normal((n, 13))
    meas[t] = 0.05 * np.abs(rng.standard_normal((n, 13)))
    ev = {}
    for r, (kf, ku) in keep.items():
        e = np.zeros((2 * n, 13))
        k = np.where(flag, kf, ku)[:, None]
        e[nk] = meas[t] + k * (meas[nk] - meas[t]) * (1.0 + 0.02 * rng.standard_normal((n, 13)))
        e[t] = meas[t]
        ev[r] = e
    return {"nk": nk, "t": t, "flag": flag, "donor": donor}, ev, meas


def test_e3_run_h3a_flags_uniform_shrinkage_as_scale_artefact():
    """End to end through the E3 evaluate code: R1 shrinks every difference by 1/2 (registered D wins, A2 does not),
    R2 recovers the gap on flagged pairs (an A2 win)."""
    s = _load("v3_e3_nkt_repair_a2h3a", "scripts/v3_e3_nkt_repair.py")
    head = (0.3, 0.8)
    pv, ev, meas = _synthetic_pairs({"head": head, "null": head, "R1": (0.15, 0.4), "R2": (0.9, 0.8)})
    gidx = [s.DECLARED["targets"].index(p) for p in ("CD56", "CD94", "CD335", "CD3")]
    out = s.run_h3a(pv, ev, meas, gidx, ["13272", "19593"], 300, 1)
    st, a2 = out["stats"], out["A2"]
    assert st["D_R1"]["verdict"] == "win" and st["D_R2"]["verdict"] == "win"
    assert a2["R1"]["registered_D_verdict"] == "win"
    assert a2["R1"]["b_flagged_recovery"]["holds"] is False
    assert a2["R1"]["verdict"] == va.A2_SCALE_ARTEFACT
    assert a2["R2"]["b_flagged_recovery"]["holds"] and a2["R2"]["c_scale_free_Dlog"]["holds"]
    assert a2["R2"]["verdict"] == "win"
    h3b = {"stats": {f"{R}-head@{c}": {"verdict": "inconclusive"} for c in s.COVERAGES for R in ("R1", "R2")}}
    v = s.verdicts({"H3a": {"sp": {"all": out}}, "H3b": {"sp": {"all": h3b}}}, "sp")
    assert v["E3.H3a_R1"] == va.A2_SCALE_ARTEFACT and v["E3.H3a_R2"] == "win"
    assert v["E3.H3a_registered_D_only"] == {"R1": "win", "R2": "win"}
    assert v["falsification"].startswith("not rejected")
    # with only the shrinking readouts, the registered reading would not reject, A2 does
    pv, ev, meas = _synthetic_pairs({"head": head, "null": head, "R1": (0.15, 0.4), "R2": (0.21, 0.56)}, seed=1)
    out = s.run_h3a(pv, ev, meas, gidx, ["13272", "19593"], 300, 1)
    v = s.verdicts({"H3a": {"sp": {"all": out}}, "H3b": {"sp": {"all": h3b}}}, "sp")
    assert v["E3.H3a_registered_D_only"] == {"R1": "win", "R2": "win"}
    assert v["E3.H3a_R1"] == v["E3.H3a_R2"] == va.A2_SCALE_ARTEFACT
    assert v["falsification"].startswith("rejected") and v["falsification_registered_D_only_reading"].startswith("not rejected")


def test_e3_a2_needs_each_primary_donor():
    """A repair present in one donor only: D can still win on the pooled pairs if the other donor's D clears
    margin/2, but (b) must be positive in each donor."""
    s = _load("v3_e3_nkt_repair_a2donor", "scripts/v3_e3_nkt_repair.py")
    head = (0.3, 0.8)
    pv, ev, meas = _synthetic_pairs({"head": head, "null": head, "R1": (0.15, 0.4), "R2": (0.9, 0.8)}, seed=2)
    # in donor 19593, R2 is a uniform 1/2 shrink of the head (no recovery on flagged pairs there)
    m = pv["donor"] == "19593"
    for i_nk, i_t, f in zip(pv["nk"][m], pv["t"][m], pv["flag"][m]):
        k = 0.15 if f else 0.4
        ev["R2"][i_nk] = meas[i_t] + k * (meas[i_nk] - meas[i_t])
    gidx = [s.DECLARED["targets"].index(p) for p in ("CD56", "CD94", "CD335", "CD3")]
    out = s.run_h3a(pv, ev, meas, gidx, ["13272", "19593"], 300, 1)
    a2 = out["A2"]["R2"]
    assert a2["b_flagged_recovery"]["per_donor"]["19593"] < 0
    assert a2["b_flagged_recovery"]["holds"] is False
    assert a2["verdict"] != "win"


# ----------------------------------------------------------------------------- A2.3 E2 decision rule
def test_a2_e2_criterion_i_rows():
    base = {"unit_a": "thin_keep0.2", "unit_b": "gene_eps1", "label": "head", "diff": 0.2, "ci95": [0.05, 0.35],
            "passes": True, "per_donor_diff": {"13272": 0.18, "19593": 0.11}}
    r = va.a2_e2_criterion_i([base])
    assert r["holds"] and r["n_passing_A2"] == 1
    for bad in ({"per_donor_diff": {"13272": 0.3, "19593": 0.05}},     # below margin / 2 in one donor
                {"per_donor_diff": {"13272": 0.3, "19593": -0.1}},     # opposite sign in one donor
                {"per_donor_diff": {"13272": 0.3, "19593": None}},     # undefined in one donor
                {"per_donor_diff": {"13272": 0.3, "19593": float("nan")}},
                {"per_donor_diff": {}},
                {"ci95": [-0.35, -0.05]},                                # interval on the other side
                {"ci95": [None, None]},
                {"passes": False}):
        assert not va.a2_e2_criterion_i([{**base, **bad}])["holds"], bad
    neg = {**base, "diff": -0.2, "ci95": [-0.35, -0.05], "per_donor_diff": {"13272": -0.1, "19593": -0.3}}
    assert va.a2_e2_criterion_i([neg])["holds"]
    assert va.a2_e2_criterion_i([{**base, "passes": False}, neg])["n_passing_A2"] == 1


def test_a2_e2_criterion_i_end_to_end_with_criterion_i():
    """er.criterion_i passes on the pooled cells (A1 reading); the donor clause decides under A2."""
    from lib import e2_response as er

    rng = np.random.default_rng(0)

    def unit(loss, n_lost, sh):
        return {"decision_loss": loss, "n_lost_cells": n_lost,
                **{f"share_{L}": v for L, v in zip(er.LABELS, sh)}}

    units = {"A": unit(0.10, 50, (0.6, 0.3, 0.1)), "B": unit(0.11, 60, (0.3, 0.5, 0.2))}
    boot = {u: {f"share_{L}": units[u][f"share_{L}"] + 0.02 * rng.standard_normal(500) for L in er.LABELS} for u in units}

    def per_donor(a_19593):
        return {"13272": {"A": unit(0.1, 25, (0.7, 0.2, 0.1)), "B": unit(0.1, 30, (0.3, 0.5, 0.2))},
                "19593": {"A": unit(0.1, 25, a_19593), "B": unit(0.1, 30, (0.3, 0.5, 0.2))}}

    # donor 19593: representation diff 0.05 and head diff -0.05, both below margin / 2 = 0.075
    weak = er.criterion_i(units, boot, per_donor=per_donor((0.35, 0.45, 0.2)))
    assert weak["holds"] is True and weak["holds_with_per_donor_clause"] is False
    a2w = va.a2_e2_criterion_i(weak["rows"])
    assert a2w["holds"] is False and a2w["n_passing_overall"] == 2 and a2w["n_passing_A2"] == 0
    strong = er.criterion_i(units, boot, per_donor=per_donor((0.5, 0.4, 0.1)))   # 19593: 0.2 and -0.1
    assert va.a2_e2_criterion_i(strong["rows"])["holds"] is True
    assert strong["holds_with_per_donor_clause"] is True
    assert va.a2_e2_verdict(False, False) == "adds nothing beyond the curve (falsified)"
    assert va.a2_e2_verdict(True, False) == va.a2_e2_verdict(False, True) == "adds information beyond the curve"


def test_a2_e2_criterion_ii():
    row = {"point": 0.03, "ci95": [0.005, 0.05], "per_donor": {"13272": 0.02, "19593": 0.012}}
    assert va.a2_e2_criterion_ii(row, 0.02)["holds"]
    assert not va.a2_e2_criterion_ii({**row, "per_donor": {"13272": 0.02, "19593": 0.009}}, 0.02)["holds"]
    assert not va.a2_e2_criterion_ii({**row, "per_donor": {"13272": 0.02, "19593": None}}, 0.02)["holds"]
    assert not va.a2_e2_criterion_ii({**row, "ci95": [0.0, 0.05]}, 0.02)["holds"]
    assert not va.a2_e2_criterion_ii({**row, "point": 0.019}, 0.02)["holds"]
    assert not va.a2_e2_criterion_ii({**row, "ci95": [None, None]}, 0.02)["holds"]
    # agrees with the registered win rule of e2_response
    from lib import e2_response as er
    for pt, lo, d1, d2 in [(0.03, 0.005, 0.02, 0.012), (0.03, 0.005, 0.02, 0.009), (0.05, -0.01, 0.03, 0.03)]:
        r = {"point": pt, "ci95": [lo, 0.1], "per_donor": {"a": d1, "b": d2}}
        assert va.a2_e2_criterion_ii(r, 0.02)["holds"] == (er.win_rule(pt, [lo, 0.1], r["per_donor"], 0.02) == "win")
