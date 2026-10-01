"""Tests for bridge_anm/lib/v3_amend.py (amendment A1 to the v3 registration).

Synthetic tests always run; tests that need the ANM engine skip without ANM_ROOT, and tests
that need the built amendment skip when registration/amendment_A1.json is absent.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "bridge_anm"))
sys.path.insert(0, str(REPO))

from lib import v3_amend as va  # noqa: E402
from lib import v3_key as vk  # noqa: E402

REG_PATH = REPO / "registration" / "registration_v3.json"
AMEND_PATH = REPO / "registration" / "amendment_A1.json"
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))


@pytest.fixture(scope="module")
def reg():
    if not REG_PATH.exists():
        pytest.skip("registration not built")
    return vk.load_registration(REG_PATH)


@pytest.fixture(scope="module")
def amended():
    if not AMEND_PATH.exists():
        pytest.skip("amendment not built")
    return va.load_registration_amended()


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
