"""Tests for the shared scoring rule in bridge_anm/lib/lineage_panels.py.

The answer key (expected_from_true), the TEDDY-alone rule
(lineage_scores_from_panel + CRITERIA readout_threshold) and the ANM event
values (event_value_under_criterion, run through the real ANM engine via
run_hard_proof.make_instance / anm_run_one) must rank lineages identically.
The legacy flag must reproduce the pre-fix functions exactly.

Run with ANM_ROOT pointing at the ANM checkout. The stored-cell test skips when
the gitignored bridge export or the pinned git commit is unavailable.
"""
from __future__ import annotations

import copy
import itertools
import json
import os
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
BRIDGE_PY = REPO / "bridge_anm"
MAIN_REPO = Path("/Users/tianchichen/Documents/GitHub/teddy_mm")
DATA = Path(os.environ.get("SCORING_BRIDGE_DIR", str(MAIN_REPO / "outputs/anm_cite_bridge")))
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))
# Last commit with the pre-fix scoring functions (reference for the legacy flag).
LEGACY_COMMIT = "ed390c2"

sys.path.insert(0, str(BRIDGE_PY))
os.environ.setdefault("ANM_ROOT", str(ANM_ROOT))

from lib import lineage_panels as lp  # noqa: E402

rhp = pytest.importorskip("run_hard_proof")

SCHEMA = json.loads((BRIDGE_PY / "schemas/cite_lineage_finite_field_v0.json").read_text())
PROTEINS = lp.all_panel_proteins()
P2L = lp.protein_to_lineage()
CRIT_IDS = ["O0", "O1", "O2"]
# Non-trivial train p95 so normalisation is exercised.
P95 = {p: 0.8 + 0.25 * i for i, p in enumerate(PROTEINS)}


@pytest.fixture(autouse=True)
def _restore_scoring():
    before = lp.get_scoring()
    lp.set_scoring("v2")
    yield
    lp.set_scoring(before)


def _events(panel_raw, p95, *, timing="simultaneous"):
    """Typed events as the exporter writes them: value = raw / train p95."""
    evs = []
    for i, p in enumerate(PROTEINS):
        lin = P2L[p]
        evs.append(
            {
                "event_id": f"syn:{p}",
                "time": 0 if timing == "simultaneous" else i,
                "protein": p,
                "action": lin,
                "modality": lp.MODALITY_FOR_LINEAGE[lin],
                "value": float(panel_raw[p]) / p95[p],
                "is_key_marker": p == lp.KEY_MARKERS[lin],
            }
        )
    return evs


def _anm(panel_raw, crit_id, p95=P95, timing="simultaneous"):
    crit = lp.CRITERIA[crit_id]
    inst = rhp.make_instance(
        "syn", {"adt_true_panel_holdout": panel_raw}, _events(panel_raw, p95, timing=timing), crit, p95
    )
    schema = copy.deepcopy(SCHEMA)
    schema["field_representation"]["readout_threshold"] = lp.anm_readout_threshold(crit, SCHEMA)
    return inst, rhp.anm_run_one(schema, inst)


def _three(crit_id, panel_raw, p95=P95):
    crit = lp.CRITERIA[crit_id]
    key = lp.expected_from_true(panel_raw, p95, crit)
    sc = lp.lineage_scores_from_panel(panel_raw, p95, crit)
    rule = rhp.decide_argmax(sc, crit["readout_threshold"])
    inst, out = _anm(panel_raw, crit_id, p95)
    return key, sc, rule, inst, out


# ---------------------------------------------------------------- basics
@pytest.mark.parametrize("cid", CRIT_IDS)
def test_event_values_bounded_and_weight_normalised(cid):
    crit = lp.CRITERIA[cid]
    w = lp.protein_weights(crit)
    assert max(w.values()) == pytest.approx(1.0)
    assert min(w.values()) >= 0.0
    for p in PROTEINS:
        for v in (0.0, 0.3, 1.0, 1.7):
            ev = lp.event_value_under_criterion(v, p == lp.KEY_MARKERS[P2L[p]], P2L[p], crit)
            assert 0.0 <= ev <= 1.0
            assert ev == pytest.approx(w[p] * min(v, 1.0))


def test_o2_priority_survives_weighting():
    crit = lp.CRITERIA["O2"]
    vals = {lin: lp.event_value_under_criterion(1.0, True, lin, crit) for lin in lp.LINEAGE_PANELS}
    assert vals["b_lineage"] == pytest.approx(1.0)
    assert vals["t_lineage"] == pytest.approx(1.3 / 1.5)
    assert vals["myeloid"] == pytest.approx(0.5 / 1.5)
    # legacy clipped all three to 1.0
    legacy = {lin: lp.event_value_under_criterion(1.0, True, lin, crit, scoring="legacy") for lin in vals}
    assert set(legacy.values()) == {1.0}


def test_thresholds_rescaled_consistently():
    assert lp.CRITERIA["O0"]["readout_threshold"] == pytest.approx(0.12)
    assert lp.CRITERIA["O1"]["readout_threshold"] == pytest.approx(0.28)
    assert lp.CRITERIA["O2"]["readout_threshold"] == pytest.approx(0.20 / 1.5)
    assert lp.CRITERIA["O2"]["expected_margin"] == pytest.approx(0.06 / 1.5)
    lp.set_scoring("legacy")
    assert lp.CRITERIA["O2"]["readout_threshold"] == pytest.approx(0.20)
    assert lp.CRITERIA["O2"]["expected_margin"] == pytest.approx(0.06)
    assert lp.rule_threshold(lp.CRITERIA["O2"], scoring="v2") == pytest.approx(0.20 / 1.5)


def test_field_gain_matches_engine():
    panel = {p: 0.0 for p in PROTEINS}
    panel["CD72"] = P95["CD72"]  # one unit event, secondary B marker
    _, out = _anm(panel, "O0")
    assert out["action_scores"]["b_lineage"] == pytest.approx(lp.field_gain(SCHEMA), rel=1e-9)


# ------------------------------------------------ CD19 = CD3 = CD16 = 1
def test_cd19_cd3_cd16_all_one_gives_b_under_o2():
    panel = {p: 0.0 for p in PROTEINS}
    for p in ("CD19", "CD3", "CD16"):
        panel[p] = P95[p]  # normalised evidence exactly 1
    key, sc, rule, inst, out = _three("O2", panel)
    assert key == "b_lineage"
    assert rule == "b_lineage"
    assert out["recommended_action"] == "b_lineage"
    assert out["Q_f"] == 1.0
    vals = {e["action"]: e["value"] for e in inst["proposed_source_events"] if e["value"] > 0}
    assert vals["b_lineage"] > vals["t_lineage"] > vals["myeloid"]
    s = out["action_scores"]
    assert s["b_lineage"] > s["t_lineage"] > s["myeloid"]


def test_cd19_cd3_cd16_legacy_bug_is_reproducible():
    """Legacy + panel-order timing: the three clipped key events are equal and ANM says myeloid."""
    lp.set_scoring("legacy")
    panel = {p: 0.0 for p in PROTEINS}
    for p in ("CD19", "CD3", "CD16"):
        panel[p] = P95[p]
    crit = lp.CRITERIA["O2"]
    inst = rhp.make_instance(
        "syn", {"adt_true_panel_holdout": panel}, _events(panel, P95, timing="panel-order"), crit, P95
    )
    out = rhp.anm_run_one(rhp.schema_for("O2", SCHEMA), inst)
    assert lp.expected_from_true(panel, P95, crit) == "b_lineage"
    assert out["recommended_action"] == "myeloid"


# ------------------------------------------------ grid agreement
def _grid_panels():
    levels = [0.0, 0.2, 0.5, 0.8, 1.0, 1.4]  # multiples of p95 (1.4 saturates)
    rng = np.random.default_rng(7)
    panels = []
    for kb, kt, km in itertools.product(levels, repeat=3):
        sec = rng.choice(levels, size=6)
        rel = {"CD19": kb, "CD3": kt, "CD16": km}
        it = iter(sec)
        for p in PROTEINS:
            if p not in rel:
                rel[p] = float(next(it))
        panels.append({p: rel[p] * P95[p] for p in PROTEINS})
    for _ in range(84):
        panels.append({p: float(rng.uniform(0, 1.3)) * P95[p] for p in PROTEINS})
    return panels  # 216 structured + 84 random


GRID = _grid_panels()


@pytest.mark.parametrize("cid", CRIT_IDS)
def test_key_rule_anm_agree_on_grid(cid):
    crit = lp.CRITERIA[cid]
    gain_c = lp.field_gain(SCHEMA) * lp.lineage_normaliser(crit)
    n_key = n_abstain_checked = 0
    for panel in GRID:
        key, sc, rule, _, out = _three(cid, panel)
        # ANM action scores are exactly G * c * S_l (same scoring, same scale)
        for lin in lp.LINEAGE_PANELS:
            assert out["action_scores"][lin] == pytest.approx(gain_c * sc[lin], abs=1e-12)
        ranked = sorted(sc.values(), reverse=True)
        top_gap = ranked[0] - ranked[1]
        near_thr = abs(ranked[0] - crit["readout_threshold"]) < 1e-9
        if top_gap > 1e-9 and not near_thr:
            assert out["recommended_action"] == rule
            n_abstain_checked += 1
        if key is not None:
            n_key += 1
            assert max(sc, key=sc.get) == key
            assert max(out["action_scores"], key=out["action_scores"].get) == key
            if rule is not None:
                assert rule == key
    assert n_key > 50 and n_abstain_checked > 250


# ------------------------------------------------ legacy on stored cells
def _load_legacy_module():
    try:
        src = subprocess.run(
            ["git", "-C", str(REPO), "show", f"{LEGACY_COMMIT}:bridge_anm/lib/lineage_panels.py"],
            check=True, capture_output=True, text=True,
        ).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    mod = types.ModuleType("lineage_panels_legacy_ref")
    exec(compile(src, "lineage_panels_legacy_ref.py", "exec"), mod.__dict__)
    return mod


def _stored_cells(n=200):
    cells = []
    with (DATA / "cite_cells_meta.jsonl").open() as f:
        for line in f:
            if line.strip():
                cells.append(json.loads(line))
            if len(cells) >= n:
                break
    want = {c["cell_id"] for c in cells}
    evs, p95 = {}, {}
    with (DATA / "cite_typed_events.jsonl").open() as f:
        for line in f:
            if not line.strip():
                continue
            ev = json.loads(line)
            if ev["cell_id"] in want:
                evs.setdefault(ev["cell_id"], []).append(ev)
                p95[ev["protein"]] = float(ev["norm_p95_train"])
            if len(evs) == len(want) and all(len(v) == len(PROTEINS) for v in evs.values()):
                break
    return cells, evs, p95


HAVE_DATA = (DATA / "cite_cells_meta.jsonl").exists() and (DATA / "cite_typed_events.jsonl").exists()


@pytest.mark.skipif(not HAVE_DATA, reason=f"bridge export not found under {DATA}")
@pytest.mark.parametrize("how", ["global_flag", "per_call"])
def test_legacy_flag_reproduces_old_values_on_200_stored_cells(how):
    old = _load_legacy_module()
    if old is None:
        pytest.skip(f"git commit {LEGACY_COMMIT} not available")
    cells, evs, p95 = _stored_cells(200)
    assert len(cells) == 200
    kw = {}
    if how == "global_flag":
        lp.set_scoring("legacy")
    else:
        kw = {"scoring": "legacy"}
    n_cmp = 0
    for cid in CRIT_IDS:
        ocrit, crit = old.CRITERIA[cid], lp.CRITERIA[cid]
        if how == "global_flag":
            assert crit["readout_threshold"] == ocrit["readout_threshold"]
            assert crit["expected_margin"] == ocrit["expected_margin"]
        for c in cells:
            t, pr = c["adt_true_panel_holdout"], c["adt_pred_panel"]
            assert lp.expected_from_true(t, p95, crit, **kw) == old.expected_from_true(t, p95, ocrit)
            assert lp.lineage_scores_from_panel(pr, p95, crit, **kw) == old.lineage_scores_from_panel(pr, p95, ocrit)
            assert lp.lineage_scores_from_panel(t, p95, crit, **kw) == old.lineage_scores_from_panel(t, p95, ocrit)
            for ev in evs[c["cell_id"]]:
                args = (float(ev["value"]), bool(ev["is_key_marker"]), ev["action"])
                assert lp.event_value_under_criterion(*args, crit, **kw) == old.event_value_under_criterion(*args, ocrit)
                n_cmp += 1
    assert n_cmp == 3 * 200 * len(PROTEINS)


@pytest.mark.skipif(not HAVE_DATA, reason=f"bridge export not found under {DATA}")
def test_v2_keeps_o0_o2_key_and_rule_decisions_on_stored_cells():
    """Rescaling margins/thresholds keeps the O0 and O2 key and rule decisions."""
    cells, _, p95 = _stored_cells(200)
    for cid in ("O0", "O2"):
        crit = lp.CRITERIA[cid]
        for c in cells:
            t, pr = c["adt_true_panel_holdout"], c["adt_pred_panel"]
            assert lp.expected_from_true(t, p95, crit) == lp.expected_from_true(t, p95, crit, scoring="legacy")
            new = rhp.decide_argmax(lp.lineage_scores_from_panel(pr, p95, crit), lp.rule_threshold(crit))
            leg = rhp.decide_argmax(
                lp.lineage_scores_from_panel(pr, p95, crit, scoring="legacy"),
                lp.rule_threshold(crit, scoring="legacy"),
            )
            assert new == leg
