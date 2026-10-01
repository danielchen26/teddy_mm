"""Tests for bridge_anm/lib/v3_key.py (v3 answer key, evidence, questions) and the frozen registration.

Synthetic tests always run. Tests that read the real data or the ANM engine skip when the
gitignored data pack (data/processed/cite) or ANM_ROOT is unavailable.
"""
from __future__ import annotations

import copy
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

from lib import v3_key as vk  # noqa: E402

REG_PATH = REPO / "registration" / "registration_v3.json"
DATA = REPO / "data" / "processed" / "cite" / "cite_arrays.npz"
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))


@pytest.fixture(scope="module")
def reg():
    if not REG_PATH.exists():
        pytest.skip("registration not built")
    return vk.load_registration(REG_PATH)


def _gate_spec(reg):
    return {k: reg["gate"][k] for k in ("b_markers", "nk_markers", "nk_min_high")}


def _H(n, **high):
    names = ["CD3", "CD19", "CD20", "CD56", "CD94", "CD335", "CD16", "CD14", "CD33", "CD11c", "CD71"]
    H = {p: np.zeros(n, dtype=bool) for p in names}
    for p, rows in high.items():
        H[p][list(rows)] = True
    return H


# ----------------------------------------------------------------------------- registration
def test_registration_hash_verified_and_tamper_detected(reg, tmp_path):
    assert reg["registration_id"] == "teddy_mm_v3"
    dst = tmp_path / "registration_v3.json"
    shutil.copy(REG_PATH, dst)
    shutil.copy(REG_PATH.with_name(REG_PATH.name + ".sha256"), tmp_path / "registration_v3.json.sha256")
    vk.load_registration(dst)  # unchanged copy verifies
    dst.write_text(dst.read_text().replace('"teddy_mm_v3"', '"teddy_mm_v3x"', 1))
    with pytest.raises(ValueError):
        vk.load_registration(dst)


def test_annotation_map_complete_and_borderline(reg):
    amap = reg["annotation_map"]
    assert len(amap) == 45
    assert set(amap.values()) == set(vk.CLASSES)
    assert amap["gdT CD158b+"] == "T" and amap["gdT TCRVD2+"] == "T" and amap["MAIT"] == "T" and amap["T reg"] == "T"
    assert amap["NK"] == "NK" and amap["NK CD158e1+"] == "NK"
    for t in ("ILC", "ILC1", "pDC", "HSC", "Lymph prog", "G/M prog", "MK/E prog", "Proerythroblast", "Erythroblast",
              "Normoblast", "Reticulocyte", "Plasma cell IGKC+", "Plasma cell IGKC-", "Plasmablast IGKC+",
              "Plasmablast IGKC-"):
        assert amap[t] == "OUT", t
    for t in ("CD14+ Mono", "CD16+ Mono", "cDC2"):
        assert amap[t] == "myeloid"
    with pytest.raises(KeyError):
        vk.annotation_class(["not a type"], reg)


def test_absent_proteins_recorded(reg):
    assert set(reg["data"]["adt_absent_checked"]) == {"CD34", "CD235a", "CD138", "TCRgd"}
    assert reg["data"]["n_adt"] == 134


# ----------------------------------------------------------------------------- gate
def test_gate_rules_synthetic(reg):
    spec = _gate_spec(reg)
    nk = spec["nk_markers"]
    # 0 T; 1 B; 2 NK (two NK markers); 3 myeloid (CD14); 4 erythroid CD71 + CD3 -> OUT; 5 nothing -> OUT;
    # 6 CD3 + CD19 -> OUT (neither rule); 7 one NK marker only -> OUT; 8 NK markers + CD33 -> myeloid
    H = _H(9, CD3=[0, 4, 6], CD19=[1, 6], CD14=[3], CD71=[4], CD33=[8])
    for p in nk[:2]:
        H[p][[2, 8]] = True
    H[nk[0]][7] = True
    g = vk.gate_class_from_high(H, spec)
    assert list(g) == ["T", "B", "NK", "myeloid", "OUT", "OUT", "OUT", "OUT", "myeloid"]


def test_gate_class_uses_strict_threshold(reg):
    names = reg["gate"]["proteins"]
    thr = reg["gate"]["thresholds"]
    adt = np.zeros((2, len(names)), dtype=np.float64)
    adt[0, names.index("CD3")] = thr["CD3"]          # equal -> not high
    adt[1, names.index("CD3")] = thr["CD3"] + 1e-3   # above -> T
    g = vk.gate_class(adt, names, reg)
    assert list(g) == ["OUT", "T"]


def test_primary_key_agreement():
    a = np.array(["T", "NK", "OUT", "B", "myeloid"])
    g = np.array(["T", "T", "OUT", "OUT", "myeloid"])
    assert list(vk.primary_key(a, g)) == ["T", vk.UNSCORED, "OUT", vk.UNSCORED, "myeloid"]


def test_q3_key(reg):
    names = reg["gate"]["proteins"]
    thr = reg["gate"]["thresholds"]["CD14"]
    adt = np.zeros((5, len(names)), dtype=np.float32)
    adt[[0, 2], names.index("CD14")] = thr + 1.0
    ct = ["CD14+ Mono", "CD16+ Mono", "CD16+ Mono", "CD14+ Mono", "NK"]
    prim = np.array(["myeloid", "myeloid", "myeloid", "myeloid", "NK"])
    q3 = vk.q3_key(prim, ct, adt, names, reg)
    assert list(q3) == ["myeloid", "OUT", vk.UNSCORED, vk.UNSCORED, "NK"]


def test_kappa():
    a = np.array(["B", "T", "NK", "OUT"] * 5)
    assert vk.cohen_kappa(a, a) == pytest.approx(1.0)
    assert vk.cohen_kappa(a, np.roll(a, 1)) < 0.0


# ----------------------------------------------------------------------------- evidence / questions
def test_evidence_independent_of_constant_size_factor(reg):
    rng = np.random.default_rng(0)
    pred = rng.gamma(2.0, 1.0, size=(200, 134))
    out = []
    for s in (1.0, 0.92663):
        r = copy.deepcopy(reg)
        p = pred * s
        r["evidence"]["teddy_head"]["q95_train_pred"] = list(np.percentile(p, 95, axis=0))
        out.append(vk.evidence(p, r))
    np.testing.assert_allclose(out[0], out[1], atol=1e-12)


def test_panels_and_questions(reg):
    pan = vk.question_panel(reg, "Q1")
    assert set(pan) == set(vk.LINEAGES)
    assert all(len(v) == reg["anm"]["panel_size"] for v in pan.values())
    flat = [p for v in pan.values() for p in v]
    assert len(flat) == len(set(flat))
    assert vk.question_panel(reg, "Q3") == {k: [p] for k, p in reg["questions"]["Q3"]["anchors"].items()}
    assert reg["questions"]["Q2"]["bar"] > reg["questions"]["Q1"]["bar"]
    e2 = reg["experiments"]["E2"]["perturbations"]["genes"]["panel_coding_genes"]
    assert set(e2) == set(flat)


def test_rule_calls_nested(reg):
    rng = np.random.default_rng(1)
    S = rng.uniform(0, 1, size=(500, 4))
    c1 = vk.rule_calls(S, reg["questions"]["Q1"]["bar"])
    c2 = vk.rule_calls(S, reg["questions"]["Q2"]["bar"])
    called2 = c2 != vk.NO_CALL
    assert np.all(c1[called2] == c2[called2])
    assert called2.sum() <= (c1 != vk.NO_CALL).sum()


def test_selective_accuracy_out_semantics():
    calls = np.array(["T", "", "B", "", "NK"])
    key = np.array(["T", "OUT", "OUT", "B", vk.UNSCORED])
    r = vk.selective_accuracy(calls, key)
    assert r["n_called_scored"] == 2 and r["accuracy_of_calls"] == pytest.approx(0.5)
    assert r["decision_accuracy"] == pytest.approx(2 / 4)
    assert r["out_decline_rate"] == pytest.approx(0.5)


# ----------------------------------------------------------------------------- ANM equivalence
def test_anm_action_score_is_fixed_multiple_of_rule(reg):
    rng = np.random.default_rng(2)
    names = reg["evidence"]["teddy_head"]["adt_names"]
    v = rng.uniform(0, 1, size=(300, len(names)))
    S = vk.class_scores(v, names, reg, "Q1")
    A = vk.anm_action_scores(v, names, reg, "Q1")
    n = reg["anm"]["panel_size"]
    G = reg["anm"]["field_gain"][str(n)]
    np.testing.assert_allclose(A, G * n * S, rtol=1e-5)
    for q in ("Q1", "Q2", "Q3"):
        Sq = vk.class_scores(v, names, reg, q)
        Aq = vk.anm_action_scores(v, names, reg, q)
        calls_rule = vk.rule_calls(Sq, reg["questions"][q]["bar"])
        calls_anm = vk.rule_calls(Aq, vk.anm_readout_threshold(reg, q))
        # exact float ties at the bar are measure-zero; allow 1e-12 slack only
        near = np.abs(Sq.max(axis=1) - reg["questions"][q]["bar"]) < 1e-12
        assert np.all(calls_rule[~near] == calls_anm[~near])


def _anm_runner():
    sys.path.insert(0, str(ANM_ROOT))
    try:
        from active_neural_matter.field import finite_field_runner as ffr  # noqa: F401
    except Exception:
        pytest.skip("ANM engine not importable (set ANM_ROOT)")
    return ffr


def test_field_star_matches_anm_engine(reg):
    ffr = _anm_runner()
    schema = json.loads((REPO / "bridge_anm/schemas/cite_lineage_finite_field_v0.json").read_text())
    schema["field_representation"].update({k: v for k, v in reg["anm"]["field_representation"].items() if k != "kind"})
    schema["source_event_schema"]["allowed_modalities"] = ["evidence"]
    rng = np.random.default_rng(3)
    vals = rng.uniform(0, 1, size=(4, 3))
    actions = [{"id": k, "label": k} for k in vk.LINEAGES]
    events = []
    for i, k in enumerate(vk.LINEAGES):
        for jj in range(3):
            events.append({"event_id": f"{k}{jj}", "time": 0, "action": k, "modality": "evidence",
                           "polarity": "support", "value": float(vals[i, jj]), "provenance": "test"})
    inst = {"instance_id": "x", "question": "q", "actions": actions, "proposed_source_events": events}
    val = ffr.validate_source_events(schema, inst)
    graph = ffr.build_graph(inst, val["field_events"])
    field = ffr.evolve_field(schema, graph, val["field_events"])
    a, e = vk.field_star(vals, reg["anm"]["field_representation"])
    for i, k in enumerate(vk.LINEAGES):
        assert field["state"][f"action:{k}"] == pytest.approx(a[i], rel=1e-9)
        for jj in range(3):
            assert field["state"][f"event:{k}{jj}"] == pytest.approx(e[i, jj], rel=1e-9)
    # closure readout through the engine's readout_coordinates
    inst["readout_coordinates"] = {k: [f"event:{k}{jj}" for jj in range(3)] for k in vk.LINEAGES}
    inst["readout_config"] = dict(reg["anm"]["closure_readout"])
    rd = ffr.readout(schema, inst, field["state"])
    rc = reg["anm"]["closure_readout"]
    for i, k in enumerate(vk.LINEAGES):
        want = rc["closure_weight"] * e[i].min() + rc["mean_coordinate_weight"] * e[i].mean() + rc["direct_action_weight"] * a[i]
        assert rd["action_scores"][k] == pytest.approx(want, rel=1e-9)


# ----------------------------------------------------------------------------- real data (skips without the pack)
@pytest.fixture(scope="module")
def pack():
    if not DATA.exists():
        pytest.skip("data pack not available")
    z = np.load(DATA, allow_pickle=False)
    return {k: z[k] for k in ("split", "sites", "donors", "cell_types", "adt_names")}


def test_registered_splits_reproduce(reg, pack):
    for name in ("train", "val", "test_primary", "test_secondary"):
        idx = vk.split_indices(pack["split"], pack["sites"], pack["donors"], reg, name)
        assert idx.size == reg["splits"][name]["n"]
    prim = vk.split_indices(pack["split"], pack["sites"], pack["donors"], reg, "test_primary")
    assert set(pack["donors"][prim].astype(str)) == {"19593", "13272"} and prim.size == 11294
    assert set(map(str, pack["cell_types"])) <= set(reg["annotation_map"])


def test_val_key_reproduces_registration(reg, pack):
    z = np.load(DATA, allow_pickle=False)
    va = np.where(pack["split"] == "val")[0]
    adt = z["adt"][va]
    keys = vk.build_keys(pack["cell_types"][va], adt, pack["adt_names"], reg)
    rep = vk.agreement_report(keys["annotation"], keys["gated"])
    assert rep["kappa5"] == pytest.approx(reg["key_quality"]["val"]["kappa5"], abs=1e-6)
    counts = {c: int(np.sum(keys["primary"] == c)) for c in list(vk.CLASSES) + [vk.UNSCORED]}
    assert counts == reg["key_quality"]["val"]["primary_key_counts"]
