"""Tests for bridge_anm/run_hard_proof.py corrections.

(a) label-cost arm: train pool and evaluation cells are disjoint and the
    abstain-calibrated thresholds are chosen on the train pool only;
(b) attribution permutation null shuffles values only within each lineage;
(c) the label-free TEDDY-alone rule reports needs_endpoint_labels=False;
plus: bootstrap vectorisation matches the old Counter loop, and --workers N
gives the same attribution as --workers 1.

Run with ANM_ROOT pointing at the ANM checkout (the module imports the engine).
"""
from __future__ import annotations

import copy
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
BRIDGE_PY = REPO / "bridge_anm"
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(REPO.parent / "ANM")))
if not (ANM_ROOT / "active_neural_matter").exists():
    pytest.skip(f"ANM checkout not found at {ANM_ROOT}; set ANM_ROOT", allow_module_level=True)
os.environ.setdefault("ANM_ROOT", str(ANM_ROOT))
sys.path.insert(0, str(BRIDGE_PY))

import run_hard_proof as hp  # noqa: E402
from lib.lineage_panels import KEY_MARKERS, LINEAGE_PANELS, MODALITY_FOR_LINEAGE  # noqa: E402

LINS = list(LINEAGE_PANELS)
P95 = {p: 1.0 for p in hp.PROTEINS}
CELL_TYPE_FOR = {"b_lineage": "Naive CD20+ B IGKC+", "t_lineage": "CD4+ T naive", "myeloid": "CD14+ Mono"}


def _panel(rng, lin, noise=0.0):
    out = {}
    for L, prots in LINEAGE_PANELS.items():
        for p in prots:
            if L == lin:
                v = rng.uniform(0.7, 1.0) if p == KEY_MARKERS[L] else rng.uniform(0.4, 0.9)
            else:
                v = rng.uniform(0.0, 0.2)
            out[p] = float(max(0.0, v + rng.normal(0.0, noise)))
    return out


def make_cells(n, seed=0):
    rng = np.random.default_rng(seed)
    cells = {}
    for i in range(n):
        lin = LINS[i % 3]
        true = _panel(rng, lin)
        pred = {p: float(max(0.0, v + rng.normal(0.0, 0.12))) for p, v in true.items()}
        cid = f"c{i:04d}"
        cells[cid] = {
            "cell_id": cid,
            "cell_type": CELL_TYPE_FOR[lin],
            "adt_pred_panel": pred,
            "adt_true_panel_holdout": true,
        }
    return cells


def make_events(cells):
    by = {}
    for cid, c in cells.items():
        evs = []
        t = 0
        for lin, prots in LINEAGE_PANELS.items():
            for p in prots:
                evs.append(
                    {
                        "event_id": f"{cid}:{p}",
                        "cell_id": cid,
                        "time": t,
                        "protein": p,
                        "action": lin,
                        "modality": MODALITY_FOR_LINEAGE[lin],
                        "value": min(1.0, max(0.0, c["adt_pred_panel"][p] / P95[p])),
                        "is_key_marker": p == KEY_MARKERS[lin],
                        "provenance": "test",
                    }
                )
                t += 1
        by[cid] = evs
    return by


# -------------------- (a) label-cost split / threshold selection --------------------
def test_split_disjoint_and_legacy_overlap():
    ids = [f"c{i:04d}" for i in range(4000)]
    perm = np.random.default_rng(17).permutation(len(ids))
    train, ev = hp.split_label_cost(ids, perm, 0.7, "disjoint")
    assert len(train) == 2800 and len(ev) == 1200
    assert not set(train) & set(ev)
    assert set(train) | set(ev) == set(ids)
    assert ev == sorted(ev)
    train_l, ev_l = hp.split_label_cost(ids, perm, 0.7, "legacy")
    assert ev_l == ids  # legacy: evaluate on every cell ...
    assert set(train_l) <= set(ev_l) and len(train_l) == 2800  # ... including the train pool
    assert train_l == train  # same pool, only the evaluation set differs
    with pytest.raises(ValueError):
        hp.split_label_cost(ids, perm, 0.7, "bogus")


def test_label_cost_thresholds_chosen_on_train_pool_only():
    cells = make_cells(420, seed=1)
    ids = sorted(cells)
    perm = np.random.default_rng(3).permutation(len(ids))
    train, ev = hp.split_label_cost(ids, perm, 0.7, "disjoint")
    targets = {"abstain_rate": 0.2, "Q_analogue": 0.9, "accuracy_strict_labeled": 0.7}

    def run(cells_):
        return hp.run_train_plus_anm_curves(
            ev, train, cells_, P95, targets, seed=17, calib_ids=train, calib_target_abstain=0.25
        )

    out = run(cells)
    sp = out["split"]
    assert sp["train_eval_overlap"] == 0
    assert sp["threshold_selection_set"] == "train_pool"
    assert sp["n_threshold_selection_cells"] == len(train)
    assert sp["threshold_selection_uses_labels"] is False

    # Scramble every evaluation cell's features: thresholds must not move.
    rng = np.random.default_rng(99)
    cells2 = copy.deepcopy(cells)
    for cid in ev:
        cells2[cid]["adt_pred_panel"] = {
            p: float(rng.uniform(0.0, 2.0)) for p in cells2[cid]["adt_pred_panel"]
        }
    out2 = run(cells2)

    def thresholds(o):
        c = o["curves"]
        return (
            [r["abstain_calibrated_to_anm_o2"]["conf_threshold"] for r in c["logistic_jev_standin"]],
            [
                r["abstain_calibrated_to_anm_o2"]["conf_threshold"]
                for r in c["mlp_small_head"]
                if "abstain_calibrated_to_anm_o2" in r
            ],
            [(r["best_grid"], r["abstain_calibrated_thr"]["thr"]) for r in c["teddy_thr_boost_grid"]],
        )

    assert out["budgets"] and out["curves"]["logistic_jev_standin"]
    assert thresholds(out) == thresholds(out2)
    # ... while the evaluation itself does see the scrambled cells
    assert [r["abstain_calibrated_to_anm_o2"]["abstain_rate"] for r in out["curves"]["logistic_jev_standin"]] != [
        r["abstain_calibrated_to_anm_o2"]["abstain_rate"] for r in out2["curves"]["logistic_jev_standin"]
    ]


# -------------------- (b) permutation null --------------------
def test_permutation_shuffles_only_within_lineage():
    cells = make_cells(30, seed=2)
    by = make_events(cells)
    rng = np.random.default_rng(5)
    moved = False
    for cid in sorted(cells):
        evs = by[cid]
        vals = hp.permute_event_values(evs, rng, "within_lineage")
        assert len(vals) == len(evs)
        for lin in LINS:
            idx = [i for i, e in enumerate(evs) if e["action"] == lin]
            assert sorted(vals[i] for i in idx) == sorted(evs[i]["value"] for i in idx)
        moved |= vals != [e["value"] for e in evs]
    assert moved  # it does shuffle
    # legacy null crosses lineages
    rng = np.random.default_rng(5)
    crossed = False
    for cid in sorted(cells):
        evs = by[cid]
        vals = hp.permute_event_values(evs, rng, "within_cell")
        assert sorted(vals) == sorted(e["value"] for e in evs)
        for lin in LINS:
            idx = [i for i, e in enumerate(evs) if e["action"] == lin]
            if sorted(vals[i] for i in idx) != sorted(evs[i]["value"] for i in idx):
                crossed = True
    assert crossed
    # deterministic given the seed
    a = [hp.permute_event_values(by[c], np.random.default_rng(7), "within_lineage") for c in sorted(cells)]
    b = [hp.permute_event_values(by[c], np.random.default_rng(7), "within_lineage") for c in sorted(cells)]
    assert a == b


def test_bootstrap_mode_matches_counter_loop():
    rng0 = np.random.default_rng(0)
    tops = np.array(rng0.choice(["CD19", "CD3", "CD16", "CD2"], size=57, p=[0.3, 0.3, 0.3, 0.1]))
    r_new, p_new = hp.bootstrap_mode(tops, np.random.default_rng(11), 300)
    rng = np.random.default_rng(11)
    r_old, p_old = [], []
    for _ in range(300):
        sample = rng.choice(tops, size=len(tops), replace=True)
        mp_, mn = Counter(sample).most_common(1)[0]
        r_old.append(mn / len(sample))
        p_old.append(str(mp_))
    assert r_new == r_old and p_new == p_old


def test_attribution_parallel_matches_serial_and_records_null():
    cells = make_cells(12, seed=4)
    by = make_events(cells)
    schema = __import__("json").loads(
        (BRIDGE_PY / "schemas/cite_lineage_finite_field_v0.json").read_text()
    )
    ids = sorted(cells)
    kw = dict(flip_grid=5, n_boot=10, n_perm=3, seed=17, perm_null="within_lineage")
    a1, c1 = hp.run_full_attribution(ids, cells, by, P95, schema, **kw)
    o0 = {cid: hp.make_instance(cid, cells[cid], by[cid], hp.CRITERIA["O0"], P95) for cid in ids}
    with hp.Runner(2, schema, o0) as runner:
        a2, c2 = hp.run_full_attribution(ids, cells, by, P95, schema, runner=runner, **kw)
    a1.pop("runtime_sec"), a2.pop("runtime_sec")
    assert a1 == a2
    for k in c1:
        assert list(c1[k]) == list(c2[k])
    assert a1["permutation_test"]["null"] == "within_lineage"
    assert a1["permutation_test"]["n_perm"] == 3


# -------------------- (c) TEDDY-alone edit cost --------------------
def test_teddy_alone_rule_is_label_free():
    cells = make_cells(30, seed=6)
    ids = sorted(cells)
    arm = hp.run_teddy_arm(ids, ids, cells, P95)
    assert arm["criterion_edit_cost"]["needs_endpoint_labels"] is False
    for c in hp.CRIT_IDS:
        assert arm["criteria"][c]["rule"]["label_fit"] is False
        ag = arm["criteria"][c]["annotation_graded"]
        assert ag["in_scope"]["n"] == 30 and ag["out_of_scope"]["n"] == 0
    legacy = hp.run_teddy_arm(ids, ids, cells, P95, legacy_edit_cost=True)
    assert legacy["criterion_edit_cost"]["needs_endpoint_labels"] is True


def test_annotation_map_groups():
    assert hp.cell_type_group("NK CD158e1+") == "nk"
    assert hp.cell_type_group("Reticulocyte") == "erythroid"
    assert hp.cell_type_group("HSC") == "progenitor"
    assert hp.cell_type_group("CD16+ Mono") == "myeloid"
    assert hp.cell_type_group("something new") == "other"
    g = hp.annotation_graded(
        ["myeloid", None, "t_lineage", "t_lineage"], ["NK", "NK", "CD8+ T naive", "B1 B IGKC+"]
    )
    assert g["out_of_scope"]["n"] == 2 and g["out_of_scope"]["n_called"] == 1
    assert g["out_of_scope"]["per_group"]["nk"]["call_counts"] == {"myeloid": 1}
    assert g["in_scope"]["n"] == 2 and g["in_scope"]["n_correct"] == 1
    assert g["in_scope"]["accuracy_strict"] == 0.5
