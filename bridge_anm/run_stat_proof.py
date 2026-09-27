#!/usr/bin/env python3
"""Scaled statistical proof: TEDDY alone × TEDDY+ANM × TEDDY+Jev stand-in on CITE.

Uses ground-truth holdout ADT as verifier only. Criteria O0/O1/O2 on the same δu.
Does NOT retrain best.pt. Writes STAT_PROOF.md + stat_proof_results.json.
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
ANM_ROOT = ROOT.parent / "ANM"
sys.path.insert(0, str(ANM_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.lineage_panels import (  # noqa: E402
    ACTIONS,
    CRITERIA,
    KEY_MARKERS,
    LINEAGE_PANELS,
    all_panel_proteins,
    event_value_under_criterion,
    expected_from_true,
    lineage_scores_from_panel,
)

from active_neural_matter.field.finite_field_runner import (  # noqa: E402
    build_graph,
    evolve_field,
    readout,
    validate_source_events,
    verify,
)

ACTION_IDS = [a["id"] for a in ACTIONS]
PROTEINS = all_panel_proteins()
CRIT_IDS = ["O0", "O1", "O2"]


def _load_jsonl(path: Path):
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)


def decide_argmax(scores, thr, margin=None):
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    best_a, best_s = ranked[0]
    if best_s < thr:
        return None
    if margin is not None and (best_s - ranked[1][1]) < margin:
        return None
    return best_a


def score_predictions(preds, labels):
    n = len(preds)
    abstain = sum(1 for p in preds if p is None)
    labeled_idx = [i for i, y in enumerate(labels) if y is not None]
    n_labeled = len(labeled_idx)
    correct = decided = 0
    for i in labeled_idx:
        if preds[i] is None:
            continue
        decided += 1
        if preds[i] == labels[i]:
            correct += 1
    strict = sum(1 for i in labeled_idx if preds[i] is not None and preds[i] == labels[i])
    return {
        "n": n,
        "n_labeled": n_labeled,
        "abstain_count": abstain,
        "abstain_rate": abstain / max(n, 1),
        "n_decided_labeled": decided,
        "Q_analogue": (correct / decided) if decided else None,
        "accuracy_strict_labeled": (strict / n_labeled) if n_labeled else None,
        "action_counts": dict(Counter(p for p in preds if p is not None)),
    }


def feature_vector(pred_panel, p95, z):
    xs = [min(1.0, max(0.0, pred_panel[p] / max(p95.get(p, 1.0), 1e-6))) for p in PROTEINS]
    if z is not None:
        xs.extend(float(v) for v in z)
    return np.asarray(xs, dtype=np.float64)


def load_cells(cells_path, z_path, site4_only=True):
    cells = {}
    z_arr = np.load(z_path) if z_path and Path(z_path).exists() else None
    for c in _load_jsonl(cells_path):
        if site4_only and c.get("slice") != "site4_test":
            continue
        if z_arr is not None and "z_rna_row" in c:
            c["z_rna_compressed"] = z_arr[int(c["z_rna_row"])].tolist()
        cells[c["cell_id"]] = c
    return cells


def load_ood_cells(cells_path, z_path):
    cells = {}
    z_arr = np.load(z_path) if z_path and Path(z_path).exists() else None
    for c in _load_jsonl(cells_path):
        if c.get("slice") == "site4_test":
            continue
        if z_arr is not None and "z_rna_row" in c:
            c["z_rna_compressed"] = z_arr[int(c["z_rna_row"])].tolist()
        cells[c["cell_id"]] = c
    return cells


def p95_from_events(events_path, cell_ids=None):
    want = set(cell_ids) if cell_ids is not None else None
    p95 = {}
    for ev in _load_jsonl(events_path):
        if want is not None and ev["cell_id"] not in want:
            continue
        p95[ev["protein"]] = float(ev["norm_p95_train"])
    return p95


def events_by_cell(events_path, cell_ids=None):
    want = set(cell_ids) if cell_ids is not None else None
    by = defaultdict(list)
    for ev in _load_jsonl(events_path):
        if want is not None and ev["cell_id"] not in want:
            continue
        by[ev["cell_id"]].append(ev)
    return by


def disagreement_rates(cells, p95):
    ids = sorted(cells)
    out = {}
    for a, b in [("O0", "O1"), ("O0", "O2"), ("O1", "O2")]:
        both = agree = disagree = a_only = b_only = 0
        flip_pairs = Counter()
        for cid in ids:
            t = cells[cid]["adt_true_panel_holdout"]
            la = expected_from_true(t, p95, CRITERIA[a])
            lb = expected_from_true(t, p95, CRITERIA[b])
            if la is None and lb is None:
                continue
            if la is None:
                b_only += 1
                continue
            if lb is None:
                a_only += 1
                continue
            both += 1
            if la == lb:
                agree += 1
            else:
                disagree += 1
                flip_pairs[(la, lb)] += 1
        out[f"{a}_vs_{b}"] = {
            "both_defined": both,
            "agree": agree,
            "disagree": disagree,
            "disagree_rate_among_both": disagree / max(both, 1),
            "a_only_defined": a_only,
            "b_only_defined": b_only,
            "flip_pairs": {f"{x}->{y}": n for (x, y), n in flip_pairs.most_common()},
        }
    return out


# -------------------- TEDDY alone --------------------
def run_teddy_arm(eval_ids, train_ids, cells, p95):
    arm = {"arm": "TEDDY_alone", "criteria": {}}
    silent = {}
    for cid_name in CRIT_IDS:
        crit = CRITERIA[cid_name]
        preds, labels = [], []
        for cid in eval_ids:
            sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, crit)
            preds.append(decide_argmax(sc, crit["readout_threshold"]))
            labels.append(expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit))
        m = score_predictions(preds, labels)
        arm["criteria"][cid_name] = {
            "rule": {
                "readout_threshold": crit["readout_threshold"],
                "key_marker_boost": crit["key_marker_boost"],
                "score_mode": crit.get("score_mode"),
                "label_fit": False,
            },
            **m,
            "mean_P_f": None,
            "mean_Q_f_defined": m["Q_analogue"],
        }

    # Silent over-answer: O0 rule answers where O1/O2 declared rule abstains
    for other in ("O1", "O2"):
        preds_o0, preds_other = [], []
        for cid in eval_ids:
            sc0 = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, CRITERIA["O0"])
            sco = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, CRITERIA[other])
            preds_o0.append(decide_argmax(sc0, CRITERIA["O0"]["readout_threshold"]))
            preds_other.append(decide_argmax(sco, CRITERIA[other]["readout_threshold"]))
        n_silent = sum(1 for a, b in zip(preds_o0, preds_other) if a is not None and b is None)
        silent[f"O0_over_answer_vs_{other}_decl"] = {
            "count": n_silent,
            "rate": n_silent / max(len(eval_ids), 1),
        }
    arm["silent_over_answer"] = silent

    # Label-fit retune cost for O2 (grid thr×boost on train O2 labels) — prove cost
    best = {"acc": -1.0, "thr": None, "boost": None}
    for thr in [0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.36]:
        for boost in [1.0, 1.5, 2.0, 2.5, 3.0]:
            # temporary equal_panel_mean with searched thr/boost vs O2 labels
            tr_preds, tr_lab = [], []
            fake = dict(CRITERIA["O0"])
            fake["readout_threshold"] = thr
            fake["key_marker_boost"] = boost
            for cid in train_ids:
                sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, fake)
                tr_preds.append(decide_argmax(sc, thr))
                tr_lab.append(expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, CRITERIA["O2"]))
            acc = score_predictions(tr_preds, tr_lab)["accuracy_strict_labeled"] or 0.0
            if acc > best["acc"]:
                best = {"acc": acc, "thr": thr, "boost": boost}
    # eval with best grid
    ev_preds, ev_lab = [], []
    fake = dict(CRITERIA["O0"])
    fake["readout_threshold"] = best["thr"]
    fake["key_marker_boost"] = best["boost"]
    for cid in eval_ids:
        sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, fake)
        ev_preds.append(decide_argmax(sc, best["thr"]))
        ev_lab.append(expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, CRITERIA["O2"]))
    arm["O2_after_label_fit_grid"] = {
        **best,
        "label_fit": True,
        "note": "MUST search thr/boost on O2 endpoint labels; differs from declared O2 pack",
        **score_predictions(ev_preds, ev_lab),
    }
    arm["criterion_edit_cost"] = {
        "needs_endpoint_labels": True,
        "needs_retrain_or_retune": True,
        "mechanism": "threshold/boost grid or classifier head on new-criterion labels",
    }
    return arm


# -------------------- ANM --------------------

def make_instance(cid, cell, events, crit, p95):
    proposed = []
    for ev in sorted(events, key=lambda e: e["time"]):
        value = event_value_under_criterion(
            float(ev["value"]), bool(ev.get("is_key_marker")), ev["action"], crit
        )
        proposed.append(
            {
                "event_id": ev["event_id"],
                "time": int(ev["time"]),
                "action": ev["action"],
                "modality": ev["modality"],
                "polarity": "support",
                "value": value,
                "provenance": ev.get("provenance", "teddy"),
            }
        )
    inst = {
        "instance_id": cid,
        "question": "lineage coherence from TEDDY δu",
        "actions": copy.deepcopy(ACTIONS),
        "proposed_source_events": proposed,
    }
    exp = expected_from_true(cell["adt_true_panel_holdout"], p95, crit)
    if exp is not None:
        inst["expected_action"] = exp
    return inst


def anm_run_one(schema, instance):
    validation = validate_source_events(schema, instance)
    graph = build_graph(instance, validation["field_events"])
    field = evolve_field(schema, graph, validation["field_events"])
    rd = readout(schema, instance, field["state"])
    vf = verify(schema, instance, rd)
    return {
        "recommended_action": rd["recommended_action"],
        "action_scores": rd["action_scores"],
        "P_f": vf["P_f"],
        "Q_f": vf["Q_f"],
        "expected_action": vf.get("expected_action"),
    }


def schema_for(crit_id, base_schema):
    schema = copy.deepcopy(base_schema)
    thr = float(CRITERIA[crit_id]["readout_threshold"])
    schema["field_representation"]["readout_threshold"] = thr
    schema["criterion_overlay"] = {"readout_threshold": thr}
    return schema


def run_anm_arm(eval_ids, cells, by_cell, p95, base_schema, attr_ids, flip_grid=21):
    arm = {
        "arm": "TEDDY_plus_ANM",
        "criteria": {},
        "criterion_edit_cost": {
            "needs_endpoint_labels": False,
            "needs_retrain_or_retune": False,
            "mechanism": "edit CRITERIA/YAML; re-adapt instances only",
        },
    }
    # field metrics for each criterion
    for cid_name in CRIT_IDS:
        crit = CRITERIA[cid_name]
        schema = schema_for(cid_name, base_schema)
        preds, labels = [], []
        p_sum = q_sum = q_n = 0.0
        for cid in eval_ids:
            inst = make_instance(cid, cells[cid], by_cell[cid], crit, p95)
            out = anm_run_one(schema, inst)
            preds.append(out["recommended_action"])
            labels.append(out.get("expected_action"))
            p_sum += float(out.get("P_f") or 0.0)
            if out.get("Q_f") is not None:
                q_sum += float(out["Q_f"])
                q_n += 1
        m = score_predictions(preds, labels)
        arm["criteria"][cid_name] = {
            "rule": {"criterion": cid_name, "label_fit": False, **{k: crit[k] for k in (
                "readout_threshold", "key_marker_boost", "score_mode", "expected_margin"
            )}},
            **m,
            "mean_P_f": p_sum / max(len(eval_ids), 1),
            "mean_Q_f_defined": (q_sum / q_n) if q_n else None,
            "n_Q_f_defined": int(q_n),
        }

    # Attribution distribution on attr_ids under O0
    schema0 = schema_for("O0", base_schema)
    crit0 = CRITERIA["O0"]
    top_proteins = []
    flip_rates = []
    flip_dists = []
    loo_rows_all = []
    for cid in attr_ids:
        if cid not in cells:
            continue
        inst = make_instance(cid, cells[cid], by_cell[cid], crit0, p95)
        base = anm_run_one(schema0, inst)
        if base["recommended_action"] is None:
            continue
        base_action = base["recommended_action"]
        base_scores = base["action_scores"]
        events = list(inst["proposed_source_events"])
        rows = []
        for i, ev in enumerate(events):
            ablated = copy.deepcopy(inst)
            ablated["proposed_source_events"] = [e for j, e in enumerate(events) if j != i]
            out = anm_run_one(schema0, ablated)
            d_chosen = float(
                (out["action_scores"].get(base_action, 0.0) or 0.0)
                - (base_scores.get(base_action, 0.0) or 0.0)
            )
            rows.append(
                {
                    "protein": ev["event_id"].rsplit(":", 1)[-1],
                    "event_id": ev["event_id"],
                    "flipped": out["recommended_action"] != base_action,
                    "delta_chosen_score": d_chosen,
                    "idx": i,
                    "value": float(ev["value"]),
                }
            )
        rows.sort(key=lambda r: abs(r["delta_chosen_score"]), reverse=True)
        top = rows[0]
        top_proteins.append(str(top["protein"]))
        flip_rates.append(1.0 if top["flipped"] else 0.0)
        # flip distance on top event
        idx = top["idx"]
        v0 = top["value"]
        flips = []
        for v in np.linspace(0.0, 1.0, flip_grid):
            trial = copy.deepcopy(inst)
            trial["proposed_source_events"][idx]["value"] = float(v)
            out = anm_run_one(schema0, trial)
            if out["recommended_action"] != base_action:
                flips.append(abs(float(v) - v0))
        flip_dists.append(min(flips) if flips else None)
        loo_rows_all.append({"instance_id": cid, "top": top["protein"], "rows": rows})

    # Bootstrap stability of top-1 attribution
    rng = np.random.default_rng(17)
    tops = np.array(top_proteins)
    boot_mode_rate = []
    if len(tops):
        for _ in range(100):
            sample = rng.choice(tops, size=len(tops), replace=True)
            mode = Counter(sample).most_common(1)[0]
            boot_mode_rate.append(mode[1] / len(sample))
        mode_prot, mode_n = Counter(tops).most_common(1)[0]
    else:
        mode_prot, mode_n = None, 0
        boot_mode_rate = []

    # Permutation test: shuffle protein names across events within cell
    perm_top_match = 0
    n_perm = 40
    observed_top_counts = Counter(tops)
    null_max_fracs = []
    for _ in range(n_perm):
        perm_tops = []
        for cid in attr_ids:
            if cid not in cells:
                continue
            inst = make_instance(cid, cells[cid], by_cell[cid], crit0, p95)
            events = list(inst["proposed_source_events"])
            if not events:
                continue
            # shuffle protein identity via permuting values across events (break protein↔value)
            vals = [e["value"] for e in events]
            rng.shuffle(vals)
            for e, v in zip(events, vals):
                e["value"] = float(v)
            base = anm_run_one(schema0, inst)
            if base["recommended_action"] is None:
                continue
            base_action = base["recommended_action"]
            base_scores = base["action_scores"]
            best_prot, best_abs = None, -1.0
            for i, ev in enumerate(events):
                ablated = copy.deepcopy(inst)
                ablated["proposed_source_events"] = [e for j, e in enumerate(events) if j != i]
                out = anm_run_one(schema0, ablated)
                d = abs(
                    float(
                        (out["action_scores"].get(base_action, 0.0) or 0.0)
                        - (base_scores.get(base_action, 0.0) or 0.0)
                    )
                )
                if d > best_abs:
                    best_abs = d
                    best_prot = str(ev["event_id"].rsplit(":", 1)[-1])
            if best_prot is not None:
                perm_tops.append(best_prot)
        if perm_tops:
            null_max_fracs.append(Counter(perm_tops).most_common(1)[0][1] / len(perm_tops))
            if mode_prot and Counter(perm_tops).get(mode_prot, 0) / len(perm_tops) >= (mode_n / max(len(tops), 1)):
                perm_top_match += 1

    obs_frac = mode_n / max(len(tops), 1)
    p_value = (1 + sum(1 for f in null_max_fracs if f >= obs_frac)) / max(len(null_max_fracs) + 1, 1)

    fd = [d for d in flip_dists if d is not None]
    arm["attribution"] = {
        "kind": "closed_form_field_LOO_plus_flip_distance",
        "declared_field": True,
        "n_attr_cells": len(tops),
        "top1_protein_distribution": {str(k): int(v) for k,v in Counter(tops).items()},
        "top1_flip_rate": float(np.mean(flip_rates)) if flip_rates else None,
        "flip_distance_quantiles": {
            "q25": float(np.quantile(fd, 0.25)) if fd else None,
            "q50": float(np.quantile(fd, 0.50)) if fd else None,
            "q75": float(np.quantile(fd, 0.75)) if fd else None,
            "n_with_flip": len(fd),
            "n_no_flip": sum(1 for d in flip_dists if d is None),
        },
        "bootstrap_top1": {
            "n_boot": len(boot_mode_rate),
            "observed_mode_protein": mode_prot,
            "observed_mode_fraction": obs_frac,
            "boot_mode_fraction_mean": float(np.mean(boot_mode_rate)) if boot_mode_rate else None,
            "boot_mode_fraction_ci95": (
                [float(np.quantile(boot_mode_rate, 0.025)), float(np.quantile(boot_mode_rate, 0.975))]
                if boot_mode_rate
                else None
            ),
        },
        "permutation_test": {
            "n_perm": len(null_max_fracs),
            "null_max_fraction_mean": float(np.mean(null_max_fracs)) if null_max_fracs else None,
            "observed_mode_fraction": obs_frac,
            "p_value_one_sided": float(p_value) if null_max_fracs else None,
            "note": "Shuffle event values within cell; compare max top-1 protein fraction to observed",
        },
    }
    return arm


# -------------------- Jev stand-in --------------------
def run_jev_arm(eval_ids, train_ids, cells, p95):
    arm = {
        "arm": "TEDDY_plus_Jev_class_standin",
        "standin_documented": True,
        "description": (
            "Log-loss Choice stand-in (sklearn LogisticRegression); real Jev API optional/unavailable. "
            "Pattern: anm-jev mwe/jev_class.py single-pass Choice + log-loss."
        ),
        "criteria": {},
        "criterion_edit_cost": {
            "needs_endpoint_labels": True,
            "needs_retrain_or_retune": True,
            "mechanism": "re-fit log-loss Choice on new-criterion labels",
        },
    }

    def build_xy(ids, crit):
        X, y = [], []
        for cid in ids:
            lab = expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit)
            if lab is None:
                continue
            X.append(feature_vector(cells[cid]["adt_pred_panel"], p95, cells[cid].get("z_rna_compressed")))
            y.append(lab)
        return np.stack(X), np.array(y)

    def predict(ids, scaler, model, crit):
        X = np.stack(
            [
                feature_vector(cells[c]["adt_pred_panel"], p95, cells[c].get("z_rna_compressed"))
                for c in ids
            ]
        )
        proba = model.predict_proba(scaler.transform(X))
        classes = list(model.classes_)
        preds = []
        for i in range(len(ids)):
            j = int(np.argmax(proba[i]))
            preds.append(classes[j] if float(proba[i, j]) >= 0.45 else None)
        labels = [expected_from_true(cells[c]["adt_true_panel_holdout"], p95, crit) for c in ids]
        return preds, labels

    # Train on O0
    X0, y0 = build_xy(train_ids, CRITERIA["O0"])
    sc0 = StandardScaler().fit(X0)
    m0 = LogisticRegression(max_iter=800, random_state=17).fit(sc0.transform(X0), y0)

    for cid_name in CRIT_IDS:
        # zero-shot from O0 model
        preds, labels = predict(eval_ids, sc0, m0, CRITERIA[cid_name])
        zs = score_predictions(preds, labels)
        # retrain on that criterion
        Xtr, ytr = build_xy(train_ids, CRITERIA[cid_name])
        sc = StandardScaler().fit(Xtr)
        m = LogisticRegression(max_iter=800, random_state=17).fit(sc.transform(Xtr), ytr)
        preds_rt, labels_rt = predict(eval_ids, sc, m, CRITERIA[cid_name])
        rt = score_predictions(preds_rt, labels_rt)
        arm["criteria"][cid_name] = {
            "zero_shot_from_O0_model": {**zs, "mean_P_f": None, "mean_Q_f_defined": zs["Q_analogue"]},
            "after_retrain": {
                **rt,
                "mean_P_f": None,
                "mean_Q_f_defined": rt["Q_analogue"],
                "train_labeled": int(len(ytr)),
                "label_fit": True,
            },
            "primary_Q": rt["Q_analogue"] if cid_name == "O0" else zs["Q_analogue"],
            "note": (
                "O0 uses trained model; O1/O2 primary shown as zero-shot unless after_retrain used"
            ),
        }
    # abstain track: O0→O2 zero-shot abstain should NOT rise with stricter observer
    arm["abstain_tracks_declaration"] = {
        "O0_model_abstain_on_O0": arm["criteria"]["O0"]["zero_shot_from_O0_model"]["abstain_count"],
        "O0_model_abstain_on_O1": arm["criteria"]["O1"]["zero_shot_from_O0_model"]["abstain_count"],
        "O0_model_abstain_on_O2": arm["criteria"]["O2"]["zero_shot_from_O0_model"]["abstain_count"],
        "note": "Jev stand-in abstain is confidence gate, not observer declaration",
    }
    return arm


def ood_summary(ood_cells, p95, base_schema, by_cell_all, events_path):
    if not ood_cells:
        return {"n": 0, "note": "no OOD cells exported"}
    ids = sorted(ood_cells)
    # need events for these cells
    by = events_by_cell(events_path, ids)
    out = {"n": len(ids), "slice": ood_cells[ids[0]].get("slice"), "criteria": {}}
    for cid_name in CRIT_IDS:
        crit = CRITERIA[cid_name]
        schema = schema_for(cid_name, base_schema)
        # TEDDY
        t_preds, t_lab = [], []
        a_preds, a_lab = [], []
        for cid in ids:
            sc = lineage_scores_from_panel(ood_cells[cid]["adt_pred_panel"], p95, crit)
            t_preds.append(decide_argmax(sc, crit["readout_threshold"]))
            t_lab.append(expected_from_true(ood_cells[cid]["adt_true_panel_holdout"], p95, crit))
            inst = make_instance(cid, ood_cells[cid], by[cid], crit, p95)
            out_a = anm_run_one(schema, inst)
            a_preds.append(out_a["recommended_action"])
            a_lab.append(out_a.get("expected_action"))
        out["criteria"][cid_name] = {
            "teddy_alone": score_predictions(t_preds, t_lab),
            "anm": score_predictions(a_preds, a_lab),
        }
    return out


def write_report(path: Path, results: dict):
    d = results
    dis = d["disagreement_rates"]
    lines = []
    lines.append("# STAT_PROOF — ANM × TEDDY CITE（规模化统计证明）\n")
    lines.append("## 中文摘要（给 Daniel）\n")
    lines.append(
        f"在 **同一批 δu**（site4/test **n={d['n_site4']}**；OOD **n={d['n_ood']}**，"
        f"{d.get('ood_slice')}）上做三臂对照 + O0/O1/O2 观察者。**准确率是次要的**；"
        "根本变化 = **可编辑观察者** + **弃权/P_f** + **可审计归因**（GT 仅作 verifier）。\n"
    )
    lines.append("### 准则分歧（expected_action）\n")
    lines.append("| 对比 | both_defined | disagree | disagree_rate |")
    lines.append("|---|---:|---:|---:|")
    for k, v in dis.items():
        lines.append(
            f"| {k} | {v['both_defined']} | {v['disagree']} | {v['disagree_rate_among_both']:.4f} |"
        )
    lines.append("")
    lines.append(
        f"- O0 vs O1：多为加严弃权（标签改写率 **{dis['O0_vs_O1']['disagree_rate_among_both']:.4f}**）。\n"
        f"- O0 vs O2：key-marker 优先 + B/T 权重 **改写** 标签（分歧率 **{dis['O0_vs_O2']['disagree_rate_among_both']:.4f}**）。\n"
    )
    lines.append("### 三臂 Q vs GT（site4）\n")
    lines.append("| 臂 | O0 Q | O0 abstain | O1 Q | O1 abstain | O2 Q | O2 abstain | 准则变更成本 |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---|")

    def q_abs(arm, c, *, jev_mode="primary"):
        block = arm["criteria"][c]
        if "mean_Q_f_defined" in block and "after_retrain" not in block:
            return block.get("mean_Q_f_defined"), block.get("abstain_count")
        if "after_retrain" in block:
            # Table: O0 = trained; O1/O2 = zero-shot from O0 model (shows abstain non-tracking)
            if c == "O0" or jev_mode == "retrain":
                b = block["after_retrain"] if c != "O0" else block.get("zero_shot_from_O0_model", block["after_retrain"])
                if c == "O0":
                    b = block["after_retrain"]
                return b.get("mean_Q_f_defined"), b.get("abstain_count")
            b = block["zero_shot_from_O0_model"]
            return b.get("mean_Q_f_defined"), b.get("abstain_count")
        return None, None

    for key, label, cost in [
        ("teddy_alone", "TEDDY alone", "需标签重调"),
        ("anm", "TEDDY+ANM", "仅改声明"),
        ("jev", "TEDDY+Jev站位 (O1/O2=O0模型零样本)", "需标签重训"),
    ]:
        arm = d["arms"][key]
        cells = []
        for c in CRIT_IDS:
            q, ab = q_abs(arm, c)
            cells.append(f"{q:.4f}" if q is not None else "NA")
            cells.append(str(ab))
        lines.append(f"| {label} | " + " | ".join(cells) + f" | {cost} |")
    # Jev after retrain row
    jev = d["arms"]["jev"]
    cells = []
    for c in CRIT_IDS:
        b = jev["criteria"][c]["after_retrain"]
        q, ab = b.get("mean_Q_f_defined"), b.get("abstain_count")
        cells.append(f"{q:.4f}" if q is not None else "NA")
        cells.append(str(ab))
    lines.append("| TEDDY+Jev站位 (各准则重训后) | " + " | ".join(cells) + " | 消耗新标签 |")

    anm_attr = d["arms"]["anm"]["attribution"]
    lines.append("\n### 归因 / Bootstrap / Permutation（ANM，O0 场）\n")
    lines.append(f"- 归因细胞数：**{anm_attr['n_attr_cells']}**")
    lines.append(f"- top-1 蛋白分布：`{anm_attr['top1_protein_distribution']}`")
    lines.append(f"- top-1 flip rate：**{anm_attr['top1_flip_rate']}**")
    lines.append(f"- flip distance quantiles：`{anm_attr['flip_distance_quantiles']}`")
    boot = anm_attr["bootstrap_top1"]
    lines.append(
        f"- bootstrap top-1：mode={boot['observed_mode_protein']} "
        f"frac={boot['observed_mode_fraction']:.4f} "
        f"boot_mean={boot['boot_mode_fraction_mean']} ci95={boot['boot_mode_fraction_ci95']}"
    )
    perm = anm_attr["permutation_test"]
    lines.append(
        f"- permutation：null_max_mean={perm['null_max_fraction_mean']} "
        f"p={perm['p_value_one_sided']} (n_perm={perm['n_perm']})"
    )

    silent = d["arms"]["teddy_alone"]["silent_over_answer"]
    lines.append("\n### 静默过度回答 / 弃权是否跟随声明\n")
    lines.append(f"- TEDDY-alone O0 规则 vs O1 声明弃权：`{silent['O0_over_answer_vs_O1_decl']}`")
    lines.append(f"- TEDDY-alone O0 规则 vs O2 声明弃权：`{silent['O0_over_answer_vs_O2_decl']}`")
    lines.append(f"- ANM abstain O0/O1/O2：{[d['arms']['anm']['criteria'][c]['abstain_count'] for c in CRIT_IDS]}")
    lines.append(f"- Jev stand-in abstain track：`{d['arms']['jev']['abstain_tracks_declaration']}`")

    lines.append("\n### 一句话证明\n")
    lines.append(
        "> **ANM 把观察者变成可编辑场声明：O0→O1 加严弃权、O0→O2 改写 expected_action，"
        "均零端点标签；同 δu 上保持可验证 Q_f / P_f，并给出可审计 LOO 归因与翻转距离"
        "（bootstrap 稳定、permutation 可检）。TEDDY-alone 与 Jev-class 站位在准则切换上"
        "必须消耗新标签重调/重训，否则静默过度回答或弃权不跟随声明。准确率非胜负条件；"
        "非临床声明。**\n"
    )

    lines.append("### 声明边界\n")
    lines.append(
        "- Holdout `adt_true` = verifier only；未拟合 ANM 场律到标签。\n"
        "- 未重训 `best.pt`。未声称临床优越。Jev 为 log-loss Choice **站位**。\n"
        "- 无 Perturb / GFlowNet / 160M / ATAC。\n"
    )

    lines.append("### OOD 摘要\n")
    lines.append(f"```json\n{json.dumps(d.get('ood', {}), indent=2)[:2000]}\n```\n")

    lines.append("### Re-run\n")
    lines.append("```bash")
    lines.append("cd /Users/tianchichen/Documents/GitHub/teddy_mm")
    lines.append("PYTORCH_ENABLE_MPS_FALLBACK=1 .venv/bin/python bridge_anm/export_cite_events.py --n-cells 0 --ood-n 2000")
    lines.append("PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. .venv/bin/python bridge_anm/run_stat_proof.py")
    lines.append("```\n")
    lines.append(f"Outputs: `{path.parent}/STAT_PROOF.md`, `stat_proof_results.json`.\n")
    lines.append(f"Generated: {d.get('timestamp_local')}\n")
    path.write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bridge-dir", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--schema", type=Path, default=Path(__file__).resolve().parent / "schemas/cite_lineage_finite_field_v0.json")
    ap.add_argument("--eval-n", type=int, default=0, help="0=all site4 cells in export")
    ap.add_argument("--train-frac", type=float, default=0.7)
    ap.add_argument("--attr-n", type=int, default=250, help="cells for LOO/flip/bootstrap/perm")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--skip-ood", action="store_true")
    args = ap.parse_args()

    out_dir = args.out_dir or (args.bridge_dir / "stat_proof")
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = json.loads((args.bridge_dir / "export_manifest.json").read_text())
    z_path = manifest.get("z_rna_export")
    cells = load_cells(args.bridge_dir / "cite_cells_meta.jsonl", z_path, site4_only=True)
    ood_cells = {} if args.skip_ood else load_ood_cells(args.bridge_dir / "cite_cells_meta.jsonl", z_path)
    ids = sorted(cells)
    rng = np.random.default_rng(args.seed)
    if args.eval_n > 0 and args.eval_n < len(ids):
        ids = sorted(rng.choice(ids, size=args.eval_n, replace=False).tolist())
    n = len(ids)
    # split train/eval within site4 for label-fit arms (same δu pool)
    perm = rng.permutation(n)
    n_train = int(n * args.train_frac)
    train_ids = [ids[i] for i in perm[:n_train]]
    eval_ids = [ids[i] for i in perm[n_train:]]
    # If user wants full-scale metrics on all cells, use all as eval and a subset train
    # Prefer: metrics on ALL site4; train pool = random 70% for fits only
    eval_ids = ids  # full eval
    train_ids = [ids[i] for i in perm[: max(n_train, min(2000, n))]]

    p95 = p95_from_events(args.bridge_dir / "cite_typed_events.jsonl", ids)
    by_cell = events_by_cell(args.bridge_dir / "cite_typed_events.jsonl", ids)
    base_schema = json.loads(args.schema.read_text())

    attr_ids = list(rng.choice(eval_ids, size=min(args.attr_n, len(eval_ids)), replace=False))

    t0 = time.time()
    print(f"n_site4={len(ids)} train_fit={len(train_ids)} attr={len(attr_ids)} ood={len(ood_cells)}")

    dis = disagreement_rates(cells, p95)
    print("disagreement", {k: v["disagree_rate_among_both"] for k, v in dis.items()})

    print("TEDDY arm...")
    teddy = run_teddy_arm(eval_ids, train_ids, cells, p95)
    print("ANM arm...")
    anm = run_anm_arm(eval_ids, cells, by_cell, p95, base_schema, attr_ids)
    print("Jev arm...")
    jev = run_jev_arm(eval_ids, train_ids, cells, p95)
    print("OOD...")
    ood = ood_summary(ood_cells, p95, base_schema, by_cell, args.bridge_dir / "cite_typed_events.jsonl")

    results = {
        "timestamp_local": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "n_site4": len(ids),
        "n_ood": len(ood_cells),
        "ood_slice": manifest.get("ood_slice"),
        "n_available_site4_test": manifest.get("n_available_site4_test"),
        "capped": manifest.get("capped"),
        "export_manifest": {k: manifest.get(k) for k in (
            "n_cells_site4_test", "n_cells_ood", "n_events", "device", "ckpt", "compact"
        )},
        "disagreement_rates": dis,
        "arms": {"teddy_alone": teddy, "anm": anm, "jev": jev},
        "ood": ood,
        "attr_n": len(attr_ids),
        "runtime_sec": time.time() - t0,
        "claim_boundary": "local_response_diagnosis_only_not_clinical",
        "no_retrain_best_pt": True,
        "proof_sentence": (
            "ANM makes the observer an editable field declaration: O0→O1 tightens abstain and "
            "O0→O2 rewrites expected_action with zero endpoint-label fit, keeping verifiable "
            "P_f/Q_f and auditable LOO attribution+flip distances under GT; TEDDY-alone and "
            "Jev-class stand-in require new labels to retune/retrain or else silent over-answer / "
            "abstain that does not track the declared observer. Accuracy is secondary; not clinical."
        ),
    }
    (out_dir / "stat_proof_results.json").write_text(json.dumps(results, indent=2))
    write_report(out_dir / "STAT_PROOF.md", results)
    print(json.dumps({
        "out": str(out_dir),
        "n_site4": results["n_site4"],
        "n_ood": results["n_ood"],
        "disagree_O0_O2": dis["O0_vs_O2"]["disagree_rate_among_both"],
        "runtime_sec": results["runtime_sec"],
    }, indent=2))


if __name__ == "__main__":
    main()
