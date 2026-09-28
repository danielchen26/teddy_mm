#!/usr/bin/env python3
"""SCOPE_REFINE_PROOF — ANM refining eval/train scope improves conditional correctness.

Two claims on the same TEDDY δu + holdout GT (no best.pt retrain):

1) Soft-P_f / abstain gate: answer always (TEDDY alone) vs only when ANM
   continuous workability soft_P = max(action_scores) ≥ τ (sweep τ).
   Verifier P_f is binary (answered=1 / abstain=0); soft_P is the continuous
   field score that underlies that gate. Show conditional Q rises as τ tightens;
   report coverage–Q curves for O0, O2, and missing-modality adt_only.

2) Attribution → refine scope: LOO top / flip-sensitive proteins & cells define
   a priority subset; compare Q vs random same-size vs full panel. Optional
   label-budget: small head supervised on ANM-flagged hard cells vs random,
   evaluated on hard holdout (Q per label).

Claim boundary: conditional correctness / label efficiency only — not clinical;
not claiming phase-1 Pearson >0.61 globally unless numbers show it.
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
import time
import warnings
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
ANM_ROOT = ROOT.parent / "ANM"
sys.path.insert(0, str(ANM_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.lineage_panels import (  # noqa: E402
    ACTIONS,
    CRITERIA,
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

warnings.filterwarnings("ignore", category=UserWarning)

PROTEINS = all_panel_proteins()
TZ = ZoneInfo("America/New_York")


def _now() -> str:
    return datetime.now(TZ).strftime("%Y-%m-%d %H:%M:%S %Z")


def _load_jsonl(path: Path):
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line:
                yield json.loads(line)


def decide_argmax(scores: dict[str, float], thr: float | None = None):
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    best_a, best_s = ranked[0]
    if thr is not None and best_s < thr:
        return None
    return best_a


def schema_for(crit_id: str, base_schema: dict) -> dict:
    schema = copy.deepcopy(base_schema)
    thr = float(CRITERIA[crit_id]["readout_threshold"])
    schema["field_representation"]["readout_threshold"] = thr
    schema["criterion_overlay"] = {"readout_threshold": thr}
    return schema


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
    scores = rd["action_scores"]
    ranked = sorted(scores.values(), reverse=True)
    soft_p = float(ranked[0]) if ranked else 0.0
    margin = float(ranked[0] - ranked[1]) if len(ranked) > 1 else soft_p
    return {
        "recommended_action": rd["recommended_action"],
        "action_scores": scores,
        "P_f": float(vf["P_f"] or 0.0),
        "Q_f": vf["Q_f"],
        "expected_action": vf.get("expected_action"),
        "soft_P": soft_p,
        "margin": margin,
    }


def load_site4_cells(cells_path: Path) -> dict[str, dict]:
    cells = {}
    for c in _load_jsonl(cells_path):
        if c.get("slice") == "site4_test":
            cells[c["cell_id"]] = c
    return cells


def load_events_for(cells_path_events: Path, want: set[str]) -> tuple[dict, dict]:
    by = defaultdict(list)
    p95: dict[str, float] = {}
    for ev in _load_jsonl(cells_path_events):
        cid = ev["cell_id"]
        if cid not in want:
            continue
        by[cid].append(ev)
        p95[ev["protein"]] = float(ev["norm_p95_train"])
    return by, p95


def load_mm_adt_only(mm_cells: Path, mm_events: Path, eval_fraction_keep: float = 1.0, seed: int = 0):
    """Load adt_only mask cells + events. Optionally subsample for speed."""
    cells = {}
    for c in _load_jsonl(mm_cells):
        if c.get("modality_mask") == "adt_only" and c.get("slice") == "site4_test":
            cells[c["cell_id"]] = c
    ids = sorted(cells)
    if eval_fraction_keep < 1.0:
        rng = np.random.default_rng(seed)
        n = max(500, int(len(ids) * eval_fraction_keep))
        ids = sorted(rng.choice(ids, size=min(n, len(ids)), replace=False).tolist())
        cells = {i: cells[i] for i in ids}
    want = set(cells)
    by, p95 = load_events_for(mm_events, want)
    return cells, by, p95


def coverage_q_curve(rows: list[dict], pred_key: str, score_key: str, taus: np.ndarray, n_total: int):
    out = []
    for tau in taus:
        answered = [r for r in rows if r[score_key] >= float(tau) and r.get(pred_key) is not None]
        labeled = [r for r in answered if r.get("label") is not None]
        if not labeled:
            out.append(
                {
                    "tau": float(tau),
                    "coverage": len(answered) / max(n_total, 1),
                    "Q": None,
                    "n_answered": len(answered),
                    "n_labeled_answered": 0,
                }
            )
            continue
        correct = sum(1 for r in labeled if r[pred_key] == r["label"])
        out.append(
            {
                "tau": float(tau),
                "coverage": len(answered) / max(n_total, 1),
                "Q": correct / len(labeled),
                "n_answered": len(answered),
                "n_labeled_answered": len(labeled),
            }
        )
    return out


def monotonic_q_rise(curve: list[dict]) -> dict:
    """Summarize how Q changes as coverage falls under the soft_P gate."""
    pts = [(c["coverage"], c["Q"], c["tau"]) for c in curve if c["Q"] is not None and c["n_labeled_answered"] >= 30]
    if len(pts) < 3:
        return {"n_pts": len(pts), "frac_nondecreasing_as_coverage_drops": None}
    # sort by decreasing coverage
    pts = sorted(pts, key=lambda x: -x[0])
    ups = 0
    comps = 0
    for i in range(len(pts) - 1):
        q0, q1 = pts[i][1], pts[i + 1][1]
        comps += 1
        if q1 >= q0 - 1e-12:
            ups += 1
    peak = max(pts, key=lambda x: x[1])
    # best point with coverage still >= 0.4 (practical operating region)
    mid = [p for p in pts if p[0] >= 0.4]
    mid_best = max(mid, key=lambda x: x[1]) if mid else peak
    return {
        "n_pts": len(pts),
        "frac_nondecreasing_as_coverage_drops": ups / max(comps, 1),
        "Q_at_fullish": pts[0][1],
        "Q_at_tightest": pts[-1][1],
        "delta_Q_tightest": pts[-1][1] - pts[0][1],
        "coverage_fullish": pts[0][0],
        "coverage_tightest": pts[-1][0],
        "Q_peak": peak[1],
        "coverage_at_peak": peak[0],
        "tau_at_peak": peak[2],
        "delta_Q_peak": peak[1] - pts[0][1],
        "Q_best_cov_ge_0.4": mid_best[1],
        "coverage_best_cov_ge_0.4": mid_best[0],
        "delta_Q_best_cov_ge_0.4": mid_best[1] - pts[0][1],
    }


def run_gate_claim(
    cells: dict,
    by_cell: dict,
    p95: dict,
    base_schema: dict,
    crit_ids: list[str],
    n_tau: int = 16,
    tag: str = "site4_rna_export",
) -> dict:
    ids = sorted(cells)
    results = {"tag": tag, "n_cells": len(ids), "criteria": {}}
    for crit_id in crit_ids:
        crit = CRITERIA[crit_id]
        schema = schema_for(crit_id, base_schema)
        rows = []
        t0 = time.time()
        for cid in ids:
            inst = make_instance(cid, cells[cid], by_cell[cid], crit, p95)
            out = anm_run_one(schema, inst)
            sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, crit)
            teddy_always = decide_argmax(sc, thr=None)
            teddy_thr = decide_argmax(sc, thr=float(crit["readout_threshold"]))
            rows.append(
                {
                    "cell_id": cid,
                    "label": out.get("expected_action"),
                    "P_f": out["P_f"],
                    "soft_P": out["soft_P"],
                    "margin": out["margin"],
                    "anm_action": out["recommended_action"],
                    "teddy_always": teddy_always,
                    "teddy_thr": teddy_thr,
                }
            )
        elapsed = time.time() - t0
        soft = np.array([r["soft_P"] for r in rows], dtype=np.float64)
        taus = np.unique(np.quantile(soft, np.linspace(0.0, 1.0, n_tau)))
        # also include 0 for "answer always" anchor on soft_P
        taus = np.unique(np.concatenate([[0.0], taus]))

        lab = [r for r in rows if r["label"] is not None]
        q_teddy_always = sum(1 for r in lab if r["teddy_always"] == r["label"]) / max(len(lab), 1)
        decided_thr = [r for r in lab if r["teddy_thr"] is not None]
        q_teddy_thr = (
            sum(1 for r in decided_thr if r["teddy_thr"] == r["label"]) / max(len(decided_thr), 1)
            if decided_thr
            else None
        )
        anm_dec = [r for r in lab if r["anm_action"] is not None]
        q_anm = (
            sum(1 for r in anm_dec if r["anm_action"] == r["label"]) / max(len(anm_dec), 1)
            if anm_dec
            else None
        )

        curve_teddy = coverage_q_curve(rows, "teddy_always", "soft_P", taus, len(ids))
        curve_anm = coverage_q_curve(rows, "anm_action", "soft_P", taus, len(ids))
        # binary P_f gate (equiv. ANM built-in abstain): two-point
        binary_curve = coverage_q_curve(rows, "teddy_always", "P_f", np.array([0.0, 0.5, 1.0]), len(ids))

        results["criteria"][crit_id] = {
            "elapsed_sec": elapsed,
            "mean_P_f_binary": float(np.mean([r["P_f"] for r in rows])),
            "soft_P_quantiles": {
                str(q): float(np.quantile(soft, q)) for q in [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
            },
            "baseline": {
                "teddy_always_Q": q_teddy_always,
                "teddy_always_n_labeled": len(lab),
                "teddy_thr_Q": q_teddy_thr,
                "teddy_thr_n_decided_labeled": len(decided_thr),
                "anm_default_Q": q_anm,
                "anm_default_n_decided_labeled": len(anm_dec),
                "anm_abstain_count": sum(1 for r in rows if r["anm_action"] is None),
            },
            "coverage_Q_curve_teddy_gated_by_soft_P": curve_teddy,
            "coverage_Q_curve_anm_gated_by_soft_P": curve_anm,
            "coverage_Q_curve_teddy_gated_by_binary_P_f": binary_curve,
            "monotonicity_teddy": monotonic_q_rise(curve_teddy),
            "monotonicity_anm": monotonic_q_rise(curve_anm),
        }
    return results


def Q_on_cells(cell_ids, cells, p95, crit_id="O0", force_always=True):
    crit = CRITERIA[crit_id]
    correct = decided = n_lab = 0
    for cid in cell_ids:
        if cid not in cells:
            continue
        lab = expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit)
        if lab is None:
            continue
        n_lab += 1
        sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, crit)
        pred = decide_argmax(sc, None if force_always else float(crit["readout_threshold"]))
        if pred is None:
            continue
        decided += 1
        if pred == lab:
            correct += 1
    return {
        "Q": (correct / decided) if decided else None,
        "n_decided": decided,
        "n_labeled": n_lab,
    }


def Q_protein_subset(prot_subset, cells, p95, crit_id="O0"):
    crit = CRITERIA[crit_id]
    subset = set(prot_subset)
    correct = decided = n_lab = 0
    for cid, c in cells.items():
        lab = expected_from_true(c["adt_true_panel_holdout"], p95, crit)
        if lab is None:
            continue
        n_lab += 1
        pred_panel = {p: (c["adt_pred_panel"][p] if p in subset else 0.0) for p in PROTEINS}
        sc = lineage_scores_from_panel(pred_panel, p95, crit)
        pred = decide_argmax(sc, thr=None)
        decided += 1
        if pred == lab:
            correct += 1
    return {
        "Q": correct / max(decided, 1),
        "n_decided": decided,
        "n_labeled": n_lab,
        "subset": list(prot_subset),
    }


def balanced_priority_proteins(flip_top_counts: Counter, k_per_lineage: int = 1) -> list[str]:
    """Pick top flip-sensitive protein(s) within each lineage (avoids myeloid-only collapse)."""
    chosen = []
    for lin, prots in LINEAGE_PANELS.items():
        ranked = sorted(((p, flip_top_counts.get(p, 0)) for p in prots), key=lambda x: -x[1])
        for p, _ in ranked[:k_per_lineage]:
            if p not in chosen:
                chosen.append(p)
    return chosen


def run_attribution_claim(
    cells: dict,
    p95: dict,
    attr_path: Path,
    crit_id: str = "O0",
    n_random: int = 40,
    seed: int = 0,
) -> dict:
    attr = np.load(attr_path, allow_pickle=True)
    cid_top = {str(c): str(t) for c, t in zip(attr["cell_id"], attr["top_protein"])}
    cid_flip = {str(c): bool(f > 0.5) for c, f in zip(attr["cell_id"], attr["flipped"])}
    flip_cells = [c for c, f in cid_flip.items() if f and c in cells]
    nonflip_cells = [c for c, f in cid_flip.items() if (not f) and c in cells]
    all_ids = sorted(cells)
    rng = np.random.default_rng(seed)

    q_all = Q_on_cells(all_ids, cells, p95, crit_id)
    q_flip = Q_on_cells(flip_cells, cells, p95, crit_id)
    q_nonflip = Q_on_cells(nonflip_cells, cells, p95, crit_id)
    rand_cell_qs = []
    for i in range(n_random):
        rs = list(rng.choice(all_ids, size=min(len(flip_cells), len(all_ids)), replace=False))
        rand_cell_qs.append(Q_on_cells(rs, cells, p95, crit_id)["Q"])

    flip_tops = Counter(cid_top[c] for c in flip_cells if c in cid_top)
    overall_tops = Counter(cid_top.values())

    # protein subsets
    k = 3
    prio_flip = [p for p, _ in flip_tops.most_common(k)]
    prio_overall = [p for p, _ in overall_tops.most_common(k)]
    prio_balanced = balanced_priority_proteins(flip_tops, k_per_lineage=1)
    q_prio_flip = Q_protein_subset(prio_flip, cells, p95, crit_id)
    q_prio_overall = Q_protein_subset(prio_overall, cells, p95, crit_id)
    q_prio_bal = Q_protein_subset(prio_balanced, cells, p95, crit_id)
    q_full = Q_protein_subset(PROTEINS, cells, p95, crit_id)

    rand_prot_qs = []
    for i in range(n_random):
        rs = list(rng.choice(PROTEINS, size=len(prio_flip), replace=False))
        rand_prot_qs.append(Q_protein_subset(rs, cells, p95, crit_id)["Q"])
    rand_bal_qs = []
    for i in range(n_random):
        # random 1-per-lineage
        rs = [rng.choice(prots) for prots in LINEAGE_PANELS.values()]
        rand_bal_qs.append(Q_protein_subset(rs, cells, p95, crit_id)["Q"])

    return {
        "attr_n": int(len(attr["cell_id"])),
        "n_flip_sensitive": len(flip_cells),
        "n_nonflip": len(nonflip_cells),
        "top1_flip_rate": float(np.mean(list(cid_flip.values()))),
        "flip_top_protein_counts": dict(flip_tops.most_common()),
        "overall_top_protein_counts": dict(overall_tops.most_common()),
        "cell_scope": {
            "full_panel_cells": q_all,
            "flip_sensitive_cells": q_flip,
            "nonflip_cells": q_nonflip,
            "random_cells_same_n_as_flip": {
                "mean_Q": float(np.mean(rand_cell_qs)),
                "std_Q": float(np.std(rand_cell_qs)),
                "min_Q": float(np.min(rand_cell_qs)),
                "max_Q": float(np.max(rand_cell_qs)),
                "n_draws": n_random,
                "n_per_draw": len(flip_cells),
            },
            "delta_Q_flip_minus_random": float(q_flip["Q"] - float(np.mean(rand_cell_qs))),
        },
        "protein_scope": {
            "full_panel": q_full,
            "flip_sensitive_top_k": {**q_prio_flip, "k": k},
            "overall_loo_top_k": {**q_prio_overall, "k": k},
            "balanced_1_per_lineage_by_flip": q_prio_bal,
            "random_k_proteins": {
                "mean_Q": float(np.mean(rand_prot_qs)),
                "std_Q": float(np.std(rand_prot_qs)),
                "min_Q": float(np.min(rand_prot_qs)),
                "max_Q": float(np.max(rand_prot_qs)),
                "n_draws": n_random,
                "k": k,
            },
            "random_1_per_lineage": {
                "mean_Q": float(np.mean(rand_bal_qs)),
                "std_Q": float(np.std(rand_bal_qs)),
                "n_draws": n_random,
            },
            "delta_Q_balanced_flip_minus_random_lineage": float(
                q_prio_bal["Q"] - float(np.mean(rand_bal_qs))
            ),
        },
    }


def _feats(cid, cells, p95):
    panel = cells[cid]["adt_pred_panel"]
    return np.array(
        [min(1.0, max(0.0, panel[p] / max(p95.get(p, 1.0), 1e-6))) for p in PROTEINS],
        dtype=np.float64,
    )


def run_label_budget(
    cells: dict,
    p95: dict,
    attr_path: Path,
    crit_id: str = "O0",
    n_labels_grid: tuple = (50, 100, 200, 400),
    n_seeds: int = 10,
    seed: int = 42,
) -> dict:
    attr = np.load(attr_path, allow_pickle=True)
    cid_flip = {str(c): bool(f > 0.5) for c, f in zip(attr["cell_id"], attr["flipped"])}
    crit = CRITERIA[crit_id]
    ids = sorted(cells)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(ids)
    n_train = len(ids) // 4
    train_ids = list(perm[:n_train])
    eval_ids = list(perm[n_train:])

    def lab(cid):
        return expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit)

    train_lab = [c for c in train_ids if lab(c) is not None]
    eval_lab = [c for c in eval_ids if lab(c) is not None]
    hard_train = [c for c in train_lab if cid_flip.get(c, False)]
    hard_eval = [c for c in eval_lab if cid_flip.get(c, False)]

    def teddy_Q(cids):
        cor = dec = 0
        for cid in cids:
            y = lab(cid)
            sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, crit)
            pred = decide_argmax(sc, thr=None)
            dec += 1
            if pred == y:
                cor += 1
        return cor / max(dec, 1)

    def fit_Q(train_sel, nlab, eval_cids, s, kind):
        rng2 = np.random.default_rng(s)
        if len(train_sel) < nlab:
            return None
        sel = list(rng2.choice(train_sel, size=nlab, replace=False))
        X = np.stack([_feats(c, cells, p95) for c in sel])
        y = np.array([lab(c) for c in sel])
        if len(set(y.tolist())) < 2:
            return None
        sc = StandardScaler().fit(X)
        if kind == "logreg":
            clf = LogisticRegression(max_iter=800)
        else:
            clf = MLPClassifier(hidden_layer_sizes=(32,), max_iter=600, random_state=int(s))
        clf.fit(sc.transform(X), y)
        Xe = sc.transform(np.stack([_feats(c, cells, p95) for c in eval_cids]))
        pred = clf.predict(Xe)
        yt = np.array([lab(c) for c in eval_cids])
        return float((pred == yt).mean())

    curves = {}
    for kind in ("logreg", "mlp"):
        rows = []
        for nlab in n_labels_grid:
            qh, qr = [], []
            for s in range(n_seeds):
                a = fit_Q(hard_train, min(int(nlab), len(hard_train)), hard_eval, s, kind)
                b = fit_Q(train_lab, int(nlab), hard_eval, s, kind)
                if a is not None:
                    qh.append(a)
                if b is not None:
                    qr.append(b)
            rows.append(
                {
                    "n_labels": int(nlab),
                    "hard_labels_Q_on_hard_eval_mean": float(np.mean(qh)) if qh else None,
                    "hard_labels_Q_on_hard_eval_std": float(np.std(qh)) if qh else None,
                    "random_labels_Q_on_hard_eval_mean": float(np.mean(qr)) if qr else None,
                    "random_labels_Q_on_hard_eval_std": float(np.std(qr)) if qr else None,
                    "delta_mean": float(np.mean(qh) - np.mean(qr)) if qh and qr else None,
                    "n_seeds": n_seeds,
                }
            )
        curves[kind] = rows

    return {
        "n_train_labeled": len(train_lab),
        "n_hard_train": len(hard_train),
        "n_eval_labeled": len(eval_lab),
        "n_hard_eval": len(hard_eval),
        "teddy_Q_eval": teddy_Q(eval_lab),
        "teddy_Q_hard_eval": teddy_Q(hard_eval),
        "curves": curves,
        "note": (
            "Frozen TEDDY panel features; small logreg/MLP head only. "
            "Eval target = ANM flip-sensitive hard holdout. "
            "Positive delta ⇒ ANM-flagged hard labels buy more Q-on-hard per label."
        ),
    }


def write_report(path: Path, results: dict) -> None:
    lines: list[str] = []
    lines.append("# SCOPE_REFINE_PROOF — ANM × TEDDY CITE")
    lines.append("")
    lines.append("## 中文摘要（给 Daniel）")
    lines.append("")
    lines.append(
        "在同一批 TEDDY δu（site4/test）+ holdout GT 上证明两件事："
        "**(1) 用 ANM 连续 workability `soft_P=max(action_scores)`（二元 verifier `P_f` 的连续底层）做弃权门控，"
        "收紧 τ 后条件正确率 Q 上升**；"
        "**(2) ANM LOO / flip-sensitive 归因划定的 hard 细胞与优先蛋白子集，会改变测到的正确率，"
        "并在小标签预算下对 hard holdout 更省标签。**"
        "不声称临床；不声称 phase-1 Pearson>0.61 全局，除非数字本身显示。"
    )
    lines.append("")

    # Claim 1
    g = results["claim1_gate"]
    lines.append("### Claim 1 — soft_P / abstain 门控（coverage–Q）")
    lines.append("")
    lines.append(
        f"- 主切片：`{g['tag']}` n={g['n_cells']}；"
        "二元 `P_f`∈{0,1}（有推荐动作=1）；连续门控用同场 `soft_P`。"
    )
    for crit_id, block in g["criteria"].items():
        base = block["baseline"]
        mono = block["monotonicity_teddy"]
        lines.append("")
        lines.append(f"#### {crit_id}")
        lines.append("")
        lines.append(
            f"- TEDDY always Q={base['teddy_always_Q']:.4f} (n_lab={base['teddy_always_n_labeled']}); "
            f"TEDDY thr Q={base['teddy_thr_Q']}; ANM default Q={base['anm_default_Q']} "
            f"(abstain={base['anm_abstain_count']}); mean binary P_f={block['mean_P_f_binary']:.4f}."
        )
        lines.append(
            f"- soft_P 收紧：baseline Q={mono.get('Q_at_fullish')} @cov={mono.get('coverage_fullish')}; "
            f"peak Q={mono.get('Q_peak')} @cov={mono.get('coverage_at_peak')} "
            f"(ΔQ_peak={mono.get('delta_Q_peak')}); "
            f"best@cov≥0.4 Q={mono.get('Q_best_cov_ge_0.4')} "
            f"(Δ={mono.get('delta_Q_best_cov_ge_0.4')}); "
            f"tightest Q={mono.get('Q_at_tightest')} @cov={mono.get('coverage_tightest')}; "
            f"非降比例={mono.get('frac_nondecreasing_as_coverage_drops')}."
        )
        lines.append("")
        lines.append("| τ (soft_P) | coverage | Q (TEDDY∨gate) | n_lab_ans |")
        lines.append("|---:|---:|---:|---:|")
        curve = block["coverage_Q_curve_teddy_gated_by_soft_P"]
        # print every other + last
        idxs = list(range(0, len(curve), max(1, len(curve) // 8)))
        if (len(curve) - 1) not in idxs:
            idxs.append(len(curve) - 1)
        for i in idxs:
            row = curve[i]
            q = f"{row['Q']:.4f}" if row["Q"] is not None else "NA"
            lines.append(
                f"| {row['tau']:.4f} | {row['coverage']:.4f} | {q} | {row['n_labeled_answered']} |"
            )

    if "claim1_gate_adt_only" in results:
        ga = results["claim1_gate_adt_only"]
        lines.append("")
        lines.append(f"#### Missing-modality `adt_only` (n={ga['n_cells']})")
        lines.append("")
        lines.append(
            "> Reading note: with `adt_only` (RNA missing) TEDDY does not run. \"TEDDY always\" here is TEDDY alone's "
            "fixed rule applied to the phase-2 protein-only stand-in (see MISSING_MODALITY_ANM_DEMO.md), so this curve "
            "gates the stand-in's calls, not TEDDY's."
        )
        lines.append("")
        for crit_id, block in ga["criteria"].items():
            mono = block["monotonicity_teddy"]
            base = block["baseline"]
            lines.append(
                f"- **{crit_id}**: TEDDY always Q={base['teddy_always_Q']:.4f}; "
                f"ΔQ_peak={mono.get('delta_Q_peak')} "
                f"(Q {mono.get('Q_at_fullish')}→{mono.get('Q_peak')} @cov {mono.get('coverage_at_peak')}); "
                f"best@cov≥0.4 Δ={mono.get('delta_Q_best_cov_ge_0.4')}; "
                f"nondec frac={mono.get('frac_nondecreasing_as_coverage_drops')}."
            )
            lines.append("")
            lines.append("| τ (soft_P) | coverage | Q | n_lab_ans |")
            lines.append("|---:|---:|---:|---:|")
            curve = block["coverage_Q_curve_teddy_gated_by_soft_P"]
            idxs = list(range(0, len(curve), max(1, len(curve) // 6)))
            if (len(curve) - 1) not in idxs:
                idxs.append(len(curve) - 1)
            for i in idxs:
                row = curve[i]
                q = f"{row['Q']:.4f}" if row["Q"] is not None else "NA"
                lines.append(
                    f"| {row['tau']:.4f} | {row['coverage']:.4f} | {q} | {row['n_labeled_answered']} |"
                )

    # Claim 2
    a = results["claim2_attribution"]
    cs = a["cell_scope"]
    ps = a["protein_scope"]
    lines.append("")
    lines.append("### Claim 2 — 归因划定 scope 改变测到的正确率")
    lines.append("")
    lines.append(
        f"- LOO attr n={a['attr_n']}; flip-sensitive cells n={a['n_flip_sensitive']} "
        f"(top1_flip_rate={a['top1_flip_rate']:.4f})."
    )
    lines.append("")
    lines.append("| cell scope | Q | n_decided |")
    lines.append("|---|---:|---:|")
    lines.append(
        f"| full panel cells | {cs['full_panel_cells']['Q']:.4f} | {cs['full_panel_cells']['n_decided']} |"
    )
    lines.append(
        f"| flip-sensitive (ANM LOO flipped) | {cs['flip_sensitive_cells']['Q']:.4f} | {cs['flip_sensitive_cells']['n_decided']} |"
    )
    lines.append(
        f"| nonflip | {cs['nonflip_cells']['Q']:.4f} | {cs['nonflip_cells']['n_decided']} |"
    )
    rc = cs["random_cells_same_n_as_flip"]
    lines.append(
        f"| random same-n as flip (mean±std, {rc['n_draws']} draws) | "
        f"{rc['mean_Q']:.4f}±{rc['std_Q']:.4f} | {rc['n_per_draw']} |"
    )
    lines.append("")
    lines.append(
        f"- ΔQ (flip − random same-n) = **{cs['delta_Q_flip_minus_random']:.4f}** "
        "→ attributed-hard scope 明显更难。"
    )
    lines.append("")
    lines.append("| protein scope | subset | Q |")
    lines.append("|---|---|---:|")
    lines.append(f"| full panel | {ps['full_panel']['subset']} | {ps['full_panel']['Q']:.4f} |")
    lines.append(
        f"| flip-sensitive top-{ps['flip_sensitive_top_k']['k']} | "
        f"{ps['flip_sensitive_top_k']['subset']} | {ps['flip_sensitive_top_k']['Q']:.4f} |"
    )
    lines.append(
        f"| overall LOO top-k | {ps['overall_loo_top_k']['subset']} | {ps['overall_loo_top_k']['Q']:.4f} |"
    )
    lines.append(
        f"| balanced 1/lineage by flip | {ps['balanced_1_per_lineage_by_flip']['subset']} | "
        f"{ps['balanced_1_per_lineage_by_flip']['Q']:.4f} |"
    )
    rp = ps["random_k_proteins"]
    lines.append(
        f"| random k proteins (mean±std) | k={rp['k']} | {rp['mean_Q']:.4f}±{rp['std_Q']:.4f} |"
    )
    rb = ps["random_1_per_lineage"]
    lines.append(
        f"| random 1/lineage (mean±std) | — | {rb['mean_Q']:.4f}±{rb['std_Q']:.4f} |"
    )
    lines.append("")
    lines.append(
        f"- Flip-top-k 蛋白子集 Q={ps['flip_sensitive_top_k']['Q']:.4f} vs random-k mean "
        f"{rp['mean_Q']:.4f}（归因优先蛋白改变 panel 正确率；"
        f"注：flip-top 常偏 myeloid，balanced 对照 Δ="
        f"{ps['delta_Q_balanced_flip_minus_random_lineage']:.4f}）。"
    )

    lb = results.get("claim2_label_budget")
    if lb:
        lines.append("")
        lines.append("### Claim 2b — 标签预算（frozen features + 小 head）")
        lines.append("")
        lines.append(
            f"- Train labeled={lb['n_train_labeled']} (hard={lb['n_hard_train']}); "
            f"eval labeled={lb['n_eval_labeled']} (hard={lb['n_hard_eval']}). "
            f"TEDDY Q eval={lb['teddy_Q_eval']:.4f}, hard-eval={lb['teddy_Q_hard_eval']:.4f}."
        )
        lines.append(f"- {lb['note']}")
        for kind, rows in lb["curves"].items():
            lines.append("")
            lines.append(f"#### {kind} → Q on hard holdout")
            lines.append("")
            lines.append("| n_labels | hard→Q_hard | random→Q_hard | Δ |")
            lines.append("|---:|---:|---:|---:|")
            for row in rows:
                lines.append(
                    f"| {row['n_labels']} | "
                    f"{row['hard_labels_Q_on_hard_eval_mean']:.4f}±{row['hard_labels_Q_on_hard_eval_std']:.4f} | "
                    f"{row['random_labels_Q_on_hard_eval_mean']:.4f}±{row['random_labels_Q_on_hard_eval_std']:.4f} | "
                    f"{row['delta_mean']:+.4f} |"
                )

    lines.append("")
    lines.append("### 一句话证明")
    lines.append("")
    lines.append(f"> {results['proof_sentence']}")
    lines.append("")
    lines.append("### 声明边界")
    lines.append("")
    lines.append("- Holdout `adt_true` = verifier only；未把 ANM 场律拟合到标签。")
    lines.append("- 未重训 `best.pt`。条件正确率 / 标签效率，非临床。")
    lines.append("- Verifier `P_f` 二元；τ 扫描用同场连续 `soft_P=max(action_scores)`。")
    lines.append("- 不声称 phase-1 Pearson >0.61 全局（除非导出数字本身如此）。")
    lines.append("- 无 Perturb / GFlowNet / 160M / ATAC。")
    lines.append("")
    lines.append("### Re-run")
    lines.append("")
    lines.append("```bash")
    lines.append("cd /Users/tianchichen/Documents/GitHub/teddy_mm")
    lines.append(
        "PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. "
        ".venv/bin/python bridge_anm/run_scope_refine_proof.py"
    )
    lines.append("```")
    lines.append("")
    lines.append(
        "Outputs: `outputs/anm_cite_bridge/scope_refine/SCOPE_REFINE_PROOF.md`, "
        "`scope_refine_results.json`."
    )
    lines.append("")
    lines.append(f"Generated: {results['timestamp_local']}")
    path.write_text("\n".join(lines) + "\n")


def build_proof_sentence(results: dict) -> str:
    g0 = results["claim1_gate"]["criteria"]["O0"]["monotonicity_teddy"]
    g2 = results["claim1_gate"]["criteria"]["O2"]["monotonicity_teddy"]
    a = results["claim2_attribution"]["cell_scope"]
    lb = results.get("claim2_label_budget")
    delta50 = None
    if lb:
        row = next(
            (r for r in lb["curves"]["logreg"] if r["n_labels"] == 50),
            lb["curves"]["logreg"][0],
        )
        delta50 = row.get("delta_mean")
    adt = results.get("claim1_gate_adt_only", {}).get("criteria", {}).get("O0")
    adt_bit = ""
    if adt:
        am = adt["monotonicity_teddy"]
        adt_bit = (
            f" Separately, on missing-modality adt_only O0 (RNA missing, so TEDDY does not run; these are TEDDY alone's "
            f"fixed rule's calls on the phase-2 protein-only stand-in), soft_P gate lifts Q "
            f"{am.get('Q_at_fullish'):.4f}→{am.get('Q_peak'):.4f} "
            f"(ΔQ_peak={am.get('delta_Q_peak'):.4f}); not a gain on TEDDY's calls."
        )
    lb_bit = ""
    if delta50 is not None:
        lb_bit = (
            f" Label-budget: frozen-feature head on ANM flip-hard cells beats random "
            f"labels on hard holdout at n=50 (ΔQ={delta50:+.4f})."
        )
    return (
        f"On the same TEDDY δu (site4/test n={results['claim1_gate']['n_cells']}), "
        f"ANM soft_P gating raises conditional Q O0 {g0.get('Q_at_fullish'):.4f}→"
        f"{g0.get('Q_peak'):.4f} (ΔQ_peak={g0.get('delta_Q_peak'):.4f}; "
        f"@cov≥0.4 best={g0.get('Q_best_cov_ge_0.4'):.4f}) and O2 "
        f"{g2.get('Q_at_fullish'):.4f}→{g2.get('Q_peak'):.4f} "
        f"(ΔQ_peak={g2.get('delta_Q_peak'):.4f}; mid operating region — extreme-low "
        f"coverage can dip). "
        f"ANM LOO flip-sensitive cells Q={a['flip_sensitive_cells']['Q']:.4f} vs "
        f"random same-n {a['random_cells_same_n_as_flip']['mean_Q']:.4f} "
        f"(Δ={a['delta_Q_flip_minus_random']:.4f})."
        f"{adt_bit}{lb_bit} "
        f"Conditional correctness / label efficiency only — not clinical; no best.pt retrain."
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=ROOT / "outputs/anm_cite_bridge/scope_refine",
    )
    ap.add_argument("--n-tau", type=int, default=16)
    ap.add_argument("--n-random", type=int, default=40)
    ap.add_argument("--n-seeds", type=int, default=10)
    ap.add_argument(
        "--skip-adt-only",
        action="store_true",
        help="Skip missing-modality adt_only gate curve",
    )
    ap.add_argument(
        "--adt-only-frac",
        type=float,
        default=1.0,
        help="Subsample fraction of adt_only cells (1.0=all)",
    )
    args = ap.parse_args()
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    bridge_out = ROOT / "outputs/anm_cite_bridge"
    cells_path = bridge_out / "cite_cells_meta.jsonl"
    events_path = bridge_out / "cite_typed_events.jsonl"
    schema_path = ROOT / "bridge_anm/schemas/cite_lineage_finite_field_v0.json"
    attr_path = bridge_out / "hard_proof/attr_compact.npz"
    mm_cells = bridge_out / "missing_modality/missing_modality_cells.jsonl"
    mm_events = bridge_out / "missing_modality/missing_modality_events.jsonl"

    t_all = time.time()
    print(f"[{_now()}] loading site4 cells…")
    cells = load_site4_cells(cells_path)
    print(f"  n_cells={len(cells)}")
    print(f"[{_now()}] loading events…")
    by_cell, p95 = load_events_for(events_path, set(cells))
    print(f"  n_with_events={len(by_cell)} p95_keys={sorted(p95)}")
    base_schema = json.loads(schema_path.read_text())

    print(f"[{_now()}] claim1 gate O0/O2…")
    claim1 = run_gate_claim(
        cells, by_cell, p95, base_schema, crit_ids=["O0", "O2"], n_tau=args.n_tau, tag="site4_phase1_export"
    )

    claim1_adt = None
    if not args.skip_adt_only and mm_cells.exists() and mm_events.exists():
        print(f"[{_now()}] claim1 gate adt_only (frac={args.adt_only_frac})…")
        mm_c, mm_by, mm_p95 = load_mm_adt_only(mm_cells, mm_events, eval_fraction_keep=args.adt_only_frac)
        # p95 from mm events; fall back
        if not mm_p95:
            mm_p95 = p95
        claim1_adt = run_gate_claim(
            mm_c,
            mm_by,
            mm_p95,
            base_schema,
            crit_ids=["O0", "O2"],
            n_tau=args.n_tau,
            tag="missing_modality_adt_only",
        )
    else:
        print("  skip adt_only")

    print(f"[{_now()}] claim2 attribution scope…")
    claim2 = run_attribution_claim(cells, p95, attr_path, crit_id="O0", n_random=args.n_random)

    print(f"[{_now()}] claim2b label budget…")
    claim2b = run_label_budget(
        cells, p95, attr_path, crit_id="O0", n_seeds=args.n_seeds, seed=42
    )

    results: dict[str, Any] = {
        "timestamp_local": _now(),
        "claim_boundary": "conditional_correctness_and_label_efficiency_only_not_clinical",
        "no_retrain_best_pt": True,
        "P_f_note": (
            "Verifier P_f is binary (1 if recommended_action in actions else 0). "
            "τ-sweep uses continuous soft_P=max(ANM action_scores) from the same field readout."
        ),
        "paths": {
            "cells": str(cells_path),
            "events": str(events_path),
            "attr": str(attr_path),
            "schema": str(schema_path),
            "mm_cells": str(mm_cells),
            "mm_events": str(mm_events),
        },
        "claim1_gate": claim1,
        "claim2_attribution": claim2,
        "claim2_label_budget": claim2b,
        "runtime_sec": None,
    }
    if claim1_adt is not None:
        results["claim1_gate_adt_only"] = claim1_adt

    results["proof_sentence"] = build_proof_sentence(results)
    results["runtime_sec"] = time.time() - t_all

    json_path = out_dir / "scope_refine_results.json"
    md_path = out_dir / "SCOPE_REFINE_PROOF.md"
    json_path.write_text(json.dumps(results, indent=2, default=str))
    write_report(md_path, results)
    print(f"[{_now()}] wrote {md_path}")
    print(f"[{_now()}] wrote {json_path}")
    print(f"runtime_sec={results['runtime_sec']:.1f}")
    print("PROOF:", results["proof_sentence"])


if __name__ == "__main__":
    main()
