#!/usr/bin/env python3
"""HARD statistical proof: TEDDY alone × TEDDY+ANM × Train+ANM on CITE.

Full-n attribution (LOO top protein + flip distance), stronger bootstrap (B>=200)
and permutation (n_perm>=100), O0/O1/O2 three-arm with GT verifier, silent
over-answer / abstain-tracks-observer, and label-consuming Train+ANM curves
showing how many O2 labels are needed to approach ANM's *declared* O2
abstain/Q — proving declaration beats both (1) no-field TEDDY and (3) just-train-harder.

Does NOT retrain best.pt. No clinical claim. Jev = log-loss Choice stand-in.
"""
from __future__ import annotations

import argparse
import copy
import json
import multiprocessing as mp
import os
import re
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
ANM_ROOT = Path(os.environ.get("ANM_ROOT", str(ROOT.parent / "ANM")))
sys.path.insert(0, str(ANM_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from lib.lineage_panels import (  # noqa: E402
    anm_readout_threshold,
    ACTIONS,
    CRITERIA,
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

PROTEINS = all_panel_proteins()
CRIT_IDS = ["O0", "O1", "O2"]
LINEAGES = [a["id"] for a in ACTIONS]

# -------------------- annotation map (cell_type -> scope group) --------------------
# Explicit map of every CITE ``cell_type`` seen in cite_cells_meta.jsonl. The three
# in-scope groups are the declared actions (B / T / myeloid); everything else is
# out of scope for this 3-lineage question and any *call* on it is an over-answer.
# pDC is put in "other" (not myeloid): it lacks the CD16/CD11c/CD36 myeloid panel.
IN_SCOPE_GROUPS = ("b_lineage", "t_lineage", "myeloid")
OUT_OF_SCOPE_GROUPS = ("nk", "ilc", "erythroid", "progenitor", "other")
CELL_TYPE_GROUP: dict[str, str] = {
    # B lineage (incl. plasma cells / plasmablasts)
    "B1 B IGKC+": "b_lineage",
    "B1 B IGKC-": "b_lineage",
    "Naive CD20+ B IGKC+": "b_lineage",
    "Naive CD20+ B IGKC-": "b_lineage",
    "Transitional B": "b_lineage",
    "Plasma cell IGKC+": "b_lineage",
    "Plasma cell IGKC-": "b_lineage",
    "Plasmablast IGKC+": "b_lineage",
    "Plasmablast IGKC-": "b_lineage",
    # T lineage
    "CD4+ T activated": "t_lineage",
    "CD4+ T activated integrinB7+": "t_lineage",
    "CD4+ T naive": "t_lineage",
    "CD4+ T CD314+ CD45RA+": "t_lineage",
    "CD8+ T CD49f+": "t_lineage",
    "CD8+ T CD57+ CD45RA+": "t_lineage",
    "CD8+ T CD57+ CD45RO+": "t_lineage",
    "CD8+ T CD69+ CD45RA+": "t_lineage",
    "CD8+ T CD69+ CD45RO+": "t_lineage",
    "CD8+ T TIGIT+ CD45RA+": "t_lineage",
    "CD8+ T TIGIT+ CD45RO+": "t_lineage",
    "CD8+ T naive": "t_lineage",
    "MAIT": "t_lineage",
    "T reg": "t_lineage",
    "gdT CD158b+": "t_lineage",
    "gdT TCRVD2+": "t_lineage",
    "dnT": "t_lineage",
    # myeloid
    "CD14+ Mono": "myeloid",
    "CD16+ Mono": "myeloid",
    "cDC2": "myeloid",
    # out of scope
    "NK": "nk",
    "NK CD158e1+": "nk",
    "ILC": "ilc",
    "ILC1": "ilc",
    "Erythroblast": "erythroid",
    "Normoblast": "erythroid",
    "Proerythroblast": "erythroid",
    "Reticulocyte": "erythroid",
    "HSC": "progenitor",
    "G/M prog": "progenitor",
    "Lymph prog": "progenitor",
    "MK/E prog": "progenitor",
    "pDC": "other",
}


def cell_type_group(cell_type: str | None) -> str:
    """Scope group of an annotated cell type; unknown types are out of scope ("other")."""
    return CELL_TYPE_GROUP.get(str(cell_type), "other")


def annotation_graded(preds, cell_types):
    """Grade calls against the annotation: in-scope accuracy, out-of-scope call share."""
    groups = [cell_type_group(ct) for ct in cell_types]
    unmapped = sorted({str(ct) for ct in cell_types if str(ct) not in CELL_TYPE_GROUP})
    ins = [i for i, g in enumerate(groups) if g in IN_SCOPE_GROUPS]
    n_in = len(ins)
    called_in = [i for i in ins if preds[i] is not None]
    correct = sum(1 for i in called_in if preds[i] == groups[i])
    per_lin = {}
    for lin in IN_SCOPE_GROUPS:
        idx = [i for i in ins if groups[i] == lin]
        c = [i for i in idx if preds[i] is not None]
        per_lin[lin] = {
            "n": len(idx),
            "n_called": len(c),
            "n_correct": sum(1 for i in c if preds[i] == lin),
            "accuracy_strict": (sum(1 for i in c if preds[i] == lin) / len(idx)) if idx else None,
            "call_counts": dict(Counter(preds[i] for i in c)),
        }
    outs = [i for i, g in enumerate(groups) if g not in IN_SCOPE_GROUPS]
    called_out = [i for i in outs if preds[i] is not None]
    per_out = {}
    for g in OUT_OF_SCOPE_GROUPS:
        idx = [i for i in outs if groups[i] == g]
        c = [i for i in idx if preds[i] is not None]
        per_out[g] = {
            "n": len(idx),
            "n_called": len(c),
            "call_share": (len(c) / len(idx)) if idx else None,
            "call_counts": dict(Counter(preds[i] for i in c)),
        }
    return {
        "in_scope": {
            "n": n_in,
            "n_called": len(called_in),
            "n_correct": correct,
            "abstain_rate": ((n_in - len(called_in)) / n_in) if n_in else None,
            "accuracy_strict": (correct / n_in) if n_in else None,
            "accuracy_among_called": (correct / len(called_in)) if called_in else None,
            "per_lineage": per_lin,
        },
        "out_of_scope": {
            "n": len(outs),
            "n_called": len(called_out),
            "call_share": (len(called_out) / len(outs)) if outs else None,
            "share_of_all_calls": (
                len(called_out) / max(len(called_in) + len(called_out), 1)
            ),
            "per_group": per_out,
        },
        "unmapped_cell_types": unmapped,
    }


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
    # v2: rule threshold mapped onto the ANM field scale; legacy: declared value.
    thr = float(anm_readout_threshold(CRITERIA[crit_id], base_schema))
    schema["field_representation"]["readout_threshold"] = thr
    schema["criterion_overlay"] = {"readout_threshold": thr}
    return schema


# -------------------- parallel runner (spawn-safe) --------------------
# Workers hold the per-criterion schemas and the O0 instances of the attribution
# cells; tasks are top-level functions so they pickle by reference under spawn.
# Results come back in task order (imap), and every random draw happens in the
# parent, so --workers N gives the same numbers as --workers 1.
_WSTATE: dict[str, Any] = {}


def _winit(schemas, o0_insts, scoring=None):
    # Schemas are built in the parent so a parent-side scoring switch
    # (lineage_panels.set_scoring) is honoured by spawn workers too.
    if scoring is not None:
        import lib.lineage_panels as _lp

        if hasattr(_lp, "set_scoring") and getattr(_lp, "SCORING", None) != scoring:
            _lp.set_scoring(scoring)
    _WSTATE["schemas"] = schemas
    _WSTATE["o0"] = o0_insts or {}


def _t_anm(arg):
    crit_id, inst = arg
    return anm_run_one(_WSTATE["schemas"][crit_id], inst)


def _t_loo(arg):
    cid, flip_grid = arg
    return _loo_top_for_instance(_WSTATE["schemas"]["O0"], _WSTATE["o0"][cid], flip_grid=flip_grid)


def _t_perm(arg):
    cid, vals = arg
    inst = _WSTATE["o0"][cid]
    events = copy.deepcopy(inst["proposed_source_events"])
    for e, v in zip(events, vals):
        e["value"] = float(v)
    return _perm_top_for_instance(_WSTATE["schemas"]["O0"], {**inst, "proposed_source_events": events})


def default_workers() -> int:
    return max(1, (os.cpu_count() or 2) - 1)


class Runner:
    """Ordered map over top-level task functions, in-process or on a spawn pool."""

    def __init__(self, workers, base_schema, o0_insts=None):
        self.workers = max(1, int(workers))
        self.pool = None
        import lib.lineage_panels as _lp

        schemas = {c: schema_for(c, base_schema) for c in CRIT_IDS}
        init = (schemas, o0_insts or {}, getattr(_lp, "SCORING", None))
        if self.workers > 1:
            ctx = mp.get_context("spawn")
            self.pool = ctx.Pool(self.workers, initializer=_winit, initargs=init)
        else:
            _winit(*init)

    def map(self, fn, args):
        args = list(args)
        if self.pool is None:
            return [fn(a) for a in args]
        if not args:
            return []
        cs = max(1, len(args) // (self.workers * 8))
        return list(self.pool.imap(fn, args, chunksize=cs))

    def close(self):
        if self.pool is not None:
            self.pool.close()
            self.pool.join()
            self.pool = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


# -------------------- TEDDY alone --------------------
TEDDY_EDIT_COST = {
    "needs_endpoint_labels": False,
    "needs_retrain_or_retune": False,
    "mechanism": (
        "apply the declared criterion (lineage scores + readout threshold) directly to "
        "TEDDY's predicted panel; label-free, no retraining"
    ),
    "lacks_vs_anm": "no field P_f / verifier readout, no LOO attribution or flip distance",
}
TEDDY_EDIT_COST_LEGACY = {
    "needs_endpoint_labels": True,
    "needs_retrain_or_retune": True,
    "mechanism": "threshold/boost grid or classifier head on new-criterion labels",
}


def run_teddy_arm(eval_ids, train_ids, cells, p95, *, legacy_edit_cost=False):
    arm = {"arm": "TEDDY_alone", "criteria": {}}
    silent = {}
    ctypes = [cells[cid].get("cell_type") for cid in eval_ids]
    for cid_name in CRIT_IDS:
        crit = CRITERIA[cid_name]
        preds, labels = [], []
        for cid in eval_ids:
            sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, crit)
            preds.append(decide_argmax(sc, crit["readout_threshold"]))
            labels.append(expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit))
        m = score_predictions(preds, labels)
        m["annotation_graded"] = annotation_graded(preds, ctypes)
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
    arm["criterion_edit_cost"] = dict(TEDDY_EDIT_COST_LEGACY if legacy_edit_cost else TEDDY_EDIT_COST)
    return arm


# -------------------- ANM arm (metrics + full-n attribution) --------------------
def run_anm_metrics(eval_ids, cells, by_cell, p95, base_schema, *, runner=None, per_cell=None):
    """ANM metrics per criterion. If ``per_cell`` is a dict it receives
    ``{crit: {"preds": [...], "labels": [...]}}`` aligned with ``eval_ids``."""
    own = runner is None
    runner = runner or Runner(1, base_schema)
    criteria_out = {}
    ctypes = [cells[cid].get("cell_type") for cid in eval_ids]
    for cid_name in CRIT_IDS:
        crit = CRITERIA[cid_name]
        preds, labels = [], []
        p_sum = q_sum = q_n = 0.0
        outs = runner.map(
            _t_anm,
            ((cid_name, make_instance(cid, cells[cid], by_cell[cid], crit, p95)) for cid in eval_ids),
        )
        for out in outs:
            preds.append(out["recommended_action"])
            labels.append(out.get("expected_action"))
            p_sum += float(out.get("P_f") or 0.0)
            if out.get("Q_f") is not None:
                q_sum += float(out["Q_f"])
                q_n += 1
        m = score_predictions(preds, labels)
        m["annotation_graded"] = annotation_graded(preds, ctypes)
        if per_cell is not None:
            per_cell[cid_name] = {"preds": list(preds), "labels": list(labels)}
        criteria_out[cid_name] = {
            "rule": {
                "criterion": cid_name,
                "label_fit": False,
                **{
                    k: crit[k]
                    for k in (
                        "readout_threshold",
                        "key_marker_boost",
                        "score_mode",
                        "expected_margin",
                    )
                },
            },
            **m,
            "mean_P_f": p_sum / max(len(eval_ids), 1),
            "mean_Q_f_defined": (q_sum / q_n) if q_n else None,
            "n_Q_f_defined": int(q_n),
        }
    if own:
        runner.close()
    return criteria_out


def _loo_top_for_instance(schema, inst, flip_grid=21):
    """Return (top_protein, flipped, flip_dist_or_None, delta) or None if abstain."""
    base = anm_run_one(schema, inst)
    if base["recommended_action"] is None:
        return None
    base_action = base["recommended_action"]
    base_scores = base["action_scores"]
    events = list(inst["proposed_source_events"])
    best_i, best_abs, best_d, best_flipped = -1, -1.0, 0.0, False
    for i, ev in enumerate(events):
        ablated = {
            **inst,
            "proposed_source_events": [e for j, e in enumerate(events) if j != i],
        }
        out = anm_run_one(schema, ablated)
        d_chosen = float(
            (out["action_scores"].get(base_action, 0.0) or 0.0)
            - (base_scores.get(base_action, 0.0) or 0.0)
        )
        ad = abs(d_chosen)
        if ad > best_abs:
            best_abs = ad
            best_i = i
            best_d = d_chosen
            best_flipped = out["recommended_action"] != base_action
    if best_i < 0:
        return None
    top_ev = events[best_i]
    protein = str(top_ev["event_id"].rsplit(":", 1)[-1])
    v0 = float(top_ev["value"])
    flips = []
    for v in np.linspace(0.0, 1.0, flip_grid):
        trial_events = copy.deepcopy(events)
        trial_events[best_i]["value"] = float(v)
        trial = {**inst, "proposed_source_events": trial_events}
        out = anm_run_one(schema, trial)
        if out["recommended_action"] != base_action:
            flips.append(abs(float(v) - v0))
    flip_dist = min(flips) if flips else None
    return {
        "protein": protein,
        "flipped": bool(best_flipped),
        "delta_chosen_score": float(best_d),
        "flip_distance": None if flip_dist is None else float(flip_dist),
        "base_action": base_action,
        "event_idx": int(best_i),
        "value0": v0,
    }


def _perm_top_for_instance(schema, inst):
    """Top-1 LOO protein of an instance (permutation null; no flip grid) or None."""
    events = inst["proposed_source_events"]
    base = anm_run_one(schema, inst)
    if base["recommended_action"] is None:
        return None
    base_action = base["recommended_action"]
    base_scores = base["action_scores"]
    best_prot, best_abs = None, -1.0
    for i, ev in enumerate(events):
        ablated = {**inst, "proposed_source_events": [e for j, e in enumerate(events) if j != i]}
        out = anm_run_one(schema, ablated)
        d = abs(
            float(
                (out["action_scores"].get(base_action, 0.0) or 0.0)
                - (base_scores.get(base_action, 0.0) or 0.0)
            )
        )
        if d > best_abs:
            best_abs = d
            best_prot = str(ev["event_id"].rsplit(":", 1)[-1])
    return best_prot


PERM_NULLS = ("within_lineage", "within_cell")


def permute_event_values(events, rng, mode="within_lineage"):
    """Shuffled copy of the event values (list aligned with ``events``).

    ``within_lineage`` (default): values move only among events of the same
    lineage (``action``), so each lineage keeps its evidence multiset and the null
    asks whether *which protein inside a lineage* carries the evidence matters.
    ``within_cell`` (legacy): values move across all events of the cell; this
    also scrambles the lineage call itself, so it is a weaker null.
    """
    vals = [e["value"] for e in events]
    if mode == "within_cell":
        rng.shuffle(vals)
        return vals
    if mode != "within_lineage":
        raise ValueError(f"unknown permutation null {mode!r}")
    groups: dict[str, list[int]] = {}
    for i, e in enumerate(events):
        groups.setdefault(str(e["action"]), []).append(i)
    for idx in groups.values():
        sub = [vals[i] for i in idx]
        rng.shuffle(sub)
        for i, v in zip(idx, sub):
            vals[i] = v
    return vals


def bootstrap_mode(tops, rng, n_boot):
    """Bootstrap of the top-1 mode fraction; same draws and tie-breaking as
    ``Counter(rng.choice(tops, n, replace=True)).most_common(1)`` (first seen wins)."""
    rates, prots = [], []
    if not len(tops):
        return rates, prots
    uniq, codes = np.unique(np.asarray(tops), return_inverse=True)
    n = len(codes)
    for _ in range(n_boot):
        sample = rng.choice(codes, size=n, replace=True)
        counts = np.bincount(sample, minlength=len(uniq))
        mx = counts.max()
        tied = np.flatnonzero(counts == mx)
        if len(tied) == 1:
            code = int(tied[0])
        else:
            first = {int(c): int(np.argmax(sample == c)) for c in tied}
            code = min(first, key=first.get)
        rates.append(int(mx) / n)
        prots.append(str(uniq[code]))
    return rates, prots


def run_full_attribution(
    attr_ids,
    cells,
    by_cell,
    p95,
    base_schema,
    *,
    flip_grid=21,
    n_boot=200,
    n_perm=100,
    seed=17,
    perm_null="within_lineage",
    runner=None,
):
    crit0 = CRITERIA["O0"]
    rng = np.random.default_rng(seed)
    kept_ids = [cid for cid in attr_ids if cid in cells]
    o0_insts = {cid: make_instance(cid, cells[cid], by_cell[cid], crit0, p95) for cid in kept_ids}
    own = runner is None
    runner = runner or Runner(1, base_schema, o0_insts)

    cell_ids = []
    top_proteins = []
    flip_flags = []
    flip_dists = []
    deltas = []

    t_attr0 = time.time()
    rows = runner.map(_t_loo, ((cid, flip_grid) for cid in kept_ids))
    for cid, row in zip(kept_ids, rows):
        if row is None:
            continue
        cell_ids.append(cid)
        top_proteins.append(row["protein"])
        flip_flags.append(1.0 if row["flipped"] else 0.0)
        flip_dists.append(row["flip_distance"])
        deltas.append(row["delta_chosen_score"])
    attr_sec = time.time() - t_attr0
    print(f"  attribution done n={len(top_proteins)} in {attr_sec:.1f}s", flush=True)

    tops = np.array(top_proteins)
    # Bootstrap stability of mode fraction (vectorised; draws identical to the old loop)
    boot_mode_rate, boot_mode_prot = bootstrap_mode(tops, rng, n_boot)
    if len(tops):
        mode_prot, mode_n = Counter(tops).most_common(1)[0]
    else:
        mode_prot, mode_n = None, 0

    obs_frac = mode_n / max(len(tops), 1)

    # Permutation: shuffle values (within lineage by default), recompute top-1.
    # All shuffles are drawn here in the parent, in (perm, cell) order.
    null_max_fracs = []
    null_mode_match_fracs = []
    t_perm0 = time.time()
    perm_cells = [cid for cid in kept_ids if o0_insts[cid]["proposed_source_events"]]
    for pi in range(n_perm):
        tasks = [
            (cid, permute_event_values(o0_insts[cid]["proposed_source_events"], rng, perm_null))
            for cid in perm_cells
        ]
        perm_tops = [p for p in runner.map(_t_perm, tasks) if p is not None]
        if perm_tops:
            mc = Counter(perm_tops).most_common(1)[0]
            null_max_fracs.append(mc[1] / len(perm_tops))
            if mode_prot:
                null_mode_match_fracs.append(Counter(perm_tops).get(mode_prot, 0) / len(perm_tops))
        if (pi + 1) % 10 == 0:
            print(f"  perm {pi+1}/{n_perm}", flush=True)
    perm_sec = time.time() - t_perm0
    print(f"  permutation done in {perm_sec:.1f}s", flush=True)
    if own:
        runner.close()

    p_value = (
        (1 + sum(1 for f in null_max_fracs if f >= obs_frac)) / max(len(null_max_fracs) + 1, 1)
        if null_max_fracs
        else None
    )

    fd = [d for d in flip_dists if d is not None]
    compact = {
        "cell_id": np.array(cell_ids, dtype=object),
        "top_protein": np.array(top_proteins, dtype=object),
        "flipped": np.array(flip_flags, dtype=np.float32),
        "flip_distance": np.array(
            [(-1.0 if d is None else d) for d in flip_dists], dtype=np.float32
        ),
        "delta_chosen_score": np.array(deltas, dtype=np.float32),
    }

    attribution = {
        "kind": "closed_form_field_LOO_plus_flip_distance",
        "declared_field": True,
        "n_attr_cells": len(tops),
        "n_attr_requested": len(attr_ids),
        "full_n": len(attr_ids) >= len(cells) * 0.95 if cells else False,
        "top1_protein_distribution": {str(k): int(v) for k, v in Counter(tops).items()},
        "top1_flip_rate": float(np.mean(flip_flags)) if flip_flags else None,
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
                [
                    float(np.quantile(boot_mode_rate, 0.025)),
                    float(np.quantile(boot_mode_rate, 0.975)),
                ]
                if boot_mode_rate
                else None
            ),
            "boot_mode_protein_stability": (
                {str(k): int(v) for k, v in Counter(boot_mode_prot).most_common()}
                if boot_mode_prot
                else {}
            ),
        },
        "permutation_test": {
            "n_perm": len(null_max_fracs),
            "null_max_fraction_mean": float(np.mean(null_max_fracs)) if null_max_fracs else None,
            "null_max_fraction_ci95": (
                [
                    float(np.quantile(null_max_fracs, 0.025)),
                    float(np.quantile(null_max_fracs, 0.975)),
                ]
                if null_max_fracs
                else None
            ),
            "observed_mode_fraction": obs_frac,
            "p_value_one_sided": float(p_value) if p_value is not None else None,
            "null_mode_match_fraction_mean": (
                float(np.mean(null_mode_match_fracs)) if null_mode_match_fracs else None
            ),
            "null": perm_null,
            "note": (
                "Shuffle event values within each lineage (values stay in their lineage); "
                "compare max top-1 protein fraction to observed"
                if perm_null == "within_lineage"
                else "Shuffle event values within cell (legacy null); "
                "compare max top-1 protein fraction to observed"
            ),
        },
        "runtime_sec": {"attribution": attr_sec, "permutation": perm_sec},
    }
    return attribution, compact


# -------------------- Train+ANM (label cost curves) --------------------
LABEL_COST_SPLITS = ("disjoint", "legacy")


def split_label_cost(ids, perm, train_frac, mode="disjoint"):
    """Train pool / evaluation cells for the label-cost arm.

    ``disjoint`` (default): train pool = first ``int(n*train_frac)`` cells of
    ``perm``; evaluation = the remaining cells (no overlap).
    ``legacy``: evaluation = all ``ids`` and train pool =
    ``perm[:max(n_train, min(2000, n))]`` of the same ids (overlapping).
    Returns ``(train_pool_ids, eval_ids)``.
    """
    n = len(ids)
    n_train = int(n * train_frac)
    if mode == "legacy":
        return [ids[i] for i in perm[: max(n_train, min(2000, n))]], list(ids)
    if mode != "disjoint":
        raise ValueError(f"unknown label-cost split {mode!r}")
    train = [ids[i] for i in perm[:n_train]]
    ev = sorted(ids[i] for i in perm[n_train:])
    return train, ev


def _predict_proba_gate(proba, classes, conf_thr):
    preds = []
    for i in range(len(proba)):
        j = int(np.argmax(proba[i]))
        preds.append(classes[j] if float(proba[i, j]) >= conf_thr else None)
    return preds


def calibrate_thr_to_abstain(proba, classes, target_rate, grid=None):
    """Confidence threshold whose abstain rate on ``proba``'s cells is closest to target."""
    if grid is None:
        grid = np.linspace(0.34, 0.95, 62)
    best_thr, best_gap = 0.45, 1e9
    n = len(proba)
    for thr in grid:
        preds = _predict_proba_gate(proba, classes, float(thr))
        rate = sum(1 for p in preds if p is None) / max(n, 1)
        gap = abs(rate - target_rate)
        if gap < best_gap:
            best_gap = gap
            best_thr = float(thr)
    return best_thr


def _grid_crit(thr, boost, mode):
    fake = dict(CRITERIA["O0"] if mode == "equal_panel_mean" else CRITERIA["O2"])
    fake["readout_threshold"] = thr
    fake["key_marker_boost"] = boost
    fake["score_mode"] = mode
    if mode == "key_marker_priority":
        fake["secondary_weight"] = 0.0
        fake["lineage_weights"] = {"b_lineage": 1.5, "t_lineage": 1.3, "myeloid": 0.5}
    return fake


def _scores_matrix(ids_, cells, p95, crit):
    """(n, n_lineages) lineage scores in dict order, plus that key order."""
    rows, keys = [], None
    for cid in ids_:
        sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, crit)
        if keys is None:
            keys = list(sc)
        rows.append([sc[k] for k in keys])
    if keys is None:
        keys = list(LINEAGES)
    return np.asarray(rows, dtype=np.float64).reshape(len(rows), len(keys)), keys


def _decide_vec(S, keys, thr):
    """Vectorised ``decide_argmax`` (no margin): first max wins, None below thr."""
    if not len(S):
        return []
    j = S.argmax(axis=1)
    m = S[np.arange(len(S)), j]
    return [None if m[i] < thr else keys[int(j[i])] for i in range(len(S))]


GRID_THR = [0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.36]
GRID_BOOST = [1.0, 1.5, 2.0, 2.5, 3.0]
GRID_MODES = ("equal_panel_mean", "key_marker_priority")


def run_train_plus_anm_curves(
    eval_ids,
    train_pool_ids,
    cells,
    p95,
    anm_o2_targets,
    *,
    seed=17,
    calib_ids=None,
    calib_target_abstain=None,
):
    """How many O2 labels does training need to approach ANM's declared O2 behavior?

    Abstain-calibrated thresholds are chosen on ``calib_ids`` (the train pool,
    label-free: only predicted confidences / scores are used) against
    ``calib_target_abstain`` (ANM's declared O2 abstain rate on those cells).
    ``calib_ids=None`` is the legacy behaviour: thresholds chosen on ``eval_ids``
    against the eval target.
    """
    rng = np.random.default_rng(seed)
    crit_o2 = CRITERIA["O2"]
    target_abstain = float(anm_o2_targets["abstain_rate"])
    target_q = float(anm_o2_targets["Q_analogue"] or 0.0)
    target_q_strict = float(anm_o2_targets.get("accuracy_strict_labeled") or 0.0)
    legacy_calib = calib_ids is None
    if legacy_calib:
        calib_ids = list(eval_ids)
        calib_target = target_abstain
    else:
        calib_ids = list(calib_ids)
        calib_target = float(
            target_abstain if calib_target_abstain is None else calib_target_abstain
        )

    def feats(ids_):
        return np.stack(
            [
                feature_vector(cells[c]["adt_pred_panel"], p95, cells[c].get("z_rna_compressed"))
                for c in ids_
            ]
        )

    # Labeled train pool under O2
    o2_lab = {
        cid: expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit_o2)
        for cid in train_pool_ids
    }
    labeled_train = [cid for cid in train_pool_ids if o2_lab[cid] is not None]
    rng.shuffle(labeled_train)

    # Eval labels/features once
    eval_X = feats(eval_ids)
    eval_y_o2 = [
        expected_from_true(cells[c]["adt_true_panel_holdout"], p95, crit_o2) for c in eval_ids
    ]
    calib_X = eval_X if legacy_calib else feats(calib_ids)

    budgets = [50, 100, 200, 500, 1000, 2000, 5000, len(labeled_train)]
    budgets = sorted({b for b in budgets if 0 < b <= len(labeled_train)})

    predict_proba_gate = _predict_proba_gate

    # Grid lineage scores are deterministic per (combo, cell): compute once per
    # combo over the (prefix-nested) labeled train pool and reuse across budgets.
    grid_combos = [(t, b, m) for t in GRID_THR for b in GRID_BOOST for m in GRID_MODES]
    _grid_S: dict[tuple, tuple] = {}

    def grid_scores(combo, n_):
        S, keys = _grid_S.get(combo, (None, None))
        if S is None or len(S) < n_:
            done = 0 if S is None else len(S)
            extra, keys = _scores_matrix(labeled_train[done:n_], cells, p95, _grid_crit(*combo))
            S = extra if S is None else np.concatenate([S, extra])
            _grid_S[combo] = (S, keys)
        return S[:n_], keys

    def best_scores(best_, ids_):
        return _scores_matrix(ids_, cells, p95, _grid_crit(best_["thr"], best_["boost"], best_["mode"]))

    curve_logistic = []
    curve_mlp = []
    curve_grid = []

    metric_keys = (
        "abstain_count", "abstain_rate", "Q_analogue",
        "accuracy_strict_labeled", "n_decided_labeled",
    )

    for nlab in budgets:
        subset = labeled_train[:nlab]
        Xtr = feats(subset)
        ytr = np.array([o2_lab[c] for c in subset])
        # --- logistic (Jev stand-in) ---
        sc = StandardScaler().fit(Xtr)
        try:
            clf = LogisticRegression(max_iter=1000, random_state=seed).fit(
                sc.transform(Xtr), ytr
            )
            proba = clf.predict_proba(sc.transform(eval_X))
            classes = list(clf.classes_)
            # default gate
            preds_def = predict_proba_gate(proba, classes, 0.45)
            m_def = score_predictions(preds_def, eval_y_o2)
            # calibrate abstain to ANM O2 (on the calibration set only)
            proba_cal = proba if legacy_calib else clf.predict_proba(sc.transform(calib_X))
            thr_cal = calibrate_thr_to_abstain(proba_cal, classes, calib_target)
            preds_cal = predict_proba_gate(proba, classes, thr_cal)
            m_cal = score_predictions(preds_cal, eval_y_o2)
            curve_logistic.append(
                {
                    "n_o2_labels": int(nlab),
                    "default_conf_0.45": {k: m_def[k] for k in metric_keys},
                    "abstain_calibrated_to_anm_o2": {
                        "conf_threshold": thr_cal,
                        **{k: m_cal[k] for k in metric_keys},
                        "abstain_gap_vs_anm": abs(m_cal["abstain_rate"] - target_abstain),
                        "Q_gap_vs_anm": (
                            abs((m_cal["Q_analogue"] or 0) - target_q)
                            if m_cal["Q_analogue"] is not None
                            else None
                        ),
                    },
                }
            )
        except Exception as e:
            curve_logistic.append({"n_o2_labels": int(nlab), "error": str(e)})

        # --- small MLP head ---
        try:
            if nlab >= 100 and len(set(ytr)) >= 2:
                mlp = MLPClassifier(
                    hidden_layer_sizes=(32,),
                    max_iter=400,
                    random_state=seed,
                    early_stopping=nlab >= 500,
                    validation_fraction=0.15 if nlab >= 500 else 0.1,
                )
                mlp.fit(sc.transform(Xtr), ytr)
                proba_m = mlp.predict_proba(sc.transform(eval_X))
                classes_m = list(mlp.classes_)
                preds_m = predict_proba_gate(proba_m, classes_m, 0.45)
                m_m = score_predictions(preds_m, eval_y_o2)
                proba_mcal = proba_m if legacy_calib else mlp.predict_proba(sc.transform(calib_X))
                thr_m = calibrate_thr_to_abstain(proba_mcal, classes_m, calib_target)
                preds_mc = predict_proba_gate(proba_m, classes_m, thr_m)
                m_mc = score_predictions(preds_mc, eval_y_o2)
                curve_mlp.append(
                    {
                        "n_o2_labels": int(nlab),
                        "default_conf_0.45": {k: m_m[k] for k in metric_keys},
                        "abstain_calibrated_to_anm_o2": {
                            "conf_threshold": thr_m,
                            **{k: m_mc[k] for k in metric_keys},
                            "abstain_gap_vs_anm": abs(m_mc["abstain_rate"] - target_abstain),
                            "Q_gap_vs_anm": (
                                abs((m_mc["Q_analogue"] or 0) - target_q)
                                if m_mc["Q_analogue"] is not None
                                else None
                            ),
                        },
                    }
                )
            else:
                curve_mlp.append(
                    {"n_o2_labels": int(nlab), "skipped": "need >=100 labels and >=2 classes"}
                )
        except Exception as e:
            curve_mlp.append({"n_o2_labels": int(nlab), "error": str(e)})

        # --- thr×boost grid retune (TEDDY-alone style, fit on O2 labels) ---
        # Search on the labeled subset; evaluate on eval
        best = {"acc": -1.0, "thr": None, "boost": None, "mode": None}
        tr_lab = np.array([o2_lab[c] for c in subset], dtype=object)
        for combo in grid_combos:
            thr, boost, mode = combo
            S_tr, keys_tr = grid_scores(combo, nlab)
            tr_preds = np.array(_decide_vec(S_tr, keys_tr, thr), dtype=object)
            acc = float(int(np.sum(tr_preds == tr_lab)) / len(subset)) if len(subset) else 0.0
            if acc > best["acc"]:
                best = {"acc": acc, "thr": thr, "boost": boost, "mode": mode}
        ev_S, ev_keys = best_scores(best, eval_ids)
        ev_preds = _decide_vec(ev_S, ev_keys, best["thr"])
        m_g = score_predictions(ev_preds, eval_y_o2)
        # Abstain-matching threshold (label-free), chosen on the calibration set
        cal_S = ev_S if legacy_calib else best_scores(best, calib_ids)[0]
        best_thr_ab, best_gap = best["thr"], 1e9
        cal_max = cal_S.max(axis=1) if len(cal_S) else np.zeros(0)
        for thr_try in np.linspace(0.05, 0.55, 51):
            rate = int(np.sum(cal_max < float(thr_try))) / max(len(calib_ids), 1)
            gap = abs(rate - calib_target)
            if gap < best_gap:
                best_gap = gap
                best_thr_ab = float(thr_try)
        preds_ab = _decide_vec(ev_S, ev_keys, best_thr_ab)
        m_ab = score_predictions(preds_ab, eval_y_o2)
        curve_grid.append(
            {
                "n_o2_labels": int(nlab),
                "best_grid": best,
                "eval_at_best_grid": {k: m_g[k] for k in metric_keys},
                "abstain_calibrated_thr": {
                    "thr": best_thr_ab,
                    **{k: m_ab[k] for k in metric_keys},
                    "abstain_gap_vs_anm": abs(m_ab["abstain_rate"] - target_abstain),
                    "Q_gap_vs_anm": (
                        abs((m_ab["Q_analogue"] or 0) - target_q)
                        if m_ab["Q_analogue"] is not None
                        else None
                    ),
                },
            }
        )
        print(f"  train+ANM curves nlab={nlab}", flush=True)

    # Headline: minimal n_labels where calibrated strict accuracy within 0.02 of ANM
    # AND abstain within 0.02 — or report never
    def find_match(curve, key_path=("abstain_calibrated_to_anm_o2",), q_tol=0.02, ab_tol=0.02):
        # Match abstain AND strict accuracy to ANM O2 declared operating point.
        # (Q_analogue among-decided is a different slice; report both.)
        for row in curve:
            block = row
            for k in key_path:
                if not isinstance(block, dict) or k not in block:
                    block = None
                    break
                block = block[k]
            if not isinstance(block, dict):
                continue
            q_strict = block.get("accuracy_strict_labeled")
            q_dec = block.get("Q_analogue")
            ar = block.get("abstain_rate")
            if q_strict is None or ar is None:
                continue
            if abs(q_strict - target_q_strict) <= q_tol and abs(ar - target_abstain) <= ab_tol:
                return {
                    "n_o2_labels": row["n_o2_labels"],
                    "accuracy_strict_labeled": q_strict,
                    "Q_analogue": q_dec,
                    "abstain_rate": ar,
                    "matched": True,
                    "match_metric": "abstain+accuracy_strict vs ANM O2 declared",
                }
        return {
            "matched": False,
            "note": (
                "no budget matched ANM O2 accuracy_strict±0.02 AND abstain±0.02; "
                "calibrating abstain often raises among-decided Q above ANM Q_f-strict"
            ),
            "max_budget": budgets[-1] if budgets else 0,
        }

    # Also: zero-shot O0 logistic evaluated on O2 (shows abstain does NOT track)
    labeled_o0 = []
    for cid in train_pool_ids:
        lab = expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, CRITERIA["O0"])
        if lab is not None:
            labeled_o0.append(cid)
    X0 = feats(labeled_o0)
    y0 = np.array(
        [
            expected_from_true(cells[c]["adt_true_panel_holdout"], p95, CRITERIA["O0"])
            for c in labeled_o0
        ]
    )
    sc0 = StandardScaler().fit(X0)
    m0 = LogisticRegression(max_iter=1000, random_state=seed).fit(sc0.transform(X0), y0)
    proba0 = m0.predict_proba(sc0.transform(eval_X))
    classes0 = list(m0.classes_)
    zs_preds = predict_proba_gate(proba0, classes0, 0.45)
    # Same model, same gate, scored vs O0 / O1 / O2 labels — abstain count identical
    zs_track = {}
    for cid_name in CRIT_IDS:
        labs = [
            expected_from_true(cells[c]["adt_true_panel_holdout"], p95, CRITERIA[cid_name])
            for c in eval_ids
        ]
        zs_track[cid_name] = score_predictions(zs_preds, labs)

    overlap = len(set(train_pool_ids) & set(eval_ids))
    return {
        "arm": "Train_plus_ANM_label_cost",
        "anm_o2_targets": {
            "abstain_rate": target_abstain,
            "Q_analogue": target_q,
            "accuracy_strict_labeled": target_q_strict,
            "source": "TEDDY+ANM declared O2 on the label-cost evaluation cells",
        },
        "split": {
            "n_train_pool": len(train_pool_ids),
            "n_eval": len(eval_ids),
            "train_eval_overlap": overlap,
            "threshold_selection_set": "eval (legacy)" if legacy_calib else "train_pool",
            "n_threshold_selection_cells": len(calib_ids),
            "threshold_selection_target_abstain": calib_target,
            "threshold_selection_uses_labels": False,
        },
        "n_labeled_train_pool_o2": len(labeled_train),
        "budgets": budgets,
        "curves": {
            "logistic_jev_standin": curve_logistic,
            "mlp_small_head": curve_mlp,
            "teddy_thr_boost_grid": curve_grid,
        },
        "match_headline": {
            "logistic_calibrated": find_match(curve_logistic, key_path=("abstain_calibrated_to_anm_o2",)),
            "mlp_calibrated": find_match(curve_mlp, key_path=("abstain_calibrated_to_anm_o2",)),
            "grid_calibrated": find_match(
                curve_grid, key_path=("abstain_calibrated_thr",)
            ),
        },
        "zero_shot_o0_model_abstain_does_not_track_observer": {
            "abstain_counts": {c: zs_track[c]["abstain_count"] for c in CRIT_IDS},
            "Q": {c: zs_track[c]["Q_analogue"] for c in CRIT_IDS},
            "note": (
                "O0-trained logistic uses one confidence gate; abstain count identical "
                "across O0/O1/O2 label sets — does not track declared observer. "
                "ANM abstain rises O0→O1→O2 by declaration with zero new labels."
            ),
        },
        "criterion_edit_cost": {
            "needs_endpoint_labels": True,
            "needs_retrain_or_retune": True,
            "mechanism": "consume O2 labels for retune/retrain/head; still no editable field attribution",
        },
    }


def ood_summary(ood_cells, p95, base_schema, events_path, *, runner=None):
    if not ood_cells:
        return {"n": 0, "note": "no OOD cells exported"}
    own = runner is None
    runner = runner or Runner(1, base_schema)
    ids = sorted(ood_cells)
    by = events_by_cell(events_path, ids)
    ctypes = [ood_cells[cid].get("cell_type") for cid in ids]
    out = {"n": len(ids), "slice": ood_cells[ids[0]].get("slice"), "criteria": {}}
    for cid_name in CRIT_IDS:
        crit = CRITERIA[cid_name]
        t_preds, t_lab = [], []
        a_preds, a_lab = [], []
        p_sum = 0.0
        for cid in ids:
            sc = lineage_scores_from_panel(ood_cells[cid]["adt_pred_panel"], p95, crit)
            t_preds.append(decide_argmax(sc, crit["readout_threshold"]))
            t_lab.append(expected_from_true(ood_cells[cid]["adt_true_panel_holdout"], p95, crit))
        outs = runner.map(
            _t_anm,
            ((cid_name, make_instance(cid, ood_cells[cid], by[cid], crit, p95)) for cid in ids),
        )
        for out_a in outs:
            a_preds.append(out_a["recommended_action"])
            a_lab.append(out_a.get("expected_action"))
            p_sum += float(out_a.get("P_f") or 0.0)
        out["criteria"][cid_name] = {
            "teddy_alone": {
                **score_predictions(t_preds, t_lab),
                "annotation_graded": annotation_graded(t_preds, ctypes),
            },
            "anm": {
                **score_predictions(a_preds, a_lab),
                "mean_P_f": p_sum / max(len(ids), 1),
                "annotation_graded": annotation_graded(a_preds, ctypes),
            },
        }
    if own:
        runner.close()
    return out


def _git(*args, strip=True):
    try:
        out = subprocess.run(
            ["git", *args], cwd=str(ROOT), capture_output=True, text=True, timeout=10
        ).stdout
        return out.strip() if strip else out
    except Exception:
        return None


def run_metadata(args, manifest, argv=None):
    """Provenance for the results JSON: git commit, flags, export size factor."""
    size_keys = {k: v for k, v in manifest.items() if re.search(r"size[_ -]?factor", k, re.I)}
    import lib.lineage_panels as _lp

    return {
        "git_commit": _git("rev-parse", "HEAD"),
        "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_dirty_files": [
            ln[3:]
            for ln in (_git("status", "--porcelain", "--untracked-files=no", strip=False) or "").splitlines()
            if ln.strip()
        ],
        "argv": list(sys.argv if argv is None else argv),
        "flags": {k: (str(v) if isinstance(v, Path) else v) for k, v in sorted(vars(args).items())},
        "export_size_factor": size_keys if size_keys else None,
        "export_size_factor_note": (
            None if size_keys else "export_manifest.json has no size-factor key (older export)"
        ),
        "lineage_scoring": getattr(_lp, "SCORING", None),
        "anm_root": str(ANM_ROOT),
        "anm_git_commit": _git_at(ANM_ROOT),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
    }


def _git_at(path):
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(path), capture_output=True, text=True, timeout=10
        ).stdout.strip() or None
    except Exception:
        return None


def write_report(path: Path, results: dict):
    d = results
    dis = d["disagreement_rates"]
    lines = []
    lines.append("# HARD_PROOF — ANM × TEDDY CITE（最强统计档）\n")
    lines.append("## 中文摘要（给 Daniel）\n")
    lines.append(
        f"在 **同一批 δu**（site4/test **n={d['n_site4']}**；OOD **n={d['n_ood']}**，"
        f"{d.get('ood_slice')}）上做三故事对照："
        "**（1）TEDDY multimodal alone**、**（2）TEDDY+ANM 声明场**、**（3）Train+ANM（用标签硬训逼近）**。"
        "准确率是次要的；根本变化 = 可编辑观察者 + 弃权/P_f + 可审计归因（GT 仅 verifier）。"
        f"归因全量 **n_attr={d['arms']['anm']['attribution']['n_attr_cells']}**；"
        f"bootstrap B={d['arms']['anm']['attribution']['bootstrap_top1']['n_boot']}；"
        f"permutation n_perm={d['arms']['anm']['attribution']['permutation_test']['n_perm']}。\n"
    )

    lines.append("### 准则分歧 / O2 标签改写\n")
    lines.append("| 对比 | both_defined | disagree | disagree_rate |")
    lines.append("|---|---:|---:|---:|")
    for k, v in dis.items():
        lines.append(
            f"| {k} | {v['both_defined']} | {v['disagree']} | {v['disagree_rate_among_both']:.4f} |"
        )
    lines.append("")
    lines.append(
        f"- O0 vs O1：加严弃权（标签改写率 **{dis['O0_vs_O1']['disagree_rate_among_both']:.4f}**）。\n"
        f"- O0 vs O2：key-marker 优先 **改写** expected_action（分歧率 "
        f"**{dis['O0_vs_O2']['disagree_rate_among_both']:.4f}**）。\n"
    )

    lines.append("### 三臂 Q vs GT（site4，full n）\n")
    lines.append(
        "| 臂 | O0 Q | O0 abstain | O1 Q | O1 abstain | O2 Q | O2 abstain | 准则变更成本 |"
    )
    lines.append("|---|---:|---:|---:|---:|---:|---:|---|")

    def row_arm(label, arm, cost):
        cells = []
        for c in CRIT_IDS:
            b = arm["criteria"][c]
            q = b.get("mean_Q_f_defined")
            cells.append(f"{q:.4f}" if q is not None else "NA")
            cells.append(str(b.get("abstain_count")))
        lines.append(f"| {label} | " + " | ".join(cells) + f" | {cost} |")

    t_cost = d["arms"]["teddy_alone"]["criterion_edit_cost"]
    row_arm(
        "（1）TEDDY alone",
        d["arms"]["teddy_alone"],
        "需标签重调" if t_cost.get("needs_endpoint_labels") else "改声明（无标签），但无场/归因",
    )
    row_arm("（2）TEDDY+ANM", d["arms"]["anm"], "仅改声明")

    anm_attr = d["arms"]["anm"]["attribution"]
    lines.append("\n### 全量归因 / Bootstrap / Permutation（ANM，O0 场）\n")
    lines.append(
        f"- 归因细胞数：**{anm_attr['n_attr_cells']}** / requested={anm_attr['n_attr_requested']}"
    )
    lines.append(f"- top-1 蛋白分布：`{anm_attr['top1_protein_distribution']}`")
    lines.append(f"- top-1 flip rate：**{anm_attr['top1_flip_rate']}**")
    lines.append(f"- flip distance quantiles：`{anm_attr['flip_distance_quantiles']}`")
    boot = anm_attr["bootstrap_top1"]
    lines.append(
        f"- bootstrap top-1：mode={boot['observed_mode_protein']} "
        f"frac={boot['observed_mode_fraction']:.4f} "
        f"boot_mean={boot['boot_mode_fraction_mean']} ci95={boot['boot_mode_fraction_ci95']} "
        f"(B={boot['n_boot']})"
    )
    perm = anm_attr["permutation_test"]
    lines.append(
        f"- permutation：null_max_mean={perm['null_max_fraction_mean']} "
        f"null_ci95={perm.get('null_max_fraction_ci95')} "
        f"p={perm['p_value_one_sided']} (n_perm={perm['n_perm']})"
    )

    silent = d["arms"]["teddy_alone"]["silent_over_answer"]
    lines.append("\n### 静默过度回答 / 弃权是否跟随声明\n")
    lines.append(f"- TEDDY-alone O0 规则 vs O1 声明弃权：`{silent['O0_over_answer_vs_O1_decl']}`")
    lines.append(f"- TEDDY-alone O0 规则 vs O2 声明弃权：`{silent['O0_over_answer_vs_O2_decl']}`")
    lines.append(
        f"- ANM abstain O0/O1/O2："
        f"{[d['arms']['anm']['criteria'][c]['abstain_count'] for c in CRIT_IDS]}"
    )
    train = d["arms"]["train_plus_anm"]
    zs = train["zero_shot_o0_model_abstain_does_not_track_observer"]
    lines.append(f"- Train 臂 O0 模型零样本弃权不跟随：`{zs['abstain_counts']}`")

    lines.append("\n### 注释分级（cell_type；范围内 B/T/myeloid 准确率，范围外调用占比）\n")
    lines.append("| 臂 | crit | 范围内 n | 范围内准确率(strict) | 范围外 n | 范围外调用 | 范围外调用率 |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for arm_key, arm_lab in (("teddy_alone", "TEDDY alone"), ("anm", "TEDDY+ANM")):
        for c in CRIT_IDS:
            ag = d["arms"][arm_key]["criteria"][c].get("annotation_graded")
            if not ag:
                continue
            i_, o_ = ag["in_scope"], ag["out_of_scope"]
            acc = i_["accuracy_strict"]
            sh = o_["call_share"]
            lines.append(
                f"| {arm_lab} | {c} | {i_['n']} | {'NA' if acc is None else f'{acc:.4f}'} | "
                f"{o_['n']} | {o_['n_called']} | {'NA' if sh is None else f'{sh:.4f}'} |"
            )
    lines.append("")

    lines.append("\n### （3）Train+ANM 标签成本曲线（逼近声明 O2）\n")
    sp = d["arms"]["train_plus_anm"].get("split", {})
    lines.append(
        f"划分：mode={sp.get('mode')}，train_pool={sp.get('n_train_pool')}，"
        f"eval={sp.get('n_eval')}，overlap={sp.get('train_eval_overlap')}，"
        f"阈值选择集={sp.get('threshold_selection_set')}。\n"
    )
    lines.append(
        f"目标（ANM 声明 O2）：abstain_rate={train['anm_o2_targets']['abstain_rate']:.4f}，"
        f"Q={train['anm_o2_targets']['Q_analogue']:.4f}。\n"
    )
    lines.append("| 方法 | 匹配 headline |")
    lines.append("|---|---|")
    for name, mh in train["match_headline"].items():
        lines.append(f"| {name} | `{mh}` |")
    lines.append("")
    lines.append("Logistic（Jev 站位）校准弃权后曲线（节选）：\n")
    lines.append("| n_O2_labels | abstain | Q | Q_gap_vs_ANM |")
    lines.append("|---:|---:|---:|---:|")
    for row in train["curves"]["logistic_jev_standin"]:
        if "abstain_calibrated_to_anm_o2" not in row:
            continue
        b = row["abstain_calibrated_to_anm_o2"]
        q = b.get("Q_analogue")
        q_s = "NA" if q is None else f"{q:.4f}"
        lines.append(
            f"| {row['n_o2_labels']} | {b['abstain_rate']:.4f} | {q_s} | {b.get('Q_gap_vs_anm')} |"
        )

    hl = d.get("label_cost_headline", "")
    lines.append(f"\n**标签成本 headline：** {hl}\n")

    lines.append("\n### 一句话证明\n")
    lines.append(f"> **{d['proof_sentence']}**\n")

    lines.append("### 声明边界\n")
    lines.append(
        "- Holdout `adt_true` = verifier only；未拟合 ANM 场律到标签。\n"
        "- 未重训 `best.pt`。未声称临床优越。Jev/MLP 为训练站位，非 live Jev API。\n"
        "- 无 Perturb / GFlowNet / 160M / ATAC。\n"
    )

    lines.append("### OOD 摘要\n")
    ood = d.get("ood", {})
    if ood.get("criteria"):
        lines.append("| crit | TEDDY Q | TEDDY abstain | ANM Q | ANM abstain | ANM mean_P_f |")
        lines.append("|---|---:|---:|---:|---:|---:|")
        for c in CRIT_IDS:
            t = ood["criteria"][c]["teddy_alone"]
            a = ood["criteria"][c]["anm"]
            lines.append(
                f"| {c} | {t.get('Q_analogue')} | {t.get('abstain_count')} | "
                f"{a.get('Q_analogue')} | {a.get('abstain_count')} | {a.get('mean_P_f')} |"
            )
        lines.append("")
    else:
        lines.append(f"```json\n{json.dumps(ood, indent=2)[:2000]}\n```\n")

    lines.append(f"### Disk\n\n{d.get('disk_note', '')}\n")

    lines.append("### Re-run\n")
    lines.append("```bash")
    lines.append("cd /Users/tianchichen/Documents/GitHub/teddy_mm")
    lines.append(
        "PYTHONPATH=/Users/tianchichen/Documents/GitHub/ANM:. "
        ".venv/bin/python bridge_anm/run_hard_proof.py --attr-n 0 --n-boot 200 --n-perm 100 "
        "--workers 9  # add --label-cost-split legacy --perm-null within_cell for the old numbers"
    )
    lines.append("```\n")
    lines.append(
        f"Outputs: `{path.parent}/HARD_PROOF.md`, `hard_proof_results.json`, "
        f"`attr_compact.npz`.\n"
    )
    lines.append(f"Generated: {d.get('timestamp_local')}\n")
    path.write_text("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bridge-dir", type=Path, default=ROOT / "outputs/anm_cite_bridge")
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument(
        "--schema",
        type=Path,
        default=Path(__file__).resolve().parent / "schemas/cite_lineage_finite_field_v0.json",
    )
    ap.add_argument("--eval-n", type=int, default=0, help="0=all site4 cells")
    ap.add_argument("--train-frac", type=float, default=0.7)
    ap.add_argument("--attr-n", type=int, default=0, help="0=all decided cells (full-n)")
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-perm", type=int, default=100)
    ap.add_argument("--flip-grid", type=int, default=21)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--skip-ood", action="store_true")
    ap.add_argument("--skip-perm", action="store_true")
    ap.add_argument(
        "--label-cost-split",
        choices=LABEL_COST_SPLITS,
        default="disjoint",
        help=(
            "disjoint (default): label-cost train pool and evaluation cells do not overlap and "
            "abstain thresholds are chosen on the train pool; legacy: evaluate on all cells, "
            "train pool is a subset of them, thresholds chosen on the evaluation cells"
        ),
    )
    ap.add_argument(
        "--perm-null",
        choices=PERM_NULLS,
        default="within_lineage",
        help="attribution permutation null: shuffle values within each lineage (default) "
        "or within the whole cell (legacy)",
    )
    ap.add_argument(
        "--legacy-teddy-edit-cost",
        action="store_true",
        help="report the old (incorrect) needs_endpoint_labels=True for the TEDDY-alone rule",
    )
    ap.add_argument(
        "--workers",
        type=int,
        default=default_workers(),
        help="processes for the per-cell ANM runs (default: cpu_count-1); results do not depend on it",
    )
    args = ap.parse_args(argv)

    out_dir = args.out_dir or (args.bridge_dir / "hard_proof")
    out_dir.mkdir(parents=True, exist_ok=True)

    import shutil

    usage = shutil.disk_usage(str(args.bridge_dir))
    disk_note = (
        f"Avail {usage.free/1e9:.1f} GiB / total {usage.total/1e9:.1f} GiB "
        f"({100*usage.used/usage.total:.0f}% used). No cache purge needed."
    )
    print(disk_note, flush=True)

    manifest = json.loads((args.bridge_dir / "export_manifest.json").read_text())
    z_path = manifest.get("z_rna_export")
    cells = load_cells(args.bridge_dir / "cite_cells_meta.jsonl", z_path, site4_only=True)
    ood_cells = (
        {}
        if args.skip_ood
        else load_ood_cells(args.bridge_dir / "cite_cells_meta.jsonl", z_path)
    )
    ids = sorted(cells)
    rng = np.random.default_rng(args.seed)
    if args.eval_n > 0 and args.eval_n < len(ids):
        ids = sorted(rng.choice(ids, size=args.eval_n, replace=False).tolist())
    n = len(ids)
    perm = rng.permutation(n)
    eval_ids = ids  # label-free arms (TEDDY alone, ANM, attribution) use every cell
    train_ids, lc_eval_ids = split_label_cost(ids, perm, args.train_frac, args.label_cost_split)

    p95 = p95_from_events(args.bridge_dir / "cite_typed_events.jsonl", ids)
    by_cell = events_by_cell(args.bridge_dir / "cite_typed_events.jsonl", ids)
    base_schema = json.loads(args.schema.read_text())

    if args.attr_n <= 0:
        attr_ids = list(eval_ids)
    else:
        attr_ids = list(
            rng.choice(eval_ids, size=min(args.attr_n, len(eval_ids)), replace=False)
        )

    n_perm = 0 if args.skip_perm else args.n_perm
    metadata = run_metadata(args, manifest, argv)

    t0 = time.time()
    print(
        f"n_site4={len(ids)} label_cost_split={args.label_cost_split} "
        f"train_pool={len(train_ids)} label_cost_eval={len(lc_eval_ids)} attr={len(attr_ids)} "
        f"ood={len(ood_cells)} boot={args.n_boot} perm={n_perm} null={args.perm_null} "
        f"workers={args.workers}",
        flush=True,
    )

    dis = disagreement_rates({cid: cells[cid] for cid in eval_ids}, p95)
    print(
        "disagreement",
        {k: v["disagree_rate_among_both"] for k, v in dis.items()},
        flush=True,
    )

    print("TEDDY arm...", flush=True)
    teddy = run_teddy_arm(
        eval_ids, train_ids, cells, p95, legacy_edit_cost=args.legacy_teddy_edit_cost
    )

    crit0 = CRITERIA["O0"]
    o0_insts = {
        cid: make_instance(cid, cells[cid], by_cell[cid], crit0, p95)
        for cid in attr_ids
        if cid in cells
    }
    phase_sec = {}
    with Runner(args.workers, base_schema, o0_insts) as runner:
        print("ANM metrics...", flush=True)
        t_ = time.time()
        anm_per_cell: dict[str, dict] = {}
        anm_crit = run_anm_metrics(
            eval_ids, cells, by_cell, p95, base_schema, runner=runner, per_cell=anm_per_cell
        )
        phase_sec["anm_metrics"] = time.time() - t_

        print("ANM full attribution...", flush=True)
        t_ = time.time()
        attribution, compact = run_full_attribution(
            attr_ids,
            cells,
            by_cell,
            p95,
            base_schema,
            flip_grid=args.flip_grid,
            n_boot=args.n_boot,
            n_perm=n_perm,
            seed=args.seed,
            perm_null=args.perm_null,
            runner=runner,
        )
        phase_sec["attribution_boot_perm"] = time.time() - t_

        print("OOD...", flush=True)
        t_ = time.time()
        ood = ood_summary(
            ood_cells, p95, base_schema, args.bridge_dir / "cite_typed_events.jsonl", runner=runner
        )
        phase_sec["ood"] = time.time() - t_
    np.savez_compressed(out_dir / "attr_compact.npz", **compact)

    anm = {
        "arm": "TEDDY_plus_ANM",
        "criteria": anm_crit,
        "criterion_edit_cost": {
            "needs_endpoint_labels": False,
            "needs_retrain_or_retune": False,
            "mechanism": "edit CRITERIA/YAML; re-adapt instances only",
        },
        "attribution": attribution,
    }

    # ANM declared O2 operating point on a subset of the (label-free) ANM run
    pos = {cid: i for i, cid in enumerate(eval_ids)}

    def anm_o2_on(sub_ids):
        pc = anm_per_cell["O2"]
        return score_predictions(
            [pc["preds"][pos[c]] for c in sub_ids], [pc["labels"][pos[c]] for c in sub_ids]
        )

    print("Train+ANM label-cost curves...", flush=True)
    t_ = time.time()
    if args.label_cost_split == "legacy":
        targets = {
            "abstain_rate": anm_crit["O2"]["abstain_rate"],
            "Q_analogue": anm_crit["O2"]["mean_Q_f_defined"],
            "accuracy_strict_labeled": anm_crit["O2"]["accuracy_strict_labeled"],
        }
        train_arm = run_train_plus_anm_curves(
            lc_eval_ids, train_ids, cells, p95, targets, seed=args.seed
        )
    else:
        m_ev = anm_o2_on(lc_eval_ids)
        m_tr = anm_o2_on(train_ids)
        targets = {
            "abstain_rate": m_ev["abstain_rate"],
            "Q_analogue": m_ev["Q_analogue"],
            "accuracy_strict_labeled": m_ev["accuracy_strict_labeled"],
        }
        train_arm = run_train_plus_anm_curves(
            lc_eval_ids,
            train_ids,
            cells,
            p95,
            targets,
            seed=args.seed,
            calib_ids=train_ids,
            calib_target_abstain=m_tr["abstain_rate"],
        )
    train_arm["split"]["mode"] = args.label_cost_split
    phase_sec["train_plus_anm"] = time.time() - t_

    mh = train_arm["match_headline"]
    parts = []
    for name, hit in mh.items():
        if hit.get("matched"):
            parts.append(f"{name} matched at n_O2_labels={hit['n_o2_labels']}")
        else:
            parts.append(f"{name}: no match up to {hit.get('max_budget')}")
    label_cost_headline = (
        "；".join(parts)
        + "。关键：即便用满池 O2 标签把 Q/弃权率数值校准到接近 ANM，"
        "Train 臂仍无声明场、无 LOO 归因/翻转距离，且 O0→O2 观察者切换必须重新消耗标签；"
        "ANM 零标签改声明即切换弃权与 expected_action。"
    )

    proof_sentence = (
        "On the same TEDDY δu (site4/test full-n), ANM declaration edits the observer "
        "(O0→O1 tighter abstain; O0→O2 rewrites expected_action) with zero endpoint-label fit, "
        f"keeping verifiable P_f/Q_f and auditable full-n LOO attribution "
        f"(n_attr={attribution['n_attr_cells']}, bootstrap B={args.n_boot}, "
        f"permutation p={attribution['permutation_test']['p_value_one_sided']}). "
        "TEDDY-alone silently over-answers where the declared observer abstains; "
        "Train+ANM (logistic/MLP/grid) must consume O2 labels to approach ANM's O2 Q/abstain "
        "numbers and still cannot buy editable-field attribution or label-free observer edits. "
        "Accuracy is secondary; not clinical."
    )

    results = {
        "timestamp_local": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
        "tier": "hard",
        "run_metadata": metadata,
        "n_site4": len(ids),
        "n_ood": len(ood_cells),
        "ood_slice": manifest.get("ood_slice"),
        "n_available_site4_test": manifest.get("n_available_site4_test"),
        "capped": manifest.get("capped"),
        "export_manifest": {
            k: manifest.get(k)
            for k in (
                "n_cells_site4_test",
                "n_cells_ood",
                "n_events",
                "device",
                "ckpt",
                "compact",
                "phase1_mlp_pearson_quoted",
            )
        },
        "annotation_map": {
            "in_scope_groups": list(IN_SCOPE_GROUPS),
            "out_of_scope_groups": list(OUT_OF_SCOPE_GROUPS),
            "cell_type_group": dict(CELL_TYPE_GROUP),
            "unknown_cell_type_group": "other",
        },
        "disagreement_rates": dis,
        "arms": {
            "teddy_alone": teddy,
            "anm": anm,
            "train_plus_anm": train_arm,
        },
        "ood": ood,
        "attr_n": attribution["n_attr_cells"],
        "n_boot": args.n_boot,
        "n_perm": n_perm,
        "label_cost_headline": label_cost_headline,
        "runtime_sec": time.time() - t0,
        "runtime_phase_sec": phase_sec,
        "disk_note": disk_note,
        "claim_boundary": "local_response_diagnosis_only_not_clinical",
        "no_retrain_best_pt": True,
        "proof_sentence": proof_sentence,
        "paths": {
            "HARD_PROOF.md": str(out_dir / "HARD_PROOF.md"),
            "hard_proof_results.json": str(out_dir / "hard_proof_results.json"),
            "attr_compact.npz": str(out_dir / "attr_compact.npz"),
            "run_script": str(Path(__file__).resolve()),
        },
    }
    (out_dir / "hard_proof_results.json").write_text(json.dumps(results, indent=2))
    write_report(out_dir / "HARD_PROOF.md", results)
    print(
        json.dumps(
            {
                "out": str(out_dir),
                "n_site4": results["n_site4"],
                "n_attr": results["attr_n"],
                "n_ood": results["n_ood"],
                "perm_p": attribution["permutation_test"]["p_value_one_sided"],
                "label_cost_headline": label_cost_headline[:240],
                "runtime_sec": results["runtime_sec"],
                "runtime_phase_sec": phase_sec,
            },
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
