#!/usr/bin/env python3
"""Missing-modality × ANM demo on TEDDY CITE.

Compares three modality_mask values on the SAME base cells:
  rna_only | adt_only | joint

Arms: (1) TEDDY alone  (2) TEDDY+ANM  (3) train-retune under missing modality.
Criteria: O0 (soft) and O2 (key-marker rewrite) — declaration edit, no retrain.
Finite-field + LOO attribution/flip per mask.

Does NOT claim win over phase-1 Pearson 0.61. No best.pt backbone retrain.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from sklearn.linear_model import LogisticRegression
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
MASKS = ("rna_only", "adt_only", "joint")
CRIT_IDS = ("O0", "O2")  # O0 soft workability; O2 key-marker rewrite

# Missing-modality declaration (no TEDDY retrain): RNA is key evidence for this family.
for _cid in CRIT_IDS:
    CRITERIA[_cid] = {
        **CRITERIA[_cid],
        "modality_reliability": {"rna_only": 1.0, "joint": 1.0, "adt_only": 0.35},
        "missing_modality_note": "adt_only scales event values x0.35 (RNA absent)",
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


def load_by_mask(cells_path, events_path):
    cells_by_mask = {m: {} for m in MASKS}
    for c in _load_jsonl(cells_path):
        m = c["modality_mask"]
        if m in cells_by_mask:
            cells_by_mask[m][c["cell_id"]] = c
    events_by_mask = {m: defaultdict(list) for m in MASKS}
    for ev in _load_jsonl(events_path):
        m = ev["modality_mask"]
        if m in events_by_mask:
            events_by_mask[m][ev["cell_id"]].append(ev)
    return cells_by_mask, events_by_mask


def p95_from_cells(cells):
    # recover from any event file values stored on cells via first cell's... 
    # actually export stores norm in events; use first events batch
    return None


def modality_reliability(mask: str, crit: dict) -> float:
    """Declared workability weight when key evidence (RNA) is absent.

    Decision family is RNA→lineage coherence. Missing RNA (adt_only) drops
    declared reliability without retraining TEDDY — O0/O2-style criterion edit.
    """
    rel = crit.get("modality_reliability") or {
        "rna_only": 1.0,
        "joint": 1.0,
        "adt_only": 0.35,  # RNA absent: key evidence missing
    }
    return float(rel.get(mask, 1.0))


def make_instance(cid, cell, events, crit, p95):
    proposed = []
    rel = modality_reliability(cell.get("modality_mask", "rna_only"), crit)
    for ev in sorted(events, key=lambda e: e["time"]):
        value = event_value_under_criterion(
            float(ev["value"]), bool(ev.get("is_key_marker")), ev["action"], crit
        )
        value = max(0.0, min(1.0, float(value) * rel))
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
        "question": "lineage coherence under missing modality",
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
        "n_admitted": len(validation["field_events"]),
    }


def schema_for(crit_id, base_schema):
    schema = copy.deepcopy(base_schema)
    # v2: rule threshold mapped onto the ANM field scale; legacy: declared value.
    thr = float(anm_readout_threshold(CRITERIA[crit_id], base_schema))
    schema["field_representation"]["readout_threshold"] = thr
    schema["criterion_overlay"] = {"readout_threshold": thr, "criterion": crit_id}
    return schema


def run_teddy_arm(eval_ids, cells, p95):
    arm = {"arm": "TEDDY_alone", "criteria": {}}
    for cid_name in CRIT_IDS:
        crit = CRITERIA[cid_name]
        preds, labels = [], []
        for cid in eval_ids:
            sc = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, crit)
            preds.append(decide_argmax(sc, crit["readout_threshold"]))
            labels.append(expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit))
        m = score_predictions(preds, labels)
        arm["criteria"][cid_name] = {**m, "mean_P_f": None, "mean_Q_f_defined": m["Q_analogue"]}
    # silent over-answer O0 vs O2
    preds_o0, preds_o2 = [], []
    for cid in eval_ids:
        sc0 = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, CRITERIA["O0"])
        sc2 = lineage_scores_from_panel(cells[cid]["adt_pred_panel"], p95, CRITERIA["O2"])
        preds_o0.append(decide_argmax(sc0, CRITERIA["O0"]["readout_threshold"]))
        preds_o2.append(decide_argmax(sc2, CRITERIA["O2"]["readout_threshold"]))
    n_silent = sum(1 for a, b in zip(preds_o0, preds_o2) if a is not None and b is None)
    arm["silent_over_answer_O0_vs_O2"] = {
        "count": n_silent,
        "rate": n_silent / max(len(eval_ids), 1),
    }
    arm["criterion_edit_cost"] = {
        "needs_endpoint_labels": True,
        "needs_retrain_or_retune": True,
    }
    return arm


def run_anm_arm(eval_ids, cells, by_cell, p95, base_schema):
    criteria_out = {}
    for cid_name in CRIT_IDS:
        crit = CRITERIA[cid_name]
        schema = schema_for(cid_name, base_schema)
        preds, labels = [], []
        p_sum = q_sum = q_n = 0.0
        n_empty = 0
        for cid in eval_ids:
            inst = make_instance(cid, cells[cid], by_cell[cid], crit, p95)
            if not inst["proposed_source_events"]:
                n_empty += 1
                preds.append(None)
                labels.append(inst.get("expected_action"))
                continue
            out = anm_run_one(schema, inst)
            preds.append(out["recommended_action"])
            labels.append(out.get("expected_action"))
            p_sum += float(out.get("P_f") or 0.0)
            if out.get("Q_f") is not None:
                q_sum += float(out["Q_f"])
                q_n += 1
        m = score_predictions(preds, labels)
        criteria_out[cid_name] = {
            **m,
            "mean_P_f": p_sum / max(len(eval_ids), 1),
            "mean_Q_f_defined": (q_sum / q_n) if q_n else None,
            "n_Q_f_defined": int(q_n),
            "n_empty_source": n_empty,
            "label_fit": False,
        }
    return {
        "arm": "TEDDY_plus_ANM",
        "criteria": criteria_out,
        "criterion_edit_cost": {
            "needs_endpoint_labels": False,
            "mechanism": "YAML/criterion declaration only (O0→O2)",
        },
    }


def _loo_top(schema, inst, flip_grid=15):
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
    }


def run_attribution(attr_ids, cells, by_cell, p95, base_schema, flip_grid=15):
    schema0 = schema_for("O0", base_schema)
    crit0 = CRITERIA["O0"]
    cell_ids, tops, flips, fds, deltas = [], [], [], [], []
    t0 = time.time()
    for k, cid in enumerate(attr_ids):
        inst = make_instance(cid, cells[cid], by_cell[cid], crit0, p95)
        if not inst["proposed_source_events"]:
            continue
        row = _loo_top(schema0, inst, flip_grid=flip_grid)
        if row is None:
            continue
        cell_ids.append(cid)
        tops.append(row["protein"])
        flips.append(1.0 if row["flipped"] else 0.0)
        fds.append(row["flip_distance"])
        deltas.append(row["delta_chosen_score"])
        if (k + 1) % 500 == 0:
            print(f"    attr {k+1}/{len(attr_ids)} kept={len(tops)}", flush=True)
    fd = [d for d in fds if d is not None]
    return {
        "n_attr_cells": len(tops),
        "n_attr_requested": len(attr_ids),
        "top1_protein_distribution": {str(k): int(v) for k, v in Counter(tops).items()},
        "top1_flip_rate": float(np.mean(flips)) if flips else None,
        "flip_distance_quantiles": {
            "q25": float(np.quantile(fd, 0.25)) if fd else None,
            "q50": float(np.quantile(fd, 0.50)) if fd else None,
            "q75": float(np.quantile(fd, 0.75)) if fd else None,
            "n_with_flip": len(fd),
            "n_no_flip": sum(1 for d in fds if d is None),
        },
        "elapsed_sec": time.time() - t0,
        "compact": {
            "cell_id": np.array(cell_ids, dtype=object),
            "top_protein": np.array(tops, dtype=object),
            "flipped": np.array(flips, dtype=np.float32),
            "flip_distance": np.array(
                [(-1.0 if d is None else d) for d in fds], dtype=np.float32
            ),
        },
    }


def feature_vector(pred_panel, p95):
    return np.asarray(
        [min(1.0, max(0.0, pred_panel[p] / max(p95.get(p, 1.0), 1e-6))) for p in PROTEINS],
        dtype=np.float64,
    )


def run_train_retune(eval_ids, train_ids, cells, p95, n_labels_grid=(50, 200, 500, 1000)):
    """Train logistic on O2 labels from train split of THIS mask; eval on eval_ids.
    Shows label cost to approach ANM-declared O2 abstain under missing modality.
    """
    crit = CRITERIA["O2"]
    # Build matrices
    def pack(ids):
        Xs, ys, keep = [], [], []
        for cid in ids:
            y = expected_from_true(cells[cid]["adt_true_panel_holdout"], p95, crit)
            if y is None:
                continue
            Xs.append(feature_vector(cells[cid]["adt_pred_panel"], p95))
            ys.append(y)
            keep.append(cid)
        if not Xs:
            return None, None, []
        return np.stack(Xs), np.array(ys), keep

    Xtr, ytr, tr_keep = pack(train_ids)
    Xev, yev, ev_keep = pack(eval_ids)
    if Xtr is None or Xev is None or len(tr_keep) < 30:
        return {"arm": "train_retune", "error": "insufficient labeled rows", "n_train_labeled": 0 if Xtr is None else len(tr_keep)}

    # Map labels to ints
    classes = sorted(set(ytr.tolist()) | set(yev.tolist()))
    cmap = {c: i for i, c in enumerate(classes)}
    ytr_i = np.array([cmap[y] for y in ytr])
    yev_i = np.array([cmap[y] for y in yev])

    scaler = StandardScaler()
    Xtr_s = scaler.fit_transform(Xtr)
    Xev_s = scaler.transform(Xev)

    rng = np.random.default_rng(17)
    curves = []
    for n_lab in n_labels_grid:
        n_use = min(int(n_lab), len(tr_keep))
        idx = rng.choice(len(tr_keep), size=n_use, replace=False)
        clf = LogisticRegression(max_iter=500)
        clf.fit(Xtr_s[idx], ytr_i[idx])
        proba = clf.predict_proba(Xev_s)
        pred_i = proba.argmax(axis=1)
        # Calibrate abstain: keep top margin quantile matching ANM-ish ~0.2 if possible
        margin = proba.max(axis=1) - np.partition(proba, -2, axis=1)[:, -2]
        # Grid thr to get abstain ~0.15–0.25
        best = None
        for thr in np.linspace(0.02, 0.5, 25):
            preds = []
            for i in range(len(ev_keep)):
                if margin[i] < thr:
                    preds.append(None)
                else:
                    preds.append(classes[pred_i[i]])
            m = score_predictions(preds, list(yev))
            score = abs((m["abstain_rate"] or 0) - 0.20)
            if best is None or score < best["gap"]:
                best = {"gap": score, "metrics": m, "thr": float(thr)}
        curves.append({
            "n_O2_labels": n_use,
            "abstain_rate": best["metrics"]["abstain_rate"],
            "Q_analogue": best["metrics"]["Q_analogue"],
            "accuracy_strict_labeled": best["metrics"]["accuracy_strict_labeled"],
            "calib_margin_thr": best["thr"],
        })
    return {
        "arm": "train_retune",
        "n_train_labeled_pool": len(tr_keep),
        "n_eval_labeled": len(ev_keep),
        "curves": curves,
        "criterion_edit_cost": {
            "needs_endpoint_labels": True,
            "note": "each observer/mask change consumes new O2 labels; no editable field",
        },
    }


def write_report(out: dict, path: Path):
    man = out.get("export_manifest", {})
    lines = []
    lines.append("# MISSING_MODALITY_ANM_DEMO — TEDDY × ANM")
    lines.append("")
    lines.append("## Framing")
    lines.append("")
    lines.append(
        f"- Source choice: **{man.get('source', out.get('source'))}** "
        "(prefer phase-2 ckpt; hybrid/simulate if phase-2 RNA pearson weak)."
    )
    lines.append(
        f"- Base site4/test cells: **n={out['n_base_cells']}** "
        f"(×3 masks → {out['n_cells_total']} cell-mask rows)."
    )
    lines.append(
        f"- Phase-1 unidirectional Pearson (quoted): **~0.61** "
        f"(measured full-ADT on export slice: {man.get('phase1_full_adt_pearson')})."
    )
    lines.append(
        f"- Panel Pearson by mask: `{json.dumps(man.get('panel_pearson_by_mask', {}))}`."
    )
    lines.append(
        "- **Do NOT claim win over 0.61.** If phase-2/hybrid predictions are weak, "
        "still show ANM abstain / attribution / declaration advantages."
    )
    lines.append(
        "- ANM semantics: missing modality = **source ablation / channel-restricted "
        "typed evidence**; declared P_f drops when key evidence absent; Q_f vs holdout GT only."
    )
    lines.append("- O0→O2 criterion edit without TEDDY retrain.")
    lines.append("- No Perturb / GFlowNet / 160M / ATAC. Not clinical.")
    lines.append("")
    lines.append("## Reading note (TEDDY's role)")
    lines.append('')
    lines.append("- **`adt_only` = RNA missing, so TEDDY does not run.** Its evidence comes from the phase-2 BidirectionalCite model run with its RNA input off: it reads the cell's 134 measured proteins and reconstructs the 9 panel proteins (panel Pearson 0.446). The same measured proteins also set the answer key.")
    lines.append('- **Arm names in the tables below.** "TEDDY alone" is the fixed rule with no ANM layer, reading each mask\'s predicted panel unscaled (O0: average each lineage\'s 3 markers, call the highest if it reaches 0.12; O2: key marker only × lineage weight B 1.5 / T 1.3 / myeloid 0.5, bar 0.20). "TEDDY+ANM" is ANM on the same evidence (×0.35 for `adt_only`). Under `adt_only` neither arm involves TEDDY, so those rows are labelled "Fixed rule (stand-in)" and "ANM (stand-in)".')
    lines.append('- **The `adt_only` abstentions follow the declared trust factor 0.35** (`modality_reliability`), which was chosen, not estimated. Under O2 (key-marker rule) no score can reach the 0.20 bar with that factor: the best case is 0.126 (`scripts/missing_modality_bounds.py`), so all `adt_only` cells are declined by construction.')
    lines.append('')
    lines.append("## Key tables by mask")
    lines.append("")
    for m in MASKS:
        block = out["by_mask"][m]
        lines.append(f"### Mask `{m}` (n={block['n']})")
        lines.append("")
        lines.append("| arm | crit | Q | abstain | mean_P_f |")
        lines.append("|---|---|---:|---:|---:|")
        # with RNA missing TEDDY does not run: both arms read the phase-2 stand-in
        stand_in = m == "adt_only"
        for arm_key, arm_name in [
            ("teddy_alone", "Fixed rule (stand-in)" if stand_in else "TEDDY alone"),
            ("anm", "ANM (stand-in)" if stand_in else "TEDDY+ANM"),
        ]:
            arm = block[arm_key]
            for c in CRIT_IDS:
                crit = arm["criteria"][c]
                lines.append(
                    f"| {arm_name} | {c} | {crit.get('Q_analogue')} | "
                    f"{crit.get('abstain_count')} | {crit.get('mean_P_f')} |"
                )
        lines.append("")
        attr = block.get("attribution", {})
        lines.append(
            f"- LOO attribution (O0): n_attr={attr.get('n_attr_cells')} / "
            f"requested={attr.get('n_attr_requested')}; "
            f"top1_flip_rate={attr.get('top1_flip_rate')}; "
            f"top1_dist=`{attr.get('top1_protein_distribution')}`."
        )
        silent = block["teddy_alone"].get("silent_over_answer_O0_vs_O2", {})
        lines.append(
            f"- {'Fixed rule (stand-in)' if stand_in else 'TEDDY-alone'} silent over-answer O0 vs O2 decl: "
            f"count={silent.get('count')} rate={silent.get('rate')}."
        )
        tr = block.get("train_retune", {})
        if tr.get("curves"):
            lines.append("- Train-retune O2 label curve (logistic, calibrated abstain≈0.2):")
            lines.append("")
            lines.append("| n_O2_labels | abstain | Q | strict |")
            lines.append("|---:|---:|---:|---:|")
            for row in tr["curves"]:
                lines.append(
                    f"| {row['n_O2_labels']} | {row['abstain_rate']:.4f} | "
                    f"{row['Q_analogue']} | {row['accuracy_strict_labeled']} |"
                )
            lines.append("")
        else:
            lines.append(f"- Train-retune: `{tr.get('error', tr)}`")
            lines.append("")
    lines.append("## Contrast summary")
    lines.append("")
    lines.append("| mask | Fixed rule O0 Q | Fixed rule O0 abstain | ANM O0 Q | ANM O0 abstain | ANM O0 mean_P_f | ANM O2 abstain |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for m in MASKS:
        b = out["by_mask"][m]
        t0 = b["teddy_alone"]["criteria"]["O0"]
        a0 = b["anm"]["criteria"]["O0"]
        a2 = b["anm"]["criteria"]["O2"]
        lines.append(
            f"| {m} | {t0.get('Q_analogue')} | {t0.get('abstain_count')} | "
            f"{a0.get('Q_analogue')} | {a0.get('abstain_count')} | "
            f"{a0.get('mean_P_f')} | {a2.get('abstain_count')} |"
        )
    lines.append("")
    lines.append("## Proof sentence")
    lines.append("")
    lines.append(f"> {out['proof_sentence']}")
    lines.append("")
    lines.append("## Blockers / honesty")
    lines.append("")
    for b in out.get("blockers", []):
        lines.append(f"- {b}")
    lines.append("")
    lines.append("## Paths")
    lines.append("")
    lines.append(f"- Report: `{path}`")
    lines.append(f"- Results JSON: `{path.with_name('missing_modality_results.json')}`")
    lines.append(f"- Export: `{man.get('events_path')}` / `{man.get('cells_path')}`")
    lines.append(f"- Attr NPZ dir: `{path.parent / 'attr'}`")
    lines.append("")
    lines.append(f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}")
    path.write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=ROOT / "outputs/anm_cite_bridge/missing_modality",
    )
    ap.add_argument("--schema", type=Path, default=ROOT / "bridge_anm/schemas/cite_lineage_finite_field_v0.json")
    ap.add_argument("--attr-n", type=int, default=1500, help="LOO cells per mask; 0=all")
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--train-frac", type=float, default=0.3, help="fraction of cells for train-retune fit")
    args = ap.parse_args()

    out_dir = args.data_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "attr").mkdir(exist_ok=True)

    cells_path = out_dir / "missing_modality_cells.jsonl"
    events_path = out_dir / "missing_modality_events.jsonl"
    manifest_path = out_dir / "export_manifest.json"
    if not cells_path.exists():
        raise SystemExit(f"missing {cells_path}; run export_missing_modality_events.py first")

    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    p95 = manifest.get("p95_train")
    if not p95:
        # recover from events
        p95 = {}
        for ev in _load_jsonl(events_path):
            p95[ev["protein"]] = float(ev["norm_p95_train"])

    base_schema = json.loads(args.schema.read_text())
    cells_by_mask, events_by_mask = load_by_mask(cells_path, events_path)

    # Align base cell order across masks
    base_ids = sorted({c["base_cell_id"] for c in cells_by_mask["rna_only"].values()})
    n_base = len(base_ids)
    rng = np.random.default_rng(args.seed)
    perm = rng.permutation(n_base)
    n_train = max(50, int(n_base * args.train_frac))
    train_base = [base_ids[i] for i in perm[:n_train]]
    eval_base = [base_ids[i] for i in perm[n_train:]]
    if len(eval_base) < 100:
        eval_base = base_ids  # tiny set fallback
        train_base = base_ids[: max(50, n_base // 3)]

    results = {
        "n_base_cells": n_base,
        "n_cells_total": sum(len(cells_by_mask[m]) for m in MASKS),
        "n_eval_base": len(eval_base),
        "n_train_base": len(train_base),
        "source": manifest.get("source"),
        "export_manifest": manifest,
        "by_mask": {},
        "blockers": [],
    }

    # Honesty blockers
    pp = manifest.get("panel_pearson_by_mask") or {}
    for m, r in pp.items():
        if r is not None and r == r and r < 0.40:
            results["blockers"].append(
                f"panel pearson under mask={m} is {r:.3f} (weak vs phase-1 ~0.61); "
                "ANM abstain/attribution still meaningful; do not inflate Pearson."
            )
    if manifest.get("source") in ("simulate_phase1_plus_adt", "hybrid_phase1_rna_phase2_adt"):
        results["blockers"].append(
            f"export source={manifest.get('source')}: phase-2 RNA path weak or absent; "
            "documented hybrid/simulate. Not a claim that phase-2 beats 0.61."
        )
    results["blockers"].append(
        "adt_only / joint may use ADT-channel evidence (phase-2 recon or observed ADT in simulate); "
        "Q_f can look strong by construction on ADT-present masks — interpret as modality-availability, not SOTA Pearson."
    )

    for m in MASKS:
        print(f"=== mask {m} ===", flush=True)
        cells = cells_by_mask[m]
        by_cell = events_by_mask[m]
        # map base → mask-specific cell_id
        def cid_for(base):
            return f"{base}__{m}"

        eval_ids = [cid_for(b) for b in eval_base if cid_for(b) in cells]
        train_ids = [cid_for(b) for b in train_base if cid_for(b) in cells]
        print(f"  eval={len(eval_ids)} train={len(train_ids)}", flush=True)

        teddy = run_teddy_arm(eval_ids, cells, p95)
        print(f"  TEDDY O0 Q={teddy['criteria']['O0'].get('Q_analogue')} abstain={teddy['criteria']['O0']['abstain_count']}", flush=True)
        anm = run_anm_arm(eval_ids, cells, by_cell, p95, base_schema)
        print(
            f"  ANM O0 Q={anm['criteria']['O0'].get('Q_analogue')} "
            f"abstain={anm['criteria']['O0']['abstain_count']} "
            f"P_f={anm['criteria']['O0'].get('mean_P_f')}",
            flush=True,
        )

        attr_ids = eval_ids if args.attr_n <= 0 else eval_ids[: min(args.attr_n, len(eval_ids))]
        print(f"  attribution n={len(attr_ids)}", flush=True)
        attr = run_attribution(attr_ids, cells, by_cell, p95, base_schema)
        compact = attr.pop("compact")
        np.savez_compressed(out_dir / "attr" / f"attr_{m}.npz", **compact)

        train = run_train_retune(eval_ids, train_ids, cells, p95)
        print(f"  train_retune curves={len(train.get('curves') or [])}", flush=True)

        results["by_mask"][m] = {
            "n": len(eval_ids),
            "teddy_alone": teddy,
            "anm": anm,
            "attribution": attr,
            "train_retune": train,
        }

    # Proof sentence
    bits = []
    for m in MASKS:
        a0 = results["by_mask"][m]["anm"]["criteria"]["O0"]
        bits.append(
            f"{m}: ANM abstain={a0['abstain_count']}/P_f={a0.get('mean_P_f')}"
        )
    n_eval = len(eval_base)
    results["proof_sentence"] = (
        f"On {n_base:,} site4/test base cells ({n_eval:,} scored per mask; {n_base - n_eval:,} held for re-tuning), "
        f"three masks act as source ablation: rna_only = TEDDY phase-1 predictions; adt_only = RNA missing, so TEDDY "
        f"does not run and the phase-2 protein-only stand-in supplies the evidence; joint = the average of the two. "
        f"ANM's P_f/abstain follow the declared per-source trust (rna_only/joint 1.0, adt_only 0.35, chosen not estimated) "
        f"({'; '.join(bits)}), O0→O2 observer edits are YAML-only "
        f"(no TEDDY retrain), and LOO/flip attribution remains auditable per mask; "
        f"the fixed rule (TEDDY alone's rule, no trust setting) answers every adt_only cell at O0; "
        f"train-retune must consume O2 labels per mask/observer and still lacks editable-field "
        f"attribution. Pearson is secondary — do not claim win over phase-1 ~0.61. Not clinical."
    )

    json_path = out_dir / "missing_modality_results.json"
    # strip numpy from nested if any
    json_path.write_text(json.dumps(results, indent=2, default=str))
    report_path = out_dir / "MISSING_MODALITY_ANM_DEMO.md"
    write_report(results, report_path)
    print("wrote", report_path)
    print("PROOF:", results["proof_sentence"][:200], "...")


if __name__ == "__main__":
    main()
