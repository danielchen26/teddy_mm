"""Declared CITE lineage panels and criterion packs (O0 / O1 / O2).

Protein names and Pearson figures are loaded from phase-1
`test_per_protein.json` at export time; this module only declares the
decision family structure. Do not invent Pearson values here.

O0->O1: stricter abstain (margin/threshold) — expected_action usually unchanged.
O2: rewrites expected_action via key-marker priority + lineage weights so
labels disagree with O0 on a measurable fraction (therapeutic B/T priority).
"""
from __future__ import annotations

from typing import Any

LINEAGE_PANELS: dict[str, list[str]] = {
    "b_lineage": ["CD19", "CD72", "CD22"],
    "t_lineage": ["CD3", "CD2", "CD5"],
    "myeloid": ["CD16", "CD11c", "CD36"],
}

KEY_MARKERS: dict[str, str] = {
    "b_lineage": "CD19",
    "t_lineage": "CD3",
    "myeloid": "CD16",
}

MODALITY_FOR_LINEAGE = {
    "b_lineage": "b_marker",
    "t_lineage": "t_marker",
    "myeloid": "myeloid_marker",
}

ACTIONS = [
    {"id": "b_lineage", "label": "B-lineage coherent"},
    {"id": "t_lineage", "label": "T-lineage coherent"},
    {"id": "myeloid", "label": "Myeloid coherent"},
]

CRITERIA: dict[str, dict[str, Any]] = {
    "O0": {
        "name": "workability_soft",
        "description": (
            "P_f workability: soft panel coherence on 3+3+3 high-Pearson lineage markers; "
            "equal protein weights; readout_threshold=0.12."
        ),
        "score_mode": "equal_panel_mean",
        "readout_threshold": 0.12,
        "protein_weight": 1.0,
        "key_marker_boost": 1.0,
        "secondary_weight": 1.0,
        "lineage_weights": {"b_lineage": 1.0, "t_lineage": 1.0, "myeloid": 1.0},
        "value_clip": [0.0, 1.0],
        "norm": "train_p95",
        "expected_from": "true_adt_panel_argmax",
        "expected_margin": 0.05,
    },
    "O1": {
        "name": "exactness_strict",
        "description": (
            "Q_f-oriented exactness: same TEDDY predictions, re-declared observer — "
            "boost key therapeutic markers (CD19/CD3/CD16)x2.0, stricter "
            "readout_threshold=0.28, larger expected margin. Mostly abstain-stricter "
            "vs O0 (expected_action usually unchanged when both defined)."
        ),
        "score_mode": "equal_panel_mean",
        "readout_threshold": 0.28,
        "protein_weight": 1.0,
        "key_marker_boost": 2.0,
        "secondary_weight": 1.0,
        "lineage_weights": {"b_lineage": 1.0, "t_lineage": 1.0, "myeloid": 1.0},
        "value_clip": [0.0, 1.0],
        "norm": "train_p95",
        "expected_from": "true_adt_panel_argmax",
        "expected_margin": 0.12,
    },
    "O2": {
        "name": "therapeutic_key_marker_priority",
        "description": (
            "Rewrites expected_action vs O0: key-marker-only panel scores "
            "(CD19/CD3/CD16) with lineage weights favoring therapeutic B/T "
            "(b=1.5, t=1.3, myeloid=0.5). Secondary panel members weight 0. "
            "Measurable label disagreement with O0 equal-panel-mean."
        ),
        "score_mode": "key_marker_priority",
        "readout_threshold": 0.20,
        "protein_weight": 1.0,
        "key_marker_boost": 3.0,
        "secondary_weight": 0.0,
        "lineage_weights": {"b_lineage": 1.5, "t_lineage": 1.3, "myeloid": 0.5},
        "value_clip": [0.0, 1.0],
        "norm": "train_p95",
        "expected_from": "true_adt_key_marker_weighted",
        "expected_margin": 0.06,
    },
}


def all_panel_proteins() -> list[str]:
    out: list[str] = []
    for prots in LINEAGE_PANELS.values():
        for p in prots:
            if p not in out:
                out.append(p)
    return out


def protein_to_lineage() -> dict[str, str]:
    return {p: lin for lin, prots in LINEAGE_PANELS.items() for p in prots}


def _norm_val(raw: float, scale: float) -> float:
    return min(1.0, max(0.0, float(raw) / max(float(scale), 1e-6)))


def lineage_scores_from_panel(
    panel: dict[str, float],
    p95: dict[str, float],
    criterion: dict[str, Any] | None = None,
    *,
    criterion_id: str | None = None,
) -> dict[str, float]:
    """Lineage scores under a declared criterion (shared by expected / pred paths)."""
    crit = criterion if criterion is not None else CRITERIA[criterion_id or "O0"]
    mode = crit.get("score_mode", "equal_panel_mean")
    lw = crit.get("lineage_weights") or {a: 1.0 for a in LINEAGE_PANELS}
    sec_w = float(crit.get("secondary_weight", 1.0))
    key_boost = float(crit.get("key_marker_boost", 1.0))
    scores: dict[str, float] = {}
    for lin, prots in LINEAGE_PANELS.items():
        key = KEY_MARKERS[lin]
        if mode == "key_marker_priority":
            num = 0.0
            den = 0.0
            for p in prots:
                v = _norm_val(panel[p], p95.get(p, 1.0))
                w = 1.0 if p == key else sec_w
                num += v * w
                den += w
            base = (num / den) if den > 0 else 0.0
            scores[lin] = float(lw.get(lin, 1.0) * base)
        else:
            vals = []
            for p in prots:
                v = _norm_val(panel[p], p95.get(p, 1.0))
                if p == key and key_boost != 1.0:
                    v = min(1.0, v * key_boost)
                vals.append(v)
            base = float(sum(vals) / len(vals)) if vals else 0.0
            scores[lin] = float(lw.get(lin, 1.0) * base)
    return scores


def expected_from_true(
    true_panel: dict[str, float],
    p95: dict[str, float],
    criterion: dict[str, Any] | None = None,
    *,
    criterion_id: str | None = None,
    margin: float | None = None,
) -> str | None:
    """Post-hoc verifier label from holdout true ADT — never used as source event."""
    crit = criterion if criterion is not None else CRITERIA[criterion_id or "O0"]
    m = float(crit["expected_margin"] if margin is None else margin)
    if crit.get("score_mode") == "key_marker_priority":
        scores = lineage_scores_from_panel(true_panel, p95, crit)
    else:
        scores = {}
        for lin, prots in LINEAGE_PANELS.items():
            vals = [_norm_val(true_panel[p], p95.get(p, 1.0)) for p in prots]
            scores[lin] = float(sum(vals) / len(vals)) if vals else 0.0
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    if ranked[0][1] - ranked[1][1] < m:
        return None
    return ranked[0][0]


def event_value_under_criterion(
    raw_norm_value: float, is_key_marker: bool, action: str, crit: dict
) -> float:
    """Map a p95-normalized TEDDY value into the declared field event value."""
    w = float(crit.get("protein_weight", 1.0))
    if crit.get("score_mode") == "key_marker_priority":
        w *= float(crit["key_marker_boost"]) if is_key_marker else float(
            crit.get("secondary_weight", 0.0)
        )
        w *= float((crit.get("lineage_weights") or {}).get(action, 1.0))
    else:
        if is_key_marker:
            w *= float(crit.get("key_marker_boost", 1.0))
    value = float(raw_norm_value) * w
    lo, hi = crit.get("value_clip", [0.0, 1.0])
    return max(float(lo), min(float(hi), value))
