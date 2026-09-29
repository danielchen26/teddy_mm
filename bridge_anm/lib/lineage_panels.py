"""Declared CITE lineage panels and criterion packs (O0 / O1 / O2).

Protein names and Pearson figures are loaded from phase-1
`test_per_protein.json` at export time; this module only declares the
decision family structure. Do not invent Pearson values here.

O0->O1: stricter abstain (margin/threshold) — expected_action usually unchanged.
O2: rewrites expected_action via key-marker priority + lineage weights so
labels disagree with O0 on a measurable fraction (therapeutic B/T priority).

Scoring (``SCORING = "v2"``, the default since fix/anm-bridge-corrections)
-------------------------------------------------------------------------
One scoring rule is shared by the answer key (``expected_from_true``), the
TEDDY-alone rule (``lineage_scores_from_panel`` + ``rule_threshold``) and the
ANM source-event values (``event_value_under_criterion``):

1. Evidence per protein ``v_p`` = value / train-p95, clipped to [0, 1]. This is
   the only clip; it is applied to the evidence, never to weighted evidence.
2. The question declares a raw weight per protein::

       w_p = protein_weight * (key_marker_boost if p is the lineage key marker
                               else secondary_weight) * lineage_weights[lineage]

   (O0: all 1; O1: key 2, secondary 1; O2: key 3 x lineage weight, secondary 0.)
3. Weights are divided by the largest weight over all panel proteins,
   ``w~_p = w_p / max_q w_q`` in [0, 1]. There is no clip after weighting, so
   the declared priority between lineages is preserved. (Old behaviour: under
   O2 the key events were multiplied by 3 x lineage weight and clipped at 1, so
   a cell with CD19 = CD3 = CD16 = 1 had three equal events and ANM called it
   myeloid because later panel events decay less.)
4. ANM event value = ``w~_p * v_p`` (still in [0, 1]).
5. Lineage score ``S_l = sum_{p in l} w~_p v_p / c`` with
   ``c = max_l sum_{p in l} w~_p`` (``lineage_normaliser``), so a saturated
   strongest lineage scores 1. ``S_l`` is proportional to the ANM field score
   of lineage ``l`` when all events enter simultaneously (the finite-graph
   field is linear and lineages are uncoupled), so the key, the rule and ANM
   rank lineages identically. For O0 this is exactly the old equal panel mean.
6. Declared margins and rule thresholds were written on the old score scale,
   whose maximum is ``max(lineage_weights)`` (1.0 for O0/O1, 1.5 for O2); they
   are divided by that scale (``declared_score_scale``), e.g. the O2 margin
   0.06 becomes 0.04 and the O2 rule threshold 0.20 becomes 0.1333. The O2 key
   labels are therefore unchanged; the O1 key now includes the O1 key-marker
   boost (the old key ignored it while the old rule applied it with a clip).

   ``CRITERIA[cid]["readout_threshold"]`` / ``["expected_margin"]`` hold these
   effective values (the declared ones stay in ``declared_readout_threshold`` /
   ``declared_expected_margin``), so callers that read CRITERIA need no change.
7. With simultaneous events (exporter default) the ANM action score is exactly
   ``field_gain(schema) * c * S_l``; ``anm_readout_threshold(crit, schema)``
   gives the field threshold that abstains exactly where the rule does.
   (With panel-order times later events decay less and ANM favours myeloid.)

``SCORING = "legacy"`` (env ``ANM_BRIDGE_SCORING=legacy``, ``set_scoring("legacy")``
or per call ``scoring="legacy"``) restores the previous functions verbatim so
old numbers can be reproduced. For whole old runs use the global flag: it also
puts the declared thresholds/margins back into CRITERIA.
"""
from __future__ import annotations

import os
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
    """Evidence = value / train-p95, clipped to [0, 1] (the only clip in v2)."""
    return min(1.0, max(0.0, float(raw) / max(float(scale), 1e-6)))


# --------------------------------------------------------------------------
# Legacy scoring (verbatim pre-fix behaviour; reachable via SCORING="legacy")
# --------------------------------------------------------------------------
def _legacy_lineage_scores_from_panel(
    panel: dict[str, float],
    p95: dict[str, float],
    criterion: dict[str, Any] | None = None,
    *,
    criterion_id: str | None = None,
) -> dict[str, float]:
    """Legacy lineage scores (TEDDY-alone rule path before the fix)."""
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


def _legacy_expected_from_true(
    true_panel: dict[str, float],
    p95: dict[str, float],
    criterion: dict[str, Any] | None = None,
    *,
    criterion_id: str | None = None,
    margin: float | None = None,
) -> str | None:
    """Legacy answer key (O1 ignored the key-marker boost)."""
    crit = criterion if criterion is not None else _declared(criterion_id or "O0")
    m = float(_declared_value(crit, "expected_margin") if margin is None else margin)
    if crit.get("score_mode") == "key_marker_priority":
        scores = _legacy_lineage_scores_from_panel(true_panel, p95, crit)
    else:
        scores = {}
        for lin, prots in LINEAGE_PANELS.items():
            vals = [_norm_val(true_panel[p], p95.get(p, 1.0)) for p in prots]
            scores[lin] = float(sum(vals) / len(vals)) if vals else 0.0
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    if ranked[0][1] - ranked[1][1] < m:
        return None
    return ranked[0][0]


def _legacy_event_value_under_criterion(
    raw_norm_value: float, is_key_marker: bool, action: str, crit: dict
) -> float:
    """Legacy event value: weight, then clip to value_clip (priority lost at 1)."""
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


# --------------------------------------------------------------------------
# Scoring switch
# --------------------------------------------------------------------------
SCORING_MODES = ("v2", "legacy")
#: Global scoring mode. Override with env ``ANM_BRIDGE_SCORING=legacy`` or
#: ``set_scoring("legacy")``; every public function also takes ``scoring=``.
SCORING: str = os.environ.get("ANM_BRIDGE_SCORING", "v2")

# Declared (question-author) thresholds/margins, on the old score scale.
_DECLARED: dict[str, dict[str, float]] = {
    cid: {"readout_threshold": float(c["readout_threshold"]), "expected_margin": float(c["expected_margin"])}
    for cid, c in CRITERIA.items()
}
_RESCALED_KEY = "_scale_rescaled_by"  # set on a criterion whose thr/margin were divided by the scale


def _mode(scoring: str | None) -> str:
    m = SCORING if scoring is None else scoring
    if m not in SCORING_MODES:
        raise ValueError(f"unknown scoring mode {m!r}; expected one of {SCORING_MODES}")
    return m


def _declared(criterion_id: str) -> dict[str, Any]:
    """A copy of CRITERIA[criterion_id] with its declared (unrescaled) thr/margin."""
    crit = dict(CRITERIA[criterion_id])
    crit.update(_DECLARED[criterion_id])
    crit.pop(_RESCALED_KEY, None)
    return crit


def _declared_value(crit: dict[str, Any], key: str) -> float:
    """Declared thr/margin of ``crit`` (undo the v2 rescale if it was applied)."""
    return float(crit[key]) * float(crit.get(_RESCALED_KEY, 1.0))


def declared_score_scale(crit: dict[str, Any]) -> float:
    """Maximum of the old score scale = max lineage weight (O0/O1: 1, O2: 1.5)."""
    lw = crit.get("lineage_weights") or {a: 1.0 for a in LINEAGE_PANELS}
    return float(max(float(lw.get(lin, 1.0)) for lin in LINEAGE_PANELS))


def _apply_scoring_to_criteria() -> None:
    """Write the effective thr/margin for the current SCORING into CRITERIA.

    v2: declared / declared_score_scale (so scores in [0, 1] give the same
    decisions the declared numbers meant on the old scale). legacy: declared.
    Callers keep reading ``CRITERIA[cid]["readout_threshold"]`` unchanged.
    """
    for cid, crit in CRITERIA.items():
        scale = declared_score_scale(crit) if SCORING == "v2" else 1.0
        crit["readout_threshold"] = _DECLARED[cid]["readout_threshold"] / scale
        crit["expected_margin"] = _DECLARED[cid]["expected_margin"] / scale
        crit["declared_readout_threshold"] = _DECLARED[cid]["readout_threshold"]
        crit["declared_expected_margin"] = _DECLARED[cid]["expected_margin"]
        crit[_RESCALED_KEY] = scale


def set_scoring(mode: str) -> None:
    """Switch the global scoring mode ("v2" default, "legacy" = old numbers)."""
    global SCORING
    SCORING = _mode(mode)
    _apply_scoring_to_criteria()


def get_scoring() -> str:
    return SCORING


def effective_margin(crit: dict[str, Any], *, scoring: str | None = None) -> float:
    """Answer-key margin on the score scale of ``scoring``."""
    declared = _declared_value(crit, "expected_margin")
    return declared if _mode(scoring) == "legacy" else declared / declared_score_scale(crit)


def rule_threshold(crit: dict[str, Any], *, scoring: str | None = None) -> float:
    """TEDDY-alone abstain threshold on the score scale of ``scoring``."""
    declared = _declared_value(crit, "readout_threshold")
    return declared if _mode(scoring) == "legacy" else declared / declared_score_scale(crit)


# --------------------------------------------------------------------------
# v2: one shared weighting
# --------------------------------------------------------------------------
def raw_protein_weight(is_key_marker: bool, lineage: str, crit: dict[str, Any]) -> float:
    """Declared weight of one panel protein under the question (before /max)."""
    w = float(crit.get("protein_weight", 1.0))
    w *= float(crit.get("key_marker_boost", 1.0)) if is_key_marker else float(
        crit.get("secondary_weight", 1.0)
    )
    w *= float((crit.get("lineage_weights") or {}).get(lineage, 1.0))
    return w


def max_raw_weight(crit: dict[str, Any]) -> float:
    return max(
        raw_protein_weight(p == KEY_MARKERS[lin], lin, crit)
        for lin, prots in LINEAGE_PANELS.items()
        for p in prots
    )


def protein_weight(is_key_marker: bool, lineage: str, crit: dict[str, Any]) -> float:
    """Normalised weight w~ = w / max_q w_q in [0, 1] (no clip after weighting)."""
    top = max_raw_weight(crit)
    return raw_protein_weight(is_key_marker, lineage, crit) / top if top > 0 else 0.0


def protein_weights(crit: dict[str, Any]) -> dict[str, float]:
    return {
        p: protein_weight(p == KEY_MARKERS[lin], lin, crit)
        for lin, prots in LINEAGE_PANELS.items()
        for p in prots
    }


def lineage_normaliser(crit: dict[str, Any]) -> float:
    """c = max_l sum_{p in l} w~_p, so a saturated strongest lineage scores 1."""
    w = protein_weights(crit)
    c = max(sum(w[p] for p in prots) for prots in LINEAGE_PANELS.values())
    return c if c > 0 else 1.0


def _v2_scores(panel: dict[str, float], p95: dict[str, float], crit: dict[str, Any]) -> dict[str, float]:
    w = protein_weights(crit)
    c = lineage_normaliser(crit)
    return {
        lin: float(sum(w[p] * _norm_val(panel[p], p95.get(p, 1.0)) for p in prots) / c)
        for lin, prots in LINEAGE_PANELS.items()
    }


# --------------------------------------------------------------------------
# Public API (dispatch on scoring mode)
# --------------------------------------------------------------------------
def lineage_scores_from_panel(
    panel: dict[str, float],
    p95: dict[str, float],
    criterion: dict[str, Any] | None = None,
    *,
    criterion_id: str | None = None,
    scoring: str | None = None,
) -> dict[str, float]:
    """Lineage scores under a declared criterion (shared by key / rule / ANM).

    v2: S_l = sum_{p in l} w~_p v_p / c in [0, 1]; compare with
    ``CRITERIA[cid]["readout_threshold"]`` (already on this scale) or
    ``rule_threshold(crit)``.
    """
    crit = criterion if criterion is not None else CRITERIA[criterion_id or "O0"]
    if _mode(scoring) == "legacy":
        return _legacy_lineage_scores_from_panel(panel, p95, crit)
    return _v2_scores(panel, p95, crit)


def expected_from_true(
    true_panel: dict[str, float],
    p95: dict[str, float],
    criterion: dict[str, Any] | None = None,
    *,
    criterion_id: str | None = None,
    margin: float | None = None,
    scoring: str | None = None,
) -> str | None:
    """Post-hoc verifier label from holdout true ADT — never used as source event.

    ``margin`` (if given) is taken on the score scale of ``scoring`` as is.
    """
    crit = criterion if criterion is not None else CRITERIA[criterion_id or "O0"]
    mode = _mode(scoring)
    if mode == "legacy":
        m = _declared_value(crit, "expected_margin") if margin is None else float(margin)
        return _legacy_expected_from_true(true_panel, p95, crit, margin=m)
    m = effective_margin(crit, scoring=mode) if margin is None else float(margin)
    scores = _v2_scores(true_panel, p95, crit)
    ranked = sorted(scores.items(), key=lambda kv: -kv[1])
    if ranked[0][1] - ranked[1][1] < m:
        return None
    return ranked[0][0]


def event_value_under_criterion(
    raw_norm_value: float,
    is_key_marker: bool,
    action: str,
    crit: dict,
    *,
    scoring: str | None = None,
) -> float:
    """Map a p95-normalised TEDDY value into the declared field event value.

    v2: w~_p * clip(v, 0, 1) with w~ = w / max w — in [0, 1] without clipping
    the weighted value, so declared lineage priority survives into the field.
    """
    if _mode(scoring) == "legacy":
        return _legacy_event_value_under_criterion(raw_norm_value, is_key_marker, action, crit)
    v = min(1.0, max(0.0, float(raw_norm_value)))
    return float(protein_weight(bool(is_key_marker), action, crit) * v)


_apply_scoring_to_criteria()


# --------------------------------------------------------------------------
# ANM readout threshold on the same scale (simultaneous events)
# --------------------------------------------------------------------------
def field_gain(schema: dict[str, Any], n_events: int | None = None) -> float:
    """Field score of an action per unit of summed event value, all events at t=0.

    Replays the ``finite_graph_scalar`` recurrence of the ANM runner on one
    lineage (action site + its event sites; lineages are uncoupled). With all
    events simultaneous every event site is symmetric, so the ANM action score
    is ``G * sum_p value_p`` exactly, i.e. ``G * c * S_l`` in v2.
    """
    cfg = schema["field_representation"]
    steps, ret = int(cfg["steps"]), float(cfg["retention"])
    dif, src = float(cfg["diffusion"]), float(cfg["source_scale"])
    n = int(n_events or len(next(iter(LINEAGE_PANELS.values()))))
    a, e = 0.0, [0.0] * n
    for tick in range(steps + 1):
        na = ret * a
        ne = [ret * x + (src if tick == 0 else 0.0) for x in e]
        # propagate from next_state (event -> action w=1, action -> event w=0.15)
        pa = na + dif * 1.0 * sum(ne)
        pe = [x + dif * 0.15 * na for x in ne]
        a, e = pa, pe
    return a / n


def anm_readout_threshold(
    crit: dict[str, Any], schema: dict[str, Any], *, scoring: str | None = None
) -> float:
    """ANM readout threshold equivalent to the TEDDY-alone rule threshold.

    v2 with simultaneous events: ANM score = G * c * S_l, so thresholding the
    field at ``G * c * rule_threshold`` abstains exactly where the rule does.
    Legacy: the declared threshold (what the old scripts wrote into the schema).
    """
    if _mode(scoring) == "legacy":
        return _declared_value(crit, "readout_threshold")
    return field_gain(schema) * lineage_normaliser(crit) * rule_threshold(crit, scoring="v2")
