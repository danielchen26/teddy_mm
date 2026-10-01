"""E4 (C6, Experiment 2 redesigned): two-channel fusion, library part.

Registered design: registration/registration_v3.json -> experiments.E4, as amended by A1
(A1.1 tie-break and bootstrap, A1.2 panel-free channel-2 inputs, A1.9 comparator). The runner is
``bridge_anm/v3_e4_fusion.py``; the numbers the registration leaves open are fixed on train/val in
``registration/addenda/E4.json``. Nothing in this module chooses a number: every function takes the
registered or addendum values as arguments.

Vocabulary
  e1, e2        channel evidence for the 12 primary-panel proteins, in ``panel_flat`` order, in [0, 1]
                (channel 1 = TEDDY + head from RNA; channel 2 = ridge from the panel-free non-panel ADT)
  m12           measured panel evidence (stored ADT / training q95 of the stored ADT, clipped); key only
  S1, S2        class scores [n, 4] (LINEAGES order): equal-weight mean over each class panel
  s1, s2        ANM per-channel source scales, trust_c / max(trust_1, trust_2) at the noise level

ANM fusion (F5). One action per lineage. Support events: every panel protein of every class, per
channel, value s_c * e_c[p]. Contradiction events: for each channel, its top class k*_c (argmax of S_c,
first in LINEAGES order on ties) contradicts every other class j with value s_c * 0.5 * S_c[k*_c].
Times: channel 1 at t = 0, channel 2 at t = 0 (time-blind, primary) or t = 2 (declared order, control).
ANM's ``finite_graph_scalar`` has a single schema ``source_scale`` (1.0, registered); the injection is
``source_scale * amplitude`` with amplitude = +value (support) or -value (contradiction), so multiplying
the event value by s_c is exactly a per-channel source scale.

Closed form. Every event site is joined only to its own action site (build_graph), so the four actions
are uncoupled stars and the field is linear: the action state is sum_e R(n_k, m_e) * amp_e, where n_k is
the number of event sites of action k, m_e the number of ticks from the event's injection to the end
(max event time + steps + 1 - t_e), and R(n, m) the action response to a unit injection at one event site
of an n-site star (``star_unit_response``). ``anm_fusion_closed`` is that formula; the runner checks it
against the engine on every cell.
"""
from __future__ import annotations

import copy
import json
import math
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from . import v3_amend as va
from . import v3_key as vk

LEVELS: tuple[str, ...] = ("L0", "L1", "L2", "L3")
# (channel-1 RNA thinned, channel-2 inputs dropped) per registered noise level
LEVEL_NOISE: dict[str, tuple[bool, bool]] = {"L0": (False, False), "L1": (True, False),
                                             "L2": (False, True), "L3": (True, True)}
CHANNELS: tuple[str, ...] = ("teddy_rna_head", "protein_regressor")
CONTRADICTION_FACTOR = 0.5
TAU_GRID: tuple[float, ...] = tuple(round(0.05 * i, 2) for i in range(1, 13))
ALPHA_GRID: tuple[float, ...] = (0.1, 1.0, 10.0, 100.0, 1000.0)
SCORE_DECIMALS = 12  # confidence scores are rounded before ranking so float noise cannot split a true tie


# --------------------------------------------------------------------------- panel and scores
def panel_dict(reg: dict[str, Any]) -> dict[str, list[str]]:
    return {k: list(v) for k, v in reg["panels"]["primary"]["proteins"].items()}


def panel_flat(reg: dict[str, Any]) -> list[str]:
    pan = panel_dict(reg)
    return [p for k in vk.LINEAGES for p in pan[k]]


def class_scores12(v12: np.ndarray, reg: dict[str, Any]) -> np.ndarray:
    """Class scores [n, 4] from a 12-column array in ``panel_flat`` order (equal-weight panel means)."""
    v12 = np.asarray(v12, dtype=np.float64)
    pan = panel_dict(reg)
    out, j = [], 0
    for k in vk.LINEAGES:
        n = len(pan[k])
        out.append(v12[:, j:j + n].mean(axis=1))
        j += n
    if j != v12.shape[1]:
        raise ValueError(f"expected {j} panel columns, got {v12.shape[1]}")
    return np.stack(out, axis=1)


def top_class(S: np.ndarray) -> np.ndarray:
    """argmax lineage index (first in LINEAGES order on ties, as numpy argmax)."""
    return np.asarray(S, dtype=np.float64).argmax(axis=1)


def calls_from_scores(S: np.ndarray, bar: float | None = None) -> np.ndarray:
    """argmax lineage; with ``bar``, no call ('') when the top score is below it."""
    S = np.asarray(S, dtype=np.float64)
    calls = np.asarray(vk.LINEAGES, dtype=object)[top_class(S)]
    if bar is not None:
        calls[S.max(axis=1) < float(bar)] = vk.NO_CALL
    return calls.astype(str)


def val_bar(conf: np.ndarray, no_call_rate: float) -> float:
    """Registered bar rule: val quantile of the confidence at the declared no-call rate, rounded to 6."""
    return float(round(float(np.quantile(np.asarray(conf, dtype=np.float64), float(no_call_rate))), 6))


# --------------------------------------------------------------------------- E4 key
def measured_evidence(adt: np.ndarray, adt_names: Iterable[str], proteins: Sequence[str],
                      q95: Sequence[float]) -> np.ndarray:
    """Stored ADT of ``proteins`` / training q95 of the stored ADT, clipped to [0, 1] (key / trust target)."""
    names = [str(x) for x in adt_names]
    cols = [names.index(p) for p in proteins]
    x = np.asarray(adt, dtype=np.float64)[:, cols]
    return np.clip(x / np.maximum(np.asarray(q95, dtype=np.float64), 1e-12)[None, :], 0.0, 1.0)


def e4_key(m12: np.ndarray, reg: dict[str, Any], tau: float) -> np.ndarray:
    """E4 key: Q1 scoring on measured panel evidence; argmax lineage if the top score >= tau, else OUT."""
    S = class_scores12(m12, reg)
    key = np.asarray(vk.LINEAGES, dtype=object)[top_class(S)]
    key[S.max(axis=1) < float(tau)] = "OUT"
    return key.astype(str)


def choose_tau(m12_val: np.ndarray, annot_val: np.ndarray, reg: dict[str, Any],
               grid: Sequence[float] = TAU_GRID) -> tuple[float, list[dict[str, float]]]:
    """tau_K: maximum Cohen kappa (5 classes) of the E4 key against the annotation class on val;
    kappa compared at 6 decimals, ties to the smaller tau."""
    rows = [{"tau": float(t), "val_kappa": round(vk.cohen_kappa(e4_key(m12_val, reg, t), annot_val), 6)} for t in grid]
    best = sorted(rows, key=lambda r: (-r["val_kappa"], r["tau"]))[0]
    return float(best["tau"]), rows


# --------------------------------------------------------------------------- channel 2
def ridge_fit_select(w_tr: np.ndarray, y_tr: np.ndarray, w_va: np.ndarray, y_va: np.ndarray,
                     alphas: Sequence[float] = ALPHA_GRID) -> tuple[float, list[dict[str, float]], Any]:
    """Ridge (sklearn, fit_intercept) from channel-2 inputs to the 12 measured panel values; alpha with the
    highest mean val R2 over the 12 targets (6 decimals, ties to the smaller alpha); returns the refit model."""
    from sklearn.linear_model import Ridge

    rows, models = [], {}
    for a in alphas:
        m = Ridge(alpha=float(a)).fit(w_tr, y_tr)
        P = m.predict(w_va)
        r2 = [1.0 - float(np.sum((y_va[:, k] - P[:, k]) ** 2)) / max(float(np.sum((y_va[:, k] - y_va[:, k].mean()) ** 2)), 1e-300)
              for k in range(y_va.shape[1])]
        rows.append({"alpha": float(a), "val_mean_r2": round(float(np.mean(r2)), 6)})
        models[float(a)] = m
    best = sorted(rows, key=lambda r: (-r["val_mean_r2"], r["alpha"]))[0]["alpha"]
    return float(best), rows, models[float(best)]


def linear_predict(w: np.ndarray, coef: np.ndarray, intercept: np.ndarray) -> np.ndarray:
    return np.asarray(w, dtype=np.float64) @ np.asarray(coef, dtype=np.float64).T + np.asarray(intercept, dtype=np.float64)[None, :]


def normalised_evidence(pred: np.ndarray, q95: Sequence[float]) -> np.ndarray:
    """prediction / training q95 of the predictions, clipped to [0, 1] (as the registered channel-1 evidence)."""
    return np.clip(np.asarray(pred, dtype=np.float64) / np.maximum(np.asarray(q95, dtype=np.float64), 1e-12)[None, :], 0.0, 1.0)


def dropout_inputs(w: np.ndarray, cell_ids: Sequence[int], fraction: float, seed: int) -> np.ndarray:
    """L2 noise: per cell, round(fraction * n_inputs) channel-2 inputs set to 0, chosen by
    numpy default_rng([seed, 2, cell_id]).permutation(n_inputs) (label-free, the same for any subset)."""
    out = np.array(w, dtype=np.float64, copy=True)
    n_in = out.shape[1]
    k = int(round(float(fraction) * n_in))
    for r, cid in enumerate(np.asarray(cell_ids, dtype=np.int64)):
        out[r, np.random.default_rng([int(seed), 2, int(cid)]).permutation(n_in)[:k]] = 0.0
    return out


# --------------------------------------------------------------------------- RNA thinning (L1)
def recover_counts(row_values: np.ndarray, tol: float = 1e-3) -> tuple[np.ndarray, float, bool]:
    """Integer counts of one stored RNA row: row / its smallest non-zero value, rounded.
    Returns (counts, max relative deviation |c - round(c)| / max(1, round(c)), within tol). The deviation is
    relative because the stored values are float32: a count of 20,000 carries an absolute error near 1e-3."""
    v = np.asarray(row_values, dtype=np.float64)
    if v.size == 0:
        return np.zeros(0, dtype=np.int64), 0.0, True
    c = v / v.min()
    r = np.round(c)
    dev = float((np.abs(c - r) / np.maximum(1.0, r)).max())
    return r.astype(np.int64), dev, dev <= tol


def thin_counts(counts: np.ndarray, fraction: float, seed: int, cell_id: int) -> np.ndarray:
    """L1 noise: Binomial(count, fraction) per gene, numpy default_rng([seed, 1, cell_id])."""
    return np.random.default_rng([int(seed), 1, int(cell_id)]).binomial(np.asarray(counts, dtype=np.int64), float(fraction))


# --------------------------------------------------------------------------- trust
def pearson(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.size == 0 or a.max() == a.min() or b.max() == b.min():
        return 0.0  # a constant column carries no linear information
    return float(np.mean((a - a.mean()) * (b - b.mean())) / (a.std() * b.std()))


def channel_trust(ev: np.ndarray, meas: np.ndarray) -> tuple[float, list[float]]:
    """Registered trust: mean over the 12 panel proteins of max(0, Pearson(channel evidence, measured))."""
    r = [pearson(ev[:, j], meas[:, j]) for j in range(ev.shape[1])]
    return float(np.mean([max(0.0, x) for x in r])), r


def source_scales(trust: dict[str, float]) -> dict[str, float]:
    """ANM source_scale per channel = trust / max trust (6 decimals)."""
    mx = max(float(trust[c]) for c in CHANNELS)
    return {c: (round(float(trust[c]) / mx, 6) if mx > 0 else 0.0) for c in CHANNELS}


# --------------------------------------------------------------------------- simple fusion rules
def fuse_average(S1: np.ndarray, S2: np.ndarray) -> np.ndarray:
    """F1: simple average of the two channels' class scores (declared rule)."""
    return 0.5 * (np.asarray(S1, dtype=np.float64) + np.asarray(S2, dtype=np.float64))


def fuse_trust_weighted(S1: np.ndarray, S2: np.ndarray, t1: float, t2: float) -> np.ndarray:
    """F2: trust-weighted average (declared rule)."""
    tot = float(t1) + float(t2)
    if tot <= 0:
        return fuse_average(S1, S2)
    return (float(t1) * np.asarray(S1, dtype=np.float64) + float(t2) * np.asarray(S2, dtype=np.float64)) / tot


# --------------------------------------------------------------------------- ANM fusion: closed form
def star_unit_response(n_sites: int, m_ticks: int, field: dict[str, Any], propagate: bool = True) -> float:
    """Action state of an n-event-site star after m ticks, for a unit amplitude injected at one event site on
    the first of those ticks (ANM finite_graph_scalar update: retention, inject * source_scale, one Jacobi
    diffusion pass with event->action weight 1.0 and action->event weight 0.15)."""
    ret, dif, src = float(field["retention"]), float(field["diffusion"]), float(field["source_scale"])
    a = 0.0
    e = np.zeros(int(n_sites))
    for tick in range(int(m_ticks)):
        na = ret * a
        ne = ret * e
        if tick == 0:
            ne[0] += src
        if propagate:
            a = na + dif * 1.0 * float(ne.sum())
            e = ne + dif * 0.15 * na
        else:
            a, e = na, ne
    return float(a)


def anm_fusion_closed(e1: np.ndarray, e2: np.ndarray, s1: float, s2: float, reg: dict[str, Any], *,
                      t2: int = 0, contradiction: bool = True, propagate: bool = True) -> np.ndarray:
    """Re-coded closed form of the F5 action scores [n, 4] (LINEAGES order); see the module docstring."""
    fld = reg["anm"]["field_representation"]
    steps = int(fld["steps"])
    S = [class_scores12(e1, reg), class_scores12(e2, reg)]
    sc = [float(s1), float(s2)]
    times = [0, int(t2)]
    total = max(times) + steps + 1
    m = [total - t for t in times]
    pan = panel_dict(reg)
    npk = np.array([len(pan[k]) for k in vk.LINEAGES], dtype=np.int64)
    n = S[0].shape[0]
    top = [top_class(S[c]) for c in range(2)]
    topv = [S[c][np.arange(n), top[c]] for c in range(2)]
    contra_to = np.zeros((n, 4), dtype=np.int64)  # number of contradiction events attached to each action
    if contradiction:
        for c in range(2):
            contra_to += (np.arange(4)[None, :] != top[c][:, None]).astype(np.int64)
    n_sites = 2 * npk[None, :] + contra_to  # [n, 4]
    cache: dict[tuple[int, int], float] = {}

    def R(nn: int, mm: int) -> float:
        if (nn, mm) not in cache:
            cache[(nn, mm)] = star_unit_response(nn, mm, fld, propagate)
        return cache[(nn, mm)]

    out = np.zeros((n, 4), dtype=np.float64)
    for nn in np.unique(n_sites):
        sel = n_sites == nn
        for c in range(2):
            # sum of the channel's support values for the class = s_c * n_k * S_c,k
            out[sel] += R(int(nn), m[c]) * (sc[c] * (npk[None, :] * S[c]))[sel]
            if contradiction:
                amp = -sc[c] * CONTRADICTION_FACTOR * topv[c][:, None] * (np.arange(4)[None, :] != top[c][:, None])
                out[sel] += R(int(nn), m[c]) * amp[sel]
    return out


# --------------------------------------------------------------------------- ANM fusion: the engine
def anm_schema(base_schema: dict[str, Any], reg: dict[str, Any], readout_threshold: float) -> dict[str, Any]:
    """The bridge's ANM schema with the registered field and the two channel modalities."""
    schema = copy.deepcopy(base_schema)
    schema["field_representation"].update({k: v for k, v in reg["anm"]["field_representation"].items() if k != "kind"})
    schema["field_representation"]["kind"] = reg["anm"]["field_representation"]["kind"]
    schema["field_representation"]["readout_threshold"] = float(readout_threshold)
    schema["source_event_schema"]["allowed_modalities"] = list(CHANNELS)
    return schema


def fusion_events(e1_row: np.ndarray, e2_row: np.ndarray, s1: float, s2: float, reg: dict[str, Any], *,
                  t2: int = 0, contradiction: bool = True) -> list[dict[str, Any]]:
    """ANM source events of one cell (see the module docstring)."""
    pan = panel_dict(reg)
    flat = panel_flat(reg)
    rows = [np.asarray(e1_row, dtype=np.float64), np.asarray(e2_row, dtype=np.float64)]
    sc = [float(s1), float(s2)]
    times = [0, int(t2)]
    events = []
    for c, ch in enumerate(CHANNELS):
        for k in vk.LINEAGES:
            for p in pan[k]:
                v = sc[c] * float(rows[c][flat.index(p)])
                events.append({"event_id": f"{ch}:{p}", "time": times[c], "action": k, "modality": ch,
                               "polarity": "support", "value": min(1.0, max(0.0, v)), "provenance": f"{ch}:{p}"})
        if contradiction:
            S = class_scores12(rows[c][None, :], reg)[0]
            kt = int(np.argmax(S))
            v = sc[c] * CONTRADICTION_FACTOR * float(S[kt])
            for j, k in enumerate(vk.LINEAGES):
                if j == kt:
                    continue
                events.append({"event_id": f"{ch}:contradicts:{k}", "time": times[c], "action": k, "modality": ch,
                               "polarity": "contradiction", "value": min(1.0, max(0.0, v)),
                               "provenance": f"{ch}: top class {vk.LINEAGES[kt]} contradicts {k}"})
    return events


def anm_fusion_engine_cell(ffr, schema: dict[str, Any], events: list[dict[str, Any]], *,
                           propagate: bool = True) -> tuple[np.ndarray, str | None, int]:
    """Run ANM's finite_field_runner on one cell: validate -> build_graph -> evolve_field -> readout.
    Returns (action scores in LINEAGES order, recommended action or None, number of rejected events)."""
    actions = [{"id": k, "label": k} for k in vk.LINEAGES]
    inst = {"instance_id": "e4", "question": "lineage from two evidence channels", "actions": actions,
            "proposed_source_events": events}
    val = ffr.validate_source_events(schema, inst)
    graph = ffr.build_graph(inst, val["field_events"])
    field = ffr.evolve_field(schema, graph, val["field_events"], disable_propagation=not propagate)
    rd = ffr.readout(schema, inst, field["state"])
    return (np.array([rd["action_scores"][k] for k in vk.LINEAGES], dtype=np.float64), rd["recommended_action"],
            len(val["rejected"]))


def import_anm(anm_root: Path | str):
    import sys

    sys.path.insert(0, str(anm_root))
    from active_neural_matter.field import finite_field_runner as ffr  # noqa: E402

    return ffr


# --------------------------------------------------------------------------- matched coverage
def ranking_order(conf: np.ndarray, cell_ids: np.ndarray, amend: dict[str, Any]) -> np.ndarray:
    """Positions sorted by confidence (desc, rounded to SCORE_DECIMALS; NaN last), ties by the A1.1 rank."""
    c = np.asarray(conf, dtype=np.float64)
    key = np.round(np.where(np.isnan(c), -np.inf, c), SCORE_DECIMALS)
    rank = va.tie_break_rank(np.asarray(cell_ids, dtype=np.int64), amend)
    return np.lexsort((rank, -key))


def select_at_coverage(conf: np.ndarray, cell_ids: np.ndarray, coverage: float, amend: dict[str, Any]) -> np.ndarray:
    """A1.1 selection (v3_amend.coverage_select) on rounded scores."""
    c = np.round(np.asarray(conf, dtype=np.float64), SCORE_DECIMALS)
    return va.coverage_select(c, np.asarray(cell_ids, dtype=np.int64), coverage, amend)


def weighted_selective_accuracy(order: np.ndarray, W: np.ndarray, correct: np.ndarray, scored: np.ndarray,
                                coverage: float) -> np.ndarray:
    """Selective accuracy at matched coverage for multiset replicates.

    W: [B, n] multiplicity of each cell in each replicate (a bootstrap draw); the replicate's N = W.sum(1);
    the ceil(coverage * N) most confident copies are called, in ``order`` (copies of one cell share its rank
    and are interchangeable, as in v3_amend.coverage_select). Returns correct / scored among the called."""
    W = np.asarray(W)
    Wo = W[:, order].astype(np.int64)
    N = Wo.sum(axis=1)
    k = np.ceil(float(coverage) * N - 1e-9).astype(np.int64)
    before = np.cumsum(Wo, axis=1) - Wo
    sel = np.clip(k[:, None] - before, 0, Wo)
    nc = (sel * np.asarray(correct, dtype=np.int64)[order][None, :]).sum(axis=1)
    ns = (sel * np.asarray(scored, dtype=np.int64)[order][None, :]).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(ns > 0, nc / np.maximum(ns, 1), np.nan)


def two_stage_weights(donor: np.ndarray, B: int, seed: int, donors: Sequence[str] | None = None) -> np.ndarray:
    """Two-stage bootstrap multiplicities [B, n]: donors with replacement, then cells with replacement within
    each drawn donor (registered; one generator, numpy default_rng(seed), donors in sorted order)."""
    donor = np.asarray(donor).astype(str)
    ds = sorted(set(donor)) if donors is None else [str(d) for d in donors]
    by = {d: np.where(donor == d)[0] for d in ds}
    rng = np.random.default_rng(int(seed))
    W = np.zeros((int(B), donor.size), dtype=np.int32)
    for b in range(int(B)):
        for d in rng.choice(np.asarray(ds, dtype=object), size=len(ds), replace=True):
            cells = by[str(d)]
            np.add.at(W[b], cells[rng.integers(0, cells.size, size=cells.size)], 1)
    return W


def decide(point: float, lo: float, hi: float, per_donor: Sequence[float], margin: float) -> str:
    """Registered decision rule: win / loss / equivalent / inconclusive."""
    pd_ = [float(x) for x in per_donor]
    if point >= margin and lo > 0 and all(d >= margin / 2 for d in pd_):
        return "win"
    if point <= -margin and hi < 0 and all(d <= -margin / 2 for d in pd_):
        return "loss"
    if lo > -margin and hi < margin:
        return "equivalent"
    return "inconclusive"


def json_dump(obj: Any) -> str:
    return json.dumps(obj, indent=1, sort_keys=True, allow_nan=False) + "\n"


def _jsonable(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return _jsonable(x.tolist())
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x


def jsonable(x: Any) -> Any:
    return _jsonable(x)
