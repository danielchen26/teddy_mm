"""E6 helpers (external confirmation, registration experiments.E6, addendum registration/addenda/E6.json).

What is here, all label-free unless it says "key":

* the external annotation map (celltype.l2 -> class) and its canonical sha256;
* the registered gate with the E6 missing-protein rule: a protein that the external panel does not measure
  removes the gate condition that uses it; a clause whose defining condition is gone is dropped
  (``reduced_gate_spec``, ``gate_class_from_spec``); with every protein present it equals
  ``v3_key.gate_class_from_high`` exactly (tests/test_v3_e6.py);
* the registered per-batch threshold estimator applied within each external donor (``per_donor_thresholds``);
* the E6 keys and the key-validation rule (annotation vs per-donor gate kappa < 0.70 -> 'key not validated');
* the rule scores on an explicit panel and the exact closed form of ANM's ``finite_graph_scalar`` field on an
  explicit panel (per class, with that class's own number of events), used for the E1.1a bridge check; the E6
  primary uses the full registered panels (addendum version 2), and panels with absent proteins removed
  (``reduced_panel``) give the reduced-panel sensitivity reading;
* the bootstrap evaluator (matched coverage, count-vector replicates, as v3_e1) and the E6 replication rule.

Nothing here reads data or chooses a number; registered numbers come from the caller.
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from . import v3_e1 as e1
from . import v3_key as vk

LINEAGES = vk.LINEAGES
CLASSES = vk.CLASSES
UNSCORED = vk.UNSCORED
NO_CALL = vk.NO_CALL
KAPPA_MIN = 0.70
KEY_NOT_VALIDATED = "key not validated"
CONTAMINATION_LABEL = "RNA possibly seen by TEDDY in pretraining"
LIN = np.asarray(LINEAGES)


# --------------------------------------------------------------------------- annotation map
def canonical_json(m: Mapping[str, str]) -> str:
    return json.dumps(dict(m), sort_keys=True, separators=(",", ":"))


def map_sha256(m: Mapping[str, str]) -> str:
    """sha256 of json.dumps(map, sort_keys=True, separators=(',', ':')) (the prep pack's convention)."""
    return hashlib.sha256(canonical_json(m).encode()).hexdigest()


def annotation_class(cell_types: Iterable[str], amap: Mapping[str, str]) -> np.ndarray:
    """Label -> class in B / T / NK / myeloid / OUT / unscored. Unknown labels raise."""
    allowed = set(CLASSES) | {UNSCORED}
    out = []
    for t in cell_types:
        t = str(t)
        if t not in amap:
            raise KeyError(f"label {t!r} has no class in the E6 annotation map")
        c = str(amap[t])
        if c not in allowed:
            raise ValueError(f"label {t!r} maps to {c!r}, not a v3 class or 'unscored'")
        out.append(c)
    return np.asarray(out, dtype=object).astype(str)


# --------------------------------------------------------------------------- gate with absent proteins
def registered_gate_clauses(gate: Mapping[str, Any]) -> list[dict[str, Any]]:
    """The registered hierarchical gate (registration gate.rules_in_order, v3_key.gate_class_from_high) as data.

    Each clause: label; all_high (every protein high); any_high (at least one high); count_high + count_min
    (at least count_min of them high); none_high (none of them high). First clause passed wins; else OUT."""
    b = list(gate["b_markers"])
    nk = list(gate["nk_markers"])
    return [
        {"name": "erythroid", "label": "OUT", "all_high": ["CD71"]},
        {"name": "T", "label": "T", "all_high": ["CD3"], "none_high": b + ["CD14"]},
        {"name": "B", "label": "B", "any_high": b, "none_high": ["CD3", "CD14"]},
        {"name": "NK", "label": "NK", "none_high": ["CD3"] + b + ["CD14", "CD33"], "count_high": nk,
         "count_min": int(gate["nk_min_high"])},
        {"name": "myeloid", "label": "myeloid", "none_high": ["CD3"] + b, "any_high": ["CD14", "CD11c", "CD33"]},
    ]


def reduced_gate_spec(gate: Mapping[str, Any], present: Iterable[str]) -> dict[str, Any]:
    """E6 missing-protein rule applied to the registered gate.

    A condition on an absent protein is removed from its clause: an absent protein leaves an any-of set or a
    count set (count_min kept) or a none-high set; a clause is dropped when an all-high protein is absent, or
    its any-of set becomes empty, or fewer count proteins remain than count_min (no cell could pass)."""
    present = set(map(str, present))
    clauses, dropped_clauses, removed = [], [], []
    for c in registered_gate_clauses(gate):
        r = {"name": c["name"], "label": c["label"]}
        drop_reason = None
        if c.get("all_high"):
            miss = [p for p in c["all_high"] if p not in present]
            if miss:
                drop_reason = f"all-high protein(s) {miss} not measured"
            r["all_high"] = [p for p in c["all_high"] if p in present]
        if c.get("any_high") is not None:
            keep = [p for p in c["any_high"] if p in present]
            removed += [(c["name"], "any_high", p) for p in c["any_high"] if p not in present]
            if not keep:
                drop_reason = drop_reason or "no any-of protein measured"
            r["any_high"] = keep
        if c.get("count_high") is not None:
            keep = [p for p in c["count_high"] if p in present]
            removed += [(c["name"], "count_high", p) for p in c["count_high"] if p not in present]
            r["count_high"], r["count_min"] = keep, int(c["count_min"])
            if len(keep) < int(c["count_min"]):
                drop_reason = drop_reason or f"fewer than {c['count_min']} count proteins measured"
        if c.get("none_high") is not None:
            r["none_high"] = [p for p in c["none_high"] if p in present]
            removed += [(c["name"], "none_high", p) for p in c["none_high"] if p not in present]
        if drop_reason:
            dropped_clauses.append({"name": c["name"], "reason": drop_reason})
        else:
            clauses.append(r)
    used = sorted({p for c in clauses for k in ("all_high", "any_high", "count_high", "none_high") for p in c.get(k, [])})
    return {"clauses": clauses, "dropped_clauses": dropped_clauses,
            "removed_conditions": [{"clause": a, "set": b, "protein": p} for a, b, p in removed],
            "proteins_used": used,
            "absent_gate_proteins": sorted(p for p in gate["proteins"] if p not in present)}


def gate_class_from_spec(H: Mapping[str, np.ndarray], spec: Mapping[str, Any], n: int | None = None) -> np.ndarray:
    """Apply a (reduced) gate spec to boolean 'high' arrays H[p]; first clause passed wins; else OUT."""
    if n is None:
        n = len(next(iter(H.values())))
    out = np.full(n, "OUT", dtype=object)
    assigned = np.zeros(n, dtype=bool)
    for c in spec["clauses"]:
        m = np.ones(n, dtype=bool)
        for p in c.get("all_high", []):
            m &= H[p]
        if c.get("any_high") is not None:
            a = np.zeros(n, dtype=bool)
            for p in c["any_high"]:
                a |= H[p]
            m &= a
        if c.get("count_high") is not None:
            cnt = np.zeros(n, dtype=np.int64)
            for p in c["count_high"]:
                cnt += H[p].astype(np.int64)
            m &= cnt >= int(c["count_min"])
        for p in c.get("none_high", []):
            m &= ~H[p]
        sel = m & ~assigned
        out[sel] = c["label"]
        assigned |= sel
    return out.astype(str)


def high_arrays(adt: np.ndarray, adt_names: Sequence[str], thresholds: Mapping[str, float]) -> dict[str, np.ndarray]:
    """'measured CLR value > threshold (strictly greater)'; NaN (not measured) is never high."""
    j = {str(n): i for i, n in enumerate(adt_names)}
    adt = np.asarray(adt, dtype=np.float64)
    out = {}
    for p, t in thresholds.items():
        x = adt[:, j[p]]
        out[p] = np.where(np.isnan(x), False, x > float(t))
    return out


def per_donor_thresholds(adt: np.ndarray, adt_names: Sequence[str], donors: np.ndarray, proteins: Sequence[str],
                         seed: int) -> dict[str, dict[str, float]]:
    """Registered estimator (gmm_nonzero_post50, v3_key.gmm_nonzero_threshold) within each donor, label-free."""
    j = {str(n): i for i, n in enumerate(adt_names)}
    donors = np.asarray(donors).astype(str)
    adt = np.asarray(adt)
    out: dict[str, dict[str, float]] = {}
    for d in sorted(set(donors.tolist())):
        m = donors == d
        out[d] = {p: vk.gmm_nonzero_threshold(np.asarray(adt[m, j[p]], dtype=np.float64), seed=int(seed)) for p in proteins}
    return out


def gate_per_donor(adt: np.ndarray, adt_names: Sequence[str], donors: np.ndarray, spec: Mapping[str, Any],
                   thresholds: Mapping[str, Mapping[str, float]]) -> np.ndarray:
    donors = np.asarray(donors).astype(str)
    out = np.empty(donors.size, dtype=object)
    for d, thr in thresholds.items():
        m = donors == d
        H = high_arrays(np.asarray(adt)[m], adt_names, {p: thr[p] for p in spec["proteins_used"]})
        out[m] = gate_class_from_spec(H, spec, n=int(m.sum()))
    return out.astype(str)


def e6_keys(annot: np.ndarray, gated: np.ndarray) -> dict[str, np.ndarray]:
    """primary = annotation and per-donor gate agree (else unscored); annotation_only = annotation class.
    An annotation 'unscored' (Doublet) stays unscored in both."""
    annot = np.asarray(annot).astype(str)
    gated = np.asarray(gated).astype(str)
    prim = np.where((annot == gated) & (annot != UNSCORED), annot, UNSCORED)
    return {"annotation": annot, "gated_per_donor": gated, "primary": prim.astype(str),
            "annotation_only": annot.copy()}


def key_validation(annot: np.ndarray, gated: np.ndarray, donors: np.ndarray) -> dict[str, Any]:
    """Pooled 5-class kappa (annotation vs per-donor gate) over cells with an annotation class (Doublet excluded);
    per donor (descriptive). Decision: kappa < 0.70 -> results on the annotation-only key, 'key not validated'."""
    annot = np.asarray(annot).astype(str)
    gated = np.asarray(gated).astype(str)
    donors = np.asarray(donors).astype(str)
    m = np.isin(annot, CLASSES)
    pooled = vk.agreement_report(annot[m], gated[m])
    per = {d: vk.agreement_report(annot[m & (donors == d)], gated[m & (donors == d)]) for d in sorted(set(donors.tolist()))}
    kap = pooled["kappa5"]
    valid = bool(kap is not None and np.isfinite(kap) and kap >= KAPPA_MIN)
    return {"n_cells_with_annotation_class": int(m.sum()), "n_annotation_unscored": int((~m).sum()),
            "pooled": pooled, "per_donor": per, "threshold": KAPPA_MIN, "kappa5": kap, "validated": valid,
            "decision_key": "primary" if valid else "annotation_only",
            "label": None if valid else KEY_NOT_VALIDATED}


# --------------------------------------------------------------------------- panels and scores
def reduced_panel(panel: Mapping[str, Sequence[str]], present: Iterable[str]) -> dict[str, list[str]]:
    present = set(map(str, present))
    out = {k: [p for p in panel[k] if p in present] for k in LINEAGES}
    empty = [k for k, v in out.items() if not v]
    if empty:
        raise ValueError(f"no measured panel protein left for {empty}")
    return out


def class_scores_panel(v134: np.ndarray, adt_names: Sequence[str], panel: Mapping[str, Sequence[str]]) -> np.ndarray:
    """Rule scores [n, 4] (LINEAGES order): equal-weight mean of the evidence over each class's panel."""
    j = {str(n): i for i, n in enumerate(adt_names)}
    v = np.asarray(v134, dtype=np.float64)
    return np.stack([v[:, [j[p] for p in panel[k]]].mean(axis=1) for k in LINEAGES], axis=1)


def anm_recoded_scores(v134: np.ndarray, adt_names: Sequence[str], panel: Mapping[str, Sequence[str]],
                       field: Mapping[str, Any], closure: Mapping[str, Any] | None = None) -> np.ndarray:
    """Exact closed form of the ANM finite_graph_scalar field on an explicit panel (events at t = 0).

    Each class is a star with its own number of events n_k (v3_key.field_star). closure=None: the action-site
    value G(n_k) * sum of its events; else the readout_coordinates readout closure_weight * min(event states)
    + mean_coordinate_weight * mean(event states) + direct_action_weight * action."""
    j = {str(n): i for i, n in enumerate(adt_names)}
    v = np.asarray(v134, dtype=np.float64)
    out = []
    for k in LINEAGES:
        a, e = vk.field_star(v[:, [j[p] for p in panel[k]]], dict(field))
        if closure is None:
            out.append(a)
        else:
            out.append(float(closure["closure_weight"]) * e.min(axis=1) + float(closure["mean_coordinate_weight"])
                       * e.mean(axis=1) + float(closure["direct_action_weight"]) * a)
    return np.stack(out, axis=1)


def calls_at(scores: np.ndarray, threshold: float) -> np.ndarray:
    """Argmax lineage (first in B, T, NK, myeloid order on ties) if the top score reaches the threshold."""
    return vk.rule_calls(scores, threshold)


def engine_calls(call_idx: np.ndarray) -> np.ndarray:
    c = np.asarray(call_idx, dtype=np.int64)
    return np.where(c >= 0, LIN[np.maximum(c, 0)], NO_CALL).astype(str)


# --------------------------------------------------------------------------- bootstrap evaluator
class Method:
    """One ranked arm: confidence order (A1.1 tie-break) and its own calls (as v3_e1_mode_b.Method)."""

    def __init__(self, name: str, conf: np.ndarray, call_idx: np.ndarray, ids: np.ndarray, amend: dict,
                 acc_covs: Sequence[str] = (), keys: Sequence[str] = (), c2_covs: Sequence[str] = ()):
        self.name = name
        self.order = e1.ranked_order(conf, ids, amend)
        self.call = LIN[np.asarray(call_idx, dtype=np.int64)]
        self.acc_covs, self.c2_covs, self.keys = list(acc_covs), list(c2_covs), list(keys)
        self.covs = list(dict.fromkeys(self.acc_covs + self.c2_covs))
        self.arr: dict = {}

    def bind(self, keys: Mapping[str, np.ndarray]) -> None:
        call = self.call[self.order]
        for kn in self.keys:
            key = np.asarray(keys[kn]).astype(str)[self.order]
            self.arr[kn] = {"scored": (key != UNSCORED).astype(np.float64),
                            "correct": (call == key).astype(np.float64),
                            "out": (key == "OUT").astype(np.float64),
                            "inscope": np.isin(key, LINEAGES).astype(np.float64)}


class Evaluator:
    """Every bootstrapped E6 quantity as a function of a count vector w over the evaluated cells.

    ``coverage_points``: name -> fixed coverage (grid) or name -> (indicator array) for a realised coverage
    (the share of cells a rule calls at its bar, recomputed per replicate as A1.1 requires)."""

    def __init__(self, methods: Sequence[Method], keys: Mapping[str, np.ndarray], grid: Sequence[float],
                 realised: Mapping[str, np.ndarray], key_names: Sequence[str]):
        self.methods = list(methods)
        self.grid = list(grid)
        self.gn = [f"g{c:.2f}" for c in self.grid]
        self.realised = {k: np.asarray(v, dtype=np.float64) for k, v in realised.items()}
        self.key_names = list(key_names)
        self.n_out = {kn: (np.asarray(keys[kn]).astype(str) == "OUT").astype(np.float64) for kn in self.key_names}
        for m in self.methods:
            m.bind(keys)
        self.N = len(next(iter(keys.values())))

    def stats(self, w: np.ndarray) -> dict[str, float]:
        w = np.asarray(w, dtype=np.float64)
        Nb = float(w.sum())
        cov = dict(zip(self.gn, self.grid))
        for k, ind in self.realised.items():
            cov[k] = float((w * ind).sum() / Nb) if Nb > 0 else float("nan")
        out: dict[str, float] = {"N": Nb}
        for k in self.realised:
            out[f"cov|{k}"] = cov[k]
        n_out_tot = {kn: float((w * a).sum()) for kn, a in self.n_out.items()}
        for m in self.methods:
            ks = [e1.coverage_k(cov[c], Nb) for c in m.covs]
            s = e1.select_ordered(w[m.order], ks)
            ci = {c: i for i, c in enumerate(m.covs)}
            for kn in m.keys:
                A = m.arr[kn]
                nsc, nco = s @ A["scored"], s @ A["correct"]
                for c in m.acc_covs:
                    out[f"{m.name}|{c}|{kn}|acc"] = _ratio(nco[ci[c]], nsc[ci[c]])
                if m.c2_covs:
                    no, ni = s @ A["out"], s @ A["inscope"]
                    for c in m.c2_covs:
                        out[f"{m.name}|{c}|{kn}|outdecl"] = 1.0 - _ratio(no[ci[c]], n_out_tot[kn])
                        out[f"{m.name}|{c}|{kn}|inscope"] = _ratio(nco[ci[c]], ni[ci[c]])
        return out


def _ratio(a: float, b: float) -> float:
    return float(a) / float(b) if b > 0 else float("nan")


def aurc_from(get, method: str, key: str, gn: Sequence[str], grid: Sequence[float]):
    """Registered AURC (trapezoid over the coverage grid / its range) for scalars or replicate arrays."""
    cols = np.stack([np.asarray(get(f"{method}|{g}|{key}|acc"), dtype=np.float64) for g in gn])
    dg = np.abs(np.diff(np.asarray(grid, dtype=np.float64))).reshape((-1,) + (1,) * (cols.ndim - 1))
    return np.sum(0.5 * (cols[1:] + cols[:-1]) * dg, axis=0) / (max(grid) - min(grid))


# --------------------------------------------------------------------------- replication rule
def _fin(x: Any) -> bool:
    return x is not None and isinstance(x, (int, float, np.integer, np.floating)) and math.isfinite(float(x))


def replication(v3_point: float | None, e6_point: float | None, ci95: Sequence[float | None],
                per_donor: Mapping[str, float | None]) -> dict[str, Any]:
    """E6 replication rule (registration experiments.E6.replication).

    The v3 direction is the sign of the v3 (site4, E1) point difference. A v3 conclusion replicates on the
    external data when (i) the E6 per-donor difference has that sign in a strict majority of the external
    donors (an undefined donor value counts against) and (ii) the E6 two-stage bootstrap 95% interval excludes
    0 on that side. A v3 point of exactly 0 or undefined gives no direction to replicate."""
    lo, hi = (list(ci95) + [None, None])[:2]
    pd = {str(k): (float(v) if _fin(v) else None) for k, v in per_donor.items()}
    out: dict[str, Any] = {"v3_point": float(v3_point) if _fin(v3_point) else None,
                           "e6_point": float(e6_point) if _fin(e6_point) else None,
                           "e6_ci95": [lo, hi], "e6_per_donor": pd, "n_donors": len(pd)}
    if not _fin(v3_point) or float(v3_point) == 0.0:
        out.update({"v3_direction": 0, "status": "no v3 direction",
                    "e6_ci_contains_0": bool(_fin(lo) and _fin(hi) and lo <= 0 <= hi), "replicates": None})
        return out
    d = 1 if float(v3_point) > 0 else -1
    n_same = sum(1 for v in pd.values() if v is not None and np.sign(v) == d)
    majority = n_same > len(pd) / 2
    ci_ok = bool((d > 0 and _fin(lo) and float(lo) > 0) or (d < 0 and _fin(hi) and float(hi) < 0))
    rep = bool(majority and ci_ok)
    out.update({"v3_direction": d, "n_donors_same_sign": n_same, "majority_same_sign": bool(majority),
                "ci_excludes_0_same_side": ci_ok, "replicates": rep,
                "status": "replicates" if rep else "does not replicate"})
    return out


def dig(d: Mapping[str, Any], path: Sequence[str]) -> Any:
    cur: Any = d
    for k in path:
        if not isinstance(cur, Mapping) or k not in cur:
            return None
        cur = cur[k]
    return cur
