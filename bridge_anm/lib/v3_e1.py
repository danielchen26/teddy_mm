"""E1 helpers (Mode B rerun on the v3 key): matched coverage on bootstrap count vectors, the
registered verdict rules, the ANM engine bridge and the E1.3 attribution helpers.

Matched coverage (registration experiments.common.matched_coverage, amended by A1.1): at coverage c a
method calls its ceil(c * N) most confident cells; ties in its score are broken by the rank of the
global cell id in the registered permutation (``v3_amend.tie_break_rank``), the same for every method.

Bootstrap replicates are count vectors: ``w[i]`` is the number of copies of cell i of the evaluated
split in the replicate, and the replicate's N is ``w.sum()``. Copies of one cell share its score and
its tie-break rank, so the k most confident copies of a replicate are the first k copies along the
method's order of the original cells (``ranked_order``). ``select_ordered`` returns how many copies
of each cell are selected; this equals ``v3_amend.coverage_select`` on the materialised replicate
(tests/test_v3_e1.py). Every endpoint, including the selection and any realised coverage that
defines a matched point, is therefore recomputed inside each replicate, as A1.1 requires.

Nothing here reads data, labels or a threshold; registered numbers come from the caller.
"""
from __future__ import annotations

import copy
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Iterable, Iterator, Sequence

import numpy as np

from . import v3_amend as va
from . import v3_key as vk

LINEAGES = vk.LINEAGES
CLASSES = vk.CLASSES


# --------------------------------------------------------------------------- matched coverage
def ranked_order(conf: np.ndarray, cell_ids: np.ndarray, amend: dict[str, Any]) -> np.ndarray:
    """Positions sorted most confident first; ties by the A1.1 permutation rank; NaN last.

    The same order ``v3_amend.coverage_select`` uses (lexsort on (rank, -confidence))."""
    conf = np.asarray(conf, dtype=np.float64)
    key = np.where(np.isnan(conf), -np.inf, conf)
    rank = va.tie_break_rank(np.asarray(cell_ids, dtype=np.int64), amend)
    return np.lexsort((rank, -key))


def coverage_k(coverage: float, n: float) -> int:
    """ceil(coverage * N), with the same 1e-9 guard as v3_amend.coverage_select."""
    return int(math.ceil(float(coverage) * float(n) - 1e-9))


def select_ordered(w_o: np.ndarray, ks: Sequence[int] | np.ndarray) -> np.ndarray:
    """Selected copies per cell, in order space, for each k in ``ks``.

    ``w_o``: copies of each cell, already permuted into the method's ranked order. Returns
    ``[len(ks), N]``: the first k copies along the order (a cell may be partly selected; its copies
    are interchangeable)."""
    w_o = np.asarray(w_o, dtype=np.float64)
    before = np.cumsum(w_o) - w_o
    ks = np.asarray(ks, dtype=np.float64).reshape(-1, 1)
    return np.clip(ks - before[None, :], 0.0, w_o[None, :])


def select_lowest(conf: np.ndarray, eligible: np.ndarray, cell_ids: np.ndarray,
                  amend: dict[str, Any]) -> np.ndarray:
    """Order of the eligible positions by ascending ``conf`` (ties by the A1.1 permutation rank).

    Used for 'the same number of lowest-margin calls' (E1.3c)."""
    pos = np.where(np.asarray(eligible, dtype=bool))[0]
    rank = va.tie_break_rank(np.asarray(cell_ids, dtype=np.int64)[pos], amend)
    c = np.asarray(conf, dtype=np.float64)[pos]
    c = np.where(np.isnan(c), np.inf, c)
    return pos[np.lexsort((rank, c))]


def aurc(acc: Sequence[float], grid: Sequence[float]) -> float:
    """Trapezoid area of selective accuracy over the coverage grid, divided by the grid's range."""
    acc = np.asarray(acc, dtype=np.float64)
    g = np.asarray(grid, dtype=np.float64)
    if acc.size != g.size or acc.size < 2:
        raise ValueError("aurc needs one accuracy per grid point (>= 2 points)")
    if np.any(~np.isfinite(acc)):
        return float("nan")
    area = float(np.sum(0.5 * (acc[1:] + acc[:-1]) * np.abs(np.diff(g))))
    return area / float(g.max() - g.min())


# --------------------------------------------------------------------------- bootstrap
def replicate_counts(groups: Sequence[np.ndarray], n: int, n_boot: int, seed: int,
                     two_stage: bool = True) -> Iterator[np.ndarray]:
    """Yield ``n_boot`` count vectors over ``n`` cells, in a fixed order from default_rng(seed).

    two_stage: draw len(groups) groups (donors) with replacement, then within each drawn group as many
    cells as it has, with replacement. One-stage (a single group, e.g. a per-donor row): cells only."""
    rng = np.random.default_rng(int(seed))
    groups = [np.asarray(g, dtype=np.int64) for g in groups]
    for _ in range(int(n_boot)):
        w = np.zeros(int(n), dtype=np.int64)
        drawn = rng.integers(0, len(groups), size=len(groups)) if two_stage else np.arange(len(groups))
        for gi in drawn:
            g = groups[int(gi)]
            w += np.bincount(g[rng.integers(0, g.size, size=g.size)], minlength=int(n))
        yield w


def percentile_ci(x: np.ndarray) -> tuple[float | None, float | None, int]:
    x = np.asarray(x, dtype=np.float64)
    f = x[np.isfinite(x)]
    if f.size == 0:
        return None, None, 0
    lo, hi = np.percentile(f, [2.5, 97.5])
    return float(lo), float(hi), int(f.size)


def verdict(point: float, boot: np.ndarray, per_donor: dict[str, float], margin: float) -> dict[str, Any]:
    """Registered decision rule (experiments.common.statistics).

    win: point >= margin, bootstrap 95% lower bound > 0, and in each primary donor the difference has the
    same sign and is >= margin / 2; loss: the mirror; equivalent: the interval inside (-margin, +margin);
    otherwise inconclusive."""
    lo, hi, nf = percentile_ci(boot)
    pd = {k: (None if v is None or not np.isfinite(v) else float(v)) for k, v in per_donor.items()}
    out = {"point": None if point is None or not np.isfinite(point) else float(point), "ci95": [lo, hi],
           "n_boot_finite": nf, "n_boot": int(np.asarray(boot).size), "per_donor": pd, "margin": float(margin)}
    if out["point"] is None or lo is None or any(v is None for v in pd.values()):
        out["verdict"] = "inconclusive"
        out["verdict_note"] = "undefined point, interval or donor value"
        return out
    p = out["point"]
    vals = list(pd.values())
    if p >= margin and lo > 0 and all(v >= margin / 2 for v in vals):
        out["verdict"] = "win"
    elif p <= -margin and hi < 0 and all(v <= -margin / 2 for v in vals):
        out["verdict"] = "loss"
    elif lo > -margin and hi < margin:
        out["verdict"] = "equivalent"
    else:
        out["verdict"] = "inconclusive"
    return out


# --------------------------------------------------------------------------- scores
def entropy_confidence(S: np.ndarray) -> np.ndarray:
    """1 - H(q) / log(K), q = S / sum(S) (registration confidence_scores.entropy); 0 when sum(S) = 0."""
    S = np.clip(np.asarray(S, dtype=np.float64), 0.0, None)
    tot = S.sum(axis=1, keepdims=True)
    q = np.divide(S, tot, out=np.full_like(S, 1.0 / S.shape[1]), where=tot > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        h = -np.sum(np.where(q > 0, q * np.log(q), 0.0), axis=1)
    return 1.0 - h / math.log(S.shape[1])


def argmax_calls(S: np.ndarray) -> np.ndarray:
    """Argmax lineage (first in LINEAGES order on ties), as v3_key.rule_calls and ANM's readout."""
    return np.asarray([LINEAGES[i] for i in np.asarray(S).argmax(axis=1)]).astype(str)


def tied_top_classes(S: np.ndarray) -> np.ndarray:
    """True where two or more classes share the top score (the call is then set by class order)."""
    S = np.asarray(S, dtype=np.float64)
    return np.sum(S == S.max(axis=1, keepdims=True), axis=1) > 1


def proba5(model, X: np.ndarray) -> np.ndarray:
    """Class probabilities in CLASSES order (zeros for classes the model never saw)."""
    out = np.zeros((X.shape[0], len(CLASSES)))
    if getattr(model, "single_class", None) is not None:
        out[:, CLASSES.index(model.single_class)] = 1.0
        return out
    P = model.predict_proba(X)
    for ci, c in enumerate(model.classes_):
        out[:, CLASSES.index(str(c))] = P[:, ci]
    return out


def lineage_conf_call(P5: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Classifier decision (registration classifier.primary.decision): the most probable lineage and its
    probability (OUT is never called; an OUT-probable cell has low lineage probabilities)."""
    L = np.asarray(P5)[:, [CLASSES.index(k) for k in LINEAGES]]
    return L.max(axis=1), argmax_calls(L)


class SingleClass:
    """A label-cost draw with one class: predicts that class (registration E1.exp5 'classifier')."""

    def __init__(self, c: str):
        self.single_class = str(c)
        self.classes_ = np.asarray([str(c)])


def fit_logreg(X: np.ndarray, y: np.ndarray, C: float, seed: int = 0):
    """Registered classifier: sklearn LogisticRegression (lbfgs, max_iter 3000), or SingleClass."""
    from sklearn.linear_model import LogisticRegression

    y = np.asarray(y).astype(str)
    cls = np.unique(y)
    if cls.size == 1:
        return SingleClass(cls[0])
    import warnings

    from sklearn.exceptions import ConvergenceWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ConvergenceWarning)
        return LogisticRegression(C=float(C), max_iter=3000, random_state=int(seed)).fit(X, y)


def nested_draw(pool: np.ndarray, n: int | str, seed: int) -> np.ndarray:
    """Label-cost draw: the first n cells of default_rng(seed).permutation(sorted pool); 'all' = the pool.

    Draws are nested across n for a given seed (declared in the E1 addendum)."""
    pool = np.sort(np.asarray(pool, dtype=np.int64))
    if n == "all":
        return pool
    perm = np.random.default_rng(int(seed)).permutation(pool.size)
    return pool[perm[: int(n)]]


# --------------------------------------------------------------------------- E1.3 attribution
def closed_form_top_sets(v_called: np.ndarray) -> np.ndarray:
    """[N, n] mask of the panel markers at the called class's maximum evidence (A1.3, exact equality)."""
    v = np.asarray(v_called, dtype=np.float64)
    return v == v.max(axis=1, keepdims=True)


def engine_top_sets(decrease: np.ndarray, rel_tol: float = 1e-12) -> np.ndarray:
    """[N, n] mask of the events whose deletion lowers the called action's score most (ties kept).

    ``decrease[i, p]`` = called action score before minus after deleting event p (engine)."""
    d = np.asarray(decrease, dtype=np.float64)
    m = d.max(axis=1, keepdims=True)
    return d >= m - rel_tol * np.maximum(1.0, np.abs(m))


def marker_shares(mask: np.ndarray, cls: np.ndarray, n_classes: int = 4) -> dict[str, np.ndarray]:
    """Fractional-credit shares (A1.3: a set of m markers gives 1/m to each).

    Returns 'all' [n_classes, n] = credit / number of calls, and 'within_class' = credit / calls of
    that class (NaN for a class with no calls)."""
    mask = np.asarray(mask, dtype=np.float64)
    cls = np.asarray(cls, dtype=np.int64)
    credit = mask / mask.sum(axis=1, keepdims=True)
    tot = np.zeros((n_classes, mask.shape[1]))
    np.add.at(tot, cls, credit)
    n_k = np.bincount(cls, minlength=n_classes).astype(np.float64)
    n = max(1, cls.size)
    with np.errstate(invalid="ignore", divide="ignore"):
        within = tot / n_k[:, None]
    return {"all": tot / n, "within_class": within, "n_calls_per_class": n_k}


def marker_share_null(mask: np.ndarray, cls: np.ndarray, n_perm: int, seed: int,
                      n_classes: int = 4, chunk: int = 50) -> dict[str, np.ndarray]:
    """Null for E1.3b: within-class permutation of marker identities per cell, same credit.

    Each permutation assigns, per cell, a uniformly random order of the called class's marker
    identities to its credit vector. Returns arrays [n_perm, n_classes, n] for 'all' and 'within_class'."""
    mask = np.asarray(mask, dtype=np.float64)
    cls = np.asarray(cls, dtype=np.int64)
    credit = mask / mask.sum(axis=1, keepdims=True)
    n, m = credit.shape
    rng = np.random.default_rng(int(seed))
    n_k = np.bincount(cls, minlength=n_classes).astype(np.float64)
    out_all = np.zeros((int(n_perm), n_classes, m))
    for s in range(0, int(n_perm), chunk):
        b = min(chunk, int(n_perm) - s)
        perm = np.argsort(rng.random((b, n, m)), axis=2)
        cp = np.take_along_axis(np.broadcast_to(credit, (b, n, m)), perm, axis=2)
        for k in range(n_classes):
            sel = cls == k
            if sel.any():
                out_all[s:s + b, k] = cp[:, sel].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        within = out_all / n_k[None, :, None]
    return {"all": out_all / max(1, n), "within_class": within}


def perm_p_two_sided(obs: float, null: np.ndarray) -> float | None:
    null = np.asarray(null, dtype=np.float64)
    null = null[np.isfinite(null)]
    if null.size == 0 or not np.isfinite(obs):
        return None
    c = null.mean()
    return float((1 + np.sum(np.abs(null - c) >= abs(obs - c) - 1e-15)) / (1 + null.size))


# --------------------------------------------------------------------------- ANM engine bridge
class AnmBridge:
    """ANM ``finite_graph_scalar`` (finite_field_runner, ANM_ROOT) on one cell's evidence events.

    One action per lineage; one support event per panel protein, value = its evidence, every event at
    t = 0, equal weights; ANM's registered default field. The engine path is that of
    ``finite_field_runner.run_instances``: preparation apparatus -> validate_source_events ->
    build_graph -> evolve_field -> readout (thresholds per question; closure readout via
    readout_coordinates over each action's event sites)."""

    def __init__(self, reg: dict[str, Any], schema_path: Path | str, anm_root: Path | str | None = None):
        root = Path(anm_root or os.environ.get("ANM_ROOT", ""))
        if not root or not (root / "active_neural_matter").exists():
            raise FileNotFoundError(f"ANM_ROOT does not hold active_neural_matter: {root!r}")
        if str(root) not in sys.path:
            sys.path.insert(0, str(root))
        from active_neural_matter.field import finite_field_runner as ffr

        self.ffr = ffr
        self.anm_root = root
        base = json.loads(Path(schema_path).read_text())
        fr = reg["anm"]["field_representation"]
        base["field_representation"].update({k: v for k, v in fr.items()})
        base["source_event_schema"]["allowed_modalities"] = ["evidence"]
        self.base = base
        self.closure_cfg = dict(reg["anm"]["closure_readout"])
        self.thresholds = {"Q1": vk.anm_readout_threshold(reg, "Q1"), "Q2": vk.anm_readout_threshold(reg, "Q2"),
                           "Q3": vk.anm_readout_threshold(reg, "Q3"),
                           "closure_Q1": float(reg["anm"]["closure_bar_Q1"]),
                           "closure_Q2": float(reg["anm"]["closure_bar_Q2"])}
        self.schemas = {}
        for name, thr in self.thresholds.items():
            s = copy.deepcopy(base)
            s["field_representation"]["readout_threshold"] = float(thr)
            self.schemas[name] = s

    def _instance(self, panel: dict[str, list[str]], values: dict[str, Sequence[float]], drop=None) -> dict:
        events = []
        for k in LINEAGES:
            for p, v in zip(panel[k], values[k]):
                if drop is not None and drop == (k, p):
                    continue
                events.append({"event_id": f"{k}:{p}", "time": 0, "action": k, "modality": "evidence",
                               "polarity": "support", "value": float(v), "provenance": "teddy_head_v3_evidence"})
        return {"instance_id": "cell", "question": "lineage call from TEDDY + head evidence",
                "actions": [{"id": k, "label": k} for k in LINEAGES], "proposed_source_events": events}

    def run(self, panel: dict[str, list[str]], values: dict[str, Sequence[float]], readouts: Iterable[str],
            closure_readouts: Iterable[str] = (), drop=None) -> dict[str, Any]:
        """Evolve the field once; return each readout's action scores [4] and call index (-1 = no call)."""
        ffr = self.ffr
        inst = self._instance(panel, values, drop)
        proposed, _ = ffr.apply_preparation_to_field_apparatus(self.base, inst, inst["proposed_source_events"],
                                                               {"source": "instance_file"})
        inst["proposed_source_events"] = proposed
        val = ffr.validate_source_events(self.base, inst)
        graph = ffr.build_graph(inst, val["field_events"])
        field = ffr.evolve_field(self.base, graph, val["field_events"])
        out: dict[str, Any] = {"n_rejected": len(val["rejected"]), "n_repaired": len(val["repaired"])}
        for r in readouts:
            rd = ffr.readout(self.schemas[r], inst, field["state"])
            out[r] = (np.asarray([rd["action_scores"][k] for k in LINEAGES]),
                      -1 if rd["recommended_action"] is None else LINEAGES.index(rd["recommended_action"]))
        closure_readouts = list(closure_readouts)
        if closure_readouts:
            ci = dict(inst)
            ci["readout_coordinates"] = {k: [f"event:{k}:{p}" for p in panel[k]
                                             if not (drop is not None and drop == (k, p))] for k in LINEAGES}
            ci["readout_config"] = dict(self.closure_cfg)
            for r in closure_readouts:
                rd = ffr.readout(self.schemas[r], ci, field["state"])
                out[r] = (np.asarray([rd["action_scores"][k] for k in LINEAGES]),
                          -1 if rd["recommended_action"] is None else LINEAGES.index(rd["recommended_action"]))
        return out

    def anm_commit(self) -> str | None:
        import subprocess

        try:
            r = subprocess.run(["git", "-C", str(self.anm_root), "rev-parse", "HEAD"], capture_output=True, text=True)
            return r.stdout.strip() or None
        except OSError:
            return None
