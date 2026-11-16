"""Scope: by what kind of readout is a response readable at a stage?

An audit (runner.py) asks at which stage of a chain a response stops being readable. This module asks a second
question at each stage: readable by whom. The same response is read by readouts that differ only in the cells they
were learned from:

  frozen   the system's own output, as it is (no fit)
  train    a readout learned from the system's training cells
  other    a readout learned from another study of the same kind of event
  units    a readout learned from the other units of the audited study (never a cell of the unit it is scored on)
  own      a readout learned with the audited unit's own cells in the fit (cross-fitted by cell)

One statistic for every scope. With A the event indicator (a treated unit against control cells), O a readout and P the
scope's out-of-sample prediction of O:

    theta = Delta(O) - a * Delta(P),    a = slope of O on P among the control cells,    u = theta / Delta(O)

Delta is the treated-minus-control difference of means. u = 0: the prediction carries the whole response; u = 1: none of
it. The slope puts a prediction made in other units (a head trained on another study's protein scale) on the scale of
the measurement, and it is what a regression of O on P with a treated offset gives when the slope is learned from the
controls alone. For several readouts at once, u is the norm of the unexplained differences over the norm of the
differences, each in units of the control spread; for one readout that is |u|.

A second target is read with the same scopes: the variation among the cells of the treated unit (what a system's
ordinary score measures). u = 1 - the mean over readouts of the squared correlation between P and O within the unit.

What is missing follows from where the scopes part (read()):

  no deficit   frozen supported
  function     train supported, frozen rejected: the system's own training cells allow it; its map does not deliver it
  coverage     units supported, train rejected: it takes cells of this kind
  unit         own supported, units rejected: it takes this unit's own cells
  upstream     own rejected at the stage the system reads, supported at the input view
  unreadable   own rejected at the input view too (or no input view was given)

The gaps between consecutive scopes are reported with paired intervals whether or not a verdict is reached, and they are
what a forecast for an intervention is made of (forecast()). Intervals: treated cells are resampled by cluster when the
treated side has several clusters, otherwise by cell; control cells by block when the controls come in blocks of many
cells (donors), otherwise by cell. The same resamples serve every stage and scope of an event.

This module fits readouts and nothing else. It never changes the audited system.
"""
from __future__ import annotations

import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .stats import ci, make_learner

SCOPES = ("frozen", "train", "other", "units", "own")
ORDER = ("frozen", "train", "units", "own")        # each scope may use everything the one before it may
DEFAULT_CLASSES = {"ridge": {"regressor": "RidgeCV", "alphas": [1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0]}}
DEFAULT_DECISION = {"margin_u": 0.10, "step_margin": 0.10, "alpha": 0.05, "bootstrap": 1000, "seed": 0, "min_clusters": 10}
MISSING = {"none": "no deficit", "function": "the system's map: its own training cells allow it",
           "coverage": "coverage: it takes cells of this kind", "unit": "the unit's own cells",
           "upstream": "lost before the stage the system reads", "unreadable": "not readable from the views given",
           "unresolved": "unresolved", "not testable": "no resolved response"}


class ScopeError(Exception):
    pass


# ------------------------------------------------------------------ readouts
def _fit(spec, X, Y):
    """A readout of one class: standardized state -> readouts. Boosting is fitted one readout at a time."""
    Y = np.asarray(Y, np.float64)
    if spec["regressor"] == "HistGradientBoostingRegressor":
        if Y.shape[1] > 4:
            raise ScopeError("the boosting class is fitted per readout; declare ridge for more than four readouts")
        return [make_learner(spec, "m").fit(X, Y[:, j]) for j in range(Y.shape[1])]
    return make_pipeline(StandardScaler(), make_learner(spec, "m")).fit(X, Y)


def _predict(model, X):
    if isinstance(model, list):
        return np.column_stack([m.predict(X) for m in model]).astype(np.float32)
    P = model.predict(X)
    return (P[:, None] if P.ndim == 1 else P).astype(np.float32)


def fold_blocks(n, k, seed, prefix="ctl"):
    """Pseudo-blocks for cells that have no unit of their own (the controls of a perturbation audit)."""
    out = np.empty(n, dtype=object)
    perm = np.random.default_rng(seed).permutation(n)
    for f in range(k):
        out[perm[f::k]] = f"{prefix}:{f}"
    return out.astype(str)


# ------------------------------------------------------------------ a study: cells, views, and the fits they share
class Study:
    """The cells of one study and what readouts may be learned from.

    Y        (n, R) measured readouts
    stages   {name: (n, d)} views of each cell, in the order of the chain
    block    (n,) the unit a cell may not be fitted with when it is scored (donor, target gene, control fold)
    train    optional {"stages": {name: (m, d)}, "Y": (m, R), "block": (m,) or None}: the system's training cells
    other    optional {"stages": ..., "Y": ...}: another study of the same kind of event
    frozen   optional (n, R): the system's own output for the readouts
    classes  readout classes, as in a declaration
    max_fit  optional cap on the cells one readout is fitted on (a seeded subsample)
    """

    def __init__(self, Y, stages, block, *, train=None, other=None, frozen=None, classes=None, folds=5, seed=0, max_fit=None):
        self.Y = np.asarray(Y, np.float64)
        if self.Y.ndim == 1:
            self.Y = self.Y[:, None]
        self.n, self.R = self.Y.shape
        self.stages = {k: np.asarray(v, np.float32) for k, v in stages.items()}
        self.block = np.asarray(block).astype(str)
        self.train, self.other = train, other
        self.frozen = None if frozen is None else np.asarray(frozen, np.float64).reshape(self.n, -1)
        self.classes = classes or DEFAULT_CLASSES
        self.folds, self.seed, self.max_fit = int(folds), int(seed), max_fit
        self._models, self._own, self._other = {}, {}, {}
        for k, v in self.stages.items():
            if len(v) != self.n:
                raise ScopeError(f"stage {k} has {len(v)} rows, the study has {self.n} cells")
        if train is not None and train.get("block") is not None:
            train["block"] = np.asarray(train["block"]).astype(str)

    # ---- the fits
    def _model(self, where, stage, cls, excluded):
        key = (where, stage, cls, frozenset(excluded))
        if key not in self._models:
            if where == "units":
                keep = ~np.isin(self.block, list(excluded))
                X, Y = self.stages[stage][keep], self.Y[keep]
            else:
                tb = self.train.get("block")
                keep = np.ones(len(self.train["Y"]), bool) if tb is None else ~np.isin(tb, list(excluded))
                X, Y = np.asarray(self.train["stages"][stage], np.float32)[keep], np.asarray(self.train["Y"], np.float64)[keep]
            if len(X) < 50:
                raise ScopeError(f"{where} readout at {stage}: {len(X)} cells left to fit on")
            if self.max_fit and len(X) > self.max_fit:             # a seeded subsample; the same for every stage and class
                pick = np.sort(np.random.default_rng(self.seed + 5).choice(len(X), self.max_fit, replace=False))
                X, Y = X[pick], Y[pick]
            self._models[key] = _fit(self.classes[cls], X, Y.reshape(len(X), -1))
        return self._models[key]

    def predict(self, scope, stage, cls, cells, treated_blocks=()):
        """Out-of-sample predictions for the given cells under one scope. treated_blocks: blocks no fit may contain
        (the treated unit of a one-unit event), on top of each cell's own block."""
        cells = np.asarray(cells)
        S = self.stages[stage]
        out = np.empty((len(cells), self.R), np.float32)
        if scope == "own":
            key = (stage, cls)
            if key not in self._own:
                P = np.empty((self.n, self.R), np.float32)
                fold = np.empty(self.n, int)
                rng = np.random.default_rng(self.seed + 11)
                for b in np.unique(self.block):                      # folds by cell, spread evenly over every block
                    idx = rng.permutation(np.where(self.block == b)[0])
                    fold[idx] = np.arange(len(idx)) % self.folds
                for f in range(self.folds):
                    P[fold == f] = _predict(_fit(self.classes[cls], S[fold != f], self.Y[fold != f]), S[fold == f])
                self._own[key] = P
            return self._own[key][cells]
        if scope == "other":
            key = (stage, cls)
            if key not in self._other:
                self._other[key] = _fit(self.classes[cls], np.asarray(self.other["stages"][stage], np.float32),
                                        np.asarray(self.other["Y"], np.float64).reshape(len(self.other["Y"]), -1))
            return _predict(self._other[key], S[cells])
        if scope not in ("units", "train"):
            raise ScopeError(f"no fitted scope {scope}")
        blocks = self.block[cells]
        train_blocks = set() if scope == "units" or self.train.get("block") is None else set(np.unique(self.train["block"]))
        for b in np.unique(blocks):
            excl = {b} | set(treated_blocks)
            if scope == "train":
                excl &= train_blocks
            m = blocks == b
            out[m] = _predict(self._model(scope, stage, cls, excl), S[cells[m]])
        return out

    def scopes_available(self):
        return [s for s in SCOPES if (s == "frozen" and self.frozen is not None) or (s == "train" and self.train is not None)
                or (s == "other" and self.other is not None) or s in ("units", "own")]


# ------------------------------------------------------------------ the statistic and its resamples
class _Side:
    """One side of an event (treated or control) as resampling elements with their sums."""

    def __init__(self, Y, cluster):
        self.Y = Y
        if cluster is None:
            self.E, self.M = len(Y), None                  # every cell is its own element
        else:
            from scipy.sparse import csr_matrix
            el, inv = np.unique(cluster, return_inverse=True)
            self.E = len(el)
            self.M = csr_matrix((np.ones(len(Y)), (inv, np.arange(len(Y)))), shape=(self.E, len(Y)))
        self.n = self._sum(np.ones((len(Y), 1)))
        self.sy, self.syy = self._sum(Y), self._sum(Y * Y)

    def _sum(self, X):
        X = np.asarray(X, np.float64)
        return X if self.M is None else np.asarray(self.M @ X)

    def moments(self, P, counts):
        """Per draw: mean of Y, mean of P, var of P, cov(P, Y), var of Y. counts: (B, E)."""
        n = counts @ self.n
        my, mp = (counts @ self.sy) / n, (counts @ self._sum(P)) / n
        vp = (counts @ self._sum(P * P)) / n - mp * mp
        cpy = (counts @ self._sum(P * self.Y)) / n - mp * my
        vy = (counts @ self.syy) / n - my * my
        return my, mp, vp, cpy, vy


def _counts(E, B, rng):
    return np.vstack([np.ones((1, E)), rng.multinomial(E, np.full(E, 1.0 / E), size=B)]).astype(np.float64)   # row 0: the data


def _u_shift(t, c, P_t, P_c, Ct, Cc, sd_c):
    """u for every draw (row 0 is the point value) and per readout at the point."""
    myt, mpt, _, _, _ = t.moments(P_t, Ct)
    myc, mpc, vpc, cpyc, _ = c.moments(P_c, Cc)
    a = np.where(vpc > 1e-12, cpyc / np.where(vpc > 1e-12, vpc, 1.0), 0.0)
    dY, dP = myt - myc, a * (mpt - mpc)
    if dY.shape[1] == 1:
        u = 1.0 - dP[:, 0] / np.where(np.abs(dY[:, 0]) > 1e-12, dY[:, 0], np.nan)
        per = u[:1]
    else:
        u = np.linalg.norm((dY - dP) / sd_c, axis=1) / np.maximum(np.linalg.norm(dY / sd_c, axis=1), 1e-12)
        per = 1.0 - dP[0] / np.where(np.abs(dY[0]) > 1e-12, dY[0], np.nan)
    return u, per, a[0]


def _u_variation(t, P_t, Ct):
    _, _, vp, cpy, vy = t.moments(P_t, Ct)
    r2 = np.where((vp > 1e-12) & (vy > 1e-12), cpy * cpy / np.where((vp > 1e-12) & (vy > 1e-12), vp * vy, 1.0), 0.0)
    return 1.0 - r2.mean(axis=1), 1.0 - r2[0]


def _verdict(iv, margin, signed):
    lo, hi = iv
    if signed:
        return "support" if (-margin < lo and hi < margin) else ("reject" if (lo > margin or hi < -margin) else "unresolved")
    return "support" if hi < margin else ("reject" if lo > margin else "unresolved")


def _union(by_class, margin, signed):
    """As verdicts.verdict_union: supported if some class supports, rejected only if every class rejects."""
    v = {c: _verdict(x["ci"], margin, signed) for c, x in by_class.items()}
    sup = [c for c in v if v[c] == "support"]
    if sup:
        return "support", min(sup, key=lambda c: max(abs(by_class[c]["ci"][0]), abs(by_class[c]["ci"][1])))
    best = min(by_class, key=lambda c: abs(by_class[c]["u"]))
    return ("reject" if all(x == "reject" for x in v.values()) else "unresolved"), best


def scope_grid(study, treated, control, *, event=None, stages=None, scopes=None, treated_cluster=None, reads=None, input_view=None,
               decision=None, variation=False, keep_draws=False):
    """The grid of u over stages and scopes for one event.

    treated, control   boolean masks or index arrays over the study's cells
    treated_cluster    labels over the study's cells: the clusters the treated side is resampled by (default: by cell)
    reads              the stage the system's map reads (default: the last stage); input_view: the first view, if any
    """
    d = {**DEFAULT_DECISION, **(decision or {})}
    ti, ci_ = (np.where(x)[0] if np.asarray(x).dtype == bool else np.asarray(x) for x in (treated, control))
    if len(ti) < 20 or len(ci_) < 20:
        raise ScopeError(f"too few cells: {len(ti)} treated, {len(ci_)} control")
    tb = np.unique(study.block[ti])
    one_unit = len(tb) == 1
    excl = tuple(tb) if one_unit else ()
    stages = list(stages or study.stages)
    scopes = [s for s in (scopes or study.scopes_available())]
    reads = reads or stages[-1]
    Yt, Yc = study.Y[ti], study.Y[ci_]
    sd_c = Yc.std(axis=0)
    sd_c[sd_c < 1e-8] = 1.0
    cb = study.block[ci_]
    by_block = len(np.unique(cb)) >= 2 and np.median(np.unique(cb, return_counts=True)[1]) > 1 and not all(str(b).startswith("ctl:") for b in np.unique(cb))
    tc = None if treated_cluster is None else np.asarray(treated_cluster).astype(str)[ti]
    by_cluster = tc is not None and len(np.unique(tc)) >= 2
    t, c = _Side(Yt, tc if by_cluster else None), _Side(Yc, cb if by_block else None)
    rng = np.random.default_rng(d["seed"] + 1)
    Ct, Cc = _counts(t.E, d["bootstrap"], rng), _counts(c.E, d["bootstrap"], rng)
    signed = study.R == 1
    # the response itself
    myt, myc = (Ct @ t.sy) / (Ct @ t.n), (Cc @ c.sy) / (Cc @ c.n)
    dY = (myt - myc) / sd_c
    dnorm = dY[:, 0] if signed else np.linalg.norm(dY, axis=1)
    div = ci(dnorm[1:], d["alpha"])
    resolved = bool(div[0] > 0 or div[1] < 0)
    out = {"event": {**(event or {}), "n_treated": int(len(ti)), "n_control": int(len(ci_)), "treated_blocks": [str(b) for b in tb],
                     "one_unit": bool(one_unit), "resampled": {"treated": f"{t.E} clusters" if by_cluster else "cells",
                                                              "control": f"{c.E} blocks" if by_block else "cells"},
                     "approximate": bool((by_cluster and t.E < d["min_clusters"]) or (by_block and c.E < d["min_clusters"]))},
           "readouts": int(study.R), "decision": d, "reads": reads, "input_view": input_view,
           "delta": {"value": float(dnorm[0]), "ci": div, "resolved": resolved, "unit": "control SD"},
           "targets": {}}
    draws = {}

    def cell(u, per=None, extra=None):
        return {"u": float(u[0]), "ci": ci(u[1:], d["alpha"]), **({"per_readout": [float(x) for x in per]} if per is not None and len(per) > 1 else {}),
                **(extra or {})}

    def one(scope, stage, cls):
        if scope == "frozen":
            return study.frozen[ti], study.frozen[ci_]
        return (study.predict(scope, stage, cls, ti, excl).astype(np.float64), study.predict(scope, stage, cls, ci_, excl).astype(np.float64))

    targets = ["shift"] + (["variation"] if variation and one_unit else [])
    for tg in targets:
        out["targets"][tg] = {"stages": {}}
    for scope in scopes:
        for stage in (["system"] if scope == "frozen" else stages):
            by = {tg: {} for tg in targets}
            for cls in (["system"] if scope == "frozen" else list(study.classes)):
                Pt, Pc = one(scope, stage, cls)
                u, per, a = _u_shift(t, c, Pt, Pc, Ct, Cc, sd_c)
                by["shift"][cls] = cell(u, per, {"slope": [float(x) for x in np.atleast_1d(a)][:8]})
                draws[("shift", stage, scope, cls)] = u
                if "variation" in by:
                    uv, perv = _u_variation(t, Pt, Ct)
                    by["variation"][cls] = cell(uv, perv)
                    draws[("variation", stage, scope, cls)] = uv
            for tg in targets:
                sg = signed and tg == "shift"
                v, best = _union(by[tg], d["margin_u"], sg)
                out["targets"][tg]["stages"].setdefault(stage, {})[scope] = {"u": by[tg][best]["u"], "ci": by[tg][best]["ci"], "verdict": v,
                                                                              "class": best, "by_class": by[tg]}
    # gaps between consecutive scopes at the stage the system reads, paired
    for tg in targets:
        T = out["targets"][tg]
        at = lambda s: ("system" if s == "frozen" else reads)
        have = [s for s in ORDER if s in T["stages"].get(at(s), {})]
        gaps = {}
        for a_, b_ in zip(have[:-1], have[1:]):
            ua = np.abs(draws[(tg, at(a_), a_, T["stages"][at(a_)][a_]["class"])])
            ub = np.abs(draws[(tg, at(b_), b_, T["stages"][at(b_)][b_]["class"])])
            g = ua - ub
            iv = ci(g[1:], d["alpha"])
            gaps[f"{a_}->{b_}"] = {"value": float(g[0]), "ci": iv, "resolved": bool(iv[0] > d["step_margin"]),
                                   "absent": bool(iv[1] < d["step_margin"])}
        if "frozen" in have:
            for b_ in have[1:]:
                g = np.abs(draws[(tg, "system", "frozen", "system")]) - np.abs(draws[(tg, reads, b_, T["stages"][reads][b_]["class"])])
                gaps.setdefault("from_frozen", {})[b_] = {"value": float(g[0]), "ci": ci(g[1:], d["alpha"])}
        if "other" in T["stages"].get(reads, {}) and "units" in T["stages"][reads]:
            g = np.abs(draws[(tg, reads, "other", T["stages"][reads]["other"]["class"])]) - np.abs(draws[(tg, reads, "units", T["stages"][reads]["units"]["class"])])
            gaps["other->units"] = {"value": float(g[0]), "ci": ci(g[1:], d["alpha"])}
        T["gaps"] = gaps
        T["reading"] = read(T, reads, input_view, resolved or tg == "variation")
    if keep_draws:
        out["_draws"] = {"|".join(k): v for k, v in draws.items()}
    return out


# ------------------------------------------------------------------ readings
def read(T, reads, input_view=None, testable=True):
    """What is missing, from the verdicts of one target. Returns {"missing", "detail", "resolved_gaps"}."""
    S = T["stages"]
    v = lambda scope, stage: S.get("system" if scope == "frozen" else stage, {}).get(scope, {}).get("verdict")
    gaps = [k for k, g in T.get("gaps", {}).items() if isinstance(g, dict) and g.get("resolved")]
    res = lambda m, why: {"missing": m, "meaning": MISSING[m], "detail": why, "resolved_gaps": gaps}
    if not testable:
        return res("not testable", "the response is not resolved: there is nothing to explain")
    fz, tr, un, ow = v("frozen", reads), v("train", reads), v("units", reads), v("own", reads)
    if fz == "support":
        return res("none", "the system's own output carries the response")
    if tr == "support":
        if fz == "reject":
            return res("function", "a readout of the system's training cells carries it; the system's own output does not")
        return res("unresolved", "a readout of the system's training cells carries it; the system's own output is "
                                 + ("unresolved" if fz else "not given"))
    if un == "support":
        if tr == "reject":
            return res("coverage", "a readout of the other units of this study carries it; no readout of the system's training cells does")
        return res("unresolved", "a readout of the other units carries it; the system's training cells are "
                                 + ("unresolved" if tr else "not given, so coverage and function cannot be told apart"))
    if ow == "support":
        if un == "reject":
            return res("unit", "only a readout fitted with the unit's own cells carries it")
        return res("unresolved", "own scope supported, units scope unresolved")
    if ow == "reject":
        iv_ = None if input_view is None else S.get(input_view, {}).get("own", {}).get("verdict")
        if iv_ == "support":
            return res("upstream", f"readable from {input_view} with the unit's own cells, not from {reads}")
        if iv_ in ("reject", None):
            return res("unreadable", "no scope reads it" + ("" if input_view is None else f", at {reads} or at {input_view}"))
        return res("unresolved", f"own scope rejected at {reads}, unresolved at {input_view}")
    return res("unresolved", "no scope reaches a verdict at the stage the system reads")


def forecast(result, target="shift"):
    """What each kind of change could remove of |u| at the system's output, read from the gaps. These are forecasts made
    with readouts: a system changed in that way has to be audited to learn what it delivers."""
    g = result["targets"][target].get("gaps", {}).get("from_frozen")
    if not g:
        return {}
    arm = {"train": ("function", "another map, the same training cells"), "units": ("coverage", "cells of this kind in the training"),
           "own": ("calibration", "the unit's own labelled cells")}
    return {arm[k][0]: {"removes": v["value"], "ci": v["ci"], "what": arm[k][1]} for k, v in g.items() if k in arm}


def family(results, target="shift", min_clusters=10, alpha=0.05, seed=0, B=2000):
    """Events of one kind read together (the donors of a study): per stage and scope, the mean u over events with an
    interval from resampling the events, and how the events were read."""
    if not results:
        return {}
    rng = np.random.default_rng(seed)
    n = len(results)
    out = {"n_events": n, "approximate": n < min_clusters, "stages": {}, "gaps": {}, "missing": {}}
    idx = rng.integers(0, n, size=(B, n))
    for stage in results[0]["targets"][target]["stages"]:
        for scope in results[0]["targets"][target]["stages"][stage]:
            u = np.array([abs(r["targets"][target]["stages"][stage][scope]["u"]) for r in results])
            vd = [r["targets"][target]["stages"][stage][scope]["verdict"] for r in results]
            out["stages"].setdefault(stage, {})[scope] = {"mean": float(u.mean()), "ci": ci(u[idx].mean(axis=1), alpha) if n >= 3 else None,
                                                           "min": float(u.min()), "max": float(u.max()),
                                                           "verdicts": {k: vd.count(k) for k in ("support", "unresolved", "reject")}}
    for k in results[0]["targets"][target].get("gaps", {}):
        if k == "from_frozen":
            continue
        g = np.array([r["targets"][target]["gaps"][k]["value"] for r in results])
        out["gaps"][k] = {"mean": float(g.mean()), "ci": ci(g[idx].mean(axis=1), alpha) if n >= 3 else None, "min": float(g.min()),
                          "max": float(g.max()), "resolved_in": int(sum(r["targets"][target]["gaps"][k]["resolved"] for r in results))}
    for r in results:
        m = r["targets"][target]["reading"]["missing"]
        out["missing"][m] = out["missing"].get(m, 0) + 1
    return out


# ------------------------------------------------------------------ a run from a file
def _load(path, keys, base):
    from pathlib import Path
    z = np.load(Path(base) / path, allow_pickle=False)
    ro = keys["readouts"]
    if isinstance(ro, str):
        Y = np.asarray(z[ro], np.float64)
        Y = Y[:, None] if Y.ndim == 1 else Y
        if keys.get("readout_columns") is not None:
            Y = Y[:, list(keys["readout_columns"])]
    else:
        Y = np.column_stack([np.asarray(z[k], np.float64) for k in ro])
    return z, Y, {name: np.asarray(z[k], np.float32) for name, k in keys["stages"].items()}


def run_spec(spec, base="."):
    """A scope run described by a file (see anm_diag/LOOP_PROTOCOL.md, version 2). One format for every kind of event:
    a perturbation, a donor, a batch, a study, a disease. Returns {"name", "kind", "checks", "events", "results", "family"}."""
    from .events import CONTRAST_KINDS, check_contrast, contrast_events
    keys, kind = spec["keys"], spec.get("kind", "perturbation")
    z, Y, stages = _load(spec["arrays"]["audit"], keys, base)
    unit = np.asarray(z[keys["unit"]]).astype(str)
    block = np.asarray(z[keys.get("block", keys["unit"])]).astype(str).astype(object)
    ctl = keys.get("control", "others")
    is_ctl = np.zeros(len(unit), bool) if ctl == "others" else (np.asarray(z[ctl["key"]]).astype(str) == str(ctl["equals"]))
    if is_ctl.any():                                         # controls without a unit of their own get pseudo-blocks
        block[is_ctl] = fold_blocks(int(is_ctl.sum()), int(spec.get("control_folds", 5)), int((spec.get("decision") or {}).get("seed", 0)))
    block = block.astype(str)
    side = {}
    for name in ("train", "other"):
        if spec["arrays"].get(name):
            zz, YY, SS = _load(spec["arrays"][name], keys, base)
            bk = keys.get(f"{name}_block")
            side[name] = {"stages": SS, "Y": YY, "block": None if not bk else np.asarray(zz[bk]).astype(str)}
    frozen = None
    if keys.get("frozen"):
        f = np.asarray(z[keys["frozen"]["array"]], np.float64)
        f = f[:, None] if f.ndim == 1 else f
        frozen = f[:, list(keys["frozen"]["columns"])] if keys["frozen"].get("columns") is not None else f
    final = set(map(str, (spec.get("roles") or {}).get("final", [])))
    if final:                                                # units set aside for a final audit are not read by anything else
        used = (set(np.unique(unit)) | set(np.unique(block))) & final
        for name in side:
            if side[name]["block"] is not None:
                used |= set(np.unique(side[name]["block"])) & final
        if used and not (spec.get("final_audit") and spec.get("hypotheses")):
            raise ScopeError(f"units set aside as final are in this run ({', '.join(sorted(used))}): a run may read them only as a final "
                             f"audit (final_audit: true) that names the hypotheses it tests")
    study = Study(Y, stages, block, train=side.get("train"), other=side.get("other"), frozen=frozen, classes=spec.get("readout_classes"),
                  folds=int(spec.get("control_folds", 5)), seed=int((spec.get("decision") or {}).get("seed", 0)), max_fit=spec.get("max_fit"))
    skip = set(map(str, spec.get("units_exclude", [])))
    cand = [u for u in np.unique(unit[~is_ctl]) if u not in skip]
    ev = spec.get("events", "each")
    wanted = cand if ev in ("each", "pooled") else [str(u) for u in ev]
    cluster = None if not keys.get("cluster") else np.asarray(z[keys["cluster"]]).astype(str)
    todo = [("pool", np.isin(unit, wanted) & ~is_ctl)] if ev == "pooled" else [(u, (unit == u) & ~is_ctl) for u in wanted]
    if spec.get("also_pooled") and ev != "pooled":
        todo.append(("pool", np.isin(unit, wanted) & ~is_ctl))
    results, checks, typed = [], [], []
    for u, t in todo:
        tb = np.unique(block[t])
        c = is_ctl if is_ctl.any() else ~np.isin(block, tb)
        e = {"id": f"{spec.get('name', 'run')}:{kind}:{u}", "kind": kind, "study": spec.get("study"), "unit": u}
        results.append(scope_grid(study, t, c, event=e, reads=spec.get("reads"), input_view=spec.get("input_view"), decision=spec.get("decision"),
                                  treated_cluster=cluster if (cluster is not None and u == "pool") else None, variation=bool(spec.get("variation"))))
        if kind in CONTRAST_KINDS:
            cs = dict(name=spec.get("name", "run"), kind=kind, treated=[u] if u != "pool" else wanted, control="others", study=spec.get("study"),
                      label_basis=spec.get("label_basis"), confounded_with=spec.get("confounded_with"), pooling=spec.get("pooling"),
                      n_control_units=int(len(np.unique(block[c]))))
            typed += contrast_events(cs)
            for f in check_contrast(cs):
                if f not in checks:
                    checks.append(f)
    if kind not in CONTRAST_KINDS and kind != "perturbation":
        checks.append({"code": "C1", "severity": "error", "message": f"unknown kind '{kind}'"})
    single = [r for r in results if r["event"]["unit"] != "pool"]
    fam = {tg: family(single, tg, (spec.get("decision") or {}).get("min_clusters", DEFAULT_DECISION["min_clusters"]))
           for tg in (single[0]["targets"] if single else [])}
    return {"name": spec.get("name", "run"), "kind": kind, "study": spec.get("study"), "checks": checks, "typed_events": typed,
            "results": results, "family": fam}


def main(argv=None):
    import argparse
    import json
    from pathlib import Path

    import yaml

    from .events import scope_events, write_jsonl
    ap = argparse.ArgumentParser(prog="anm_diag scope", description="the grid of u over stages and scopes for the events of a file")
    ap.add_argument("spec")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    spec = yaml.safe_load(Path(a.spec).read_text())
    res = run_spec(spec, Path(a.spec).parent)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for c in res["checks"]:
        print(f"[{c['severity']}] {c['code']}: {c['message']}")
    if any(c["severity"] == "error" for c in res["checks"]):
        raise SystemExit("the declaration has errors; nothing was written")
    events = list(res["typed_events"])
    for r in res["results"]:
        events += scope_events(r, f"{res['name']}:{r['event']['unit']}", context=res.get("study"))
        print(f"{r['event']['id']}: " + "; ".join(f"{tg}: {T['reading']['missing']}" for tg, T in r["targets"].items()))
    (out / "scope.json").write_text(json.dumps({k: v for k, v in res.items() if k != "typed_events"}, indent=1, default=float))
    write_jsonl(events, out / "scope_events.jsonl")
    return 0
