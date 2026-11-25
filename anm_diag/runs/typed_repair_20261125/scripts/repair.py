#!/usr/bin/env python3
"""The repair of this run: `Full`, the rule that derives a row, and the halves. DESIGN.md, Section 3.

`Full` is anm_diag.cycle.LinearRows with two changes. Nothing in the package is changed.

 1. One scale. Every column of the stage is divided by its spread among the control cells of the fit. The readout of
    LinearRows centres the columns and penalises all of them alike, so a block of small columns was not used
    (runs/typed_state_20261125).
 2. The rows are applied. LinearRows updates the readout of the control cells by a gain K = M D' (D M D' + S)^-1,
    with S the error of the row against the error of the control fit. That gain trusts the control fit: for IRF1 it
    closed 27% of the row on the very cells it was fitted on. Here the uncertainty of the control fit is widened by one
    factor t2 >= 1 until the rows are no longer surprising:

        t2 = max(1, (r'r - trace(R)) / trace(P)),      P = (cv / n) D M D',   R = diag(scale^2 se^2),   r = the residuals of the rows
        readout = controls + (cv / n) t2 M D' (t2 P + R)^-1 r

    With rows that the control fit already explains, t2 = 1 and this is the update of LinearRows. With one row far
    from it, the share of the row that is closed is about 1 - (error of the row)^2 / (residual)^2.

Arms: controls (no record), rows (the rule above), exact (the rows closed on the fitted cells, R -> 0), kalman (the
weight of LinearRows on the scaled stage: t2 = 1), cells (fitted on every cell given).
"""
import sys
from pathlib import Path

import numpy as np

RUNS = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
from anm_diag import cycle  # noqa: E402
from _common import SEED, unit  # noqa: E402,F401

STEP = cycle.STEP
step = lambda iv: "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")


def one_scale(x, C):
    sd = x[C].std(0)
    return np.where(sd > 1e-12, sd, 1.0)


_KEPT = {}


class Full(cycle.LinearRows):
    """derived: {event: (change, squared standard error)}: the change of an event that has no protein record."""

    def __init__(self, x, y, o, C, events, other=None, derived=None, cells=True):
        self.xc, self.sd = x[C].mean(0), one_scale(x, C)
        F = self.feat(x)
        scale = o[C].std() / y[C].std()
        t = o[C].mean() + scale * (y - y[C].mean())
        pen = np.r_[np.ones(x.shape[1]), 0.0]
        key = (x.shape[1], len(C), float(x[C].sum()), hash(np.ascontiguousarray(x[C[:16]]).tobytes()), hash(t[C].tobytes()))    # by content: the same control cells, stage and readout
        if key not in _KEPT:                                    # the readout of the control cells is the same for every event read against them
            lam, cv = cycle._cv_weight(F[C], t[C], pen)
            Mi = np.linalg.inv(F[C].T @ F[C] / len(C) + lam * np.diag(pen))
            if len(_KEPT) >= 10:
                _KEPT.pop(next(iter(_KEPT)))
            _KEPT[key] = (lam, cv, Mi, Mi @ (F[C].T @ t[C] / len(C)))
        _KEPT[key] = _KEPT.pop(key)                             # the one used last is dropped last
        lam, cv, Mi, thc = _KEPT[key]
        self.th = {"controls": thc}
        self.log = {"n_controls": int(len(C)), "ridge_weight": lam, "cv_error": cv, "rows": []}
        vc = y[C].var(ddof=1) / len(C)
        if not events:
            return
        derived = derived or {}
        D = np.array([F[T].mean(0) - F[C].mean(0) for _n, T in events])
        dy = np.array([derived[n][0] if n in derived else y[T].mean() - y[C].mean() for n, T in events])
        se2 = np.array([derived[n][1] if n in derived else y[T].var(ddof=1) / len(T) + vc for n, T in events])
        k = cv / len(C)
        P, R, r = k * (D @ Mi @ D.T), np.diag(scale ** 2 * se2), scale * dy - D @ thc
        t2 = max(1.0, float((r @ r - np.trace(R)) / np.trace(P)))
        up = lambda w, Rn: thc + k * w * (Mi @ D.T) @ np.linalg.solve(w * P + Rn, r)
        self.th["rows"], self.th["kalman"] = up(t2, R), up(1.0, R)
        self.th["exact"] = thc + k * (Mi @ D.T) @ (np.linalg.pinv(P) @ r)
        closed = 1 - np.abs(scale * dy - D @ self.th["rows"]) / np.maximum(np.abs(r), 1e-12)
        self.log.update(widening=t2, rows=[{"event": n, "n_cells": int(len(T)), "change": float(a), "derived": bool(n in derived), "standard_error": float(np.sqrt(s)), "share_of_the_residual_closed": float(c_)}
                                           for (n, T), a, s, c_ in zip(events, dy, se2, closed)])
        if cells and not derived:
            tr = np.r_[C, np.concatenate([T for _n, T in events])]
            lam2, _cv = cycle._cv_weight(F[tr], t[tr], pen)
            self.th["cells"] = np.linalg.solve(F[tr].T @ F[tr] / len(tr) + lam2 * np.diag(pen), F[tr].T @ t[tr] / len(tr))

    def feat(self, x):
        return np.hstack([(x - self.xc) / self.sd, np.ones((len(x), 1))])


def derive(V, y, C, T, recorded):
    """The change of an event that has no protein record, from its cells' RNA and the records of the events of a type.
    All cells are of the fit. V: the input view. C: control cells; T: the cells of the event; recorded: [(event, cells)].
    In the input view at one scale: the event's shift, projected on the mean shift of the recorded events; times their
    mean measured change. No protein value of the event enters. Returns the change, its squared error (that of a measured
    row of as many cells, from the spread of the control cells), and how the rule does on the records, one left out."""
    X = V / one_scale(V, C)
    xc, yc = X[C].mean(0), y[C].mean()
    Ds = {n: X[Tn].mean(0) - xc for n, Tn in recorded}
    dys = {n: float(y[Tn].mean() - yc) for n, Tn in recorded}

    def rule(Dt, names):
        Dm, dym = np.mean([Ds[n] for n in names], 0), float(np.mean([dys[n] for n in names]))
        ratio = float(Dt @ Dm / (Dm @ Dm))
        return ratio * dym, ratio, float(Dt @ Dm / (np.linalg.norm(Dt) * np.linalg.norm(Dm) + 1e-300))

    Dt = X[T].mean(0) - xc
    change, ratio, cos = rule(Dt, list(Ds))
    on_records = {}
    if len(Ds) > 1:
        for n in Ds:
            dv, _r, _c = rule(Ds[n], [m for m in Ds if m != n])
            on_records[n] = {"derived": dv, "measured": dys[n], "derived_over_measured": dv / dys[n] if dys[n] else float("nan")}
    se2 = float(y[C].var(ddof=1) * (1 / len(T) + 1 / len(C)))
    return change, se2, {"change": change, "share_of_the_types_shift": ratio, "cosine_with_the_types_shift": cos, "the_rule_on_the_records": on_records}


def folds(first):
    return ((0, first), (1, ~first))


def cross(xs, y, o, treated, unit_, first, sites, repair=Full, derived=None, cells=True):
    """cycle.cross_fit with a stage that may differ by fold. xs: one array, or (the array of the fold fitted on `first`,
    the array of the other fold): a stage whose types or scale were chosen on the fitted half only.
    derived: None, or a function fold -> {event: (change, squared error)}."""
    n = len(y)
    out, inside, fits = {}, np.full(n, np.nan), []
    for f, fit in folds(first):
        x = xs[f] if isinstance(xs, (tuple, list)) else xs
        C = np.where(fit & ~treated)[0]
        ev = [(s, np.where(fit & treated & (unit_ == s))[0]) for s in sites]
        kw = {} if repair is not Full else {"derived": None if derived is None else derived(f), "cells": cells}
        R = repair(x, y, o, C, [(s, T) for s, T in ev if len(T)], None, **kw)
        fits.append(getattr(R, "log", None))
        for a in R.th:
            out.setdefault(a, np.full(n, np.nan))[~fit] = R.out(a, x[~fit])
        inside[fit] = R.out("rows", x[fit])
    return out, inside, fits


def gain(a, b):
    """|a| - |b| on paired resamples: what the second removes that the first leaves."""
    d = np.abs(a) - np.abs(b)
    iv = cycle.interval(d[1:])
    return {"gain": float(d[0]), "ci": iv, "reading": step(iv), "no_harm": bool(iv[0] > -STEP)}


def u_of(y, p, treated, cluster=None):
    u, slope = cycle.output_u(y, p, treated, cluster, SEED)
    return u, slope


def val(u, slope=None):
    v = cycle.val(u)
    if slope is not None:
        v["slope"] = float(slope)
    return v
