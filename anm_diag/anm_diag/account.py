"""The account of a change (GENERAL.md, Section 8). Nothing here names a model, a modality, a readout or a stage.

An event changes a readout. Which part of the system holds that change? A stage is any array a head could read: a view
of the model's input, a step inside the model, the state its head reads. For one event, one readout and one stage, three
questions, each answered by its own kind of evidence, on halves of the cells (fitted on one, read on the other):

  read      what a readout fitted on control cells reads of the change.                       Needs no record of the event.
  carried   what the stage holds of it. Among cells whose stage is alike, the event's cells and the control cells still
            differ by some amount; `carried` is one minus that amount over the whole change. It is the partial effect of
            the arm given the stage, by residuals on residuals (the readout given the stage, the arm given the stage).
            The evidence is inside the arms. It never fits the difference between the two means.   Needs the event's cells.
  others    what a readout fitted on the cells of other events forecasts for this one.         Needs the other events only.

and one number that is NOT evidence that a stage holds a change:

  found_again   what a readout made to follow the event's mean change on one half reads on the other half. It is near
                one whenever the event's shift of the stage can be found again, whatever the stage holds of the change
                (toy_account.py plants an event whose change no stage holds; this number calls it closed).

`chain` puts the readings of several stages in the order of the model and splits the change into four parts:

  not in the input       1 - carried(input)                  no stage holds it cell by cell (see below)
  lost before the state  carried(input) - carried(state)     the model, between its input and the state its head reads
  held, not read         carried(state) - read(state)        the head: a readout with the records of events reads it
  read                   read(state)                         nothing to do

These parts are about single cells: what a head could read from a cell. A readout that is measured in another way than
the input (a protein beside RNA) is often not in a cell's input in full, and its mean change can still follow from the
mean shift of a stage by a law that holds across events. That is the other level, `others` for a few events and `Law`
for many: what the system can forecast of an event's mean change without its record. A part that is not in the input
needs an edge of the event graph only if the other events do not forecast the change either.

`Law` does the same for many events at once from their rows (the mean shift of a stage, the mean change of every
readout): a readout with a prior from the control cells, the rows of the other events as evidence, and for what no row
through the stage explains a direct part whose size depends on the edge between the event and the readout. Every event
is forecast with itself left out; an event whose direct part stands out and has no edge is a candidate for a missing edge.
"""
import numpy as np

from . import cycle

STEP, MARGIN, ALPHA = cycle.STEP, cycle.MARGIN, cycle.ALPHA
DRAWS, MIN_SNR, MIN_OVERLAP = 300, 4.0, 0.10
LAMBDAS = 10.0 ** np.arange(-3, 3.01, 0.5)
interval = lambda d: [float(v) for v in np.percentile(d, [100 * ALPHA / 2, 100 * (1 - ALPHA / 2)])]
step = lambda iv: "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")


def one_scale(x, C):
    """The spread of every column among the control cells of the fit (a column that does not vary is left as it is)."""
    sd = x[C].std(0)
    return np.where(sd > 1e-12, sd, 1.0)


def ridge(F, t, free=0, seed=0, folds=5):
    """A ridge readout of centred columns F. The first `free` columns are not penalised. The weight is chosen by folds."""
    n = len(t)
    F = F.copy()
    F[:, :free] *= 1e3
    fold = np.random.default_rng(seed).permutation(n) % folds
    G, b, err = F.T @ F, F.T @ t, np.zeros(len(LAMBDAS))
    for k in range(folds):
        m = fold == k
        w, V = np.linalg.eigh((G - F[m].T @ F[m]) / (n - m.sum()))
        c, Fv = V.T @ ((b - F[m].T @ t[m]) / (n - m.sum())), F[m] @ V
        for i, lam in enumerate(LAMBDAS):
            err[i] += ((t[m] - Fv @ (c / (w + lam))) ** 2).sum()
    w, V = np.linalg.eigh(G / n)
    th = V @ ((V.T @ (b / n)) / (w + LAMBDAS[int(err.argmin())]))
    th[:free] *= 1e3
    return th


def _stage(xs, f):
    return xs[f] if isinstance(xs, (tuple, list)) else xs


_KEPT = {}


def _control_fit(x, y, C, seed=0):
    """The ridge readout of control cells of centred columns, with what it leaves uncertain: (readout, Mi, cv / n, root).
    The same control cells, stage and readout give the same fit for every event read against them: kept by content."""
    key = (x.shape[1], len(C), float(x[C].sum()), hash(np.ascontiguousarray(x[C[:16]]).tobytes()), hash(y[C].tobytes()), seed)
    if key not in _KEPT:
        if len(_KEPT) >= 12:
            _KEPT.pop(next(iter(_KEPT)))
        _KEPT[key] = _control_fit_(x, y, C, seed)
    _KEPT[key] = _KEPT.pop(key)                                 # the one used last is dropped last
    return _KEPT[key]


def _control_fit_(x, y, C, seed=0):
    Xc, yc, n = x[C] - x[C].mean(0), y[C] - y[C].mean(), len(C)
    G, b = Xc.T @ Xc, Xc.T @ yc
    fold, err = np.random.default_rng(seed).permutation(n) % 5, np.zeros(len(LAMBDAS))
    for k in range(5):
        m = fold == k
        w, V = np.linalg.eigh((G - Xc[m].T @ Xc[m]) / (n - m.sum()))
        c, Fv = V.T @ ((b - Xc[m].T @ yc[m]) / (n - m.sum())), Xc[m] @ V
        for i, lam in enumerate(LAMBDAS):
            err[i] += ((yc[m] - Fv @ (c / (w + lam))) ** 2).sum()
    w, V = np.linalg.eigh(G / n)
    lam, cv = LAMBDAS[int(err.argmin())], float(err.min() / n)
    return V @ ((V.T @ (b / n)) / (w + lam)), (V / (w + lam)) @ V.T, cv / n, (V / np.sqrt(w + lam))


def _closer(Mi, D):
    """For rows with shifts D (events by columns): the matrix A such that readout + A (rows - D readout) is the smallest
    change of a readout, in the measure of what the control cells leave uncertain, that closes the rows."""
    D = np.atleast_2d(D)
    return Mi @ D.T @ np.linalg.pinv(D @ Mi @ D.T), D


def reading(xs, y, unit, first, event, others=(), control="control", o=None, draws=DRAWS, seed=0):
    """The account of one event's change of one readout at one stage.

    xs: the stage, cells by columns; or (the stage for the half fitted on `first`, the stage for the other half) when the
    stage itself was chosen on the fitted half. unit: the event of every cell. first: the half of every cell. others:
    the other events whose rows may teach a readout. o: the system's own output as trained, if it has one.

    Returns the change, how far the event's shift of the stage stands out, and each share of the change with its
    interval. `overlap`: the share of the difference between the arms that the stage does not tell (one: the stage does
    not tell the event's cells from control cells; zero: it tells them apart completely). Below MIN_OVERLAP no cell of
    one arm is like a cell of the other, `carried` cannot be read from cells, and the stage is taken to hold the change
    as the identity of the event (`carried_by`: "identity"): a head with the record reads it; whether it is more than
    the identity only other events can say.
    `_draws` holds the resamples (row 0 is the value). A resample draws the cells that are read, the control readout
    from what its fit leaves uncertain, and the rows from their errors. The same seed gives the same resamples at every
    stage, so differences between stages are paired."""
    rng = np.random.default_rng(seed + 7)
    others = [h for h in others if h != event]
    keys = ["read", "carried", "found_again"] + (["others"] if others else []) + (["as_trained"] if o is not None else [])
    num, den = {k: np.zeros(draws + 1) for k in keys if k != "carried"}, np.zeros(draws + 1)
    sxy, sxx, see = np.zeros(draws + 1), np.zeros(draws + 1), 0.0
    vis = []
    for f, fit in ((0, first), (1, ~first)):
        x = _stage(xs, f)
        Cf, Tf = np.where(fit & (unit == control))[0], np.where(fit & (unit == event))[0]
        Cs, Ts = np.where(~fit & (unit == control))[0], np.where(~fit & (unit == event))[0]
        idx = [(np.arange(len(Ts)), np.arange(len(Cs)))] + [(rng.integers(0, len(Ts), len(Ts)), rng.integers(0, len(Cs), len(Cs))) for _ in range(draws)]
        zth = np.vstack([np.zeros(x.shape[1]), rng.normal(size=(draws, x.shape[1]))])        # for the control readout: what its fit leaves uncertain
        zrow = np.vstack([np.zeros(1 + len(others)), rng.normal(size=(draws, 1 + len(others)))])
        x = (x - x[Cf].mean(0)) / one_scale(x, Cf)
        xc, yc = x[Cf].mean(0), y[Cf].mean()
        th_c, Mi, k, root = _control_fit(x, y, Cf, seed)
        TH = th_c + np.sqrt(k) * zth @ root.T                                               # a readout of control cells for every resample
        shift = lambda a: x[a].mean(0) - xc
        row = lambda a: (y[a].mean() - yc, np.sqrt(y[a].var(ddof=1) / len(a) + y[Cf].var(ddof=1) / len(Cf)))
        Df, (dyf, sef) = shift(Tf), row(Tf)
        arms = [np.where(fit & (unit == h))[0] for h in others]
        keep = [j for j, a in enumerate(arms) if len(a) > 1]
        ro = [row(arms[j]) for j in keep]
        A_own, D_own = _closer(Mi, Df)
        A_oth, D_oth = _closer(Mi, np.array([shift(arms[j]) for j in keep])) if keep else (None, None)
        r_oth, s_oth = np.array([r[0] for r in ro]), np.array([r[1] for r in ro])
        F, R = np.r_[Cf, Tf], np.r_[Cs, Ts]
        arm = (unit == event).astype(float)
        mu = x[F].mean(0)
        e_res = arm[R] - arm[F].mean() - (x[R] - mu) @ ridge(x[F] - mu, arm[F] - arm[F].mean(), seed=seed)
        y_res = y[R] - y[F].mean() - (x[R] - mu) @ ridge(x[F] - mu, y[F] - y[F].mean(), seed=seed)
        ey, ee = e_res * y_res, e_res * e_res
        see += float(((arm[R] - arm[R].mean()) ** 2).sum())
        nC = len(Cs)
        xT, xC = x[Ts], x[Cs]
        a_tr = np.polyfit(o[Cf], y[Cf], 1)[0] if o is not None else None                     # the output as trained, on the scale of the readout among control cells
        for b, (it, ic) in enumerate(idx):
            Ds = xT[it].mean(0) - xC[ic].mean(0)
            den[b] += y[Ts[it]].mean() - y[Cs[ic]].mean()
            th = TH[b]
            num["read"][b] += Ds @ th
            num["found_again"][b] += Ds @ (th + A_own @ (dyf + sef * zrow[b, :1] - D_own @ th))
            if others:
                num["others"][b] += Ds @ (th + A_oth @ (r_oth + s_oth * zrow[b, 1:1 + len(keep)] - D_oth @ th) if keep else th)
            if o is not None:
                num["as_trained"][b] += a_tr * (o[Ts[it]].mean() - o[Cs[ic]].mean())
            sxy[b] += ey[nC + it].sum() + ey[ic].sum()
            sxx[b] += ee[nC + it].sum() + ee[ic].sum()
        Ds = xT.mean(0) - xC.mean(0)
        vis.append((float(Df @ Ds / (np.linalg.norm(Df) * np.linalg.norm(Ds) + 1e-300)), float((Df @ Df) / (x.shape[1] * (1 / len(Tf) + 1 / len(Cf))))))
    D = {k: num[k] / den for k in num}
    overlap = float(sxx[0] / see)
    by = "cells" if overlap >= MIN_OVERLAP else "identity"
    D["carried"] = 1 - (sxy / sxx) / (den / 2) if by == "cells" else np.ones(draws + 1)
    T, C = unit == event, unit == control
    change = float(y[T].mean() - y[C].mean())
    se = float(np.sqrt(y[T].var(ddof=1) / T.sum() + y[C].var(ddof=1) / C.sum()))
    return {"event": str(event), "change": change, "snr": abs(change) / se, "cells": {"event": int(T.sum()), "control": int(C.sum())},
            "shift_of_the_stage": {"cosine_between_the_halves": float(np.mean([v[0] for v in vis])), "times_its_noise": float(np.mean([v[1] for v in vis]))},
            "overlap": overlap, "carried_by": by,
            "shares": {k: {"share": float(D[k][0]), "ci": interval(D[k][1:])} for k in keys}, "_draws": D}


def gain(a, b, key="carried"):
    """What the second reading holds that the first does not: the difference of a share on paired resamples."""
    d = b["_draws"][key] - a["_draws"][key]
    iv = interval(d[1:])
    return {"gain": float(d[0]), "ci": iv, "reading": step(iv)}


PARTS = ("not_in_the_input", "lost_before_the_state", "held_not_read", "read")
LEVER = {"not_in_the_input": "no stage holds this part cell by cell, so no head can read it from a cell. What the event graph must supply is read across events: "
                             "the records of other events if they forecast the change through the stage, an edge if they do not",
         "lost_before_the_state": "the model, between its input and the state its head reads",
         "held_not_read": "the head: a readout with the cells of recorded events",
         "read": "nothing: a readout of control cells reads this part"}
FORECAST, NOT_FORECAST = 0.75, 0.50                           # across events: the other events forecast a change if the interval of `others` lies above the first, and do not if it lies below the second


def across_events(r):
    """What the rows of the other events say of this event's change through a stage (None if there are none)."""
    if "others" not in r["shares"]:
        return None
    o = r["shares"]["others"]
    says = "the other events forecast it" if o["ci"][0] >= FORECAST else ("the other events do not forecast it" if o["ci"][1] < NOT_FORECAST else "unresolved")
    return {"others": o["share"], "ci": o["ci"], "reading": says}


def says(ch):
    """One line for an account: the open parts, largest first, then the unresolved ones; and for a part that no stage
    holds cell by cell, what the other events say."""
    if ch["snr"] < MIN_SNR:
        return "no change to account for (the change does not stand out from its error)"
    ae = ch.get("across_events")
    def name(p, is_open=False):
        t = f"{p['part'].replace('_', ' ')} {p['share']:.0%}"
        if is_open and p["part"] == "not_in_the_input" and ae:
            t += f" [across events: {ae['reading']}, {ae['others']:.0%}]"
        return t
    if ch["open"]:
        return "; ".join(name(p, True) for p in ch["open"]) + ("" if not ch["unresolved"] else " (unresolved: " + "; ".join(name(p) for p in ch["unresolved"]) + ")")
    if ch["unresolved"]:
        return "unresolved: " + "; ".join(name(p) for p in ch["unresolved"])
    return "read: a readout of control cells reads the change, and no other part is as large as a step"


def chain(by_stage, inputs, state):
    """The four parts of a change, cell by cell, from the readings of its stages. inputs: the names of the views of the
    model's input; state: the stage the head reads. The input is taken to hold what the view that carries most holds,
    and at least what the state holds (if the state carries more than every view, nothing is lost before it and
    `input` names the state). Every part has the interval of its own resamples. `open`: the parts resolved as more than a step, largest first, each with the place it names.
    `unresolved`: the parts whose interval holds a step. `across_events`: what the rows of the other events forecast of
    the change through the state; a part that is not in the input cell by cell needs an edge only if they do not."""
    inp = max([*inputs, state], key=lambda s: by_stage[s]["shares"]["carried"]["share"])     # a view of the input is only a part of it: the input holds at least what the state holds
    ci_, cs, rd = by_stage[inp]["_draws"]["carried"], by_stage[state]["_draws"]["carried"], by_stage[state]["_draws"]["read"]
    d = {"not_in_the_input": 1 - ci_, "lost_before_the_state": ci_ - cs, "held_not_read": cs - rd, "read": rd}
    parts = {k: {"share": float(v[0]), "ci": interval(v[1:]), "reading": step(interval(v[1:]))} for k, v in d.items()}
    r = by_stage[state]
    out = {"event": r["event"], "change": r["change"], "snr": r["snr"], "input": inp, "state": state, "parts": parts,
           "carried_by": {inp: by_stage[inp]["carried_by"], state: by_stage[state]["carried_by"]}, "across_events": across_events(r), "open": [], "unresolved": []}
    if r["snr"] >= MIN_SNR:
        for key, want in (("open", "a step"), ("unresolved", "unresolved")):
            out[key] = [{"part": k, "share": parts[k]["share"], "ci": parts[k]["ci"], "place": LEVER[k]} for k in sorted(PARTS[:3], key=lambda k: -parts[k]["share"]) if parts[k]["reading"] == want]
    out["says"] = says(out)
    return out


def must_fail(xs, y, unit, first, event, control="control", draws=DRAWS, seed=0):
    """Two controls of the measure. A stage with its cells permuted holds nothing: it must carry nothing. An event made
    of control cells has no change: it must have no change to account for."""
    rng = np.random.default_rng(seed + 11)
    perm = rng.permutation(len(y))
    xp = tuple(x[perm] for x in xs) if isinstance(xs, (tuple, list)) else xs[perm]
    p = reading(xp, y, unit, first, event, control=control, draws=draws, seed=seed)
    C = np.where(unit == control)[0]
    fake = rng.choice(C, min(int((unit == event).sum()), len(C) // 3), replace=False)
    un = unit.astype(object).copy()
    un[fake] = "made-up"
    keep = (un == control) | (un == "made-up")
    sub = (lambda x: x[keep])
    m = reading(tuple(sub(x) for x in xs) if isinstance(xs, (tuple, list)) else sub(xs), y[keep], un[keep].astype(str), first[keep], "made-up", control=control, draws=draws, seed=seed)
    return {"permuted_stage": {k: p["shares"][k] for k in ("read", "carried", "found_again")},
            "made_up_event": {"change": m["change"], "snr": m["snr"], "no_change_to_account_for": bool(m["snr"] < MIN_SNR)}}


# --- many events at once: the law of their rows, and the edges ---------------------------------------------------------
class Law:
    """Every event of a library and every readout, from rows. In spreads of a readout among control cells:

        change(event, readout) = shift(event) . readout  +  direct(event, readout)  +  error

    readout: a prior from the control cells (their ridge readout and its uncertainty, widened by one factor that the rows
    choose), updated by the rows of the other events. direct: zero with a small spread and heavy tails for a pair with no
    edge; for a pair with an edge of class k a mean and a spread of its own, learned from the other pairs of that class.
    edges: {(event, readout): class}. A class with fewer than `min_class` recorded pairs gives no forecast of its size."""

    def __init__(self, x, Y, unit, control="control", edges=None, nu=4.0, min_cells=30, min_class=3, seed=0):
        C = np.where(unit == control)[0]
        X = (x - x[C].mean(0)) / one_scale(x, C)
        self.events = [str(e) for e in sorted(set(unit.tolist())) if e != control and (unit == e).sum() >= min_cells]
        self.readouts, self.nu, self.min_class = list(Y), nu, min_class
        self.n = np.array([(unit == e).sum() for e in self.events])
        self.D = np.array([X[unit == e].mean(0) for e in self.events])
        nc, p = len(C), X.shape[1]
        self.times_noise = (self.D ** 2).sum(1) / (p * (1 / self.n + 1 / nc))
        Xc = X[C]
        G = Xc.T @ Xc / nc
        w, V = np.linalg.eigh(G)
        fold = np.random.default_rng(seed).permutation(nc) % 5
        self.r = {}
        for name, y in Y.items():
            sd, yc = y[C].std(), y[C] - y[C].mean()
            err = np.zeros(len(LAMBDAS))
            for k in range(5):
                m = fold == k
                wk, Vk = np.linalg.eigh((G * nc - Xc[m].T @ Xc[m]) / (nc - m.sum()))
                c, Fv = Vk.T @ (Xc[~m].T @ yc[~m] / (nc - m.sum())), Xc[m] @ Vk
                for i, lam in enumerate(LAMBDAS):
                    err[i] += ((yc[m] - Fv @ (c / (wk + lam))) ** 2).sum()
            lam, cv = LAMBDAS[int(err.argmin())], float(err.min() / nc)
            th = V @ ((V.T @ (Xc.T @ yc / nc)) / (w + lam))
            DV = self.D @ V
            K0 = (cv / nc) * (DV / (w + lam)) @ DV.T / sd ** 2
            dy = np.array([y[unit == e].mean() - y[C].mean() for e in self.events])
            se2 = np.array([y[unit == e].var(ddof=1) / (unit == e).sum() + y[C].var(ddof=1) / nc for e in self.events])
            noise = (se2 + float(th @ G @ th) * (1 / self.n + 1 / nc)) / sd ** 2          # the row's own error, and that of the event's shift of the stage
            cls = np.array([(edges or {}).get((e, name)) for e in self.events], dtype=object)
            self.r[name] = {"sd": float(sd), "cv": cv, "lam": float(lam), "K0": K0, "change": dy / sd, "se": np.sqrt(se2) / sd, "noise": noise,
                            "read": self.D @ th / sd, "cls": cls, "tau2": 1.0, "w": np.ones(len(dy))}
        self.classes = sorted({c for q in self.r.values() for c in q["cls"] if c is not None})
        self.sigma0, self.cl = None, {}

    def _var_direct(self, q):
        s2 = np.array([self.sigma0 ** 2 if c is None else self.cl[c]["sigma"] ** 2 for c in q["cls"]])
        return s2 / q["w"]

    def _mean_direct(self, q, without=None):
        return np.array([0.0 if c is None else self.cl[c]["mu"] for c in q["cls"]])

    def fit(self, iters=40):
        res = np.concatenate([q["change"] - q["read"] for q in self.r.values()])
        cls = np.concatenate([q["cls"] for q in self.r.values()])
        none = np.array([c is None for c in cls])
        self.sigma0 = float(1.4826 * np.median(np.abs(res[none] - np.median(res[none])))) + 1e-6
        edge_res = res[~none]
        pooled = (float(edge_res.mean()), float(edge_res.std(ddof=1))) if len(edge_res) >= self.min_class else (0.0, float("nan"))
        for c in self.classes:
            v = res[cls == c]
            self.cl[c] = {"n": int(len(v)), "mu": float(v.mean()) if len(v) >= self.min_class else pooled[0],
                          "sigma": float(v.std(ddof=1)) if len(v) >= self.min_class else pooled[1], "own": bool(len(v) >= self.min_class)}
        taus = 10.0 ** np.arange(0, 4.01, 0.1)
        for _ in range(iters):
            Ed2, Ed, W = [], [], []
            for q in self.r.values():
                Vd, rho = self._var_direct(q), q["change"] - q["read"] - self._mean_direct(q)
                Vd = np.where(np.isfinite(Vd), Vd, 1e6)
                best = None
                for t2 in taus:
                    L = np.linalg.cholesky(t2 * q["K0"] + np.diag(Vd + q["noise"]))
                    a = np.linalg.solve(L, rho)
                    nll = 0.5 * a @ a + np.log(np.diag(L)).sum()
                    if best is None or nll < best[0]:
                        best = (nll, t2)
                q["tau2"] = float(best[1])
                Si = np.linalg.inv(q["tau2"] * q["K0"] + np.diag(Vd + q["noise"]))
                m = Vd * (Si @ rho)                                   # the direct part beyond its class mean, given every row
                v = Vd - Vd * np.diag(Si) * Vd
                q["direct_dev"], q["direct_var"] = m, v
                Ed.append(m + self._mean_direct(q)); Ed2.append(m * m + v); W.append(q["w"])
            Ed, Ed2, W = np.concatenate(Ed), np.concatenate(Ed2), np.concatenate(W)
            s2 = np.array([self.sigma0 ** 2 if c is None else self.cl[c]["sigma"] ** 2 for c in cls])
            s2 = np.where(np.isfinite(s2), s2, 1e6)
            Wn = (self.nu + 1) / (self.nu + Ed2 / s2)
            self.sigma0 = float(np.sqrt((Wn[none] * Ed2[none]).mean())) + 1e-9
            for c in self.classes:
                m = cls == c
                if self.cl[c]["own"]:
                    self.cl[c]["mu"] = float((Wn[m] * Ed[m]).sum() / Wn[m].sum())
                    self.cl[c]["sigma"] = float(np.sqrt(max((Wn[m] * Ed2[m]).mean(), 1e-12)))
            k = 0
            for q in self.r.values():
                q["w"] = Wn[k:k + len(q["w"])]
                k += len(q["w"])
        return self

    def forecast(self):
        """Every pair with its event left out: the part a readout of control cells reads, what the rows of the other events
        add through the stage, the part by its edge, the forecast with its spread, and how far the measured change is."""
        rows = []
        dev = {}                                                    # the direct parts of the pairs of a class, as fitted: (event, value)
        for name, q in self.r.items():
            for i, c in enumerate(q["cls"]):
                if c is not None:
                    dev.setdefault(c, []).append((self.events[i], float(q["direct_dev"][i] + self.cl[c]["mu"])))
        t95 = 1.96
        for name, q in self.r.items():
            Vd = np.where(np.isfinite(self._var_direct(q)), self._var_direct(q), 1e6)
            V = Vd + q["noise"]
            rho = q["change"] - q["read"] - self._mean_direct(q)
            Si = np.linalg.inv(q["tau2"] * q["K0"] + np.diag(V))
            m_loo = rho - (Si @ rho) / np.diag(Si)                  # what the rows of the other events add through the stage
            v_loo = np.maximum(1 / np.diag(Si) - V, 0.0)
            for i, e in enumerate(self.events):
                c = q["cls"][i]
                if c is None:
                    edge, sd_d, why = 0.0, self.sigma0, None
                else:
                    oth = [v for h, v in dev[c] if h != e]
                    if len(oth) >= self.min_class - 1 and len(oth) >= 2:
                        edge, sd_d, why = float(np.mean(oth)), float(np.std(oth, ddof=1) * np.sqrt(1 + 1 / len(oth))), None
                    else:
                        edge, sd_d, why = float("nan"), float("nan"), f"an edge of class '{c}' with {len(oth)} other recorded pairs: its size cannot be forecast"
                f = q["read"][i] + m_loo[i] + (edge if edge == edge else 0.0)
                sd = float(np.sqrt(v_loo[i] + (sd_d ** 2 if sd_d == sd_d else 0.0) + q["noise"][i]))
                z = float((q["change"][i] - f) / sd) if why is None else float("nan")
                rows.append({"event": e, "readout": name, "cells": int(self.n[i]), "edge": c, "change": float(q["change"][i]), "se": float(q["se"][i]),
                             "shift_times_noise": float(self.times_noise[i]), "read": float(q["read"][i]), "from_the_other_events": float(m_loo[i]),
                             "by_its_edge": edge, "forecast": float(f), "sd": sd, "interval": [float(f - t95 * sd), float(f + t95 * sd)], "z": z, "no_forecast": why,
                             "direct": float(q["change"][i] - q["read"][i] - m_loo[i]), "direct_sd": float(np.sqrt(v_loo[i] + q["noise"][i]))})
        return rows

    def calibration(self, rows=None):
        """How the forecasts did: the share of measured changes inside the 95% interval, the mean squared z, and the error
        against two forecasts that use less: no change, and a readout of control cells alone."""
        rows = rows or self.forecast()
        ok = [r for r in rows if r["no_forecast"] is None]
        z = np.array([r["z"] for r in ok])
        rm = lambda v: float(np.sqrt(np.mean(np.square(v))))
        out = {"pairs": len(ok), "inside_95": float(np.mean(np.abs(z) < 1.96)), "mean_z2": float(np.mean(z ** 2)), "median_abs_z": float(np.median(np.abs(z))),
               "error": {"no_change": rm([r["change"] for r in ok]), "control_readout": rm([r["change"] - r["read"] for r in ok]), "forecast": rm([r["change"] - r["forecast"] for r in ok])}}
        big = [r for r in ok if abs(r["change"]) / r["se"] >= MIN_SNR]
        out["pairs_with_a_change"] = {"pairs": len(big), "inside_95": float(np.mean([abs(r["z"]) < 1.96 for r in big])) if big else None,
                                      "error": {"no_change": rm([r["change"] for r in big]), "control_readout": rm([r["change"] - r["read"] for r in big]),
                                                "forecast": rm([r["change"] - r["forecast"] for r in big])} if big else None}
        for c in self.classes:
            sel = [r for r in ok if r["edge"] == c]
            out.setdefault("by_class", {})[c] = {"pairs": len(sel), "inside_95": float(np.mean([abs(r["z"]) < 1.96 for r in sel])) if sel else None,
                                                 "error_forecast": rm([r["change"] - r["forecast"] for r in sel]) if sel else None,
                                                 "error_without_edges": rm([r["change"] - r["read"] - r["from_the_other_events"] for r in sel]) if sel else None}
        return out

    def flags(self, rows=None, q=0.05):
        """Pairs with no edge whose direct part stands out (false discovery rate q over all such pairs). An event whose
        shift of the stage stands out: a direction the law does not cover, and its record would teach the readout. An
        event that does not move the stage: nothing a stage holds explains it, so the graph lacks an edge."""
        from math import erfc, sqrt
        rows = rows or self.forecast()
        cand = [r for r in rows if r["edge"] is None]
        pv = np.array([erfc(abs(r["direct"]) / np.hypot(r["direct_sd"], self.sigma0) / sqrt(2)) for r in cand])
        order = np.argsort(pv)
        thr = q * (np.arange(len(pv)) + 1) / len(pv)
        passed = np.where(pv[order] <= thr)[0]
        keep = order[: passed.max() + 1] if len(passed) else []
        out = []
        for i in keep:
            r = cand[i]
            moved = r["shift_times_noise"] >= 4.0
            out.append({"event": r["event"], "readout": r["readout"], "direct": r["direct"], "errors": r["direct"] / np.hypot(r["direct_sd"], self.sigma0), "shift_times_noise": r["shift_times_noise"],
                        "candidate": "a direction of the stage the law does not cover: its record would teach the readout" if moved else "a missing edge: the event does not move the stage"})
        return sorted(out, key=lambda r: -abs(r["errors"]))

    def summary(self):
        return {"events": len(self.events), "readouts": self.readouts, "spread_of_a_direct_part_without_an_edge": self.sigma0,
                "widening_of_the_control_prior": {k: q["tau2"] for k, q in self.r.items()}, "classes": self.cl,
                "events_whose_shift_stands_out": {"beyond_2x_noise": int((self.times_noise > 2).sum()), "beyond_4x_noise": int((self.times_noise > 4).sum())}}
