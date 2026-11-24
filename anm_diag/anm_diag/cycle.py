"""One cycle of the loop for any system (GENERAL.md). Nothing in this module names a model, a modality, a readout or a
stage: a system is arrays and names that a run declares.

  select        which events are admitted, on a selection share of their cells (A1, A2; protocol 3d)
  cross_fit     a repair, fitted on one half of the cells and read on the other, for every arm
  read          the reading of a round: A3, every arm against the reference stage, gate 3, the two controls (3a, 3b)
  grid, load    the result of the audit tool, in memory or from its files
  place         where a loss is, across any stages: the model's own input, the stage the model reads, other measured inputs
  target        the result as one target of `outer.open_items` (target_by_event: a type with too few events to pool)
  one_type      the four steps above for one type, from its cells and its grid to its target

  one_units     the same for a type whose events have no control cells: a unit (a donor, a batch) against the other units

The audit itself (the grid of u over stages and scopes) is `scope.run_spec`, which is general already.
The repair is a plug-in. `LinearRows` is the one used so far: a linear readout of a stage, fitted on a type's control
cells and updated once by the rows of its events. `UnitRows` is the same readout where no cell is a control cell. Any
object with the same two methods can take its place.

Two shapes of event have no control cells of their own (GENERAL.md, Section 6):
  a unit against the other units      a donor, a batch, a study: `select(control=None)`, `UnitRows`, `one_units`
  a group of units against another    a disease: patients against healthy donors. The healthy donors' cells are the control
                                      cells; `control_cluster` makes that side count by unit, as the audit tool does
"""
import hashlib

import numpy as np

from . import scope as SC
from .stats import ci as _ci

MARGIN, STEP, ALPHA, DRAWS, MIN_SNR = 0.10, 0.10, 0.05, 1000, 3.0
MIN_POOL = 3                                                 # fewer events than this are read one by one: resampling two events says nothing about a family
interval = lambda d: [float(v) for v in _ci(d, ALPHA)]
verdict = lambda iv: SC._verdict(iv, MARGIN, True)
step_reading = lambda iv: "loss" if iv[0] > STEP else ("no added loss" if iv[1] < STEP else "unresolved")
val = lambda u: {"u": float(u[0]), "ci": interval(u[1:]), "verdict": verdict(interval(u[1:]))}


def diff(a, b):
    """|a| - |b| on paired resamples, read against the step margin."""
    d = np.abs(a) - np.abs(b)
    iv = interval(d[1:])
    return {"added": float(d[0]), "ci": iv, "reading": step_reading(iv)}


# --- selection ---------------------------------------------------------------------------------------------------------
def in_split(cell, seed, of=3):
    """A cell is in the selection share if the SHA-256 of '<seed>|<cell name>' is divisible by `of`."""
    return int(hashlib.sha256(f"{seed}|{cell}".encode()).hexdigest(), 16) % of == 0


def screen(y, yc, min_snr=MIN_SNR, units=None):
    """A1 on the cells given: the change against its interval (Welch, cells as units). units: the unit of every reference
    cell where the reference side is made of units (other donors, healthy donors); that side then counts by unit."""
    d = y.mean() - yc.mean()
    vc = yc.var(ddof=1) / len(yc)
    if units is not None:                                    # the variance of a mean of cells that come in units
        el, inv = np.unique(units, return_inverse=True)
        dev = np.bincount(inv, weights=yc - yc.mean())
        vc = len(el) / max(len(el) - 1, 1) * float((dev ** 2).sum()) / len(yc) ** 2
    half = 1.96 * np.sqrt(y.var(ddof=1) / len(y) + vc)
    return {"cells": int(len(y)), "change_in_control_sd": float(d / yc.std(ddof=1)), "snr": float(abs(d) / half), "A1": bool(abs(d) / half >= min_snr)}


def select(ids, unit, y, seed, control="control", own=(), min_snr=MIN_SNR, reference_unit=None):
    """Every unit screened on its cells in the selection share against the control cells in it. `own`: units that act on
    the readout itself (A2), refused. The cells of the share enter no reading of u.
    control=None: no cell is a control cell; every unit is screened against the other units, which count as units.
    reference_unit: the unit of every cell, where the control cells come in units of their own (healthy donors)."""
    ids, unit = np.asarray(ids).astype(str), np.asarray(unit).astype(str)
    share = np.array([in_split(i, seed) for i in ids])
    if control is None:
        ref = lambda g: ((unit != g) & share, unit)
    else:
        ref = lambda g: ((unit == control) & share, reference_unit)
    one = lambda g, m, ru: screen(y[(unit == g) & share], y[m], min_snr, None if ru is None else np.asarray(ru).astype(str)[m])
    rows = sorted(({"event": str(g), **one(g, *ref(g)), "A2_own": str(g) in own} for g in sorted(set(unit) - {control})), key=lambda x: -x["snr"])
    return {"share": share, "pass_A1": [x["event"] for x in rows if x["A1"]], "admitted": [x["event"] for x in rows if x["A1"] and not x["A2_own"]], "table": rows}


# --- the reading of an output --------------------------------------------------------------------------------------------
def output_u(y, pred, treated, cluster=None, seed=0, control_cluster=None):
    """u of an output for one event, as the tool computes it for the system's own output: 1 - a * shift(pred) / shift(y),
    a the slope of y on pred among the control cells. Returns the draws (row 0 is the value on the data) and the slope.
    Treated cells are resampled by `cluster` when given, otherwise by cell; the cells of the other side by
    `control_cluster` when given (the other units of a donor event; healthy donors), otherwise by cell. The same seed
    gives the same resamples for every output, so differences between outputs are paired."""
    ti, ci_ = np.where(treated)[0], np.where(~treated)[0]
    t, c = SC._Side(y[ti, None], None if cluster is None else cluster[ti]), SC._Side(y[ci_, None], None if control_cluster is None else control_cluster[ci_])
    rng = np.random.default_rng(seed + 1)
    Ct, Cc = SC._counts(t.E, DRAWS, rng), SC._counts(c.E, DRAWS, rng)
    u, _, a = SC._u_shift(t, c, pred[ti, None], pred[ci_, None], Ct, Cc, np.ones(1))
    return u, float(np.atleast_1d(a)[0])


def halves(groups, seed=0, share=0.5):
    rng, first = np.random.default_rng(seed), np.zeros(len(groups), bool)
    for g in np.unique(groups):
        idx = rng.permutation(np.where(groups == g)[0])
        first[idx[: int(round(share * len(idx)))]] = True
    return first


# --- a repair ------------------------------------------------------------------------------------------------------------
GRID, FOLDS = [10.0 ** (k / 2) for k in range(-6, 5)], 5


def _cv_weight(F, t, pen, seed=0):
    fold = np.random.default_rng(seed).permutation(len(F)) % FOLDS
    err = np.zeros(len(GRID))
    for k in range(FOLDS):
        tr, te = fold != k, fold == k
        M, v = F[tr].T @ F[tr] / tr.sum(), F[tr].T @ t[tr] / tr.sum()
        for i, lam in enumerate(GRID):
            err[i] += ((F[te] @ np.linalg.solve(M + lam * np.diag(pen), v) - t[te]) ** 2).sum()
    err /= len(F)
    return GRID[int(err.argmin())], float(err.min())


class LinearRows:
    """A linear readout of a stage for one type: fitted on its control cells C (arm `controls`), updated once by the rows
    of its events (`rows`), by the rows of another type (`rows_other`), or fitted on every cell given (`cells`).
    x: the stage; y: the measured readout; o: the system's own output, whose scale the readout keeps."""

    def __init__(self, x, y, o, C, events, other=None):
        self.xc = x[C].mean(0)
        F = self.feat(x)
        scale = o[C].std() / y[C].std()
        t = o[C].mean() + scale * (y - y[C].mean())
        pen = np.r_[np.ones(x.shape[1]), 0.0] * x[C].var(0).mean()
        lam, cv = _cv_weight(F[C], t[C], pen)
        Minv = np.linalg.inv(F[C].T @ F[C] / len(C) + lam * np.diag(pen))
        self.th = {"controls": Minv @ (F[C].T @ t[C] / len(C))}
        self.log = {"n_controls": int(len(C)), "ridge_weight": lam, "cv_error": cv, "rows": []}

        def update(D, dy, se2):
            S = len(C) * scale ** 2 * se2 / cv
            K = Minv @ D.T @ np.linalg.inv(D @ Minv @ D.T + np.diag(S))
            before = scale * dy - D @ self.th["controls"]
            th = self.th["controls"] + K @ before
            return th, 1 - np.abs(scale * dy - D @ th) / np.maximum(np.abs(before), 1e-12)

        vc = y[C].var(ddof=1) / len(C)
        if events:
            D = np.array([F[T].mean(0) - F[C].mean(0) for _n, T in events])
            dy = np.array([y[T].mean() - y[C].mean() for _n, T in events])
            se2 = np.array([y[T].var(ddof=1) / len(T) + vc for _n, T in events])
            self.th["rows"], closed = update(D, dy, se2)
            self.log["rows"] = [{"event": n, "n_cells": int(len(T)), "measured_change": float(a), "standard_error": float(np.sqrt(s)), "share_of_the_residual_closed": float(c)}
                                for (n, T), a, s, c in zip(events, dy, se2, closed)]
            tr = np.r_[C, np.concatenate([T for _n, T in events])]
            lam2, _cv = _cv_weight(F[tr], t[tr], pen)
            self.th["cells"] = np.linalg.solve(F[tr].T @ F[tr] / len(tr) + lam2 * np.diag(pen), F[tr].T @ t[tr] / len(tr))
            self.log["ridge_weight_cells"] = lam2
        if other is not None:                               # the rows of another type, in units of the control spread
            xo, yo, Ao, to, sites = other
            Fo, Co = self.feat(xo), np.where(Ao == 0)[0]
            conv = y[C].std() / yo[Co].std()
            ev = [(g, np.where((Ao == 1) & (to == g))[0]) for g in sites]
            D = np.array([Fo[T].mean(0) - Fo[Co].mean(0) for _n, T in ev])
            dy = conv * np.array([yo[T].mean() - yo[Co].mean() for _n, T in ev])
            se2 = conv ** 2 * np.array([yo[T].var(ddof=1) / len(T) + yo[Co].var(ddof=1) / len(Co) for _n, T in ev])
            self.th["rows_other"], _c = update(D, dy, se2)

    def feat(self, x):
        return np.hstack([x - self.xc, np.ones((len(x), 1))])

    def out(self, arm, x):
        return self.feat(x) @ self.th[arm]


class UnitRows:
    """The readout of `LinearRows` for events that have no control cells: a unit (a donor, a batch) read against the other
    units. events: [(unit, its cells)], every unit of the fit; C is not used.
      controls   a readout of the variation among the cells of a unit (x and y centred by unit): what can be learned with
                 no record of how any unit differs from the others. It takes the place of the readout of control cells
      rows       the same, updated once by the row of every unit (the unit against the other units of the fit)
      cells      fitted on every cell, not centred
      unseen|U   for every unit U: the arm `rows` of a fit that holds no cell of U (its readout, its rows, its centre)
    The ridge weight is chosen once, on all units of the fit."""

    def __init__(self, x, y, o, C, events, other=None):
        events = [(str(n), np.asarray(T)) for n, T in events if len(T)]
        A = np.concatenate([T for _n, T in events])
        scale = o[A].std() / y[A].std()
        t = o[A].mean() + scale * (y - y[A].mean())
        self.xc, p = x[A].mean(0), x.shape[1]
        st = []                                                 # per unit: n, sums of x, t, y, y*y, and the scatter of x and of x with t about the unit's own mean
        for _n, T in events:
            xm, tm = x[T].mean(0), t[T].mean()
            st.append((len(T), x[T].sum(0), t[T].sum(), y[T].sum(), (y[T] ** 2).sum(), (x[T] - xm).T @ (x[T] - xm), (x[T] - xm).T @ (t[T] - tm)))
        xw = np.vstack([x[T] - x[T].mean(0) for _n, T in events])
        tw = np.concatenate([t[T] - t[T].mean() for _n, T in events])
        pen = np.ones(p) * xw.var(0).mean()
        lam, cv = _cv_weight(xw, tw, pen)
        self.th, self.log = {}, {"n_cells": int(len(A)), "units": len(events), "ridge_weight": lam, "cv_error": cv, "rows": []}

        def fit(keep):
            """The readout of the units `keep`: within-unit weights, the centre, and the rows of those units among themselves."""
            n = np.array([st[k][0] for k in keep], np.float64)
            N, sx, sy = n.sum(), np.array([st[k][1] for k in keep]), np.array([st[k][3] for k in keep])
            syy, tbar = np.array([st[k][4] for k in keep]), sum(st[k][2] for k in keep) / n.sum()
            Sw = sum(st[k][5] for k in keep) / N
            Minv = np.linalg.inv(Sw + lam * np.diag(pen))
            w = Minv @ (sum(st[k][6] for k in keep) / N)
            xc = sx.sum(0) / N
            D = sx / n[:, None] - (sx.sum(0) - sx) / (N - n)[:, None]                 # a unit against the other units of the fit
            dy = sy / n - (sy.sum() - sy) / (N - n)
            var = lambda s1, s2, m: (s2 - s1 ** 2 / m) / np.maximum(m - 1, 1)
            se2 = var(sy, syy, n) / n + var(sy.sum() - sy, syy.sum() - syy, N - n) / (N - n)
            S = N * scale ** 2 * se2 / cv
            K = Minv @ D.T @ np.linalg.inv(D @ Minv @ D.T + np.diag(S))
            before = scale * dy - D @ w
            w2 = w + K @ before
            closed = 1 - np.abs(scale * dy - D @ w2) / np.maximum(np.abs(before), 1e-12)
            full = lambda v: np.r_[v, tbar + (self.xc - xc) @ v]                     # in the features of this object, which are centred at self.xc
            return full(w), full(w2), dy, se2, closed

        every = list(range(len(events)))
        self.th["controls"], self.th["rows"], dy, se2, closed = fit(every)
        self.log["rows"] = [{"event": n, "n_cells": int(len(T)), "measured_change": float(a), "standard_error": float(np.sqrt(s)), "share_of_the_residual_closed": float(c)}
                            for (n, T), a, s, c in zip(events, dy, se2, closed)]
        F = self.feat(x[A])
        pen1 = np.r_[pen, 0.0]
        lam2, _cv = _cv_weight(F, t[A], pen1)
        self.th["cells"] = np.linalg.solve(F.T @ F / len(A) + lam2 * np.diag(pen1), F.T @ t[A] / len(A))
        self.log["ridge_weight_cells"] = lam2
        if len(events) > 2:
            for k, (name, _T) in enumerate(events):
                self.th[f"unseen|{name}"] = fit([j for j in every if j != k])[1]

    def feat(self, x):
        return np.hstack([x - self.xc, np.ones((len(x), 1))])

    def out(self, arm, x):
        return self.feat(x) @ self.th[arm]


def cross_fit(x, y, o, treated, unit, sites, first, repair=LinearRows, other=None, fixed=None):
    """Every cell is read by a repair fitted without it (two halves). Returns the output of every arm on the cells it was
    not fitted on, the output of the arm `rows` on the cells it was fitted on (gate 3), and the logs of the fits.
    fixed: {arm: a fitted repair applied as it is} (the repair of another type)."""
    n = len(y)
    out, inside, fits = {}, np.full(n, np.nan), []
    for fold in (0, 1):
        fit = first if fold == 0 else ~first
        rd = ~fit
        C = np.where(fit & ~treated)[0]
        ev = [(s, np.where(fit & treated & (unit == s))[0]) for s in sites]
        R = repair(x, y, o, C, [(s, T) for s, T in ev if len(T)], other)
        fits.append(getattr(R, "log", None))
        for a in R.th:
            out.setdefault(a, np.full(n, np.nan))[rd] = R.out(a, x[rd])
        for a, Rf in (fixed or {}).items():
            out.setdefault(a, np.full(n, np.nan))[rd] = Rf.out("rows", x[rd])
        inside[fit] = R.out("rows", x[fit])
    return out, inside, fits


def read(y, o, treated, unit, sites, out, inside, reference, event_reference=None, events=None, seed=0, control_cluster=None, references=None):
    """The reading of a round for one type. sites: the events whose rows repair it, pooled; events: the events also read
    one by one (default: the sites). reference: the draws of the reference stage for the pool (from the grid, paired by
    seed); event_reference: {event: draws}. Cells of events outside `sites` enter no pooled reading.
    control_cluster: the unit of every cell where the side an event is read against is made of units (`output_u`).
    references: {arm: draws} for an arm that is read against another reference than the round's (the arm `unseen` of a
    unit event, against a reader that never saw the unit)."""
    pool = ~treated | np.isin(unit, sites)
    cl = unit[pool] if len(sites) > 1 else None
    cc = None if control_cluster is None else np.asarray(control_cluster)
    output = lambda yy, pp, tt, c_, m: output_u(yy, pp, tt, c_, seed, None if cc is None else cc[m])
    u0, s0 = output(y[pool], o[pool], treated[pool], cl, pool)
    ref = np.asarray(reference, np.float64)
    e = {"events": list(sites), "resampled_by": "event" if cl is not None else "cell", "resampling_units": len(sites) if cl is not None else None,
         "approximate": bool(cl is not None and len(sites) < 10), "reference": val(ref),            # an interval from fewer than ten events is approximate, as in the tool
         "as_trained": {**val(u0), "slope": s0, "beyond_reference": diff(u0, ref)}, "arms": {}, "by_event": {}}
    U = {}
    for a, p in out.items():
        u, slope = output(y[pool], p[pool], treated[pool], cl, pool)
        U[a] = u
        T_, C_ = treated[pool], ~treated[pool]
        gain = np.abs(u0) - np.abs(u)
        ra = None if not references or a not in references else np.asarray(references[a], np.float64)
        e["arms"][a] = {**val(u), "slope": slope, "valid": bool(slope > 0 and np.sign(p[pool][T_].mean() - p[pool][C_].mean()) == np.sign(y[pool][T_].mean() - y[pool][C_].mean())),
                        "gain": float(gain[0]), "gain_ci": interval(gain[1:]), "beyond_reference": diff(u, ref if ra is None else ra), **({"reference": val(ra)} if ra is not None else {})}
    for s in (events if events is not None else (sites if len(sites) > 1 else [])):
        m = ~treated | (unit == s)
        ue0, _s = output(y[m], o[m], treated[m], None, m)
        r = None if not event_reference or s not in event_reference else np.asarray(event_reference[s], np.float64)
        e["by_event"][s] = {**({"reference": val(r)} if r is not None else {}), "as_trained": {**val(ue0), **({"beyond_reference": diff(ue0, r)} if r is not None else {})}, "arms": {}}
        for a, p in out.items():
            ue, _s = output(y[m], p[m], treated[m], None, m)
            e["by_event"][s]["arms"][a] = {**val(ue), **({"beyond_reference": diff(ue, r)} if r is not None else {})}
    def step_between(a, b):                                          # what arm b removes beyond arm a, on paired resamples
        d = np.abs(U[a]) - np.abs(U[b])
        iv = interval(d[1:])
        return {"gain_difference": float(d[0]), "ci": iv, "reading": "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")}
    e["rows_minus_controls"] = step_between("controls", "rows")      # what the rows add over a readout of control cells (3b)
    if "unseen" in U:
        e["rows_minus_unseen"] = step_between("unseen", "rows")      # what the unit's own row adds over the rows of the other units
    ui, _s = output(y[pool], inside[pool], treated[pool], cl, pool)
    g3, g0 = diff(ui, ref), e["as_trained"]["beyond_reference"]
    e["A3"] = {**g0, "holds": bool(g0["ci"][0] > STEP)}             # a repair has something to remove: the output as trained loses a step more than the reference stage
    e["gate3"] = {"on_the_cells_it_was_fitted_on": val(ui), "against_the_reference": g3, "passes": bool(g3["ci"][1] < STEP),
                  "control_the_output_as_trained": {**g0, "passes": bool(g0["ci"][1] < STEP)}}
    e["closure"] = grade(e)
    return e


def grade(e, attainable=None):
    r = e["arms"]["rows"]
    attainable = attainable or ("absolute" if e["reference"]["verdict"] == "support" else "relative")
    if e["as_trained"]["verdict"] == "support":
        g = "nothing to close"
    elif not e["A3"]["holds"]:
        g = "not a target for a repair"
    elif not e["gate3"]["passes"]:
        g = "not tested (gate 3)"
    elif r["verdict"] == "support" and r["valid"] and e["as_trained"]["verdict"] == "reject":
        g = "absolute"
    elif r["beyond_reference"]["reading"] == "no added loss" and r["valid"] and e["as_trained"]["beyond_reference"]["ci"][0] > STEP:
        g = "relative"
    elif r["beyond_reference"]["reading"] == "unresolved":
        g = "loss no longer resolved"
    else:
        g = "not closed"
    order = ["not closed", "loss no longer resolved", "relative", "absolute"]
    return {"grade": g, "attainable": attainable, "reached_what_was_attainable": bool(g == "nothing to close" or (g in order and order.index(g) >= order.index(attainable)))}


# --- where a loss is -------------------------------------------------------------------------------------------------------
def grid(res, target="shift"):
    """A result of scope.run_spec as {event: {stage: {scope: cell}}} and its draws as {'event|target|stage|scope|class': draws}."""
    by = {r["event"]["unit"]: r["targets"][target]["stages"] for r in res["results"]}
    draws = {f"{r['event']['unit']}|{k}": np.asarray(v, np.float64) for r in res["results"] for k, v in (r.get("_draws") or {}).items()}
    return by, draws


def load(folder, target="shift", readout=None):
    """The files the tool writes (scope.json, scope_draws.npz), in the form of `grid`. readout: the index of one readout of
    a run made with `each_readout`; its grid and its draws are returned under the name of the target."""
    import json
    from pathlib import Path
    res = json.loads((Path(folder) / "scope.json").read_text()) if not isinstance(folder, dict) else folder["json"]
    z = np.load(Path(folder) / "scope_draws.npz") if not isinstance(folder, dict) else folder["draws"]
    if readout is None:
        return {r["event"]["unit"]: r["targets"][target]["stages"] for r in res["results"]}, {k: np.asarray(z[k], np.float64) for k in z.files}
    tag = f"|{target}@{int(readout)}|"
    return ({r["event"]["unit"]: r["targets"][target]["each_readout"][int(readout)]["stages"] for r in res["results"]},
            {k.replace(tag, f"|{target}|"): np.asarray(z[k], np.float64) for k in z.files if tag in k})


def merge(*grids):
    """Grids of the same events read at different stages (one run of the tool per view), as one."""
    by, draws = {}, {}
    for b, d in grids:
        for u, S in b.items():
            by.setdefault(u, {}).update(S)
        draws.update(d)
    return by, draws


def _classes(cell):
    by = cell.get("by_class") or {cell["class"]: cell}
    return {c: {"u": x["u"], "ci": x["ci"], "verdict": verdict(x["ci"])} for c, x in by.items()}


def _reader(e, stages, scope, rule="union"):
    """The reader of a view that may have several stages and reader classes, and what the view reads. rule "union" is the
    tool's own (scope._union): the view is supported if any declared reader supports it, rejected only if every one
    rejects, and otherwise unresolved; the reader is the supporting one with the tightest interval, else the one with the
    smallest |u|. rule "smallest": the reader with the smallest |u| and its own verdict."""
    R = [(st, c, x) for st in stages for c, x in e["by_stage"][st][scope].items()]
    sup = [r for r in R if r[2]["verdict"] == "support"]
    if rule == "union" and sup:
        st, c, x = min(sup, key=lambda r: max(abs(r[2]["ci"][0]), abs(r[2]["ci"][1])))
        return st, c, {**x, "verdict": "support"}
    st, c, x = min(R, key=lambda r: abs(r[2]["u"]))
    if rule == "union":
        x = {**x, "verdict": "reject" if all(r[2]["verdict"] == "reject" for r in R) else "unresolved"}
    return st, c, x


def place(by, draws, model_stage, inputs, others=(), control=None, scope="own", target="shift", rule="union", linear=("ridge",), between=()):
    """Where a loss before the model is, for every event of a grid. model_stage: the stage the model reads. inputs: stages
    that are views of the model's own input. others: stages of measured inputs the model does not receive. control: a
    stage that must read nothing (a permuted view). rule: how the reader of a view is chosen among its stages and reader
    classes (`_reader`). linear: the names of the reader classes that are linear (a view read by linear readers only is
    not yet a basis for a conclusion about the input). between: stages inside the model that lie between its input and
    `model_stage`, in order (the output of the encoder of one modality, a fusion step); what is lost from the input to
    the model's stage is then read step by step along them.
    For an event: the reading at every stage by reader class, the gap from the input to the model's stage on paired
    resamples (below zero throughout: the input view is dominated and cannot be read), and what the other inputs read.
    For the family: the positive control (some reader at some stage supports some event) and the control view."""
    stages = [*inputs, *between, model_stage, *others]
    res, seen = {"rule": rule, "inputs": list(inputs), "between": list(between), "linear": list(linear), "events": {}}, []
    for u, S in by.items():
        e = {"by_stage": {st: {sc: _classes(S[st][sc]) for sc in ("units", "own") if sc in S[st]} for st in stages + ([control] if control else [])}}
        seen += [x["verdict"] for st in stages for sc in e["by_stage"][st].values() for x in sc.values()]
        dr = lambda st, c: draws[f"{u}|{target}|{st}|{scope}|{c}"]
        ms, mc, mx = _reader(e, [model_stage], scope, rule)
        st_in, ic, ix = _reader(e, inputs, scope, rule)
        g = diff(dr(ms, mc), dr(st_in, ic))                          # what is added from the input to the stage the model reads
        e["model_stage"] = {"stage": ms, "reader": mc, **mx}
        e["input"] = {"stage": st_in, "reader": ic, **ix}
        e["from_the_input_to_the_model_stage"] = {**g, "dominated": bool(g["ci"][1] < 0)}
        path = [("the input", st_in, ic)] + [(st, *_reader(e, [st], scope, rule)[:2]) for st in between] + [(ms, ms, mc)]
        e["along_the_model"] = [{"from": a[0], "to": b[0], **diff(dr(b[1], b[2]), dr(a[1], a[2]))} for a, b in zip(path, path[1:])]
        e["lost_between"] = [[x["from"], x["to"]] for x in e["along_the_model"] if x["reading"] == "loss"]
        e["other_inputs"] = {}
        for st in others:
            _s, oc, ox = _reader(e, [st], scope, rule)
            e["other_inputs"][st] = {"reader": oc, **ox, "reads_more_than_the_input": diff(dr(st_in, ic), dr(st, oc))}
        if control:
            e["control"] = min(e["by_stage"][control][scope].values(), key=lambda x: abs(x["u"]))
        res["events"][u] = e
    res["positive_control"] = {"holds": "support" in seen, "readings": len(seen), "supported": seen.count("support"), "unresolved": seen.count("unresolved"), "rejected": seen.count("reject")}
    if control:
        res["control_view"] = {"holds": bool(all(abs(e["control"]["u"]) >= 0.9 for e in res["events"].values())), "values": {u: e["control"]["u"] for u, e in res["events"].items()}}
    return res


def before(placed, event, evidence=None):
    """The part of a placing that `outer.open_items` reads for one event or for the pool."""
    e, pc = placed["events"][event], placed["positive_control"]
    v = lambda x: {k: x[k] for k in ("u", "ci", "verdict")}
    classes = sorted({c for st in placed["inputs"] for sc in e["by_stage"][st].values() for c in sc})       # every reader class of every stage of the input view
    iv = {"own": v(e["input"]), "gap_to_later": {k: e["from_the_input_to_the_model_stage"][k] for k in ("added", "ci", "reading")}, "classes": classes,
          "nonlinear": any(c not in placed["linear"] for c in classes),
          "positive_control": pc["holds"], "readout_control": None if "control_view" not in placed else placed["control_view"]["holds"], "evidence": evidence,
          "positive_control_note": f"of {pc['readings']} readings of the family {pc['supported']} are supported, {pc['unresolved']} unresolved, {pc['rejected']} rejected",
          "note": f"reader of the input: {e['input']['stage']}, {e['input']['reader']}", "along_the_model": e["along_the_model"], "lost_between": e["lost_between"],
          "other_inputs": {st: {**v(x), "reads_more_than_the_input": x["reads_more_than_the_input"]} for st, x in e["other_inputs"].items()}}
    return {"reference": {"units": v(e["model_stage"]), "evidence": evidence}, "input_view": iv}


def target(readout, group, type_, e, placed=None, events_alone=(), attainable=None, evidence=None, pair=None, **more):
    """A round's reading as one target of outer.open_items. e: the result of `read`; placed: the result of `place` on the
    same events; events_alone: events placed one by one; pair: the name of this type and the type its arms `other` and
    `rows_other` come from. A repair is read only where it had something to remove (A3) and the output as trained is not
    supported already; elsewhere the target is described by the output as trained."""
    v = lambda x: {k: x[k] for k in ("u", "ci", "verdict")}
    none = {"added": None, "ci": None, "reading": None}
    sound = e["as_trained"]["verdict"] == "support"
    repaired = bool(e["A3"]["holds"] and not sound)
    arm = (lambda x: x["arms"]["rows"]) if repaired else (lambda x: x["as_trained"])
    own = arm(e)
    t = {"readout": readout, "group": group, "type": type_, "events_n": len(e["events"]), "admitted": True, "response": {"snr": None}, "attainable": attainable or e["closure"]["attainable"],
         "evidence": evidence, "reference": {"units": v(e["reference"]), "evidence": evidence},
         "repair": "read" if repaired else ("none: the output as trained is supported" if sound else "none: the output as trained loses no step more than the reference stage (A3)"),
         "model": {"own": v(own), "added": own["beyond_reference"], "reading": None, "rows_of_type": len(e["events"]), "evidence": evidence,
                   "events": {s: {**v(arm(x)), "added": arm(x).get("beyond_reference", none), "in_pool": s in e["events"], **({"reference": x["reference"]} if "reference" in x else {})}
                              for s, x in e["by_event"].items()}},
         "gates": {"evaluable": bool(e["gate3"]["passes"]) if repaired else True, "why": [], "guard_kept": True}, "arms": {}, **more}
    if not repaired and not sound and e["A3"]["reading"] == "no added loss":
        t["control_adds_a_loss"] = False                        # nothing in the model can be closed; whether what is before it is open or undecided is read from the reference stage
    if repaired and not e["gate3"]["passes"]:
        g3 = e["gate3"]["against_the_reference"]
        t["gates"]["why"].append({"question": "gate 3 does not pass: on the cells the repair was fitted on, the output loses more than the reference stage", "value": g3["added"], "evidence": evidence})
    if repaired and pair:
        t["arms"] = {a: {**e["arms"][a]["beyond_reference"], "pair": pair, "evidence": evidence} for a in ("other", "rows_other") if a in e["arms"]}
    if placed:
        if "pool" in placed["events"]:
            t["input_view"] = before(placed, "pool", evidence)["input_view"]
        for s in events_alone:
            if s in t["model"]["events"] and s in placed["events"]:
                t["model"]["events"][s]["before"] = before(placed, s, evidence)
    return t


def target_by_event(readout, group, type_, reads, placed=None, attainable=None, evidence=None, **more):
    """A type whose events are too few to be read as a pool: every event on its own, resampled by cell, and no statement
    about a family. reads: {event: the result of `read` with that event as its only site}. The target carries its least
    favourable event; where no event leaves the model anything to remove (A3) and not all are read already, what is lost
    is lost before the model."""
    v = lambda x: {k: x[k] for k in ("u", "ci", "verdict")}
    ev, sound, repaired = {}, {}, {}
    for s, e in reads.items():
        sound[s] = e["as_trained"]["verdict"] == "support"
        repaired[s] = bool(e["A3"]["holds"] and not sound[s])
        a = e["arms"]["rows"] if repaired[s] else e["as_trained"]
        ev[s] = {**v(a), "added": a["beyond_reference"], "in_pool": False, "reference": v(e["reference"])}
        if repaired[s] and "unseen" in e["arms"]:                # a unit event: the same repair without the unit's own row and cells
            ev[s]["unseen"] = {**v(e["arms"]["unseen"]), "added": e["arms"]["unseen"]["beyond_reference"], "own_row_adds": e["rows_minus_unseen"]}
        if placed and s in placed["events"]:
            ev[s]["before"] = before(placed, s, evidence)
    order = ["loss", "unresolved", "no added loss"]
    worst = min(ev, key=lambda s: (order.index(ev[s]["added"]["reading"]), -abs(ev[s]["u"])))
    failed = [s for s in ev if repaired[s] and not reads[s]["gate3"]["passes"]]
    t = {"readout": readout, "group": group, "type": type_, "events_n": len(ev), "admitted": True, "response": {"snr": None}, "evidence": evidence, "reference": None,
         "attainable": attainable or ("relative" if any(e["closure"]["attainable"] == "relative" for e in reads.values()) else "absolute"),
         "repair": {s: "read" if repaired[s] else "none" for s in ev},
         "model": {"own": v(ev[worst]), "added": ev[worst]["added"], "reading": None, "rows_of_type": len(ev), "evidence": evidence, "events": ev},
         "gates": {"evaluable": not failed, "why": [{"question": "gate 3 does not pass for an event read on its own: " + ", ".join(failed), "evidence": evidence}] if failed else [], "guard_kept": True},
         "arms": {}, **more}
    open_ = [s for s in ev if not sound[s]]                 # no event leaves the model a step to remove: nothing in the model can be closed
    if open_ and all(reads[s]["A3"]["reading"] == "no added loss" for s in open_):
        t["control_adds_a_loss"] = False
    return t


def one_type(readout, group, type_, x, y, o, treated, unit, sites, by, draws, model_stage, inputs=(), others=(), control=None, alone=(), first=None, repair=LinearRows,
             other=None, fixed=None, pair=None, evidence=None, seed=0, target_name="shift", between=(), control_cluster=None, **more):
    """One round for one type: from its cells and its grid to its target.
    x: the stage the model reads; y: the readout, measured outside; o: the system's own output; treated, unit: the record.
    sites: the admitted events whose rows repair the type; alone: admitted events of the type read on their own.
    by, draws: the grid of the audit tool for these events (and "pool"), at `model_stage`, at the stages of the model's
    input (`inputs`), at stages inside the model (`between`), of other measured inputs (`others`) and of a control view.
    other, fixed: the cells and the fitted repair of another type, for the arms that test what carries across (cross_fit).
    control_cluster: the unit of every cell where the control cells come in units of their own (a disease: the events are
    the patient donors, the control cells are those of the healthy donors); that side is then resampled by unit."""
    if not sites:
        raise ValueError("no admitted event: the type is not a target (protocol 3c)")
    first = halves(np.where(treated, unit, "control"), seed) if first is None else first
    out, inside, fits = cross_fit(x, y, o, treated, unit, sites, first, repair, other, fixed)
    ref = lambda u, sc: draws[f"{u}|{target_name}|{model_stage}|{sc}|{by[u][model_stage][sc]['class']}"]
    placed = place(by, draws, model_stage, list(inputs), others, control, target=target_name, between=between) if inputs else None
    if len(sites) >= MIN_POOL:
        refs = {**{s: ref(s, "units") for s in sites}, **{s: ref(s, "own") for s in alone}}
        e = read(y, o, treated, unit, sites, out, inside, ref("pool", "units"), refs, events=[*sites, *alone], seed=seed, control_cluster=control_cluster)
        k = f"pool|{target_name}|system|frozen|system"
        if k in draws and abs(e["as_trained"]["u"] - draws[k][0]) > 2e-4:
            raise ValueError("the resamples of the cycle and of the tool do not pair: the seeds or the cells differ")
        t = target(readout, group, type_, e, placed, events_alone=[*sites, *alone], evidence=evidence, pair=pair, **more)
        return {"read_as": "a pool", "reading": e, "placed": placed, "target": t, "fits": fits}
    reads = {s: read(y, o, treated, unit, [s], out, inside, ref(s, "units" if s in sites else "own"), seed=seed, control_cluster=control_cluster) for s in [*sites, *alone]}
    t = target_by_event(readout, group, type_, reads, placed, evidence=evidence, **more)
    return {"read_as": "event by event", "reading": reads, "placed": placed, "target": t, "fits": fits}


def one_units(readout, group, type_, x, y, o, unit, events, by, draws, model_stage, inputs=(), others=(), control=None, first=None, repair=UnitRows, evidence=None, seed=0,
              target_name="shift", between=(), **more):
    """One round for a type whose events have no control cells: every admitted unit is read against the other units.
    unit: the unit of every cell (a donor, a batch); events: the admitted units. Every unit of the type is in the fit; no
    cell is a control cell. The fit of a half holds half the cells of every unit, so the arm `rows` has the unit's own
    row, as an event with control cells has its own; the arm `unseen` is the same repair without any cell of the unit. Each event is read
    on its own: its cells by cell, the other units by unit, as the audit tool reads it. `rows`, `controls` and `cells`
    are read against a reader of the stage with the unit's own cells, `unseen` against one that never saw the unit."""
    if not len(events):
        raise ValueError("no admitted event: the type is not a target (protocol 3c)")
    unit = np.asarray(unit).astype(str)
    first = halves(unit, seed) if first is None else first
    out, inside, fits = cross_fit(x, y, o, np.ones(len(y), bool), unit, sorted(set(unit)), first, repair)
    placed = place({s: by[s] for s in events}, draws, model_stage, list(inputs), others, control, target=target_name, between=between) if inputs else None
    reads = {}
    for s in events:
        ref = lambda sc: draws[f"{s}|{target_name}|{model_stage}|{sc}|{by[s][model_stage][sc]['class']}"]
        arms = {a: p for a, p in out.items() if "|" not in a}
        refs = None
        if f"unseen|{s}" in out:
            arms["unseen"], refs = out[f"unseen|{s}"], {"unseen": ref("units")}
        e = read(y, o, unit == s, unit, [s], arms, inside, ref("own"), seed=seed, control_cluster=unit, references=refs)
        k = f"{s}|{target_name}|system|frozen|system"
        if k in draws and abs(e["as_trained"]["u"] - draws[k][0]) > 2e-4:
            raise ValueError("the resamples of the cycle and of the tool do not pair: the seeds or the cells differ")
        reads[s] = e
    t = target_by_event(readout, group, type_, reads, placed, evidence=evidence, **more)
    return {"read_as": "unit by unit", "reading": reads, "placed": placed, "target": t, "fits": fits}

