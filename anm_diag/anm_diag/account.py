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
MIN_ROWS = 8                                                 # rows of other events below which a part of a row's own cannot be told from what a readout would learn
LAMBDAS = 10.0 ** np.arange(-3, 3.01, 0.5)
interval = lambda d: [float(v) for v in np.percentile(d, [100 * ALPHA / 2, 100 * (1 - ALPHA / 2)])]
step = lambda iv: "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")


def _units(cluster, cells):
    """The units an arm is resampled by: its clusters if it is made of two or more, else None (its cells)."""
    if cluster is None:
        return None
    u, inv = np.unique(cluster[cells], return_inverse=True)
    return inv if len(u) >= 2 else None


def _within(v, inv):
    """v centred: by unit where the cells come in units, else as a whole."""
    if inv is None:
        return v - v.mean(0)
    n = np.bincount(inv).astype(float)
    if v.ndim == 1:
        return v - (np.bincount(inv, weights=v) / n)[inv]
    S = np.zeros((len(n), v.shape[1]))
    np.add.at(S, inv, v)
    return v - (S / n[:, None])[inv]


def _var_of_mean(v, inv):
    """The variance of the mean of cells: by cell, or by unit where the cells come in units."""
    if inv is None:
        return float(v.var(ddof=1) / len(v))
    dev = np.bincount(inv, weights=v - v.mean())
    return float(len(dev) / (len(dev) - 1) * (dev ** 2).sum() / len(v) ** 2)


def one_scale(x, C, inv=None):
    """The spread of every column among the control cells of the fit, within their unit where they come in units (a
    column that does not vary is left as it is)."""
    sd = _within(x[C], inv).std(0)
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


def _control_fit(x, y, C, seed=0, inv=None):
    """The ridge readout of the reference cells of centred columns, with what it leaves uncertain: (readout, Mi, cv / n,
    root). Where the reference cells come in units they are centred by unit: the readout is then one of the variation
    inside a unit, which is what can be learned with no record of how any unit differs from another.
    The same cells, stage and readout give the same fit for every event read against them: kept by content."""
    key = (x.shape[1], len(C), float(x[C].sum()), hash(np.ascontiguousarray(x[C[:16]]).tobytes()), hash(y[C].tobytes()), seed, None if inv is None else hash(inv.tobytes()))
    if key not in _KEPT:
        if len(_KEPT) >= 12:
            _KEPT.pop(next(iter(_KEPT)))
        _KEPT[key] = _control_fit_(x, y, C, seed, inv)
    _KEPT[key] = _KEPT.pop(key)                                 # the one used last is dropped last
    return _KEPT[key]


def _control_fit_(x, y, C, seed=0, inv=None):
    Xc, yc, n = _within(x[C], inv), _within(y[C], inv), len(C)
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


WIDENINGS = 10.0 ** np.arange(0, 8.01, 0.25)


def _teacher(Mi, k, D, rho, s, nu=4.0, iters=30):
    """What rows with errors teach a readout. The readout of the reference cells is known up to k Mi; its error among
    events may be wider by a factor that the rows themselves choose (their likelihood). With MIN_ROWS rows or more, a row
    may also hold a part of its own that no readout explains (small for most rows, heavy-tailed, its spread from the
    rows, as `Law`): such a row teaches less. With fewer rows that part cannot be told from what a readout would learn,
    and every row teaches. Returns the matrix A such that readout + A (rows - D readout) is the readout given the rows,
    the rows' shifts, and the factor. With exact rows and a factor without bound this is `_closer`: the smallest change
    of the readout that closes the rows. A row inside its error teaches nothing."""
    D = np.atleast_2d(D)
    K0, s2, m = k * D @ Mi @ D.T, np.maximum(s, 1e-12) ** 2, len(rho)

    def widen(extra):
        best = (np.inf, WIDENINGS[0])
        for t2 in WIDENINGS:
            C = t2 * K0 + np.diag(s2 + extra)
            sign, logdet = np.linalg.slogdet(C)
            nll = 0.5 * (logdet + rho @ np.linalg.solve(C, rho))
            if nll < best[0] - 1e-9:
                best = (nll, t2)
        return best[1]
    Vd = np.zeros(m)
    t2 = widen(Vd)
    if m >= MIN_ROWS:
        sigma0, w = float(1.4826 * np.median(np.abs(rho - np.median(rho)))) + 1e-9, np.ones(m)
        for _ in range(iters):
            Vd = sigma0 ** 2 / w
            t2 = widen(Vd)
            Si = np.linalg.inv(t2 * K0 + np.diag(s2 + Vd))
            d = Vd * (Si @ rho)
            Ed2 = d * d + Vd - Vd * np.diag(Si) * Vd
            w = (nu + 1) / (nu + Ed2 / sigma0 ** 2)
            sigma0 = float(np.sqrt((w * Ed2).mean())) + 1e-12
        Vd = sigma0 ** 2 / w
    return t2 * k * Mi @ D.T @ np.linalg.inv(t2 * K0 + np.diag(s2 + Vd)), D, float(t2)


def reference(unit, event, control="control"):
    """The cells an event is read against: its control cells, or, where no cell is a control cell (control=None: a donor,
    a batch, a study), the cells of every other unit."""
    return (unit == control) if control is not None else (unit != event)


def reading(xs, y, unit, first, event, others=(), control="control", o=None, draws=DRAWS, seed=0, cluster=None):
    """The account of one event's change of one readout at one stage.

    xs: the stage, cells by columns; or (the stage for the half fitted on `first`, the stage for the other half) when the
    stage itself was chosen on the fitted half. unit: the event of every cell. first: the half of every cell. others:
    the other events whose rows may teach a readout. o: the system's own output as trained, if it has one.

    Three shapes of event are read by the same three readings (GENERAL.md, Section 6):
      an event with control cells            control: their label. Cells are the units of every resample.
      a unit against the other units         control=None: no cell is a control cell; a donor, a batch or a study is read
                                             against the cells of every other unit. cluster: the unit of every cell.
      a group of units against another       control: the label of the reference cells (healthy donors); cluster: the
                                             donor of every cell.
    cluster: an arm made of two or more units is resampled by unit, and the reference readout is one of the variation
    inside a unit (the cells centred by unit): with no record of an event, nothing says how one unit differs from another.

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
    ref = reference(unit, event, control)
    cl = None if cluster is None else np.asarray(cluster)
    keys = ["read", "carried", "found_again"] + (["others"] if others else []) + (["as_trained"] if o is not None else [])
    num, den = {k: np.zeros(draws + 1) for k in keys if k != "carried"}, np.zeros(draws + 1)
    sxy, sxx, see = np.zeros(draws + 1), np.zeros(draws + 1), 0.0
    vis, wide = [], []
    for f, fit in ((0, first), (1, ~first)):
        x = _stage(xs, f)
        Cf, Tf = np.where(fit & ref)[0], np.where(fit & (unit == event))[0]
        Cs, Ts = np.where(~fit & ref)[0], np.where(~fit & (unit == event))[0]
        uT, uC, uF = _units(cl, Ts), _units(cl, Cs), _units(cl, Cf)
        if uT is None and uC is None:                                                       # cells are the units of both arms
            idx = [(np.arange(len(Ts)), np.arange(len(Cs)))] + [(rng.integers(0, len(Ts), len(Ts)), rng.integers(0, len(Cs), len(Cs))) for _ in range(draws)]
            WT = np.array([np.bincount(it, minlength=len(Ts)) for it, _ic in idx], np.float64)
            WC = np.array([np.bincount(ic, minlength=len(Cs)) for _it, ic in idx], np.float64)
        else:                                                                               # an arm made of units is drawn by unit
            def weights(inv, n):
                inv = np.arange(n) if inv is None else inv
                g = int(inv.max()) + 1
                return np.vstack([np.ones(n), rng.multinomial(g, np.full(g, 1.0 / g), size=draws)[:, inv].astype(np.float64)])
            WT, WC = weights(uT, len(Ts)), weights(uC, len(Cs))
        zth = np.vstack([np.zeros(x.shape[1]), rng.normal(size=(draws, x.shape[1]))])        # for the control readout: what its fit leaves uncertain
        zrow = np.vstack([np.zeros(1 + len(others)), rng.normal(size=(draws, 1 + len(others)))])
        x = (x - x[Cf].mean(0)) / one_scale(x, Cf, uF)
        xc, yc = x[Cf].mean(0), y[Cf].mean()
        th_c, Mi, k, root = _control_fit(x, y, Cf, seed, uF)
        TH = th_c + np.sqrt(k) * zth @ root.T                                               # a readout of control cells for every resample
        v_ref = _var_of_mean(y[Cf], uF)
        shift = lambda a: x[a].mean(0) - xc
        row = lambda a: (y[a].mean() - yc, np.sqrt(_var_of_mean(y[a], _units(cl, a)) + v_ref))
        left = y - x @ th_c                                                                 # what the reference readout leaves in a cell
        v_left = _var_of_mean(left[Cf], uF)
        beside = lambda a: np.sqrt(_var_of_mean(left[a], _units(cl, a)) + v_left)            # the error of a row beside its readout: the shift's error is in both and cancels
        Df, (dyf, sef) = shift(Tf), row(Tf)
        arms = [np.where(fit & (unit == h))[0] for h in others]
        keep = [j for j, a in enumerate(arms) if len(a) > 1]
        ro = [row(arms[j]) for j in keep]
        A_own, D_own = _closer(Mi, Df)
        r_oth, s_oth = np.array([r[0] for r in ro]), np.array([beside(arms[j]) for j in keep])
        if keep:
            D_oth = np.array([shift(arms[j]) for j in keep])
            A_oth, D_oth, t2 = _teacher(Mi, k, D_oth, r_oth - D_oth @ th_c, s_oth)
            wide.append(t2)
        F, R = np.r_[Cf, Tf], np.r_[Cs, Ts]
        arm = (unit == event).astype(float)
        mu = x[F].mean(0)
        e_res = arm[R] - arm[F].mean() - (x[R] - mu) @ ridge(x[F] - mu, arm[F] - arm[F].mean(), seed=seed)
        y_res = y[R] - y[F].mean() - (x[R] - mu) @ ridge(x[F] - mu, y[F] - y[F].mean(), seed=seed)
        ey, ee = e_res * y_res, e_res * e_res
        see += float(((arm[R] - arm[R].mean()) ** 2).sum())
        nC = len(Cs)
        sT, sC = WT.sum(1), WC.sum(1)
        DS = (WT @ x[Ts]) / sT[:, None] - (WC @ x[Cs]) / sC[:, None]                         # the shift of the stage on the cells read, for every resample
        den += (WT @ y[Ts]) / sT - (WC @ y[Cs]) / sC
        num["read"] += np.einsum("bp,bp->b", DS, TH)
        num["found_again"] += np.einsum("bp,bp->b", DS, TH + (dyf + sef * zrow[:, 0] - TH @ D_own[0])[:, None] * A_own[:, 0][None, :])
        if others:
            num["others"] += np.einsum("bp,bp->b", DS, TH + (r_oth[None, :] + s_oth[None, :] * zrow[:, 1:1 + len(keep)] - TH @ D_oth.T) @ A_oth.T if keep else TH)
        if o is not None:                                                                    # the output as trained, on the scale of the readout among the reference cells
            a_tr = np.polyfit(_within(o[Cf], uF), _within(y[Cf], uF), 1)[0]
            num["as_trained"] += a_tr * ((WT @ o[Ts]) / sT - (WC @ o[Cs]) / sC)
        sxy += WT @ ey[nC:] + WC @ ey[:nC]
        sxx += WT @ ee[nC:] + WC @ ee[:nC]
        vis.append((float(Df @ DS[0] / (np.linalg.norm(Df) * np.linalg.norm(DS[0]) + 1e-300)), float((Df @ Df) / (x.shape[1] * (1 / len(Tf) + 1 / len(Cf))))))
    D = {k: num[k] / den for k in num}
    overlap = float(sxx[0] / see)
    by = "cells" if overlap >= MIN_OVERLAP else "identity"
    D["carried"] = 1 - (sxy / sxx) / (den / 2) if by == "cells" else np.ones(draws + 1)
    T, C = np.where(unit == event)[0], np.where(ref)[0]
    change = float(y[T].mean() - y[C].mean())
    se = float(np.sqrt(_var_of_mean(y[T], _units(cl, T)) + _var_of_mean(y[C], _units(cl, C))))
    cells = {"event": int(len(T)), "control": int(len(C))}
    if cl is not None:
        cells.update(units_of_the_event=int(len(set(cl[T].tolist()))), units_of_the_reference=int(len(set(cl[C].tolist()))))
    return {"event": str(event), "change": change, "snr": abs(change) / se, "cells": cells, "read_against": "control cells" if control is not None else "the other units",
            "shift_of_the_stage": {"cosine_between_the_halves": float(np.mean([v[0] for v in vis])), "times_its_noise": float(np.mean([v[1] for v in vis]))},
            "overlap": overlap, "carried_by": by, **({"widening_by_the_rows_of_the_others": wide} if wide else {}),
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
FORECAST, NOT_FORECAST = 0.75, 0.50                           # across events: what the other events forecast misses the change by |1 - others|. They forecast it if the interval of
                                                               # `others` lies within a quarter of the change on either side, and do not if it lies beyond a half on one side


def across_events(r):
    """What the rows of the other events say of this event's change through a stage (None if there are none)."""
    if "others" not in r["shares"]:
        return None
    o = r["shares"]["others"]
    lo, hi = o["ci"]
    says = ("the other events forecast it" if (lo >= FORECAST and hi <= 2 - FORECAST) else
            "the other events do not forecast it" if (hi < NOT_FORECAST or lo > 2 - NOT_FORECAST) else "unresolved")
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


def must_fail(xs, y, unit, first, event, control="control", draws=DRAWS, seed=0, cluster=None):
    """Two controls of the measure. A stage with its cells permuted holds nothing: it must carry nothing. An event made
    of reference cells has no change: it must have no change to account for."""
    rng = np.random.default_rng(seed + 11)
    perm = rng.permutation(len(y))
    xp = tuple(x[perm] for x in xs) if isinstance(xs, (tuple, list)) else xs[perm]
    p = reading(xp, y, unit, first, event, control=control, draws=draws, seed=seed, cluster=cluster)
    ref = reference(unit, event, control)
    C = np.where(ref)[0]
    fake = rng.choice(C, min(int((unit == event).sum()), len(C) // 3), replace=False)
    un = unit.astype(object).copy()
    un[fake] = "made-up"
    keep = ref | (un == "made-up")
    cl = None
    if cluster is not None:
        cl = np.asarray(cluster).astype(object).copy()
        cl[fake] = "made-up"
        cl = cl[keep].astype(str)
    sub = (lambda x: x[keep])
    m = reading(tuple(sub(x) for x in xs) if isinstance(xs, (tuple, list)) else sub(xs), y[keep], un[keep].astype(str), first[keep], "made-up", control=control, draws=draws, seed=seed, cluster=cl)
    return {"permuted_stage": {k: p["shares"][k] for k in ("read", "carried", "found_again")},
            "made_up_event": {"change": m["change"], "snr": m["snr"], "no_change_to_account_for": bool(m["snr"] < MIN_SNR)}}


def credit(r, key="found_again"):
    """A repair that reads a stage is credited up to what the stage carries (LOOP_PROTOCOL.md, 4b). r: a reading of the
    stage; key: the share the repair closed on cells it was not fitted on. Returns the share that is credited and what
    the repair claims beyond what the stage carries; a step beyond means the repair has followed the mean."""
    d = r["_draws"]
    c, x = np.minimum(d[key], d["carried"]), d[key] - d["carried"]
    return {"credited": {"share": float(c[0]), "ci": interval(c[1:])},
            "beyond_what_the_stage_carries": {"share": float(x[0]), "ci": interval(x[1:]), "reading": step(interval(x[1:]))}}


ROUTE = {   # the open part of an account -> (lever of the outer step, what is proposed, how the proposal is judged)
    "read": (None, "nothing: a readout with no record of the event reads the change", None),
    "held_not_read": ("M: the head", "a readout taught by the records of the events of the same type",
                      "what it reads on cells it was not fitted on, credited up to what the state carries (`credit`)"),
    "held_as_identity": ("E: records", "the state holds the change only as the identity of the event: its own record answers it, which is the graph's memory and not the model",
                         "whether it is more than the identity is asked of the other events (`others`)"),
    "lost_before_the_state": ("M: the model before the state", "a state kept by type; if that does not bring it back, training of what comes before the state (the owner's decision)",
                              "what the candidate state carries against the mean state and against a matched control (`gain`)"),
    "not_in_the_input, forecast": ("E: records", "the records of events of the same type teach the readout: ask for them by type",
                                   "the forecast of the event left out, through the state (`others`, `Law.forecast`)"),
    "not_in_the_input, not forecast": ("G: the graph", "an edge between the event and the readout: declared if it follows from the declarations, else proposed (`events.edge`)",
                                       "rows the proposal did not see, forecast without the edge (`events.test_edges`); its size from its class and its scale (`Law`)"),
    "not_in_the_input, undecided": ("E: records", "more recorded events of the type, so that it can be said whether they forecast the change",
                                    "the interval of `others` at the state"),
}


def plan(acc, state=None):
    """What a round does next for one account: for every open part (and every part still unresolved), the lever, the
    proposal, how it is judged, and what the readings already made credit it with. This is the step after a rejection:
    the account names the place, the place names the lever, and the proposal is judged on cells or rows that did not
    make it. Nothing here fits anything."""
    ch = acc["chain"]
    state = state or ch["state"]
    out = []
    taught, own = acc.get("a_readout_taught_by_the_other_events"), acc.get("a_repair_that_follows_the_mean")
    for group in ("open", "unresolved"):
        for it in ch[group]:
            part, key = it["part"], it["part"]
            now = None
            if part == "held_not_read":
                if ch["carried_by"].get(state) == "identity":
                    key = "held_as_identity"
                    now = {"its_own_record": own["credited"]} if own else None
                else:
                    now = {"a_readout_taught_by_the_other_events": taught["credited"]} if taught else None
            elif part == "lost_before_the_state":
                g = {k: v["carried"] for k, v in acc.get("gains", {}).items() if k.endswith(f" over {state}") and "carried" in v}
                now = {"candidate_states": g} if g else None
            elif part == "not_in_the_input":
                a = (ch.get("across_events") or {}).get("reading")
                key = "not_in_the_input, " + ("forecast" if a == "the other events forecast it" else "not forecast" if a == "the other events do not forecast it" else "undecided")
                now = {"forecast_by_the_other_events": ch["across_events"]} if ch.get("across_events") else None
            lever, proposal, judged = ROUTE[key]
            out.append({"part": part, "resolved": group == "open", "share": it["share"], "ci": it["ci"], "reading": key, "lever": lever, "proposal": proposal, "judged_by": judged, "credited_now": now})
    return out


def account(stages, y, unit, first, event, others=(), inputs=(), state=None, o=None, pairs=(), control="control", draws=DRAWS, seed=0, controls=True, cluster=None):
    """The account of one event's change of one readout: the readings at every stage, the chain, the gains of the named
    pairs of stages (b over a), what a repair that follows the mean is credited with at the state, and the two controls
    of the measure at the state. stages: {name: array, or the two arrays of a stage chosen on the fitted half}.
    control=None and cluster: an event with no control cells, and arms made of units (see `reading`).
    The resamples are left out of what is returned."""
    state = state or list(stages)[-1]
    R = {n: reading(xs, y, unit, first, event, others=others, control=control, o=o if n == state else None, draws=draws, seed=seed, cluster=cluster) for n, xs in stages.items()}
    out = {"stages": {n: {k: v for k, v in r.items() if k != "_draws"} for n, r in R.items()}, "chain": chain(R, [i for i in inputs if i in R], state),
           "gains": {f"{b} over {a}": {k: gain(R[a], R[b], k) for k in ("carried", "read", "others", "found_again") if k in R[a]["_draws"] and k in R[b]["_draws"]}
                     for a, b in pairs if a in R and b in R},
           "a_repair_that_follows_the_mean": credit(R[state])}
    if "others" in R[state]["_draws"]:
        out["a_readout_taught_by_the_other_events"] = credit(R[state], "others")
    out["plan"] = plan(out, state)
    if controls:
        out["must_fail"] = must_fail(stages[state], y, unit, first, event, control=control, draws=draws, seed=seed, cluster=cluster)
    return out


def run_spec(spec, base=".", draws=DRAWS):
    """The account of every event of a SPEC (the file format of `anm_diag scope`). Read from the SPEC: keys.stages,
    keys.readouts (the first; of an array with one column per readout, the column `account.column`), keys.unit,
    keys.control (left out for a study of units in which no cell is a control cell: every unit is then read against the
    other units), keys.cluster (the donor of every cell, where an arm is made of donors), keys.frozen (the system's own
    output, if any), events, and under
    `account`: state (default Z), inputs (default V), pairs, others (default: the events of the SPEC), seed, and
    unit_length (the stages put at unit length, as a head reads them; default Z)."""
    from pathlib import Path
    keys, acc = spec["keys"], spec.get("account") or {}
    z = np.load(Path(base) / spec["arrays"]["audit"], allow_pickle=False)
    unit_ = np.asarray(z[keys["unit"]]).astype(str)
    ctl = keys.get("control")
    if ctl:
        is_ctl = np.asarray(z[ctl["key"]]).astype(str) == str(ctl["equals"])
        un, control = np.where(is_ctl, "control", unit_), "control"
    else:
        is_ctl, un, control = np.zeros(len(unit_), bool), unit_, None
    cluster = np.asarray(z[keys["cluster"]]).astype(str) if keys.get("cluster") else (un if control is None else None)
    ro = keys["readouts"][0] if isinstance(keys["readouts"], (list, tuple)) else keys["readouts"]
    col = lambda v: v[:, int(acc["column"])] if v.ndim == 2 else v
    y = col(np.asarray(z[ro], np.float64))
    o = col(np.asarray(z[keys["frozen"]["array"]], np.float64)) if keys.get("frozen") else None
    at_unit_length = lambda x: x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-12)
    stages = {n: at_unit_length(np.asarray(z[k], np.float64)) if n in acc.get("unit_length", ["Z"]) else np.asarray(z[k], np.float64) for n, k in keys["stages"].items()}
    seed = int(acc.get("seed", (spec.get("decision") or {}).get("seed", 0)))
    first = cycle.halves(un, seed)
    ev = spec.get("events", "each")
    events = [str(e) for e in (sorted(set(un.tolist()) - {"control"} - set(map(str, spec.get("units_exclude", [])))) if ev in ("each", "pooled") else ev)]
    state = acc.get("state", "Z" if "Z" in stages else list(stages)[-1])
    inputs = acc.get("inputs", [n for n in ("V",) if n in stages])
    oth = [str(e) for e in acc.get("others", events)]
    res = {"name": spec.get("name"), "readout": spec.get("readout_name", ro), "state": state, "inputs": inputs, "cells": int(len(un)), "control_cells": int(is_ctl.sum()),
           "read_against": "control cells" if control else "the other units", "events": {}}
    for e in events:
        res["events"][e] = account(stages, y, un, first, e, others=[h for h in oth if h != e], inputs=inputs, state=state, o=o,
                                   pairs=[tuple(p) for p in acc.get("pairs", [])], control=control, draws=draws, seed=seed, cluster=cluster)
    return res


def main(argv=None):
    """python -m anm_diag account SPEC.yaml --out DIR [--draws N]"""
    import argparse
    import json
    from pathlib import Path

    import yaml
    ap = argparse.ArgumentParser(prog="anm_diag account")
    ap.add_argument("spec")
    ap.add_argument("--out", required=True)
    ap.add_argument("--draws", type=int, default=DRAWS)
    a = ap.parse_args(argv)
    spec = yaml.safe_load(open(a.spec))
    res = run_spec(spec, Path(a.spec).resolve().parent, a.draws)
    Path(a.out).mkdir(parents=True, exist_ok=True)
    (Path(a.out) / "account.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    for e, r in res["events"].items():
        st = r["stages"][res["state"]]["shares"]
        print(f"{e:10s} change {r['chain']['change']:+.3f} | at {res['state']}: read {st['read']['share']:.0%}, carried {st['carried']['share']:.0%}"
              + (f", others {st['others']['share']:.0%}" if "others" in st else "") + f" | {r['chain']['says']}")
    return 0


# --- many events at once: the law of their rows, and the edges ---------------------------------------------------------
MEDIAN_Z2 = 0.4549364                                          # the median of the square of a standard normal
OWN_PART, OWN_SHARE = 1.0, 0.05                                # a direct part of a pair's own: of the order of one spread of the readout; the share of pairs that have one, to begin with


MIN_CONFORMAL, MIN_CONFORMAL_EVENTS = 40, 8                    # the recorded pairs, and the events they come from, that a kind needs before its own pairs set the width of its interval


def _mixture_half(bulk, wide, share, alpha=ALPHA):
    """The half-width c of the central (1 - alpha) interval of a mixture of two normals with mean zero: variance `bulk`
    with weight 1 - share and variance `wide` with weight share. (1 - share) P(|N(0, bulk)| > c) + share P(|N(0, wide)| > c)
    = alpha."""
    from math import erfc, sqrt
    out = lambda c: (1 - share) * erfc(c / sqrt(2 * bulk)) + share * erfc(c / sqrt(2 * wide)) - alpha
    lo, hi = 0.0, 10.0 * sqrt(wide)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if out(mid) > 0 else (lo, mid)
    return 0.5 * (lo + hi)


def _calibrated_spread(r2, base):
    """The spread s at which the median of r2 / (base + s^2) is that of the square of a standard normal: how wide a direct
    part has to be for rows left out to stand as far as they do. 0 where they are already as wide as they should be.
    The median does not move with a few rows that stand far out."""
    f = lambda s2: np.median(r2 / (base + s2 + 1e-300)) - MEDIAN_Z2
    if len(r2) == 0 or f(0.0) <= 0:
        return 0.0
    lo, hi = 0.0, float(r2.max() / MEDIAN_Z2)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return float(np.sqrt(0.5 * (lo + hi)))


class Law:
    """Every event of a library and every readout, from rows. In spreads of a readout among control cells:

        change(event, readout) = shift(event) . readout  +  direct(event, readout)  +  error

    readout: a prior from the control cells (their ridge readout and its uncertainty, widened by one factor that the rows
    choose), updated by the rows of the other events. direct: zero with a small spread and heavy tails for a pair with no
    edge; for a pair with an edge of class k a mean and a spread of its own, learned from the other pairs of that class.
    edges: {(event, readout): class} (`events.edge_table`). A class with fewer than `min_class` recorded pairs gives no
    forecast of its size.
    types: {event: type} from the event graph. The error of a readout is then a property of a type of event: the rows of
    an event update the readout for the events of its type and for no other, and an event with no type is taught by no
    row. How far the control prior is widened is chosen by type as well, by the rows of the type. For an event with no
    type the control prior stands as the control cells give it: what such an event shows beyond it is its direct part,
    and is flagged; a widening that the events of no type chose together would hide it, and would be the teaching the
    rule excludes. None: every event is of one type (one readout error, one widening for all).
    scale: {(event, readout): s}. The direct part of an edge is its class's share times s, so that one class holds in
    conditions where the readout has more or less to lose (for the knockout of a protein's own gene: how far the
    readout stands above its background among control cells). None: s = 1, one size per class.
    control=None: no cell is a control cell. Every unit is an event, read against the mean of the units: the readout
    "with no record" is one of the variation inside a unit, the spread of a readout is its spread inside a unit, and
    the rows of the units teach each other (one type, unless `types` says otherwise).
    The spread of a direct part with no edge is the one at which the recorded pairs with no edge and no type, each left
    out, are as wide as they should be (the median of their squared z is that of a standard normal): a few pairs with
    a large direct part of their own do not widen it. Such a pair is an edge the graph does not hold yet: left out, it
    stands where a direct part of the order of the readout's own spread would put it and not where a pair with none
    would be, and by that much its row teaches the readout less (`OWN_PART`).
    known: {(event, readout): value}. A part of the direct part that is known by declaration and not asked of the rows
    (what the arithmetic of a readout gives: `events.derived`). It is taken out of a recorded row before the row
    teaches anything, and put back into its forecast. None: nothing is known."""

    def __init__(self, x, Y, unit, control="control", edges=None, nu=4.0, min_cells=30, min_class=3, seed=0, types=None, scale=None, known=None):
        self.events = [str(e) for e in sorted(set(unit.tolist())) if e != control and (unit == e).sum() >= min_cells]
        self.units_only = control is None
        if self.units_only:                                         # no control cells: the reference is the inside of every unit, and the mean of the units
            keep = np.isin(unit, self.events)
            x, unit, Y = x[keep], unit[keep], {k: v[keep] for k, v in Y.items()}
            inv = np.unique(unit, return_inverse=True)[1]
            C = np.arange(len(unit))
            mean_of_units = lambda v: np.mean([v[unit == e].mean(0) for e in self.events], axis=0)
            self._centre, self._scale = mean_of_units(x), one_scale(x, C, inv)
        else:
            C, inv = np.where(unit == control)[0], None
            self._centre, self._scale = x[C].mean(0), one_scale(x, C)
        X = (x - self._centre) / self._scale
        self.readouts, self.nu, self.min_class = list(Y), nu, min_class
        self.n = np.array([(unit == e).sum() for e in self.events])
        self.D = np.array([X[unit == e].mean(0) for e in self.events])
        nc, p = len(C), X.shape[1]
        self.times_noise = (self.D ** 2).sum(1) / (p * (1 / self.n + 1 / nc))
        self.typed = types is not None
        self.type = ["(all)" if types is None else str(types.get(e, f"({e})")) for e in self.events]
        self.declared = np.array([types is None or e in types for e in self.events])       # the events whose rows teach others: those of a declared type (all, in a law with one type)
        same = np.array([[a == b for b in self.type] for a in self.type])
        by = {}
        for i, (e, t) in enumerate(zip(self.events, self.type)):
            by.setdefault(t if (types is None or e in types) else None, []).append(i)
        # the events whose rows choose a widening together: {type: its events}; an event in no group keeps the control prior as it is
        self.groups = {t: np.array(idx) for t, idx in by.items() if t is not None}
        Xc = _within(X, inv) if self.units_only else X[C]
        G = Xc.T @ Xc / nc
        w, V = np.linalg.eigh(G)
        self._nc, self._w, self._V = nc, w, V
        fold = np.random.default_rng(seed).permutation(nc) % 5
        self.r, parts = {}, []
        for k in range(5):                                          # the folds of the control cells are the same for every readout
            m = fold == k
            wk, Vk = np.linalg.eigh((G * nc - Xc[m].T @ Xc[m]) / (nc - m.sum()))
            parts.append((m, wk, Vk, Xc[m] @ Vk, Xc[~m]))
        for name, y in Y.items():
            centre = float(mean_of_units(y)) if self.units_only else float(y[C].mean())
            yc = _within(y, inv) if self.units_only else y[C] - centre
            sd = yc.std()
            err = np.zeros(len(LAMBDAS))
            for m, wk, Vk, Fv, Xr in parts:
                c = Vk.T @ (Xr.T @ yc[~m] / (nc - m.sum()))
                for i, lam in enumerate(LAMBDAS):
                    err[i] += ((yc[m] - Fv @ (c / (wk + lam))) ** 2).sum()
            lam, cv = LAMBDAS[int(err.argmin())], float(err.min() / nc)
            th = V @ ((V.T @ (Xc.T @ yc / nc)) / (w + lam))
            DV = self.D @ V
            K0 = same * ((cv / nc) * (DV / (w + lam)) @ DV.T / sd ** 2)                    # the rows of another type do not bear on this type's readout
            dy = np.array([y[unit == e].mean() - centre for e in self.events])
            left = yc - Xc @ th if self.units_only else y - centre - X @ th                 # what the readout leaves in a cell
            own = lambda v: np.array([v[unit == e].var(ddof=1) / (unit == e).sum() for e in self.events])
            # the error of the reference: of the mean of the control cells, or of the mean of the units
            ref_y, ref_left = (own(y).sum() / len(self.events) ** 2, own(left).sum() / len(self.events) ** 2) if self.units_only else (y[C].var(ddof=1) / nc, cv / nc)
            se2 = own(y) + ref_y
            # the error of a row beside its readout: the change and the shift are means over the same cells, so the error of the
            # shift is in both and cancels; what is left is the mean of what the readout leaves in the cells of the two arms
            noise = (own(left) + ref_left) / sd ** 2
            cls = np.array([(edges or {}).get((e, name)) for e in self.events], dtype=object)
            sc = np.array([float((scale or {}).get((e, name), 1.0)) for e in self.events])
            kn = np.array([float((known or {}).get((e, name), 0.0)) for e in self.events])
            self.r[name] = {"sd": float(sd), "cv": cv, "lam": float(lam), "K0": K0, "change": dy / sd, "se": np.sqrt(se2) / sd, "noise": noise, "scale": sc, "known": kn,
                            "read": self.D @ th / sd, "cls": cls, "tau2": 1.0, "tau": np.ones(len(dy)), "w": np.ones(len(dy)), "own": np.zeros(len(dy)), "th": th, "thGth": float(th @ G @ th), "centre": centre}
        self.classes = sorted({c for q in self.r.values() for c in q["cls"] if c is not None})
        self.sigma0, self.cl, self.fixed, self.own_share = None, {}, False, OWN_SHARE

    def _res(self, q):
        """What a recorded row shows beyond the control readout and beyond what is known by declaration."""
        return q["change"] - q["read"] - q["known"]

    def _prior(self, q):
        """The readout's prior among the recorded rows: the control prior, widened by the factor the rows of each type chose."""
        return q["tau2"] * q["K0"] if not self.typed else np.sqrt(np.outer(q["tau"], q["tau"])) * q["K0"]

    def _widening(self, q, t):
        """The widening for an event of type `t`: its type's, where the rows of that type chose one; else 1, the control
        prior as it is."""
        return q["tau2"] if not self.typed else q["tau2"].get(t, 1.0)

    def _var_direct(self, q):
        """The spread of the direct part of every recorded pair, as its row is weighed: a pair with an edge by its class; a
        pair with none by the spread with no edge, widened by how likely it is to hold a direct part of its own."""
        s2 = np.array([self.sigma0 ** 2 if c is None else (self.cl[c]["sigma"] * s) ** 2 for c, s in zip(q["cls"], q["scale"])]) / q["w"]
        n, g = q["noise"], q["own"]
        mixed = 1 / ((1 - g) / (self.sigma0 ** 2 + n) + g / (OWN_PART ** 2 + n)) - n
        return np.where([c is None for c in q["cls"]], mixed, s2)

    def _mean_direct(self, q, without=None):
        return np.array([0.0 if c is None else self.cl[c]["mu"] * s for c, s in zip(q["cls"], q["scale"])])

    def fit(self, iters=40, sigma0=None, fallback=None):
        """The spread of a direct part with no edge, the classes and the widening of the control prior, from the rows.
        `sigma0`: hold the spread at this value. `fallback`: the value to hold it at where this law has fewer than
        MIN_ROWS recorded pairs with no edge and no type (a stratum with too few such rows takes it from a larger
        library). The widening is still chosen by the rows."""
        res = np.concatenate([self._res(q) for q in self.r.values()])
        cls = np.concatenate([q["cls"] for q in self.r.values()])
        sc = np.concatenate([q["scale"] for q in self.r.values()])
        none = np.array([c is None for c in cls])
        res = np.where(none, res, res / sc)                          # a pair with an edge: its direct part as a share of its scale
        plain = none & np.tile(~self.declared if self.typed else np.ones(len(self.events), bool), len(self.r))     # no edge, and taught by no row of a type
        if sigma0 is None and fallback is not None and plain.sum() < MIN_ROWS:
            sigma0 = fallback
        if plain.sum() < MIN_ROWS:
            plain = none
        self.fixed = sigma0 is not None
        self.spread_from = "given" if self.fixed else f"{int(plain.sum())} recorded pairs with no edge" + (" and no type" if self.typed and plain.sum() < none.sum() else "")
        self.sigma0 = float(sigma0) if self.fixed else float(1.4826 * np.median(np.abs(res[plain] - np.median(res[plain])))) + 1e-6
        edge_res = res[~none]
        pooled = (float(edge_res.mean()), float(edge_res.std(ddof=1))) if len(edge_res) >= self.min_class else (0.0, float("nan"))
        for c in self.classes:
            v = res[cls == c]
            self.cl[c] = {"n": int(len(v)), "mu": float(v.mean()) if len(v) >= self.min_class else pooled[0],
                          "sigma": float(v.std(ddof=1)) if len(v) >= self.min_class else pooled[1], "own": bool(len(v) >= self.min_class)}
        taus = 10.0 ** np.arange(0, 4.01, 0.1)
        for _ in range(iters):
            Ed2, Ed, W, R2, BASE = [], [], [], [], []
            for q in self.r.values():
                Vd, rho = self._var_direct(q), self._res(q) - self._mean_direct(q)
                Vd = np.where(np.isfinite(Vd), Vd, 1e6)
                chosen = {}
                for g, b in self.groups.items():
                    best = None
                    for t2 in taus:
                        L = np.linalg.cholesky(t2 * q["K0"][np.ix_(b, b)] + np.diag((Vd + q["noise"])[b]))
                        a = np.linalg.solve(L, rho[b])
                        nll = 0.5 * a @ a + np.log(np.diag(L)).sum()
                        if best is None or nll < best[0]:
                            best = (nll, t2)
                    chosen[g] = float(best[1])
                    q["tau"][b] = best[1]
                q["tau2"] = chosen if self.typed else chosen["(all)"]
                Si = np.linalg.inv(self._prior(q) + np.diag(Vd + q["noise"]))
                m = Vd * (Si @ rho)                                   # the direct part beyond its class mean, given every row
                v = Vd - Vd * np.diag(Si) * Vd
                q["direct_dev"], q["direct_var"] = m, v
                Ed.append(m + self._mean_direct(q)); Ed2.append(m * m + v); W.append(q["w"])
                R2.append(((Si @ rho) / np.diag(Si)) ** 2); BASE.append(np.maximum(1 / np.diag(Si) - Vd, 0.0))     # each pair left out: how far it stands, and its spread before its own direct part
            Ed, Ed2, W, R2, BASE = np.concatenate(Ed), np.concatenate(Ed2), np.concatenate(W), np.concatenate(R2), np.concatenate(BASE)
            Ed, Ed2 = np.where(none, Ed, Ed / sc), np.where(none, Ed2, Ed2 / sc ** 2)
            s2 = np.array([self.sigma0 ** 2 if c is None else self.cl[c]["sigma"] ** 2 for c in cls])
            s2 = np.maximum(np.where(np.isfinite(s2), s2, 1e6), 1e-18)
            Wn = (self.nu + 1) / (self.nu + Ed2 / s2)
            if not self.fixed:
                self.sigma0 = _calibrated_spread(R2[plain], BASE[plain]) + 1e-9
            # how likely a pair with no edge holds a direct part of its own: where it stands, left out, under the two spreads
            dens = lambda var: np.exp(-0.5 * R2 / var) / np.sqrt(var)
            a, b = self.own_share * dens(BASE + OWN_PART ** 2), (1 - self.own_share) * dens(BASE + self.sigma0 ** 2)
            own = np.where(none, a / (a + b + 1e-300), 0.0)
            self.own_share = float(np.clip(own[plain].mean(), 0.005, 0.5))
            for c in self.classes:
                m = cls == c
                if self.cl[c]["own"]:
                    self.cl[c]["mu"] = float((Wn[m] * Ed[m]).sum() / Wn[m].sum())
                    self.cl[c]["sigma"] = float(np.sqrt(max((Wn[m] * Ed2[m]).mean(), 1e-12)))
            k = 0
            for q in self.r.values():
                q["w"] = np.where(none[k:k + len(q["w"])], 1.0, Wn[k:k + len(q["w"])])
                q["own"] = own[k:k + len(q["w"])]
                k += len(q["w"])
        return self

    def forecast(self):
        """Every pair with its event left out: the part a readout of control cells reads, what the rows of the other events
        add through the stage, the part by its edge, the forecast with its spread, and how far the measured change is."""
        rows = []
        dev = {}                                                    # the direct parts of the pairs of a class, as fitted, as shares of their scale: (event, value)
        for name, q in self.r.items():
            for i, c in enumerate(q["cls"]):
                if c is not None:
                    dev.setdefault(c, []).append((self.events[i], float(q["direct_dev"][i] / q["scale"][i] + self.cl[c]["mu"])))
        t95 = 1.96
        for name, q in self.r.items():
            Vd = np.where(np.isfinite(self._var_direct(q)), self._var_direct(q), 1e6)
            V = Vd + q["noise"]
            rho = self._res(q) - self._mean_direct(q)
            Si = np.linalg.inv(self._prior(q) + np.diag(V))
            m_loo = rho - (Si @ rho) / np.diag(Si)                  # what the rows of the other events add through the stage
            v_loo = np.maximum(1 / np.diag(Si) - V, 0.0)
            for i, e in enumerate(self.events):
                c = q["cls"][i]
                if c is None:
                    edge, sd_d, why = 0.0, self.sigma0, None
                else:
                    oth = [v for h, v in dev[c] if h != e]
                    if len(oth) >= self.min_class - 1 and len(oth) >= 2:
                        edge, sd_d, why = float(np.mean(oth) * q["scale"][i]), float(np.std(oth, ddof=1) * np.sqrt(1 + 1 / len(oth)) * q["scale"][i]), None
                    else:
                        edge, sd_d, why = float("nan"), float("nan"), f"an edge of class '{c}' with {len(oth)} other recorded pairs: its size cannot be forecast"
                f = q["read"][i] + q["known"][i] + m_loo[i] + (edge if edge == edge else 0.0)
                sd = float(np.sqrt(v_loo[i] + (sd_d ** 2 if sd_d == sd_d else 0.0) + q["noise"][i]))
                z = float((q["change"][i] - f) / sd) if why is None else float("nan")
                rows.append({"event": e, "readout": name, "cells": int(self.n[i]), "edge": c, "type": self.type[i] if self.typed else None, "change": float(q["change"][i]), "se": float(q["se"][i]),
                             "shift_times_noise": float(self.times_noise[i]), "read": float(q["read"][i]), "from_the_other_events": float(m_loo[i]),
                             "by_its_edge": edge, "forecast": float(f), "sd": sd, "interval": [float(f - t95 * sd), float(f + t95 * sd)], "z": z, "no_forecast": why,
                             "known_by_declaration": float(q["known"][i]),
                             "direct": float(q["change"][i] - q["read"][i] - q["known"][i] - m_loo[i]), "direct_sd": float(np.sqrt(v_loo[i] + q["noise"][i]))})
        return rows

    def _new_noise(self, q, n):
        """The error of a row not yet measured, of an event with n cells: what the readout leaves among control cells."""
        return (1 / n + (1 / (n * len(self.events)) if self.units_only else 1 / self._nc)) * q["cv"] / q["sd"] ** 2

    def conformal(self, alpha=ALPHA, minimum=MIN_CONFORMAL, events=MIN_CONFORMAL_EVENTS):
        """How wide the interval of a pair with no edge has to be, from the recorded pairs with no edge themselves. Each is
        left out and forecast as a new event is (its row's error taken at what the readout leaves among control cells);
        how far it then stands, in spreads of that forecast, is its score. A new pair that is exchangeable with the
        recorded pairs of its kind exceeds the returned quantile with probability `alpha` at most, whatever the shape
        of their distribution. Nothing about that shape is assumed. The pairs of one event are not independent of each
        other (they share its cells and its place in the stage), so a kind also needs pairs from `events` events.
        {"no type" | "of a type" | "any": {"quantile", "pairs", "events"}}; the quantile is None where a kind has too few."""
        score, kind, who = [], [], []
        for q in self.r.values():
            Vd = np.where(np.isfinite(self._var_direct(q)), self._var_direct(q), 1e6)
            V = Vd + q["noise"]
            rho = self._res(q) - self._mean_direct(q)
            Si = np.linalg.inv(self._prior(q) + np.diag(V))
            left_out = (Si @ rho) / np.diag(Si)
            v_loo = np.maximum(1 / np.diag(Si) - V, 0.0)
            as_new = np.array([self._new_noise(q, n) for n in self.n])
            none = np.array([c is None for c in q["cls"]])
            score += list((np.abs(left_out) / np.sqrt(v_loo + self.sigma0 ** 2 + as_new))[none])
            kind += list(np.where(self.declared, "of a type", "no type")[none])
            who += list(np.arange(len(self.events))[none])
        score, kind, who, out = np.array(score), np.array(kind), np.array(who), {}
        for k in ("no type", "of a type", "any"):
            m = np.ones(len(score), bool) if k == "any" else kind == k
            v, n_events = np.sort(score[m]), int(len(set(who[m].tolist())))
            rank = int(np.ceil((len(v) + 1) * (1 - alpha)))
            out[k] = {"quantile": float(v[rank - 1]) if len(v) >= minimum and n_events >= events and rank <= len(v) else None, "pairs": int(len(v)), "events": n_events}
        return out

    def class_values(self):
        """The direct parts of the recorded pairs of every class, as fitted, each as a share of its scale:
        {class: [(event, readout, value)]}."""
        out = {}
        for name, q in self.r.items():
            for i, c in enumerate(q["cls"]):
                if c is not None:
                    out.setdefault(c, []).append((self.events[i], name, float(q["direct_dev"][i] / q["scale"][i] + self.cl[c]["mu"])))
        return out

    def new(self, x, unit, edges=None, classes=None, types=None, scale=None, shares=None, known=None, width=None):
        """Events with no measured change, from the stage of their cells alone: for every pair the part a readout of
        control cells reads, what the rows of the recorded events add through the stage, the part by its edge, and the
        forecast with its spread. `x`, `unit`: the stage and the event of the new cells. `classes`: {class: [values]}, the
        direct parts of recorded pairs of a class as shares of their scale (this law's own when not given; a library
        with no pair of a class takes them from one that has). types: {event: type} of the new events (in a law with
        types, the rows that teach are those of the event's type; an event with no type is taught by none).
        scale: {(event, readout): s} of the new pairs. shares: {(event, readout): (share, its spread)}, the pair's own
        share of its scale where the pair was measured elsewhere (the same knockout in another condition): the edge is
        then forecast at that share of this condition's scale, and its class is not asked. The error of a row not yet
        measured is taken at what the readout leaves among control cells.
        known: {(event, readout): (value, its spread)}, a part of the direct part known by declaration (`events.derived`).
        width: how the 95% interval of a pair with no edge is formed. None: 1.96 spreads of the forecast, as up to
        version 5. Else the result of `conformal()` (this law's, or a library's where this law has too few pairs): the
        quantile of the pair's kind, else of any kind; where neither exists, the central 95% of the mixture the law
        itself assumes (most pairs with the spread with no edge, a share `own_share` with a direct part of their own),
        and the row says that its width is assumed. A pair with an edge keeps 1.96 spreads."""
        edges = edges or {}
        vals = {c: [v for _e, _r, v in rows] for c, rows in self.class_values().items()}
        vals.update(classes or {})
        X = (x - self._centre) / self._scale
        out = []
        for e in sorted(set(unit.tolist())):
            m = unit == e
            n, d = int(m.sum()), X[m].mean(0)
            dV = d @ self._V
            tn = "(all)" if not self.typed else str((types or {}).get(str(e), f"({e})"))
            of_type = np.array([t == tn for t in self.type], float)
            for name, q in self.r.items():
                DV = self.D @ self._V
                k = of_type * ((q["cv"] / self._nc) * (dV / (self._w + q["lam"])) @ DV.T / q["sd"] ** 2)
                knn = (q["cv"] / self._nc) * float((dV / (self._w + q["lam"])) @ dV) / q["sd"] ** 2
                Vd = np.where(np.isfinite(self._var_direct(q)), self._var_direct(q), 1e6)
                Si = np.linalg.inv(self._prior(q) + np.diag(Vd + q["noise"]))
                rho = self._res(q) - self._mean_direct(q)
                t2 = self._widening(q, tn)
                add = float(t2 * k @ Si @ rho)
                var = max(float(t2 * knn - t2 ** 2 * k @ Si @ k), 0.0)
                c = edges.get((str(e), name))
                why = None
                sn = float((scale or {}).get((str(e), name), 1.0))
                own = (shares or {}).get((str(e), name))
                if c is None:
                    edge, sd_d = 0.0, self.sigma0
                elif own is not None:
                    edge, sd_d = float(own[0] * sn), float(abs(own[1] * sn))
                elif len(vals.get(c, [])) >= max(self.min_class - 1, 2):
                    edge, sd_d = float(np.mean(vals[c]) * sn), float(np.std(vals[c], ddof=1) * np.sqrt(1 + 1 / len(vals[c])) * abs(sn))
                else:
                    edge, sd_d, why = float("nan"), float("nan"), f"an edge of class '{c}' with {len(vals.get(c, []))} recorded pairs: its size cannot be forecast"
                read = float(d @ q["th"] / q["sd"])
                noise = self._new_noise(q, n)
                kv, ks = (known or {}).get((str(e), name), (0.0, 0.0))
                f = read + add + kv + (edge if edge == edge else 0.0)
                sd = float(np.sqrt(var + (sd_d ** 2 if sd_d == sd_d else 0.0) + noise + ks ** 2))
                half, how = 1.96 * sd, "1.96 spreads of the forecast"
                if width is not None and c is None:
                    k = "of a type" if (self.typed and types is not None and str(e) in types) or not self.typed else "no type"
                    use = next((j for j in (k, "any") if width.get(j, {}).get("quantile") is not None), None)
                    if use is not None:
                        half, how = width[use]["quantile"] * sd, width[use].get("from") or f"the {width[use]['pairs']} recorded pairs with no edge ({use}), each left out"
                    else:
                        half, how = _mixture_half(sd ** 2, sd ** 2 - sd_d ** 2 + OWN_PART ** 2, self.own_share), "assumed: the mixture of the law, for want of recorded pairs"
                out.append({"event": str(e), "readout": name, "cells": n, "edge": c, "type": tn if self.typed else None, "taught_by": int(of_type.sum()), "widening": float(t2), "scale": sn,
                            "size_of_the_edge_from": None if c is None else ("the pair's own share, measured elsewhere" if own is not None else "its class"),
                            "shift_times_noise": float((d ** 2).sum() / (len(d) * (1 / n + 1 / self._nc))),
                            "read": read, "from_the_recorded_events": add, "by_its_edge": edge, "known_by_declaration": float(kv), "forecast": float(f), "sd": sd,
                            "interval": [float(f - half), float(f + half)], "width_from": how, "no_forecast": why,
                            "parts_of_the_spread": {"readout": float(np.sqrt(var)), "direct": float(sd_d), "row": float(np.sqrt(noise)), "known": float(ks)}})
        return out

    def measure(self, y, unit, readout):
        """The change of new events once measured, in the law's units: {event: (change, its error)}."""
        q = self.r[readout]
        out = {}
        for e in sorted(set(unit.tolist())):
            v = y[unit == e]
            ref = v.var(ddof=1) / (len(v) * len(self.events)) if self.units_only else q["sd"] ** 2 / self._nc
            out[str(e)] = (float((v.mean() - q["centre"]) / q["sd"]), float(np.sqrt(v.var(ddof=1) / len(v) + ref) / q["sd"]))
        return out

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
                        "change": r["change"], "forecast": r["forecast"], "z": r["z"], "flag": "a direction the law does not cover" if moved else "a missing edge",
                        "candidate": "a direction of the stage the law does not cover: its record would teach the readout" if moved else "a missing edge: the event does not move the stage"})
        return sorted(out, key=lambda r: -abs(r["errors"]))

    def summary(self):
        return {"events": len(self.events), "readouts": self.readouts, "spread_of_a_direct_part_without_an_edge": self.sigma0, "spread_from": self.spread_from,
                "widening_of_the_control_prior": {k: q["tau2"] for k, q in self.r.items()}, "classes": self.cl,
                "events_whose_shift_stands_out": {"beyond_2x_noise": int((self.times_noise > 2).sum()), "beyond_4x_noise": int((self.times_noise > 4).sum())}}
