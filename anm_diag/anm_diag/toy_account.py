"""A planted system for the account of a change (account.py), where the part of every change is known.

Five hidden factors of a cell. The readout is measured outside the model.

  factor  among control cells   the readout follows it   the input view holds it   the state holds it
  1       varies                yes                      yes                       yes
  2       does not vary         yes                      yes                       yes        (a pathway that is off in control cells)
  3       varies                no                       yes                       yes
  4       varies                yes                      yes                       no         (the model loses it before its state)

  event      acts on                      so its change is                                          and the account should say
  c1 c2      factor 1                     read by a readout of control cells                        read
  r1 r2      factor 2                     held by the state, and no control cell can teach it       held, not read (the head); r1 is forecast from r2
  l          factor 4                     held by the input, lost before the state                  lost before the state (the model)
  m          factor 3, and the readout    direct, though the event moves the state                  direct (an edge). A readout made to follow its
             directly                                                                               mean calls it closed: `found_again` is near one.
  n1..n4     the readout directly         direct; the event graph has an edge for them              direct; forecast by the class of their edge
  x          the readout directly         direct; the event graph has no edge                       direct; flagged as a missing edge
  null       factor 3                     none                                                      no change to account for
  q*         a little of factors 1, 3     none or small, through the state                          (they fill the library the law is fitted on)

A share of the cells of an event is not edited (`edited`), as in a screen. Truth is computed on the same system with a very
large number of cells and exact least squares, so the estimates are compared with what the planted system itself gives.

  python -m anm_diag.toy_account [seeds]            (that many systems as planted; prints the first in full)
  python -m anm_diag.toy_account suite [file]       (twenty as planted and five of every variant)
"""
import sys

import numpy as np

from . import account as AC
from . import cycle

A = np.array([0.6, 0.5, 0.0, 0.6])                             # how the readout follows the four factors
SD = np.array([1.0, 0.0, 1.0, 1.0])                            # the spread of a factor among control cells
EVENTS = {"c1": ((0, -1.5), 0.0), "c2": ((0, -1.0), 0.0), "r1": ((1, -1.5), 0.0), "r2": ((1, -1.2), 0.0), "l": ((3, -1.5), 0.0), "m": ((2, 1.5), -0.75),
          "n1": (None, -1.2), "n2": (None, -0.8), "n3": (None, -0.5), "n4": (None, -1.0), "x": (None, -0.9), "null": ((2, 1.5), 0.0)}
EDGES = ("n1", "n2", "n3", "n4")
TOL = 0.15


def planted(seed=0, n_control=4000, n_event=800, edited=0.75, noise=0.1, fill=30, big=False):
    """Cells, the two stages (V: the input view, 40 columns; Z: the state, 30 columns), the readout, the event of a cell."""
    rng = np.random.default_rng(seed)
    Lv, Lz = np.linalg.qr(rng.normal(size=(40, 4)))[0].T * 2.0, np.linalg.qr(rng.normal(size=(30, 3)))[0].T * 2.0     # the factors lie along directions at right angles
    ev = dict(EVENTS)
    frng = np.random.default_rng(1000 + seed)
    for k in range(fill):
        ev[f"q{k:02d}"] = ((int(frng.choice([0, 2])), float(frng.normal(scale=0.3))), 0.0)
    unit, F, direct = [], [], []
    for name, (shift, d) in [("control", (None, 0.0)), *ev.items()]:
        n = n_control if name == "control" else n_event
        if big:                                               # the same share of control cells as in the system that is read: what a stage carries is defined for the mixture of cells that is read
            n = (30000 * n_control // n_event if name == "control" else 30000) if name in EVENTS or name == "control" else 0
        if not n:
            continue
        f = rng.normal(size=(n, 4)) * SD
        ed = (rng.random(n) < edited).astype(float) if name != "control" else np.zeros(n)
        if shift is not None:
            f[:, shift[0]] += ed * shift[1]
        F.append(f); direct.append(ed * d); unit += [name] * n
    F, direct, unit = np.vstack(F), np.concatenate(direct), np.array(unit)
    y = F @ A + direct + rng.normal(scale=0.7, size=len(F))
    V = F @ Lv + rng.normal(scale=noise, size=(len(F), 40))
    Z = F[:, :3] @ Lz + rng.normal(scale=noise, size=(len(F), 30))
    return {"V": V, "Z": Z, "y": y, "unit": unit, "edges": {(e, "y"): "own" for e in EDGES}, "loadings": {"V": (Lv, [0, 1, 2, 3]), "Z": (Lz, [0, 1, 2])}, "noise": noise, "edited": edited}


_TRUTH = {}


def truth(**kw):
    """The parts of every planted change, on the same system with very many cells and exact least squares. They are
    properties of the planted system and not of a seed; computed once for a set of sizes and kept."""
    key = tuple(sorted(kw.items()))
    if key not in _TRUTH:
        _TRUTH[key] = _truth(**kw)
    return _TRUTH[key]


def _truth(**kw):
    d = planted(0, big=True, **kw)
    y, un, C = d["y"], d["unit"], d["unit"] == "control"
    out = {}
    for e in EVENTS:
        T = un == e
        change = y[T].mean() - y[C].mean()
        row = {"change": float(change)}
        for st in ("V", "Z"):
            x = d[st]
            L, idx = d["loadings"][st]                                # read: exact, from the planted moments of control cells (a factor that does not vary among them teaches them nothing)
            S2 = np.diag(SD[idx] ** 2)
            th = np.linalg.solve(L.T @ S2 @ L + d["noise"] ** 2 * np.eye(L.shape[1]), L.T @ S2 @ A[idx])
            shift, dr = EVENTS[e]
            D = d["edited"] * shift[1] * L[idx.index(shift[0])] if shift is not None and shift[0] in idx else np.zeros(L.shape[1])
            pop = d["edited"] * ((shift[1] * A[shift[0]] if shift is not None else 0.0) + dr)
            read = float(D @ th / pop) if abs(pop) > 0.05 else float("nan")
            m = C | T
            X = np.c_[x[m], T[m].astype(float), np.ones(m.sum())]
            b = np.linalg.lstsq(X, y[m], rcond=None)[0]
            a = T[m].astype(float)
            ov = float(((a - np.c_[x[m], np.ones(m.sum())] @ np.linalg.lstsq(np.c_[x[m], np.ones(m.sum())], a, rcond=None)[0]) ** 2).sum() / ((a - a.mean()) ** 2).sum())
            row[st] = {"read": float(read), "overlap": ov, "carried": (float(1 - b[-2] / change) if ov >= AC.MIN_OVERLAP else 1.0) if abs(change) > 0.05 else float("nan")}
        if abs(change) > 0.05:
            row["parts"] = {"direct": 1 - row["V"]["carried"], "lost_before_the_state": row["V"]["carried"] - row["Z"]["carried"],
                            "held_not_read": row["Z"]["carried"] - row["Z"]["read"], "read": row["Z"]["read"]}
        out[e] = row
    return out


def run(seed=0, draws=200, verbose=False, **kw):
    """One planted system: every check of the suite, each True or False."""
    d, tr = planted(seed, **kw), truth(**kw)
    y, un = d["y"], d["unit"]
    first = cycle.halves(un, seed)
    checks, rep, resolved = {}, {}, {}
    for e in EVENTS:
        m = np.isin(un, ["control", e, "r1", "r2"])
        oth = [h for h in ("r1", "r2") if h != e] if e in ("r1", "r2") else ()
        R = {st: AC.reading(d[st][m], y[m], un[m], first[m], e, others=oth, draws=draws, seed=seed) for st in ("V", "Z")}
        ch = AC.chain(R, ["V"], "Z")
        rep[e] = {"says": ch["says"], "parts": {k: round(v["share"], 2) for k, v in ch["parts"].items()}, "found_again_at_Z": round(R["Z"]["shares"]["found_again"]["share"], 2)}
        if e == "null":
            checks["null: no change to account for"] = ch["says"].startswith("no change")
            continue
        for k in AC.PARTS:
            lo, hi = ch["parts"][k]["ci"]
            checks[f"{e}: the planted {k} is inside the interval (or within {TOL})"] = bool(lo - 0.05 <= tr[e]["parts"][k] <= hi + 0.05 or abs(ch["parts"][k]["share"] - tr[e]["parts"][k]) <= TOL)
        tp = {k: tr[e]["parts"][k] for k in AC.PARTS[:3]}
        top = max(tp.values())
        if top < 0.25:
            checks[f"{e}: nothing open"] = not ch["open"]
        else:                                                   # the account points at the planted place: the first part that is open, or if none is, the first that is unresolved
            ok = [k for k, v in tp.items() if v >= top - TOL]
            first_part = (ch["open"] or ch["unresolved"] or [{"part": None}])[0]["part"]
            checks[f"{e}: the account points at {' or '.join(ok)}"] = bool(first_part in ok)
            resolved[e] = bool(ch["open"] and ch["open"][0]["part"] in ok)
        if e == "r1":
            lo, hi = R["Z"]["shares"]["others"]["ci"]
            checks["r1: its forecast from r2 holds its whole change (one inside the interval, or within 0.15)"] = bool(lo - 0.05 <= 1 <= hi + 0.05 or abs(R["Z"]["shares"]["others"]["share"] - 1) <= TOL)
        if e == "m":
            checks["m: a readout made to follow its mean calls it closed (found again >= 0.5) while the state carries < 0.2"] = bool(
                R["Z"]["shares"]["found_again"]["share"] >= 0.5 and R["Z"]["shares"]["carried"]["share"] < 0.2)
    mf = AC.must_fail(d["Z"][np.isin(un, ["control", "c1"])], y[np.isin(un, ["control", "c1"])], un[np.isin(un, ["control", "c1"])], first[np.isin(un, ["control", "c1"])], "c1", draws=draws, seed=seed)
    checks["a permuted stage carries nothing (|carried| < 0.15)"] = bool(abs(mf["permuted_stage"]["carried"]["share"]) < 0.15)
    checks["an event made of control cells has no change to account for"] = mf["made_up_event"]["no_change_to_account_for"]
    for st in ("V", "Z"):
        law = AC.Law(d[st], {"y": y}, un, edges=d["edges"], seed=seed).fit()
        rows = law.forecast()
        cal, fl = law.calibration(rows), law.flags(rows)
        by = {r["event"]: r for r in rows}
        flagged = {f["event"]: f["candidate"] for f in fl}
        rep[f"law at {st}"] = {"inside_95": round(cal["inside_95"], 2), "flagged": flagged, "tau2": law.r["y"]["tau2"], "sigma0": round(law.sigma0, 3),
                               "edge": {e: (round(by[e]["forecast"], 2), round(by[e]["change"], 2)) for e in EDGES}}
        checks[f"law at {st}: at least 85% of the measured changes inside the 95% interval"] = bool(cal["inside_95"] >= 0.85)
        checks[f"law at {st}: x is flagged as a missing edge"] = bool("x" in flagged and flagged["x"].startswith("a missing edge"))
        checks[f"law at {st}: at least 3 of the 4 events with an edge inside their interval"] = bool(sum(abs(by[e]["z"]) < 1.96 for e in EDGES) >= 3)
        checks[f"law at {st}: with their edge the 4 events are forecast better than without"] = bool(
            cal["by_class"]["own"]["error_forecast"] < cal["by_class"]["own"]["error_without_edges"])
        checks[f"law at {st}: no event that fills the library is flagged"] = not any(e.startswith("q") for e in flagged)
        if st == "Z":
            checks["law at Z: l is flagged (the state does not hold what moves it)"] = "l" in flagged
        else:
            checks["law at V: l is not flagged (the input view holds it)"] = "l" not in flagged
            checks["law at V: c1, c2 are not flagged"] = not ({"c1", "c2"} & set(flagged))
    rep["resolved"] = resolved
    if verbose:
        for k, v in rep.items():
            print(f"  {k:10s} {v}")
        print("  planted:", {e: {k: round(v, 2) for k, v in tr[e]["parts"].items()} for e in EVENTS if "parts" in tr[e]})
    return checks, rep


VARIANTS = {"as planted": {}, "every cell edited": {"edited": 1.0}, "half of the cells edited": {"edited": 0.5}, "noisy cells": {"noise": 0.5},
            "few cells": {"n_control": 1500, "n_event": 300}}


def suite(seeds=20, other=5, out=None, say=print):
    """The suite: `seeds` systems as planted and `other` of every variant. Returns the counts and every wrong check."""
    res = {}
    for name, kw in VARIANTS.items():
        n = seeds if name == "as planted" else other
        wrong, total, reps = [], 0, {}
        for s in range(n):
            checks, rep = run(s, **kw)
            total += len(checks)
            wrong += [f"seed {s}: {k}" for k, v in checks.items() if not v]
            reps[s] = rep
        res_ = [v for r in reps.values() for v in r["resolved"].values()]
        by_event = {e: [int(sum(r["resolved"].get(e, False) for r in reps.values())), n] for e in reps[0]["resolved"]}
        res[name] = {"systems": n, "checks": total, "wrong": len(wrong), "which": wrong, "the_planted_place_resolved_as_a_step": {"of_all": [int(sum(res_)), len(res_)], "by_event": by_event},
                     "first_system": reps.get(0), "planted": {e: v.get("parts") for e, v in truth(**kw).items()}}
        say(f"{name}: {total} checks on {n} planted systems; {len(wrong)} wrong; the planted place resolved as a step in {sum(res_)} of {len(res_)} accounts")
        for w in wrong:
            say("  wrong: " + w)
        if out:
            import json
            from pathlib import Path
            Path(out).write_text(json.dumps(res, indent=1, default=float) + "\n")
    return res


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "suite":
        suite(out=sys.argv[2] if len(sys.argv) > 2 else None)
        return
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    wrong, total = [], 0
    for s in range(seeds):
        checks, _ = run(s, verbose=(s == 0))
        total += len(checks)
        wrong += [f"seed {s}: {k}" for k, v in checks.items() if not v]
    print(f"{total} checks on {seeds} planted systems; {len(wrong)} wrong")
    for w in wrong:
        print("  wrong:", w)


if __name__ == "__main__":
    main()
