"""The candidate for version 6: what changes against version 5 (runs/freeze_v5_20261126/scripts/_v5.py), and nothing else.

6a  The arithmetic of a readout is on the graph. For every event with an edge, the direct part of the edge's readout is
    taken as measured where the pair is recorded and as the edge forecasts it where it is not; `events.derived` gives
    what that moves in every other readout of the same set of columns. A recorded row is freed of that part before it
    teaches anything; a forecast carries it.
6d  The 95% interval of a pair with no edge is as wide as the recorded pairs with no edge stand when each is left out
    and forecast as a new event (`Law.conformal`): from the stratum's own pairs of the same kind; else the library's of
    that kind; else its own of any kind; else the library's; else the mixture the law assumes. The row says which.
6e  How far a share carries between conditions is known from a few pairs measured twice. The spread of a carried share
    is widened by what so few pairs leave open: Student's t for their number in place of 1.96.

The declarations, the classes, the scales, the shares and the spread with no edge are those of version 5.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.stats import t as student

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "freeze_v5_20261126" / "scripts"))
import _v5 as V  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402

OVER = {r: ("the 24 columns of the melanoma panel", 24) for r in V.MELANOMA}       # the readouts that are a ratio, and over how many columns


def spreads(*laws):
    """One spread of every readout in its own units, among the control cells of the laws' condition."""
    return {name: q["sd"] for law in laws for name, q in law.r.items()}


def edge_parts(rows):
    """{(event, readout): (the direct part of a pair with an edge, its spread)}: as measured for the rows of recorded events
    (`Law.forecast`), as the edge forecasts it for the rows of new ones (`Law.new`)."""
    out = {}
    for r in rows:
        if r["edge"] is None or r["no_forecast"] is not None:
            continue
        out[(r["event"], r["readout"])] = (r["direct"], r["direct_sd"]) if "direct" in r else (r["by_its_edge"], r["parts_of_the_spread"]["direct"])
    return out


def values(known):
    return {k: v for k, (v, _sd) in known.items()}


def library_law(x, Y, unit, table, scales, over=OVER, types=V.PATHWAY):
    """As `_v5.library_law`, fitted twice: the first fit gives the direct parts of the pairs with an edge, the second holds
    what their arithmetic moves in the other readouts as known."""
    first = V.library_law(x, Y, unit, table, scales, types)
    known = EV.derived(edge_parts(first["law"].forecast()), spreads(first["law"]), over)
    law = AC.Law(x, Y, unit, edges=first["edges"], types=types, scale=first["scales"], known=values(known)).fit()
    vals = law.class_values()
    return {"law": law, "edges": first["edges"], "scales": first["scales"], "classes": {c: [v for _e, _r, v in rows] for c, rows in vals.items()},
            "shares": {(e, r): v for c, rows in vals.items() if not c.endswith(EV.AT_THE_BACKGROUND) for e, r, v in rows if (e, r) in first["scales"]},
            "known": known, "width": law.conformal(), "before": first["law"]}


def stratum_laws(specs, lib, over=OVER, types=V.PATHWAY, forecast_edges=None):
    """The laws of one stratum, fitted twice. specs: {key: {"x", "Y", "unit", "table", "scales", "sigma0"}}; several where
    the readouts were recorded on different units. forecast_edges: {(event, readout): (value, spread)}, the direct
    part of an edge on a readout that is not recorded for its event, as the edge forecasts it.
    Returns {key: (law, edges, scales)} and the known parts {(event, readout): (value, spread)} of the stratum."""
    def fit(sp, known=None):
        ed, sc = EV.by_scale(V.restrict(sp["table"], list(sp["Y"])), sp["scales"])
        law = AC.Law(sp["x"], sp["Y"], sp["unit"], edges=ed, types=types, scale=sc, known=known)
        return (law.fit(sigma0=sp["sigma0"]) if sp.get("sigma0") is not None else law.fit(fallback=lib["law"].sigma0)), ed, sc
    first = {k: fit(sp) for k, sp in specs.items()}
    direct = dict(forecast_edges or {})
    for law, _ed, _sc in first.values():
        direct.update(edge_parts(law.forecast()))                   # a recorded pair is taken as measured
    sd = {**(spreads(lib["law"]) if lib.get("law") is not None else {}), **spreads(*[v[0] for v in first.values()])}
    known = EV.derived(direct, sd, over)
    return {k: fit(sp, values(known)) for k, sp in specs.items()}, known, {k: v[0] for k, v in first.items()}


def widths(law, lib):
    """How wide the interval of a pair with no edge is, by kind: the stratum's own recorded pairs of the kind, else the
    library's, else the stratum's of any kind, else the library's of any kind; else nothing (the law's mixture)."""
    own, of_lib = law.conformal(), lib.get("width") or {}
    out = {}
    for kind in ("no type", "of a type"):
        for src, name, k in ((own, "the stratum's own", kind), (of_lib, "the library's", kind), (own, "the stratum's own", "any"), (of_lib, "the library's", "any")):
            if src.get(k, {}).get("quantile") is not None:
                out[kind] = {"quantile": src[k]["quantile"], "pairs": src[k]["pairs"], "from": f"{name} {src[k]['pairs']} recorded pairs with no edge ({k}), each left out"}
                break
    return out


def forecast(law, ed, sc, x, unit, lib, twice, known=None, types=V.PATHWAY):
    """As `_v5.forecast`, with the parts known by declaration and the widths of version 6."""
    own = lib["classes"].get("own", [])
    spread = EV.carried_spread(twice)
    if spread is not None:
        spread *= float(student.ppf(0.975, len(twice))) / 1.96      # the root mean square of n differences about zero: n degrees of freedom
    elif len(own) >= 2:
        spread = float(np.std(own, ddof=1))
    shares = {k: (v, spread) for k, v in lib["shares"].items() if spread is not None}
    rows = law.new(x, unit, edges=ed, classes=lib["classes"], types=types, scale=sc, shares=shares, known=known, width=widths(law, lib))
    for r in rows:
        r["spread_of_a_carried_share"] = spread if r["size_of_the_edge_from"] and r["size_of_the_edge_from"].startswith("the pair") else None
    return rows
