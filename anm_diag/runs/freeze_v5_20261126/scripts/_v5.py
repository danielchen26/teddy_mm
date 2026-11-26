"""What version 5 of the protocol holds fixed for the two screens: the declarations, and the system that forecasts.

A readout of the melanoma screen is the centred log-ratio of one antibody over the 24 columns of the panel; of the THP-1
screen, the log of one plus the counts of one antibody. The readouts read before keep their names (PDL1, MHC1, HLAE).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402

ISOTYPES = ["Rat_IgG2a", "Mouse_IgG1", "Mouse_IgG2a", "Mouse_IgG2b"]
# the antibody of a readout, and the genes of the protein it measures (rule S3 and rule 4j read these)
MELANOMA = {"PDL1": ("CD274", ["CD274"]), "MHC1": ("HLA_A", ["HLA-A", "HLA-B", "HLA-C"]), "HLAE": ("HLA_E", ["HLA-E"]),
            "CD117": ("CD117", ["KIT"]), "CD119": ("CD119", ["IFNGR1"]), "CD140a": ("CD140a", ["PDGFRA"]), "CD140b": ("CD140b", ["PDGFRB"]), "CD172a": ("CD172a", ["SIRPA"]),
            "CD184": ("CD184", ["CXCR4"]), "CD202b": ("CD202b", ["TEK"]), "CD29": ("CD29", ["ITGB1"]), "CD309": ("CD309", ["KDR"]), "CD44": ("CD44", ["CD44"]), "CD47": ("CD47", ["CD47"]),
            "CD49f": ("CD49f", ["ITGA6"]), "CD58": ("CD58", ["CD58"]), "CD59": ("CD59", ["CD59"]), "CD61": ("CD61", ["ITGB3"]), "CD9": ("CD9", ["CD9"]), "CD279": ("CD279", ["PDCD1"])}
THP1 = {"PDL1": ("PDL1", ["CD274"]), "PDL2": ("PDL2", ["PDCD1LG2"]), "CD86": ("CD86", ["CD86"]), "CD366": ("CD366", ["HAVCR2"])}
PARTNERS = {"MHC1": ["B2M"], "HLAE": ["B2M"]}                       # the second chain of both molecules: without it neither is made
PATHWAY = {g: "interferon pathway" for g in ("IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1")}     # the one declared type (run 40)
USED = {"melanoma": ["PDL1", "MHC1", "HLAE"], "thp1": ["PDL1"]}
SCALE_NAME = "how far the readout stands above the background of the antibody panel among control cells, in spreads of the readout"
unit_length = lambda z: z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-12)


def declared(sites, panel=MELANOMA, partners=PARTNERS):
    """The edges that follow from the declarations, as records and as a table {(site, readout): class}."""
    decl = {"readouts": [{"id": r, "features": genes} for r, (_ab, genes) in panel.items()]}
    edges = EV.own_edges(decl, sites, partners=partners, scale=SCALE_NAME)
    return edges, EV.edge_table(edges)


def scale(log, background, clr, panel=24):
    """scale(readout, condition) of rule 4k from control cells: how far the readout stands above the background, in its
    spreads. log, background, clr: the log count of the readout's antibody, the mean log count of the isotype
    antibodies, and the readout, over the control cells of one condition."""
    return float((log - background).mean() * (1 - 1 / panel) / clr.std())


def restrict(table, readouts, sites=None):
    return {k: v for k, v in table.items() if k[1] in readouts and (sites is None or k[0] in sites)}


def library_law(x, Y, unit, table, scales, types=PATHWAY):
    """The law of a library of recorded events, with what it gives a stratum that has too few rows of its own: the
    spread with no edge, the shares of the recorded pairs of every class, and the share of every pair with a scale."""
    ed, sc = EV.by_scale(restrict(table, list(Y)), scales)
    law = AC.Law(x, Y, unit, edges=ed, types=types, scale=sc).fit()
    values = law.class_values()
    return {"law": law, "edges": ed, "scales": sc, "classes": {c: [v for _e, _r, v in rows] for c, rows in values.items()},
            "shares": {(e, r): v for c, rows in values.items() if not c.endswith(EV.AT_THE_BACKGROUND) for e, r, v in rows if (e, r) in sc}}


def stratum_law(x, Y, unit, table, scales, lib, types=PATHWAY):
    """The law of the recorded rows of one stratum. The spread with no edge is its own where it has MIN_ROWS recorded
    pairs with no edge and no type, else the library's."""
    ed, sc = EV.by_scale(restrict(table, list(Y)), scales)
    law = AC.Law(x, Y, unit, edges=ed, types=types, scale=sc).fit(fallback=lib["law"].sigma0)
    return law, ed, sc


def measured_twice(lib, strata_laws):
    """The recorded pairs with a share in the library and in a stratum: [(share in the library, share in the stratum)]."""
    out = []
    for law, _ed, sc in strata_laws:
        for c, rows in law.class_values().items():
            out += [(lib["shares"][(e, r)], v) for e, r, v in rows if (e, r) in sc and (e, r) in lib["shares"]]
    return out


def forecast(law, ed, sc, x, unit, lib, twice, types=PATHWAY):
    """Units not yet measured, by the frozen system: the control readout of the state, the rows of the unit's type, and
    for a pair with an edge its size. The size: the pair's own share of the scale where the library holds one (its
    spread: how far the shares of pairs measured twice moved, or the spread of the class while fewer than three were);
    else the class's share; for an edge at the background the recorded pairs of that class, in spreads."""
    own = lib["classes"].get("own", [])
    spread = EV.carried_spread(twice) or (float(np.std(own, ddof=1)) if len(own) >= 2 else None)
    shares = {k: (v, spread) for k, v in lib["shares"].items() if spread is not None}
    rows = law.new(x, unit, edges=ed, classes=lib["classes"], types=types, scale=sc, shares=shares)
    for r in rows:
        r["spread_of_a_carried_share"] = spread if r["size_of_the_edge_from"] and r["size_of_the_edge_from"].startswith("the pair") else None
    return rows


def scored(rows, law, Y, units):
    """The rows of a forecast with the measured change."""
    for r in rows:
        ch, se = law.measure(Y[r["readout"]], units, r["readout"])[r["event"]]
        ok = r["no_forecast"] is None
        r.update(change=ch, se=se, snr=abs(ch) / se, z=(ch - r["forecast"]) / r["sd"] if ok else None, inside=bool(r["interval"][0] <= ch <= r["interval"][1]) if ok else None)
    return rows


rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None


def calibration(rows):
    ok = [r for r in rows if r["no_forecast"] is None]
    z = np.array([r["z"] for r in ok])
    return {"pairs": len(rows), "with_a_forecast": len(ok), "inside": int(sum(r["inside"] for r in ok)), "mean_z2": float(np.mean(z ** 2)) if len(ok) else None,
            "median_z2": float(np.median(z ** 2)) if len(ok) else None,
            "error": {"no_change": rms([r["change"] for r in ok]), "control_readout": rms([r["change"] - r["read"] for r in ok]), "forecast": rms([r["change"] - r["forecast"] for r in ok])}}
