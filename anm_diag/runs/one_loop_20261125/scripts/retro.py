#!/usr/bin/env python3
"""Batch 2 again, with the edges as records, types, and a scale. RETROSPECTIVE: every value was read before. results/retro.json.

  python retro.py --cells batch2_cells.tsv --states STATES_DIR --read read_cells.npz --pack batch2_pack.npz
                  --library LIBRARY.npz --scale scale_pack.npz

The inputs are those of runs/batch2_forecast_20261125 (the cells, their states, the readouts read before, the readouts of
the batch) and scale_pack.npz (scripts/scale.py). Nothing is opened.

B. The laws of batch 2 built again under rules 4j and 4k: the edges come from the event graph (events.own_edges), the
   events of the interferon pathway are one type, and the direct part of an edge is its class's share of a scale that is
   measured among the control cells of each condition. The class's share comes from the library under interferon alone.
   R1, R2, R3 of the design are scored.
C. The pairs the law of the library flags are proposed as edges (events.propose_edges) and asked of the rows of the other
   two conditions, forecast without them (events.test_edges, rule G2).
D. The accounts of batch 2 (its results/reading.json) are sent through the step after a rejection (account.plan,
   outer.open_parts) with the edges of the graph.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
B2 = RUN.parent / "batch2_forecast_20261125"
sys.path.insert(0, str(B2 / "scripts"))
from _b2 import AC, EDGES, MIN_ROWS, STRATA, stratum, table, unit  # noqa: E402
from anm_diag import events as EV  # noqa: E402
from anm_diag import outer  # noqa: E402

DECL = {"readouts": [{"id": "PDL1", "features": ["CD274"]}, {"id": "MHC1", "features": ["HLA-A", "HLA-B", "HLA-C"]}, {"id": "HLAE", "features": ["HLA-E"]}]}
PARTNERS = {"MHC1": ["B2M"], "HLAE": ["B2M"]}                       # the second chain of both molecules: without it neither is made
PATHWAY = {g: "interferon pathway" for g in ("IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1")}     # the type declared in run 40; no other event has a type
CONDITION = {"library": "IFNγ", "mel_co": "Co-culture", "mel_ctl": "Control"}
SCALE_NAME = "how far the readout stands above the background of the antibody panel among control cells, in spreads of the readout"
rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else None


def scales(pack):
    """scale(readout, condition) of rule 4k, from control cells."""
    P = np.load(pack)
    cond, panel = P["condition"].astype(str), int(P["panel"])
    out = {}
    for r in ("PDL1", "MHC1", "HLAE"):
        out[r] = {}
        for c in sorted(set(cond.tolist())):
            m = cond == c
            above, spread = float((P[f"log__{r}"][m] - P["background"][m]).mean()), float((P[f"log__{r}"][m] - P["centre"][m]).std())
            out[r][c] = {"control_cells": int(m.sum()), "above_background": above, "spread_of_the_readout": spread, "scale": above * (1 - 1 / panel) / spread}
    return out


def scored(rows, law, y, units):
    """The rows of Law.new with their measured change."""
    for r in rows:
        ch, se = law.measure(y[r["readout"]], units, r["readout"])[r["event"]]
        r.update(change=ch, se=se, snr=abs(ch) / se, z=(ch - r["forecast"]) / r["sd"], inside=bool(r["interval"][0] <= ch <= r["interval"][1]))
    return rows


def calib(rows):
    z = np.array([r["z"] for r in rows])
    return {"pairs": len(rows), "inside": int(sum(r["inside"] for r in rows)), "mean_z2": float(np.mean(z ** 2)),
            "error": {"no_change": rms([r["change"] for r in rows]), "control_readout": rms([r["change"] - r["read"] for r in rows]), "forecast": rms([r["change"] - r["forecast"] for r in rows])}}


def main():
    ap = argparse.ArgumentParser()
    for k in ("cells", "states", "read", "pack", "library", "scale"):
        ap.add_argument("--" + k, type=Path, required=True)
    a = ap.parse_args()
    cells = table(a.cells)
    sc = scales(a.scale)
    S = lambda g, edges: {(e, r): sc[r][CONDITION[g]]["scale"] for (e, r) in edges}
    res = {"note": "retrospective: every value here was read before (the library in development, batch 1 on 23 November, batch 2 this morning); nothing is a forecast",
           "scale": sc, "types": PATHWAY}

    # --- the library under interferon, and the edges of the graph
    d = np.load(a.library)
    un = np.where(d["A"] == 0, "control", d["target"].astype(str))
    sites = sorted(set(un.tolist()) - {"control"})
    strata = {g: stratum(g, cells, a.states, read=a.read, new=a.pack) for g in STRATA}
    for c in strata.values():
        sites = sorted(set(sites) | {u.replace("new:", "") for u in c["unit"].tolist() if u != "control"})
    declared = EV.own_edges(DECL, sites, partners=PARTNERS, scale=SCALE_NAME)
    edges = EV.edge_table(declared)
    assert edges == {k: v for k, v in EDGES.items() if k[0] in sites}, "the declared edges are not the table of batch 2"
    ZL, YL = unit(d["Z"]).astype(np.float64), {k: d["readout__" + k].astype(np.float64) for k in ("PDL1", "MHC1", "HLAE")}
    lib = AC.Law(ZL, YL, un, edges=edges, types=PATHWAY, scale=S("library", edges)).fit()
    lib_one = AC.Law(ZL, YL, un, edges=edges).fit()                 # as batch 2 had it: one type, one size per class
    shares = lib.class_values()["own"]
    sizes = lib_one.class_values()["own"]
    lib_rows = lib.forecast()
    flags = lib.flags(lib_rows)
    res["library"] = {"events": len(lib.events), "spread_of_a_direct_part_without_an_edge": lib.sigma0, "widening_of_the_control_prior": lib.summary()["widening_of_the_control_prior"],
                      "calibration": {k: v for k, v in lib.calibration(lib_rows).items() if k != "by_class"},
                      "class_own": [{"event": e, "readout": r, "share_of_the_scale": v, "scale": sc[r]["IFNγ"]["scale"], "direct_part_with_one_size": w} for (e, r, v), (_e, _r, w) in zip(shares, sizes)],
                      "flags": flags}
    share_vals, size_vals = [v for _e, _r, v in shares], [v for _e, _r, v in sizes]
    for e in declared:                                              # the edge carries the share its pair showed under interferon
        v = [x for h, r, x in shares if (h, r) == (e["site"], e["readout"])]
        if v:
            e["share"] = {"value": v[0], "measured_in": "the library under interferon"}
    print(f"library: {len(lib.events)} events; spread without an edge {lib.sigma0:.3f}; widening {res['library']['widening_of_the_control_prior']}")
    print("  class own, shares of the scale:", ", ".join(f"{e}->{r} {v:+.2f}" for e, r, v in shares), f"| mean {np.mean(share_vals):+.2f}, spread {np.std(share_vals, ddof=1):.2f}")
    print("  class own, one size:           ", ", ".join(f"{e}->{r} {v:+.2f}" for e, r, v in sizes), f"| mean {np.mean(size_vals):+.2f}, spread {np.std(size_vals, ddof=1):.2f}")
    print("  flags:", [(f["event"], f["readout"], round(f["direct"], 2), f["flag"]) for f in flags])

    # --- B. the strata, with types and a scale
    res["strata"], laws, new_rows, own10 = {}, {}, {}, []
    for g, (study, cond, readouts) in STRATA.items():
        c = strata[g]
        rec, new = c["role"] != "event", c["role"] == "event"
        ed = edges if study == "frangieh2021" else {}
        s = S(g, ed) if ed else None
        Y = {r: c[f"y_{r}"] for r in readouts}
        law = AC.Law(c["Z"][rec], {r: Y[r][rec] for r in readouts}, c["unit"][rec], edges=ed, types=PATHWAY, scale=s)
        free = sum(1 for e in law.events for r in readouts if (e, r) not in ed)
        law.fit(sigma0=None if free >= MIN_ROWS else lib.sigma0)
        genes = np.char.replace(c["unit"][new], "new:", "")
        Yn = {r: Y[r][new] for r in readouts}
        rows = scored(law.new(c["Z"][new], genes, edges=ed, classes={"own": share_vals}, types=PATHWAY, scale=s), law, Yn, genes)
        one = law.new(c["Z"][new], genes, edges=ed, classes={"own": size_vals}, types=PATHWAY)        # the same law, one size per class
        for r, o in zip(rows, one):
            r["forecast_with_one_size"], r["interval_with_one_size"] = o["forecast"], o["interval"]
        loo = law.forecast()
        laws[g], new_rows[g] = (law, c, loo), rows
        res["strata"][g] = {"study": study, "condition": cond, "recorded_events": law.events, "types": {e: t for e, t in zip(law.events, law.type)},
                            "spread_of_a_direct_part_without_an_edge": {"value": law.sigma0, "from": "the library under interferon" if law.fixed else "the stratum's own rows"},
                            "widening_of_the_control_prior": law.summary()["widening_of_the_control_prior"],
                            "recorded_pairs_each_left_out": [{k: r[k] for k in ("event", "readout", "type", "edge", "change", "se", "read", "from_the_other_events", "by_its_edge", "forecast", "sd", "interval", "z")} for r in loo],
                            "pairs_of_the_batch": rows, "calibration": calib(rows)}
        print(f"{g}: spread without an edge {law.sigma0:.3f} ({'library' if law.fixed else 'own rows'}); widening {res['strata'][g]['widening_of_the_control_prior']}; batch: {res['strata'][g]['calibration']}")
        own10 += [dict(r, stratum=g) for r in rows if r["edge"] == "own"]
    # B2M with no signal given (batch 1): forecast by a law that has not seen it
    law, c, _loo = laws["mel_ctl"]
    rec = (c["role"] != "event") & (c["unit"] != "B2M")
    ed, s = edges, S("mel_ctl", edges)
    part = AC.Law(c["Z"][rec], {r: c[f"y_{r}"][rec] for r in ("PDL1", "MHC1", "HLAE")}, c["unit"][rec], edges=ed, types=PATHWAY, scale=s).fit(sigma0=lib.sigma0)
    m = c["unit"] == "B2M"
    Yb = {r: c[f"y_{r}"][m] for r in ("PDL1", "MHC1", "HLAE")}
    b2m = scored(part.new(c["Z"][m], c["unit"][m], edges=ed, classes={"own": share_vals}, types=PATHWAY, scale=s), part, Yb, c["unit"][m])
    for r, o in zip(b2m, part.new(c["Z"][m], c["unit"][m], edges=ed, classes={"own": size_vals}, types=PATHWAY)):
        r["forecast_with_one_size"], r["interval_with_one_size"] = o["forecast"], o["interval"]
    own10 += [dict(r, stratum="mel_ctl", batch=1) for r in b2m if r["edge"] == "own"]

    allrows = [dict(r, stratum=g) for g, rows in new_rows.items() for r in rows]
    free = [r for r in allrows if r["edge"] is None]
    committed = json.loads((B2 / "results" / "outcome.json").read_text())
    res["calibration_of_the_batch"] = {"as_committed_before_the_batch_was_read": {"pairs": committed["pairs"], "inside": committed["inside"], "error": committed["error"], "mean_z2": 0.30},
                                       "with_types_a_scale_and_the_error_of_a_row_beside_its_readout": calib(allrows), "pairs_without_an_edge": calib(free)}

    # R1: with types the pathway events still forecast each other, and no event without a type is moved by them
    r1 = {}
    for g in ("mel_co", "thp1"):
        rows = [r for r in laws[g][2] if r["type"] == "interferon pathway"]
        r1[g] = {"pathway_pairs": len(rows), "inside": int(sum(abs(r["z"]) < 1.96 for r in rows)), "error": {"control_readout": rms([r["change"] - r["read"] for r in rows]), "forecast": rms([r["change"] - r["forecast"] for r in rows])}}
    moved = max(abs(r["from_the_recorded_events"]) for r in free)
    r1["pairs_without_an_edge"] = {"pairs": len(free), "largest_part_from_the_recorded_events": moved, "error_forecast": rms([r["change"] - r["forecast"] for r in free]),
                                   "error_control_readout": rms([r["change"] - r["read"] for r in free])}
    r1["outcome"] = "holds" if all(v["error"]["forecast"] < v["error"]["control_readout"] for v in (r1["mel_co"], r1["thp1"])) and moved < 1e-9 else "does not hold"
    res["R1"] = r1

    # R2: the share of the class times the scale, against one size per class, on the ten own-gene pairs outside the library
    assert len(own10) == 10
    lib_share = {k: v[0] for k, v in EV.edge_shares(declared, 0.0).items()}
    for r in own10:
        g = r["stratum"]
        r["scale"] = sc[r["readout"]][CONDITION[g]]["scale"]
        r["inside_with_one_size"] = bool(r["interval_with_one_size"][0] <= r["change"] <= r["interval_with_one_size"][1])
        r["share_measured"] = (r["change"] - r["read"]) / r["scale"]
        r["share_of_the_same_pair_under_interferon"] = lib_share.get((r["event"], r["readout"]))
        r["forecast_from_the_same_pair_under_interferon"] = r["read"] + lib_share[(r["event"], r["readout"])] * r["scale"]       # explored: the edge carries its own share
    e_scale, e_one = rms([r["change"] - r["forecast"] for r in own10]), rms([r["change"] - r["forecast_with_one_size"] for r in own10])
    res["R2"] = {"pairs": [{k: r.get(k) for k in ("stratum", "batch", "event", "readout", "cells", "change", "se", "read", "scale", "forecast", "interval", "inside", "forecast_with_one_size", "interval_with_one_size",
                                                   "inside_with_one_size", "share_measured", "share_of_the_same_pair_under_interferon", "forecast_from_the_same_pair_under_interferon")} for r in own10],
                 "error": {"share_of_the_class_times_the_scale": e_scale, "one_size_per_class": e_one, "no_edge": rms([r["change"] - r["read"] for r in own10]),
                           "explored_the_pair's_own_share_under_interferon_times_the_scale": rms([r["change"] - r["forecast_from_the_same_pair_under_interferon"] for r in own10])},
                 "inside": {"share_of_the_class_times_the_scale": int(sum(r["inside"] for r in own10)), "one_size_per_class": int(sum(r["inside_with_one_size"] for r in own10))},
                 "error_by_readout": {ro: {"pairs": sum(r["readout"] == ro for r in own10), "no_edge": rms([r["change"] - r["read"] for r in own10 if r["readout"] == ro]),
                                           "one_size_per_class": rms([r["change"] - r["forecast_with_one_size"] for r in own10 if r["readout"] == ro]),
                                           "share_of_the_class_times_the_scale": rms([r["change"] - r["forecast"] for r in own10 if r["readout"] == ro]),
                                           "explored_the_pair's_own_share_under_interferon_times_the_scale": rms([r["change"] - r["forecast_from_the_same_pair_under_interferon"] for r in own10 if r["readout"] == ro])}
                                      for ro in ("PDL1", "MHC1", "HLAE")},
                 "half_width_of_the_interval": {"share_of_the_class_times_the_scale": float(np.mean([np.diff(r["interval"])[0] / 2 for r in own10])), "one_size_per_class": float(np.mean([np.diff(r["interval_with_one_size"])[0] / 2 for r in own10]))}}
    res["R2"]["outcome"] = "holds" if e_scale < e_one and res["R2"]["inside"]["share_of_the_class_times_the_scale"] >= 8 else "does not hold"

    # R3: for the same pair, the change divided by the scale varies less across conditions than the change
    by = {}
    for e, r, v in shares:
        row = next(x for x in lib_rows if x["event"] == e and x["readout"] == r)
        by.setdefault((e, r), {})["IFNγ"] = {"change": row["change"], "direct": row["change"] - row["read"], "scale": sc[r]["IFNγ"]["scale"]}
    for r in own10:
        by.setdefault((r["event"], r["readout"]), {})[CONDITION[r["stratum"]]] = {"change": r["change"], "direct": r["change"] - r["read"], "scale": r["scale"]}
    cv = lambda v: float(np.std(v) / abs(np.mean(v)))
    r3 = []
    for (e, r), cs in by.items():
        if len(cs) < 2:
            continue
        ch, sh = [v["change"] for v in cs.values()], [v["change"] / v["scale"] for v in cs.values()]
        r3.append({"event": e, "readout": r, "conditions": cs, "variation_of_the_change": cv(ch), "variation_of_the_change_over_the_scale": cv(sh), "less_with_the_scale": bool(cv(sh) < cv(ch))})
    res["R3"] = {"how": "variation: the spread over the conditions divided by the size of the mean", "pairs": r3, "less_with_the_scale": int(sum(x["less_with_the_scale"] for x in r3)), "of": len(r3)}
    res["R3"]["outcome"] = "holds" if res["R3"]["less_with_the_scale"] == len(r3) else ("holds for PD-L1 and MHC-I, fails for HLA-E" if all(x["less_with_the_scale"] for x in r3 if x["readout"] != "HLAE") else "does not hold")

    # --- C. the flags of the library, proposed as edges, asked of the rows of the other two conditions
    tested = []
    for g in ("mel_co", "mel_ctl"):
        law, c, loo = laws[g]
        rows = [r for r in loo] + [r for r in new_rows[g]] + (b2m if g == "mel_ctl" else [])
        prop = EV.propose_edges(flags, known=declared, stratum=g)
        have = {(r["event"], r["readout"]) for r in rows}
        tested += [dict(e, proposed_in="the library under interferon") for e in EV.test_edges([p for p in prop if (p["site"], p["readout"]) in have], rows)]
    seen = {(e["site"], e["readout"]) for e in tested}
    untested = [dict(e, proposed_in="the library under interferon") for e in EV.propose_edges(flags, known=declared) if (e["site"], e["readout"]) not in seen]     # no row elsewhere: stays proposed
    graph = declared + tested + untested
    res["edges"] = graph
    res["edge_checks"] = EV.check_edges(graph, sites=sites, readouts=["PDL1", "MHC1", "HLAE"])
    for e in tested:                                                # explored: is the direct part of a flagged pair the same share of the scale in both conditions?
        e["explored_share_of_the_scale"] = {"under_interferon": e["evidence"]["direct"] / sc[e["readout"]]["IFNγ"]["scale"],
                                            "where_tested": (e["test"]["change"] - e["test"]["forecast_without_the_edge"]) / sc[e["readout"]][CONDITION[e["stratum"]]]["scale"]}
    res["C"] = {"proposed": sorted({(e["site"], e["readout"]) for e in EV.propose_edges(flags, known=declared)}),
                "tested": [{"edge": e["edge_id"], "status": e["status"], "proposed_with": e["evidence"], "test": e["test"], "explored_share_of_the_scale": e["explored_share_of_the_scale"]} for e in tested],
                "not_tested": [e["edge_id"] for e in untested],
                "counts": {**{s: sum(e["status"] == s for e in tested) for s in ("admitted", "refused", "undecided")}, "proposed, no row elsewhere": len(untested)}}
    print("C:", [(e["edge_id"], e["status"], round(e["test"]["z"], 1), {k: round(v, 3) for k, v in e["explored_share_of_the_scale"].items()}) for e in tested], "| checks:", res["edge_checks"])

    # --- D. the accounts of batch 2 through the step after a rejection
    reading = json.loads((B2 / "results" / "reading.json").read_text())
    accounts = []
    for g, v in reading["strata"].items():
        for key, acc in v["accounts"].items():
            ev, r = key.split("|")
            acc = dict(acc, plan=AC.plan(acc))
            accounts.append({"readout": r, "group": g, "type": PATHWAY.get(ev), "event": ev, "account": acc})
    items = outer.open_parts(accounts, edges=[e for e in graph if e["status"] in ("declared", "admitted")])
    res["step_after_a_rejection"] = {"accounts": len(accounts), "items": items,
                                     "by_reading_and_state": {f"{k[0]} | {k[1]}": n for k, n in sorted({(i["reading"], i["state"]): sum(1 for j in items if (j["reading"], j["state"]) == (i["reading"], i["state"])) for i in items}.items())}}
    print("D:", res["step_after_a_rejection"]["by_reading_and_state"])

    print("calibration:", json.dumps(res["calibration_of_the_batch"], default=float)[:900])
    print("R1:", json.dumps(res["R1"], default=float))
    print("R2:", json.dumps({k: res["R2"][k] for k in ("error", "inside", "half_width_of_the_interval", "outcome")}, default=float))
    for r in res["R2"]["pairs"]:
        print(f"   {r['stratum']:7s} {r['event']:6s}->{r['readout']:5s} change {r['change']:+.2f} | scale {r['scale']:+.2f} | share x scale {r['forecast']:+.2f} [{r['interval'][0]:+.2f}, {r['interval'][1]:+.2f}] {'in ' if r['inside'] else 'OUT'}"
              f" | one size {r['forecast_with_one_size']:+.2f} {'in ' if r['inside_with_one_size'] else 'OUT'} | measured share {r['share_measured']:+.2f}, under interferon {r['share_of_the_same_pair_under_interferon']:+.2f} -> {r['forecast_from_the_same_pair_under_interferon']:+.2f}")
    print("R2 by readout:", json.dumps(res["R2"]["error_by_readout"], default=float))
    print("R3:", res["R3"]["outcome"], [(x["event"], x["readout"], round(x["variation_of_the_change"], 2), round(x["variation_of_the_change_over_the_scale"], 2)) for x in r3])
    (RUN / "results" / "retro.json").write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
