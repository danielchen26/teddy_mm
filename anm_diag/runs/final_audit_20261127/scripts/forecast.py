#!/usr/bin/env python3
"""The forecasts of the final audit by three systems, before any value of a final unit or of a final control cell is
opened. results/forecasts.json, written once.

  python forecast.py --cells final_cells.tsv --states MEL_STATES... --states-thp THP_STATES... --library LIBRARY.npz
                     --packs batch3_pack_A.npz batch3_pack_B.npz batch4_pack.npz

S0  the model alone: the state of the unit's cells, read by a readout of the control cells of its stratum
S2  version 5 (`_v5.py`, frozen at d777070)
S3  version 6 (`_v6.py`, frozen at b4f6196)

Every law stands on what was read before: the control cells and every unit of batches 1 to 4. The melanoma strata are
forecast as batch 4 was (one law per stratum over its 20 readouts; the graph: the declared edges and the edges admitted
so far in that stratum). THP-1 as the frozen system has always forecast it: one law over its four readouts, the one
declared type, no edge (its panel has no isotype antibodies, so no scale is defined there).

A final unit is read against the final control cells of its stratum (`_fa.relative`, `_fa.against`): the forecast is
of the change against that reference, and its spread carries the error of that reference. THP-1 has no final control
cell and is read against its control cells.

The script asserts that no cell of a final unit and no final control cell has a value in the files it reads.
`core` is also what `preflight.py` runs on cells read before.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fa as F  # noqa: E402
from _fa import V, V6  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
key = lambda r: f"{r['stratum']}|{r['event']}|{r['readout']}"


def core(cells, st, packs, library, new_units, on_record, check=True):
    """Everything forecast for the cells of kind `new`. new_units: {stratum: [unit]}. Returns the content of forecasts.json
    without the hashes of its inputs."""
    S = {g: F.stratum(g, cells, st, new_units.get(g)) for g in F.STRATA}
    for g, c in S.items():
        study = F.STRATA[g][0]
        old = np.isin(c["kind"], ["control", "recorded"])
        assert packs.has(study, c["ids"][old]).all(), "a cell read before has no antibody columns in the files given"
        assert not check or not packs.has(study, c["ids"][~old]).any(), "a cell of the final audit has values in the files given"
    d = np.load(library)
    lib_ids, un = d["ids"].astype(str), np.where(d["A"] == 0, "control", d["target"].astype(str))
    YL, logL, bgL = packs.readouts("frangieh2021", lib_ids)
    CL = un == "control"
    scale = {r: {"library": V.scale(logL[r][CL], bgL[CL], YL[r][CL])} for r in V.MELANOMA}
    for g in F.MEL:
        c = S[g]
        C = c["kind"] == "control"
        _Y, log, bg = packs.readouts("frangieh2021", c["ids"][C])
        for r in V.MELANOMA:
            scale[r][g] = V.scale(log[r], bg, _Y[r])
    sc_of = lambda g: {r: scale[r][g] for r in V.MELANOMA}
    sites = sorted((set(un.tolist()) | {u for g in F.MEL for u in S[g]["unit"].tolist()}) - {"control", "reference"})
    declared, tab = V.declared(sites)
    admitted = [e for e in on_record if e["status"] == "admitted"]
    tab_of = {g: {**tab, **{(e["site"], e["readout"]): e["class"] for e in admitted if e["stratum"] == g}} for g in F.MEL}

    # --- the library under interferon, by both versions; what it flags is asked of the strata that have not been asked
    ZL = V.unit_length(d["Z"].astype(np.float64))
    lib5, lib6 = V.library_law(ZL, YL, un, tab, sc_of("library")), V6.library_law(ZL, YL, un, tab, sc_of("library"))
    flagged = EV.propose_edges(lib6["law"].flags(), known=declared)
    asked = {(e["site"], e["readout"], e["stratum"]) for e in on_record if e["status"] in ("admitted", "refused")}
    to_test = {g: [e for e in flagged if (e["site"], e["readout"], g) not in asked and e["site"] in set(new_units.get(g, []))] for g in F.MEL}

    res = {"units": sum(len(v) for v in new_units.values()), "code": {"S2": "d777070, runs/freeze_v5_20261126/scripts/_v5.py", "S3": "b4f6196, runs/v6_candidate_20261126/scripts/_v6.py"},
           "scale": scale, "edges": {"declared": declared, "on_record": on_record, "flagged_by_the_library": [e["edge_id"] for e in flagged], "to_test_on_the_rows_of_the_audit": to_test},
           "library": {"events": len(lib6["law"].events), "spread_with_no_edge": {"S2": lib5["law"].sigma0, "S3": lib6["law"].sigma0}, "width": lib6["width"],
                       "shares": [{"event": e, "readout": r, "share": v} for (e, r), v in lib6["shares"].items()]},
           "strata": {}, "rows": []}

    def add(g, c, r5, r6, l5, l6, extra):
        by5 = {(r["event"], r["readout"]): r for r in r5}
        assert len(by5) == len(r5) == len(r6)
        for r in r6:
            o = by5[(r["event"], r["readout"])]
            assert abs(o["read"] - r["read"]) < 1e-9 and r["no_forecast"] is None and o["no_forecast"] is None, (g, r["event"], r["readout"], r["no_forecast"])
            untyped = r["type"] is None or r["type"].startswith("(")
            res["rows"].append({"stratum": g, "event": r["event"], "readout": r["readout"], "cells": r["cells"], "reference_cells": int(c["ref"].sum()), "edge": r["edge"], "type": None if untyped else r["type"],
                                "scale": r["scale"] if r["edge"] is not None and not r["edge"].endswith(EV.AT_THE_BACKGROUND) else None, "size_of_the_edge_from": r["size_of_the_edge_from"],
                                "shift_times_noise": r["shift_times_noise"], "S0": {"forecast": r["read"]},
                                "S2": {"forecast": o["forecast"], "sd": o["sd"], "interval": o["interval"], "from_the_recorded_events": o["from_the_recorded_events"], "by_its_edge": o["by_its_edge"]},
                                "S3": {"forecast": r["forecast"], "sd": r["sd"], "interval": r["interval"], "width_from": r["width_from"], "from_the_recorded_events": r["from_the_recorded_events"],
                                       "by_its_edge": r["by_its_edge"], "known_by_declaration": r["known_by_declaration"], "parts_of_the_spread": r["parts_of_the_spread"]}})
        res["strata"][g] = {"recorded_events": l6.events, "control_cells": int((c["kind"] == "control").sum()), "reference_cells": int(c["ref"].sum()), "own_reference": c["own_reference"],
                            "new_units": sorted(set(c["unit"][c["kind"] == "new"].tolist())), "spread_with_no_edge": {"S2": l5.sigma0, "S3": l6.sigma0},
                            "spread_from": {"S2": l5.spread_from, "S3": l6.spread_from}, "width": l6.conformal(), **extra}

    # --- the melanoma strata: every unit read so far, every readout
    law5, law6 = {}, {}
    for g in F.MEL:
        c = S[g]
        rec = np.isin(c["kind"], ["control", "recorded"])
        Y, _l, _b = packs.readouts("frangieh2021", c["ids"][rec])
        law5[g] = V.stratum_law(c["Z"][rec], Y, c["unit"][rec], tab_of[g], sc_of(g), lib5)
        laws, _known, _first = V6.stratum_laws({"all": {"x": c["Z"][rec], "Y": Y, "unit": c["unit"][rec], "table": tab_of[g], "scales": sc_of(g)}}, lib6)
        law6[g] = laws["all"]
    twice5, twice6 = V.measured_twice(lib5, list(law5.values())), V.measured_twice(lib6, list(law6.values()))
    res["pairs_measured_twice"] = {"S2": {"pairs": len(twice5), "spread_of_a_carried_share": EV.carried_spread(twice5)}, "S3": {"pairs": len(twice6), "spread_of_a_carried_share": EV.carried_spread(twice6)}}
    for g in F.MEL:
        c = S[g]
        new, n_ref = c["kind"] == "new", int(c["ref"].sum())
        if not new.any():
            continue
        x5, x6, u = F.relative(c, law5[g][0], new), F.relative(c, law6[g][0], new), c["unit"][new]
        r5 = F.against(V.forecast(*law5[g], x5, u, lib5, twice5), law5[g][0], n_ref)
        first = V6.forecast(*law6[g], x6, u, lib6, twice6)
        known = EV.derived(V6.edge_parts(first), V6.spreads(law6[g][0]), V6.OVER)
        r6 = F.against(V6.forecast(*law6[g], x6, u, lib6, twice6, known=known), law6[g][0], n_ref)
        add(g, c, r5, r6, law5[g][0], law6[g][0], {"new_pairs_with_a_known_part": sum(1 for v in known.values() if abs(v[0]) > 0)})

    # --- THP-1: one law over its four readouts, the declared type, no edge
    c = S["thp1"]
    new = c["kind"] == "new"
    if new.any():
        rec = np.isin(c["kind"], ["control", "recorded"])
        Y, _l, _b = packs.readouts("papalexi2021", c["ids"][rec])
        thp = AC.Law(c["Z"][rec], Y, c["unit"][rec], types=V.PATHWAY).fit()
        thp_lib = {"law": thp, "classes": {}, "shares": {}, "width": thp.conformal()}
        l5, ed, scl = V.stratum_law(c["Z"][rec], Y, c["unit"][rec], {}, {}, thp_lib)
        n_ref = int(c["ref"].sum())
        r5 = F.against(V.forecast(l5, ed, scl, F.relative(c, l5, new), c["unit"][new], thp_lib, []), l5, n_ref)
        r6 = F.against(V6.forecast(l5, ed, scl, F.relative(c, l5, new), c["unit"][new], thp_lib, []), l5, n_ref)
        add("thp1", c, r5, r6, l5, l5, {})

    rows = res["rows"]
    res["pairs"] = len(rows)
    with_edge = {(r["stratum"], r["event"]) for r in rows if r["edge"] is not None}
    removes = {(r["stratum"], r["event"]) for r in rows if r["edge"] is not None and r["scale"] is not None and r["stratum"] in F.MEL}
    res["sets"] = {"the_systems_can_differ": [key(r) for r in rows if (r["stratum"], r["event"]) in with_edge or r["type"] is not None],
                   "F3": [key(r) for r in rows if r["edge"] is not None and r["scale"] is not None],
                   "F4": [key(r) for r in rows if (r["stratum"], r["event"]) in removes and r["edge"] is None],
                   "with_an_edge": [key(r) for r in rows if r["edge"] is not None]}
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--states-thp", type=Path, nargs="+", required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    a = ap.parse_args()
    out_file = F.RUN / "results" / "forecasts.json"
    if out_file.exists():
        raise SystemExit("forecasts.json exists: the forecasts were made once and are not made again")
    st = {"frangieh2021": F.States(a.states), "papalexi2021": F.States(a.states_thp)}
    res = core(F.table(a.cells), st, F.Packs(a.packs), a.library, {g: F.final_units(g) for g in F.STRATA}, F.edges_on_record())
    res = {"cells_file_sha256": sha(a.cells), "states_sha256": [sha(p) for p in a.states + a.states_thp], "packs_sha256": [sha(p) for p in a.packs], **res}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    rows = res["rows"]
    print("pairs forecast:", res["pairs"], {g: sum(r["stratum"] == g for r in rows) for g in F.STRATA}, "| sets", {k: len(v) for k, v in res["sets"].items()})
    for g, v in res["strata"].items():
        hw = lambda s: np.median([(r[s]["interval"][1] - r[s]["interval"][0]) / 2 for r in rows if r["stratum"] == g and r["edge"] is None])
        print(f"  {g:7s} recorded {len(v['recorded_events'])} | reference {v['reference_cells']} ({'final control cells' if v['own_reference'] else 'its control cells'}) | spread S2 {v['spread_with_no_edge']['S2']:.3f} S3 {v['spread_with_no_edge']['S3']:.3f} "
              f"| width {v['width']['no type']} | median half-width S2 {hw('S2'):.3f} S3 {hw('S3'):.3f}")
    print("edges to test:", {g: [e["edge_id"] for e in v] for g, v in res["edges"]["to_test_on_the_rows_of_the_audit"].items()})
    for r in rows:
        if r["edge"] is not None:
            print(f"  {r['stratum']:7s} {r['event']:6s}->{r['readout']:6s} {r['edge']:26s} scale {r['scale']} | S0 {r['S0']['forecast']:+.2f} | S2 {r['S2']['forecast']:+.2f} | S3 {r['S3']['forecast']:+.2f} [{r['S3']['interval'][0]:+.2f}, {r['S3']['interval'][1]:+.2f}] {r['size_of_the_edge_from']}")


if __name__ == "__main__":
    main()
