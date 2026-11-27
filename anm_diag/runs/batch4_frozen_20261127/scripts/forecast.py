#!/usr/bin/env python3
"""The forecasts of batch 4 by three systems, before its values are opened. results/forecasts.json, once.

  python forecast.py --cells batch4_cells.tsv --states B2_STATES.npz B3_STATES.npz B4_STATES.npz --library LIBRARY.npz --packs batch3_pack_A.npz batch3_pack_B.npz

S0, the model alone: the part of a forecast that a readout of control cells reads from the state.
S2, version 5: runs/freeze_v5_20261126/scripts/_v5.py.
S3, version 6: runs/v6_candidate_20261126/scripts/_v6.py, the code of commit b4f6196.

Each stratum has one law over the 20 readouts of the panel. Its recorded rows are all its units read so far, on every
readout (the antibody columns opened for batch 3). The graph: the declared edges, and the edges admitted in run 43 and in
batch 3, each in the stratum that admitted it. The edges proposed earlier with no row to test them are listed here, so
that their test on the rows of the batch is fixed before those rows are opened; they enter no forecast.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _b4 as B  # noqa: E402
from _b4 import V, V6  # noqa: E402
from anm_diag import events as EV  # noqa: E402

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    a = ap.parse_args()
    out_file = B.RUN / "results" / "forecasts.json"
    if out_file.exists():
        raise SystemExit("forecasts.json exists: the forecasts were made once and are not made again")
    cells = B.table(a.cells)
    st, packs = B.States(a.states), B.Packs(a.packs)
    S = {g: B.stratum(g, cells, st) for g in B.STRATA}
    how = B.chosen()
    for g, c in S.items():
        new = c["kind"] == "new"
        assert packs.has(B.STUDY, c["ids"][~new]).all() and not packs.has(B.STUDY, c["ids"][new]).any(), "a value of the batch is at hand, or one the forecasts stand on is missing"
    d = np.load(a.library)
    lib_ids, un = d["ids"].astype(str), np.where(d["A"] == 0, "control", d["target"].astype(str))

    # --- the scale of every readout in every condition, from control cells
    scale = {r: {} for r in B.READOUTS}
    for g, ids in (("library", lib_ids[un == "control"]), *((g, S[g]["ids"][S[g]["kind"] == "control"]) for g in S)):
        clr, log, bg = packs.readouts(B.STUDY, ids)
        for r in B.READOUTS:
            scale[r][g] = V.scale(log[r], bg, clr[r])
    sc_of = lambda g: {r: scale[r][g] for r in B.READOUTS}

    # --- the event graph
    sites = sorted((set(un.tolist()) | {u for g in S for u in S[g]["unit"].tolist()}) - {"control"})
    declared, tab = V.declared(sites)
    on_record, waiting = B.edges_on_record()
    admitted = [e for e in on_record if e["status"] == "admitted"]
    tab_of = {g: {**tab, **{(e["site"], e["readout"]): e["class"] for e in admitted if e["stratum"] == g}} for g in S}

    # --- the library under interferon, by both versions
    YL, _log, _bg = packs.readouts(B.STUDY, lib_ids)
    ZL = V.unit_length(d["Z"].astype(np.float64))
    lib5, lib6 = V.library_law(ZL, YL, un, tab, sc_of("library")), V6.library_law(ZL, YL, un, tab, sc_of("library"))
    flags = lib6["law"].flags()
    proposed = EV.propose_edges(flags, known=declared + on_record + waiting)
    to_test = waiting + proposed
    print(f"library: spread with no edge {lib5['law'].sigma0:.3f} (version 5), {lib6['law'].sigma0:.3f} (version 6); width {lib6['width']['no type']}", flush=True)
    print("  flags of version 6:", [(f["event"], f["readout"]) for f in flags], "| proposed here:", [e["edge_id"] for e in proposed], flush=True)

    # --- the laws of the strata: every unit read so far, every readout
    law5, law6, known_rec = {}, {}, {}
    for g, c in S.items():
        rec = c["kind"] != "new"
        Y, _l, _b = packs.readouts(B.STUDY, c["ids"][rec])
        law5[g] = V.stratum_law(c["Z"][rec], Y, c["unit"][rec], tab_of[g], sc_of(g), lib5)
        laws, known_rec[g], _first = V6.stratum_laws({"all": {"x": c["Z"][rec], "Y": Y, "unit": c["unit"][rec], "table": tab_of[g], "scales": sc_of(g)}}, lib6)
        law6[g] = laws["all"]
    twice5, twice6 = V.measured_twice(lib5, list(law5.values())), V.measured_twice(lib6, list(law6.values()))

    res = {"cells_file_sha256": sha(a.cells), "states_sha256": [sha(p) for p in a.states], "packs_sha256": [sha(p) for p in a.packs], "units": sum(len(v) for v in how.values()),
           "code": {"S2": "d777070, runs/freeze_v5_20261126/scripts/_v5.py", "S3": "b4f6196, runs/v6_candidate_20261126/scripts/_v6.py"},
           "scale": scale, "edges": {"declared": declared, "on_record": on_record, "proposed_here": proposed, "to_test_on_the_rows_of_the_batch": to_test},
           "library": {"events": len(lib6["law"].events), "spread_with_no_edge": {"S2": lib5["law"].sigma0, "S3": lib6["law"].sigma0}, "width": lib6["width"], "flags": [f"{f['event']}->{f['readout']}" for f in flags],
                       "shares": [{"event": e, "readout": r, "share": v} for (e, r), v in lib6["shares"].items()]},
           "pairs_measured_twice": {"S2": {"pairs": len(twice5), "spread_of_a_carried_share": EV.carried_spread(twice5)}, "S3": {"pairs": len(twice6), "spread_of_a_carried_share": EV.carried_spread(twice6)}},
           "strata": {}, "rows": []}
    for g, c in S.items():
        new = c["kind"] == "new"
        r5 = V.forecast(*law5[g], c["Z"][new], c["unit"][new], lib5, twice5)
        first = V6.forecast(*law6[g], c["Z"][new], c["unit"][new], lib6, twice6)
        known = EV.derived(V6.edge_parts(first), V6.spreads(law6[g][0]), V6.OVER)
        r6 = V6.forecast(*law6[g], c["Z"][new], c["unit"][new], lib6, twice6, known=known)
        by5 = {(r["event"], r["readout"]): r for r in r5}
        assert len(by5) == len(r5) == len(r6)
        for r in r6:
            o = by5[(r["event"], r["readout"])]
            assert abs(o["read"] - r["read"]) < 1e-9 and r["no_forecast"] is None and o["no_forecast"] is None
            untyped = r["type"] is None or r["type"].startswith("(")
            res["rows"].append({"stratum": g, "event": r["event"], "chosen": how[g][r["event"]], "readout": r["readout"], "cells": r["cells"], "edge": r["edge"], "type": None if untyped else r["type"],
                                "scale": r["scale"] if r["edge"] is not None and not r["edge"].endswith(EV.AT_THE_BACKGROUND) else None, "size_of_the_edge_from": r["size_of_the_edge_from"],
                                "shift_times_noise": r["shift_times_noise"], "S0": {"forecast": r["read"]},
                                "S2": {"forecast": o["forecast"], "sd": o["sd"], "interval": o["interval"], "from_the_recorded_events": o["from_the_recorded_events"], "by_its_edge": o["by_its_edge"]},
                                "S3": {"forecast": r["forecast"], "sd": r["sd"], "interval": r["interval"], "width_from": r["width_from"], "from_the_recorded_events": r["from_the_recorded_events"],
                                       "by_its_edge": r["by_its_edge"], "known_by_declaration": r["known_by_declaration"], "parts_of_the_spread": r["parts_of_the_spread"]}})
        l5, l6 = law5[g][0], law6[g][0]
        res["strata"][g] = {"recorded_events": l6.events, "control_cells": int((c["kind"] == "control").sum()), "new_units": sorted(how[g]),
                            "spread_with_no_edge": {"S2": l5.sigma0, "S3": l6.sigma0}, "spread_from": {"S2": l5.spread_from, "S3": l6.spread_from}, "width": l6.conformal(),
                            "recorded_pairs_with_a_known_part": sum(1 for (e, _r) in known_rec[g] if e in l6.events), "new_pairs_with_a_known_part": sum(1 for v in known.values() if abs(v[0]) > 0)}
    res["pairs"] = len(res["rows"])
    key = lambda r: f"{r['stratum']}|{r['event']}|{r['readout']}"
    res["sets"] = {"by lot": [key(r) for r in res["rows"] if r["chosen"] == "by lot"], "by an item on record": [key(r) for r in res["rows"] if r["chosen"] != "by lot"],
                   "with an edge": [key(r) for r in res["rows"] if r["edge"] is not None]}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    rows = res["rows"]
    print("pairs forecast:", res["pairs"], {g: sum(r["stratum"] == g for r in rows) for g in S}, "| sets", {k: len(v) for k, v in res["sets"].items()})
    for g, v in res["strata"].items():
        print(f"  {g:7s} recorded {len(v['recorded_events'])} | spread S2 {v['spread_with_no_edge']['S2']:.3f} ({v['spread_from']['S2']}) | S3 {v['spread_with_no_edge']['S3']:.3f} | width {v['width']['no type']} | known parts: recorded {v['recorded_pairs_with_a_known_part']}, new {v['new_pairs_with_a_known_part']}")
        hw = lambda s: np.median([(r[s]["interval"][1] - r[s]["interval"][0]) / 2 for r in rows if r["stratum"] == g and r["edge"] is None])
        print(f"          median half-width of a pair with no edge: S2 {hw('S2'):.3f}, S3 {hw('S3'):.3f} | largest difference of the two forecasts {max(abs(r['S2']['forecast'] - r['S3']['forecast']) for r in rows if r['stratum'] == g):.3f}")
    print("pairs measured twice:", res["pairs_measured_twice"])
    print("edges to test on the rows of the batch:", [e["edge_id"] for e in to_test])
    for r in rows:
        if r["edge"] is not None:
            print(f"  {r['stratum']:7s} {r['event']:6s}->{r['readout']:6s} {r['edge']:26s} scale {r['scale']} | S0 {r['S0']['forecast']:+.2f} | S2 {r['S2']['forecast']:+.2f} | S3 {r['S3']['forecast']:+.2f} [{r['S3']['interval'][0]:+.2f}, {r['S3']['interval'][1]:+.2f}] {r['size_of_the_edge_from']}")


if __name__ == "__main__":
    main()
