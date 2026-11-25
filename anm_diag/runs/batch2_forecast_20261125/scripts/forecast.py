#!/usr/bin/env python3
"""The forecasts of batch 2, from RNA alone. results/forecasts.json. Made once, before any protein value of the batch is opened.

  python forecast.py --cells batch2_cells.tsv --states STATES_DIR --read read_cells.npz --library LIBRARY.npz

For every stratum the law of its recorded rows (anm_diag.account.Law at the state): the readout of its control cells, the
recorded events (development cells and batch 1) with their measured changes, and the edges known beforehand. Then, for every
unit of the batch and every readout, Law.new on the states of its cells: the part the control readout reads, what the
recorded events add through the state, the part by its edge, the forecast and its 95% interval.

What a stratum cannot tell from its own rows it takes from the library under interferon (206 events, all development): the
spread of a direct part without an edge when it has fewer than 8 recorded pairs without an edge, and the direct parts of the
class of edges 'own' (pooled with the recorded pairs of that class in the strata of the same study).
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from _b2 import AC, EDGES, MIN_ROWS, OWN, ROLES, RUN, STRATA, stratum, table, unit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, required=True)
    ap.add_argument("--read", type=Path, required=True)
    ap.add_argument("--library", type=Path, required=True)
    a = ap.parse_args()
    out_file = RUN / "results" / "forecasts.json"
    if out_file.exists():
        raise SystemExit("forecasts.json exists: the forecasts were made once and are not made again")
    cells = table(a.cells)
    # --- the library under interferon: the spread without an edge, and the class of edges
    d = np.load(a.library)
    un = np.where(d["A"] == 0, "control", d["target"].astype(str))
    lib = AC.Law(unit(d["Z"]).astype(np.float64), {k: d["readout__" + k].astype(np.float64) for k in OWN}, un, edges=EDGES).fit()
    res = {"cells_file_sha256": hashlib.sha256(a.cells.read_bytes()).hexdigest(), "edges_known_beforehand": OWN,
           "library": {"events": len(lib.events), "spread_of_a_direct_part_without_an_edge": lib.sigma0, "class_own": [{"event": e, "readout": r, "direct": v} for e, r, v in lib.class_values()["own"]]},
           "strata": {}}
    laws = {}
    for g, (study, cond, readouts) in STRATA.items():
        c = stratum(g, cells, a.states, read=a.read)
        rec = c["role"] != "event"
        edges = EDGES if study == "frangieh2021" else {}
        law = AC.Law(c["Z"][rec], {r: c[f"y_{r}"][rec] for r in readouts}, c["unit"][rec], edges=edges)
        free = sum(1 for e in law.events for r in readouts if (e, r) not in edges)
        law.fit(sigma0=None if free >= MIN_ROWS else lib.sigma0)
        laws[g] = (law, c)
    pooled = [v for _e, _r, v in lib.class_values()["own"]]
    for g, (law, c) in laws.items():
        if STRATA[g][0] == "frangieh2021":
            pooled += [v for _e, _r, v in law.class_values().get("own", [])]
    for g, (law, c) in laws.items():
        study, cond, readouts = STRATA[g]
        edges = EDGES if study == "frangieh2021" else {}
        new = c["role"] == "event"
        genes = np.char.replace(c["unit"][new], "new:", "")
        rows = law.new(c["Z"][new], genes, edges=edges, classes={"own": pooled})
        for r in rows:
            r["unit"] = f"{r['event']}@{study} | {cond}"
            r["role"] = ROLES["units"][r["unit"]]
        rec_rows = law.forecast()
        s = law.summary()
        res["strata"][g] = {"study": study, "condition": cond, "control_cells": int((c["unit"] == "control").sum()), "recorded_events": law.events,
                            "pairs_without_an_edge_among_the_recorded": sum(1 for e in law.events for r in readouts if (e, r) not in edges),
                            "spread_of_a_direct_part_without_an_edge": {"value": law.sigma0, "from": "the library under interferon" if law.fixed else "the stratum's own rows"},
                            "widening_of_the_control_prior": s["widening_of_the_control_prior"],
                            "recorded_pairs_each_left_out": [{k: r[k] for k in ("event", "readout", "cells", "edge", "change", "se", "read", "from_the_other_events", "by_its_edge", "forecast", "sd", "z")} for r in rec_rows],
                            "forecasts": rows}
        print(f"{g}: {len(law.events)} recorded events; spread without an edge {law.sigma0:.3f} ({'library' if law.fixed else 'own rows'}); widening {s['widening_of_the_control_prior']}")
        for r in rows:
            print(f"   {r['event']:9s} -> {r['readout']:5s} {r['role']:11s} cells {r['cells']:4d} shift {r['shift_times_noise']:5.1f}x | read {r['read']:+.2f} | recorded {r['from_the_recorded_events']:+.2f} | edge {r['by_its_edge']:+.2f} "
                  f"| forecast {r['forecast']:+.2f} [{r['interval'][0]:+.2f}, {r['interval'][1]:+.2f}]")
    res["class_own_pooled"] = {"values": pooled, "mean": float(np.mean(pooled)), "spread_of_a_new_pair": float(np.std(pooled, ddof=1) * np.sqrt(1 + 1 / len(pooled)))}
    res["pairs"] = sum(len(v["forecasts"]) for v in res["strata"].values())
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("pairs forecast:", res["pairs"], "| class own pooled:", len(pooled), "values, mean", round(res["class_own_pooled"]["mean"], 2))


if __name__ == "__main__":
    main()
