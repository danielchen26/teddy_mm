"""Is the weak response to in-silico events a matter of input-token coverage? From records already in the repository.

An in-silico event on a gene can only be applied to a cell that carries that gene among its input tokens. For each target and
channel of the forward run (runs/forward_typed_events): on what share of the cells the event could be applied, how far it
moved Z there, and how that compares with the mean shift of Z under the real knockouts."""
import json

from _paths import RES, RUNS

fw = json.load(open(RUNS / "forward_typed_events" / "forward_record.json"))["datasets"]
J = json.load(open(RES / "join_closure_admission.json"))["datasets"]
out = {}
for ds, key in (("frangieh2021", next(k for k in fw if "frangieh" in k)), ("papalexi2021", next(k for k in fw if "papalexi" in k))):
    n_cells = fw[key]["n_confirmatory"]
    r = J[ds]["readouts"]["PDL1"]
    rows = {}
    for s in r["sites"]:
        rows.setdefault(s["site"], {})[s["channel"]] = {"n_applied": s["n"], "share_applied": s["n"] / n_cells if n_cells else None,
                                                         "z_move": s["z_move"], "resolved": s["resolved"]}
    out[ds] = {"n_cells": n_cells, "real_knockout_z_shift": r["ko_direction_norm"], "insilico_share_of_measured": r["insilico_share_of_measured"], "sites": rows}
    print(ds, "cells", n_cells, "| real knockouts move Z by", round(r["ko_direction_norm"], 2), "| in-silico at most", f"{100 * r['insilico_share_of_measured']:.1f}% of the measured PD-L1 change")
    for g, ch in rows.items():
        print("   ", g, " | ".join(f"{c}: applied on {v['n_applied']} cells, Z moves {v['z_move']:.3f} ({100 * v['z_move'] / r['ko_direction_norm']:.1f}% of the real shift)" for c, v in ch.items()))
json.dump(out, open(RES / "token_coverage.json", "w"), indent=1)
