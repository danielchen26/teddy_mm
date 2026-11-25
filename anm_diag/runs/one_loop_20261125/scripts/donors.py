#!/usr/bin/env python3
"""The account of eight donors, each against the other seven: events with no control cells. EXPLORED. results/donors.json.

  python donors.py --donors DONORS.npz [--draws 200]

DONORS: the cells of run 32 (runs/unit_events_20261123) that were read there, two thirds of eight donors of one study:
unit (the donor of a cell), RNA512 (the input view), Z (the state), Y (the measured proteins), frozen (the model's own
output). The pairs of a donor and a protein are the 26 that run 32 admitted on the third it kept apart for selection
(its results/selection.json); that third is not here.

Every pair gets the account of rule 4i: the donor read against the cells of the other seven (control=None), the arms
resampled by donor, the rows of the other donors as `others`; then the step after a rejection of rule 4l
(account.plan, outer.open_parts), and for a part that lies outside the input and that the other donors do not forecast,
a proposed edge between the donor and the protein (events.edge, basis "account").

This was the first run of the code on real donors and was seen before the design was written: it has no hypothesis and
is counted nowhere. A donor is also a batch of the measurement in this study: its change cannot be told from its
batch's.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN.parents[1]))
from anm_diag import account as AC  # noqa: E402
from anm_diag import cycle  # noqa: E402
from anm_diag import events as EV  # noqa: E402
from anm_diag import outer  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", type=Path, required=True)
    ap.add_argument("--draws", type=int, default=200)
    a = ap.parse_args()
    sel = json.loads((RUN.parent / "unit_events_20261123" / "results" / "selection.json").read_text())
    d = np.load(a.donors)
    un = d["unit"].astype(str)
    donors = sorted(set(un.tolist()))
    first = cycle.halves(un, 0)
    unit = lambda x: x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-12)
    stages = {"V": d["RNA512"].astype(np.float64), "Z": unit(d["Z"].astype(np.float64))}
    res = {"note": "explored: seen before the design was written; no hypothesis; counted nowhere", "cells": int(len(un)), "donors": donors, "draws": a.draws, "pairs": {}}
    accounts, t0 = [], time.time()
    for name, v in sel["proteins"].items():
        for e in v["admitted"]:
            j = v["index"]
            acc = AC.account(stages, d["Y"][:, j].astype(np.float64), un, first, e, others=[u for u in donors if u != e], inputs=["V"], state="Z", o=d["frozen"][:, j].astype(np.float64),
                             control=None, cluster=un, draws=a.draws, pairs=[("Z", "V")], controls=False)
            z, ch = acc["stages"]["Z"], acc["chain"]
            res["pairs"][f"{name}|{e}"] = {"protein": name, "donor": e, "change": z["change"], "snr": z["snr"], "cells": z["cells"], "read_against": z["read_against"], "overlap": {s: acc["stages"][s]["overlap"] for s in stages},
                                           "shares": {s: {k: {"share": x["share"], "ci": x["ci"]} for k, x in acc["stages"][s]["shares"].items()} for s in stages},
                                           "parts": ch["parts"], "open": ch["open"], "unresolved": ch["unresolved"], "across_events": ch.get("across_events"), "says": ch["says"], "plan": acc["plan"]}
            accounts.append({"readout": name, "group": "hao2021", "type": "donor", "event": e, "account": acc})
            print(f"{name:10s} {e[-2:]} change {z['change']:+.2f} ({z['snr']:.0f} errors) | read {z['shares']['read']['share']:4.0%} | carried V {acc['stages']['V']['shares']['carried']['share']:4.0%} Z {z['shares']['carried']['share']:4.0%}"
                  f" | others {z['shares']['others']['share']:5.0%} | {ch['says']}", flush=True)
    items = outer.open_parts(accounts)
    proposed = [EV.edge(i["event"], i["readout"], "unit", basis="account", stratum="hao2021", evidence={"share_outside_the_input": i["value"], "ci": i["ci"]},
                        provenance="the account places the change outside the input, and the other donors do not forecast it")
                for i in items if i["reading"] == "outside the input, not forecast" and i["state"] == "open"]
    res["step_after_a_rejection"] = {"items": items, "by_reading_and_state": {}, "edges_proposed": proposed, "edge_checks": EV.check_edges(proposed),
                                     "note": "a proposed edge enters no forecast; it is tested by rows it did not see (the same donor measured again, or in another batch), and this study has none"}
    for i in items:
        k = f"{i['reading']} | {i['state']}"
        res["step_after_a_rejection"]["by_reading_and_state"][k] = res["step_after_a_rejection"]["by_reading_and_state"].get(k, 0) + 1
    kinds = {}
    for p in res["pairs"].values():
        k = "nothing open" if not p["open"] else " + ".join(sorted(x["part"] for x in p["open"]))
        kinds[k] = kinds.get(k, 0) + 1
    res["pairs_by_open_part"] = kinds
    res["seconds"] = round(time.time() - t0)
    print("pairs by open part:", kinds)
    print("the step after a rejection:", res["step_after_a_rejection"]["by_reading_and_state"], "| edges proposed:", len(proposed))
    (RUN / "results" / "donors.json").write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
