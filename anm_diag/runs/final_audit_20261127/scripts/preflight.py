#!/usr/bin/env python3
"""Before anything of the final audit is opened: the code of its forecasts, on cells read before. results/preflight.json.

  python preflight.py --cells batch4_cells.tsv --states MEL_STATES... --states-thp THP_STATES... --library LIBRARY.npz
                      --packs batch3_pack_A.npz batch3_pack_B.npz --pack4 batch4_pack.npz

1. The same forecasts. `forecast.core`, given the cells of batch 4 as they were when batch 4 was forecast, has to give
   the committed forecasts of batch 4 again, number for number: the code that forecasts the final audit is the code
   that was frozen, and nothing changed it since.
2. A reference that is not the law's control cells. The final units are read against control cells the laws never
   saw. Batch 4 again, with one control cell in three taken out of everything and used as the reference (by a hash of
   its name): the share of pairs inside and the median of squared z should be what batch 4 gave.
3. THP-1. Its six final units are forecast by a law over four readouts. Units of THP-1 read in batch 2 are forecast
   the same way from the rest and compared with what was measured.
All of it on values read before: it finds faults of the code and is no evidence.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _fa as F  # noqa: E402
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("final_forecast", Path(__file__).resolve().parent / "forecast.py")      # by its path: the scripts of earlier runs have a forecast.py too
FC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(FC)

B4 = F.RUNS / "batch4_frozen_20261127" / "results"
third = lambda n: int(hashlib.sha256(("preflight|" + n).encode()).hexdigest(), 16) % 3 == 0


def measured(c, packs, study, e, r, sd):
    """As read.py: the change of the unit against the reference, in spreads of the readout among the law's control cells."""
    Y, _l, _b = packs.readouts(study, c["ids"], [r])
    y, m, ref = Y[r], c["unit"] == e, c["ref"]
    return float((y[m].mean() - y[ref].mean()) / sd)


def score(rows, S, packs, spread):
    out = {}
    for s in ("S2", "S3"):
        z2, inside = [], 0
        for r in rows:
            ch = measured(S[r["stratum"]], packs, F.STRATA[r["stratum"]][0], r["event"], r["readout"], spread[r["stratum"]][r["readout"]])
            z2.append(((ch - r[s]["forecast"]) / r[s]["sd"]) ** 2)
            inside += r[s]["interval"][0] <= ch <= r[s]["interval"][1]
        out[s] = {"pairs": len(rows), "inside": int(inside), "share_inside": inside / len(rows), "median_z2": float(np.median(z2)), "mean_z2": float(np.mean(z2))}
    return out


def spreads(S, packs):
    out = {}
    for g, c in S.items():
        C = c["kind"] == "control"
        Y, _l, _b = packs.readouts(F.STRATA[g][0], c["ids"][C])
        out[g] = {r: float((y - y.mean()).std()) for r, y in Y.items()}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--states-thp", type=Path, nargs="+", required=True)
    ap.add_argument("--library", type=Path, required=True)
    ap.add_argument("--packs", type=Path, nargs="+", required=True)
    ap.add_argument("--pack4", type=Path, required=True)
    a = ap.parse_args()
    st = {"frangieh2021": F.States(a.states), "papalexi2021": F.States(a.states_thp)}
    cells = F.table(a.cells)
    before, after = F.Packs(a.packs), F.Packs(a.packs + [a.pack4])
    committed = json.loads((B4 / "forecasts.json").read_text())
    units4 = {g: sorted({r["event"] for r in committed["rows"] if r["stratum"] == g}) for g in F.MEL}
    record4 = [e for e in F.edges_on_record() if not e["edge_id"].startswith("TXNDC17") and not (e["edge_id"] == "DAG1->PDL1@mel_co")]   # the edges on record when batch 4 was forecast
    res = {}
    # --- 1. the committed forecasts of batch 4 again
    r1 = FC.core(cells, st, before, a.library, units4, record4)
    by = {FC.key(r): r for r in r1["rows"]}
    worst = 0.0
    for r in committed["rows"]:
        o = by[FC.key(r)]
        for s, k in (("S0", "forecast"), ("S2", "forecast"), ("S2", "sd"), ("S3", "forecast"), ("S3", "sd")):
            worst = max(worst, abs(o[s][k] - r[s][k]))
        worst = max(worst, abs(o["S3"]["interval"][0] - r["S3"]["interval"][0]), abs(o["S3"]["interval"][1] - r["S3"]["interval"][1]))
    res["1_the_forecasts_of_batch_4_again"] = {"pairs": len(committed["rows"]), "largest_difference": worst, "the_same": bool(worst < 1e-9 and len(by) == len(committed["rows"]))}
    print("1. batch 4 again: largest difference", worst, flush=True)
    assert res["1_the_forecasts_of_batch_4_again"]["the_same"], "the code does not give the committed forecasts of batch 4"
    # --- 2. a reference the laws never saw
    S = {g: F.stratum(g, cells, st) for g in F.STRATA}
    sp = spreads(S, after)
    res["2_a_reference_the_laws_never_saw"] = {"as_batch_4_was_read": score(r1["rows"], S, after, sp)}
    alt = {k: v.astype(object) for k, v in cells.items()}
    move = (alt["role"] == "control") & (alt["study"] == "frangieh2021") & np.array([third(n) for n in alt["cell"]])
    alt["role"][move] = "reference"
    r2 = FC.core(alt, st, before, a.library, units4, record4, check=False)
    S2 = {g: F.stratum(g, alt, st) for g in F.STRATA}
    assert all(S2[g]["own_reference"] and abs(S2[g]["ref"].sum() / (S[g]["kind"] == "control").sum() - 1 / 3) < 0.03 for g in F.MEL), "the reference is not one control cell in three"
    res["2_a_reference_the_laws_never_saw"]["one_control_cell_in_three_as_the_reference"] = {**score(r2["rows"], S2, after, spreads(S2, after)), "reference_cells": {g: int(S2[g]["ref"].sum()) for g in F.MEL}}
    print("2.", json.dumps(res["2_a_reference_the_laws_never_saw"]), flush=True)
    # --- 3. THP-1: units read in batch 2, forecast again from the rest
    c = S["thp1"]
    b2 = sorted(set(c["unit"][c["role"] == "batch2"].tolist()))
    pick = b2[::max(len(b2) // 6, 1)][:6]
    alt = {k: v.astype(object) for k, v in cells.items()}
    alt["role"][(alt["study"] == "papalexi2021") & np.isin(alt["target"], pick)] = "event"
    alt["role"][(alt["study"] == "frangieh2021") & (alt["role"] == "event")] = "batch4"
    r3 = FC.core(alt, st, after, a.library, {"thp1": pick}, F.edges_on_record(), check=False)
    S3 = {g: F.stratum(g, alt, st) for g in F.STRATA}
    res["3_thp1_units_read_before_forecast_again"] = {"units": pick, **score(r3["rows"], S3, after, spreads(S3, after)), "width_from": sorted({r["S3"]["width_from"] for r in r3["rows"]})}
    print("3.", json.dumps(res["3_thp1_units_read_before_forecast_again"]), flush=True)
    (F.RUN / "results" / "preflight.json").write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
