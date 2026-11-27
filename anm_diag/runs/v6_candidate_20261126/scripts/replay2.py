#!/usr/bin/env python3
"""Batch 2 again, by the candidate for version 6. RETROSPECTIVE: every value was read before. results/replay_batch2.json.

  python replay2.py --cells batch2_cells.tsv --states STATES_DIR --read read_cells.npz --pack batch2_pack.npz --library LIBRARY.npz --scale scale_pack.npz

The 79 pairs of runs/batch2_forecast_20261125, assembled as runs/freeze_v5_20261126/scripts/check.py assembles them,
with the candidate in place of version 5. Only the three readouts read then are used, so the arithmetic of the
readout knows of the edges to those three alone.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN.parent / "batch2_forecast_20261125" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _v6 as V6  # noqa: E402
from _v6 import V  # noqa: E402
from _b2 import STRATA, stratum, table  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402

CONDITION = {"library": "IFNγ", "mel_co": "Co-culture", "mel_ctl": "Control"}


def main():
    ap = argparse.ArgumentParser()
    for k in ("cells", "states", "read", "pack", "library", "scale"):
        ap.add_argument("--" + k, type=Path, required=True)
    a = ap.parse_args()
    P = np.load(a.scale)
    cond = P["condition"].astype(str)
    sc = {r: {c: V.scale(P[f"log__{r}"][cond == c], P["background"][cond == c], P[f"log__{r}"][cond == c] - P["centre"][cond == c], int(P["panel"])) for c in sorted(set(cond.tolist()))}
          for r in V.USED["melanoma"]}
    d = np.load(a.library)
    un = np.where(d["A"] == 0, "control", d["target"].astype(str))
    cells = table(a.cells)
    strata = {g: stratum(g, cells, a.states, read=a.read, new=a.pack) for g in STRATA}
    sites = sorted((set(un.tolist()) | {u.replace("new:", "") for c in strata.values() for u in c["unit"].tolist()}) - {"control"})
    _edges, tab = V.declared(sites)
    over = {r: V6.OVER[r] for r in V.USED["melanoma"]}
    ZL, YL = V.unit_length(d["Z"].astype(np.float64)), {k: d["readout__" + k].astype(np.float64) for k in V.USED["melanoma"]}
    lib5 = V.library_law(ZL, YL, un, tab, {r: sc[r]["IFNγ"] for r in sc})
    lib6 = V6.library_law(ZL, YL, un, tab, {r: sc[r]["IFNγ"] for r in sc}, over=over)
    out = {"note": "retrospective: every value here was read before; nothing is a forecast", "strata": {}}
    rows5, rows6, law5, law6 = [], [], {}, {}
    for g, (study, _c, readouts) in STRATA.items():
        s = strata[g]
        rec, mel = s["role"] != "event", study == "frangieh2021"
        Y = {r: s[f"y_{r}"][rec] for r in readouts}
        scales = {r: sc[r][CONDITION[g]] for r in sc} if mel else {}
        law5[g] = V.stratum_law(s["Z"][rec], Y, s["unit"][rec], tab if mel else {}, scales, lib5)
        laws, _known, _first = V6.stratum_laws({"all": {"x": s["Z"][rec], "Y": Y, "unit": s["unit"][rec], "table": tab if mel else {}, "scales": scales}}, lib6, over=over if mel else {})
        law6[g] = laws["all"]
    twice5, twice6 = V.measured_twice(lib5, [law5["mel_co"], law5["mel_ctl"]]), V.measured_twice(lib6, [law6["mel_co"], law6["mel_ctl"]])
    for g, (study, _c, readouts) in STRATA.items():
        s = strata[g]
        new, mel = s["role"] == "event", study == "frangieh2021"
        genes = np.char.replace(s["unit"][new], "new:", "")
        Yn = {r: s[f"y_{r}"][new] for r in readouts}
        r5 = V.scored(V.forecast(*law5[g], s["Z"][new], genes, lib5, twice5), law5[g][0], Yn, genes)
        first = V6.forecast(*law6[g], s["Z"][new], genes, lib6, twice6)
        known = EV.derived(V6.edge_parts(first), V6.spreads(law6[g][0]), over if mel else {})
        r6 = V.scored(V6.forecast(*law6[g], s["Z"][new], genes, lib6, twice6, known=known), law6[g][0], Yn, genes)
        rows5 += [dict(r, stratum=g) for r in r5]; rows6 += [dict(r, stratum=g) for r in r6]
        out["strata"][g] = {"version 5": V.calibration(r5), "candidate": V.calibration(r6), "spread_with_no_edge": {"version 5": law5[g][0].sigma0, "candidate": law6[g][0].sigma0},
                            "width_from": sorted({r["width_from"] for r in r6 if r["edge"] is None})}
    out["all"] = {"version 5": V.calibration(rows5), "candidate": V.calibration(rows6), "as_committed_for_the_batch": {"pairs": 79, "inside": 77, "mean_z2": 0.30, "error": 0.231}}
    out["pairs_with_an_edge"] = [{k: r.get(k) for k in ("stratum", "event", "readout", "edge", "change", "forecast", "interval", "inside")} for r in rows6 if r["edge"] is not None]
    (RUN / "results" / "replay_batch2.json").write_text(json.dumps(out, indent=1, default=float) + "\n")
    for k in ("version 5", "candidate"):
        c = out["all"][k]
        print(f"{k:10s} pairs {c['pairs']} inside {c['inside']} mean z2 {c['mean_z2']:.2f} median z2 {c['median_z2']:.2f} error {c['error']['forecast']:.3f} (model alone {c['error']['control_readout']:.3f})")
    for g, v in out["strata"].items():
        print(f"  {g:7s} inside {v['version 5']['inside']} -> {v['candidate']['inside']} of {v['candidate']['pairs']} | spread {v['spread_with_no_edge']['version 5']:.3f} -> {v['spread_with_no_edge']['candidate']:.3f} | {v['width_from']}")


if __name__ == "__main__":
    main()
