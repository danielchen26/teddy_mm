#!/usr/bin/env python3
"""The frozen system on cells read before. RETROSPECTIVE: nothing here is a forecast. results/check.json.

  python check.py --cells batch2_cells.tsv --states STATES_DIR --read read_cells.npz --pack batch2_pack.npz
                  --library LIBRARY.npz --scale scale_pack.npz --donors DONORS.npz

1. Batch 2 again (runs/batch2_forecast_20261125), by the system of version 5 exactly as batch 3 will use it (_v5.py):
   the calibration of its 79 pairs by stratum, and its eight own-gene pairs.
2. Eight donors with no control cells (runs/unit_events_20261123): each donor left out, forecast from its state and the
   rows of the other seven (`Law(control=None)`), for the proteins of the 26 pairs run 32 admitted.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN.parent / "batch2_forecast_20261125" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _v5 as V  # noqa: E402
from _b2 import STRATA, stratum, table  # noqa: E402
from anm_diag import account as AC  # noqa: E402

CONDITION = {"library": "IFNγ", "mel_co": "Co-culture", "mel_ctl": "Control"}


def main():
    ap = argparse.ArgumentParser()
    for k in ("cells", "states", "read", "pack", "library", "scale", "donors"):
        ap.add_argument("--" + k, type=Path, required=True)
    a = ap.parse_args()
    res = {"note": "retrospective: every value here was read before; nothing is a forecast"}
    P = np.load(a.scale)
    cond = P["condition"].astype(str)
    sc = {r: {c: V.scale(P[f"log__{r}"][cond == c], P["background"][cond == c], P[f"log__{r}"][cond == c] - P["centre"][cond == c], int(P["panel"])) for c in sorted(set(cond.tolist()))}
          for r in V.USED["melanoma"]}
    res["scale"] = sc
    d = np.load(a.library)
    un = np.where(d["A"] == 0, "control", d["target"].astype(str))
    cells = table(a.cells)
    strata = {g: stratum(g, cells, a.states, read=a.read, new=a.pack) for g in STRATA}
    sites = sorted((set(un.tolist()) | {u.replace("new:", "") for c in strata.values() for u in c["unit"].tolist()}) - {"control"})
    _edges, tab = V.declared(sites)
    lib = V.library_law(V.unit_length(d["Z"].astype(np.float64)), {k: d["readout__" + k].astype(np.float64) for k in V.USED["melanoma"]}, un, tab, {r: sc[r]["IFNγ"] for r in sc})
    res["library"] = {"events": len(lib["law"].events), "spread_of_a_direct_part_without_an_edge": lib["law"].sigma0, "spread_from": lib["law"].spread_from, "share_of_pairs_with_a_part_of_their_own": lib["law"].own_share,
                      "classes": {c: len(v) for c, v in lib["classes"].items()}, "shares": {f"{e}->{r}": v for (e, r), v in lib["shares"].items()},
                      "calibration": {k: v for k, v in lib["law"].calibration().items() if k in ("pairs", "inside_95", "mean_z2")}, "flags": [(f["event"], f["readout"], f["flag"]) for f in lib["law"].flags()]}
    print("library:", {k: res["library"][k] for k in ("spread_of_a_direct_part_without_an_edge", "spread_from", "share_of_pairs_with_a_part_of_their_own", "calibration")})
    print("  classes", res["library"]["classes"], "| shares", {k: round(v, 2) for k, v in res["library"]["shares"].items()}, "| flags", res["library"]["flags"])
    laws, allrows = {}, []
    res["strata"] = {}
    for g, (study, c, readouts) in STRATA.items():
        s = strata[g]
        rec, new = s["role"] != "event", s["role"] == "event"
        mel = study == "frangieh2021"
        law, ed, scl = V.stratum_law(s["Z"][rec], {r: s[f"y_{r}"][rec] for r in readouts}, s["unit"][rec], tab if mel else {}, {r: sc[r][CONDITION[g]] for r in sc} if mel else {}, lib)
        laws[g] = (law, ed, scl)
    twice_before = V.measured_twice(lib, [laws["mel_co"], laws["mel_ctl"]])
    for g, (study, c, readouts) in STRATA.items():
        s, (law, ed, scl) = strata[g], laws[g]
        new = s["role"] == "event"
        genes = np.char.replace(s["unit"][new], "new:", "")
        rows = V.scored(V.forecast(law, ed, scl, s["Z"][new], genes, lib, twice_before), law, {r: s[f"y_{r}"][new] for r in readouts}, genes)
        allrows += [dict(r, stratum=g) for r in rows]
        res["strata"][g] = {"spread_of_a_direct_part_without_an_edge": law.sigma0, "spread_from": law.spread_from, "widening_of_the_control_prior": law.summary()["widening_of_the_control_prior"],
                            "calibration": V.calibration(rows), "pairs_without_an_edge": V.calibration([r for r in rows if r["edge"] is None]), "pairs_of_the_batch": rows}
        print(g, "spread", round(law.sigma0, 3), law.spread_from, "|", {k: (round(v, 2) if isinstance(v, float) else v) for k, v in res["strata"][g]["calibration"].items() if k != "error"}, res["strata"][g]["calibration"]["error"])
    res["batch_2"] = {"as_committed": {"pairs": 79, "inside": 77, "mean_z2": 0.30}, "all": V.calibration(allrows), "pairs_without_an_edge": V.calibration([r for r in allrows if r["edge"] is None]),
                      "pairs_measured_twice_before_the_batch": len(twice_before)}
    own = [r for r in allrows if r["edge"] is not None]
    res["batch_2"]["pairs_with_an_edge"] = [{k: r.get(k) for k in ("stratum", "event", "readout", "edge", "scale", "size_of_the_edge_from", "change", "forecast", "interval", "inside", "no_forecast")} for r in own]
    print("batch 2, all:", res["batch_2"]["all"]); print("   no edge:", res["batch_2"]["pairs_without_an_edge"])
    for r in own:
        print(f"   {r['stratum']:7s} {r['event']:6s}->{r['readout']:5s} {r['edge']:24s} change {r['change']:+.2f} forecast {r['forecast']:+.2f} [{r['interval'][0]:+.2f}, {r['interval'][1]:+.2f}] {r['inside']} | {r['size_of_the_edge_from']}")
    # what the system holds for batch 3: the units of batch 2 are recorded rows now
    after = []
    for g in ("mel_co", "mel_ctl"):
        s = strata[g]
        law, ed, scl = V.stratum_law(s["Z"], {r: s[f"y_{r}"] for r in STRATA[g][2]}, np.char.replace(s["unit"], "new:", ""), tab, {r: sc[r][CONDITION[g]] for r in sc}, lib)
        after.append((law, ed, scl))
    twice = V.measured_twice(lib, after)
    from anm_diag import events as EV
    res["after_batch_2"] = {"pairs_measured_twice": len(twice), "spread_of_a_carried_share": EV.carried_spread(twice), "pairs": [[round(x, 3), round(y, 3)] for x, y in twice]}
    print("pairs measured twice now:", len(twice), "spread of a carried share:", res["after_batch_2"]["spread_of_a_carried_share"])

    # --- donors: each left out
    sel = json.loads((RUN.parent / "unit_events_20261123" / "results" / "selection.json").read_text())
    D = np.load(a.donors)
    du = D["unit"].astype(str)
    donors = sorted(set(du.tolist()))
    Z = V.unit_length(D["Z"].astype(np.float64))
    prot = {name: v["index"] for name, v in sel["proteins"].items() if v["admitted"]}
    admitted = {(name, e) for name, v in sel["proteins"].items() for e in v["admitted"]}
    rows = []
    for u in donors:
        keep = du != u
        law = AC.Law(Z[keep], {n: D["Y"][keep, j].astype(np.float64) for n, j in prot.items()}, du[keep], control=None).fit()
        rows += V.scored(law.new(Z[~keep], du[~keep]), law, {n: D["Y"][~keep, j].astype(np.float64) for n, j in prot.items()}, du[~keep])
        print("donor", u[-2:], "left out: spread", round(law.sigma0, 3), "share with a part of their own", round(law.own_share, 2), flush=True)
    adm = [r for r in rows if (r["readout"], r["event"]) in admitted]
    res["donors"] = {"donors": len(donors), "proteins": len(prot), "all_pairs": V.calibration(rows), "the_26_admitted_pairs": V.calibration(adm),
                     "pairs": [{k: r[k] for k in ("event", "readout", "change", "se", "read", "from_the_recorded_events", "forecast", "sd", "inside")} for r in adm]}
    print("donors, all pairs:", res["donors"]["all_pairs"]); print("donors, the 26:", res["donors"]["the_26_admitted_pairs"])
    (RUN / "results" / "check.json").write_text(json.dumps(res, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
