#!/usr/bin/env python3
"""The forecasts of batch 3 by three systems, before stage B is opened. results/forecasts.json, once.

  python forecast.py --cells batch3_cells.tsv --states B2_STATES_frangieh.npz B3_STATES_frangieh.npz --states-thp B2_STATES_papalexi.npz
                     --read read_cells.npz --pack2 batch2_pack.npz --library LIBRARY.npz --pack-a batch3_pack_A.npz --old OLD --work DIR

S0, the model alone: the part of S2's forecast that a readout of control cells reads from the state.
S1, the system as it forecast batch 2: scripts/s1.py, with the package of commit 903532a in OLD.
S2, version 5: runs/freeze_v5_20261126/scripts/_v5.py with the package of this checkout.

For a readout read before, a stratum's recorded rows are all its units read so far, and the seven new units are
forecast. For a readout never read, they are the teachers of roles.json alone, and every other unit is forecast. With
no signal given there is no teacher: the law holds a placeholder event (a copy of the control cells; `_b3.py`) and the
library's spread with no edge. THP-1 has no library of its own on other cells: where its law has too few rows, the
spread with no edge is that of the THP-1 law on PD-L1 over the units read so far.

The flags of the library's law are proposed as edges here, so that their test on the rows of stage B is fixed before
those rows are opened. A pair that run 43 already tested is not proposed again; the proposal of run 43 that had no row
elsewhere is tested now. A proposed edge enters no forecast.
"""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _b3 as B  # noqa: E402
from _b3 import V  # noqa: E402
from anm_diag import account as AC  # noqa: E402
from anm_diag import events as EV  # noqa: E402

sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
PATHWAY = "interferon pathway"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=Path, required=True)
    ap.add_argument("--states", type=Path, nargs="+", required=True)
    ap.add_argument("--states-thp", type=Path, nargs="+", required=True)
    for k in ("read", "pack2", "library", "pack-a", "old", "work"):
        ap.add_argument("--" + k, type=Path, required=True)
    a = ap.parse_args()
    out_file = B.RUN / "results" / "forecasts.json"
    if out_file.exists():
        raise SystemExit("forecasts.json exists: the forecasts were made once and are not made again")
    cells = B.table(a.cells)
    st = {"frangieh2021": B.States(a.states), "papalexi2021": B.States(a.states_thp)}
    packs = B.Packs([a.pack_a])
    used = B.used_values(a.read, a.pack2)
    S = {g: B.stratum(g, cells, st[B.STRATA[g][0]]) for g in B.STRATA}
    d = np.load(a.library)
    lib_ids, un = d["ids"].astype(str), np.where(d["A"] == 0, "control", d["target"].astype(str))
    mel_all, mel_new = list(B.PANEL["frangieh2021"]), B.NEW["frangieh2021"]

    # --- the scale of every readout in every condition, from control cells
    scale = {r: {} for r in mel_all}
    for g, ids in (("library", lib_ids[un == "control"]), ("mel_co", S["mel_co"]["ids"][S["mel_co"]["kind"] == "control"]), ("mel_ctl", S["mel_ctl"]["ids"][S["mel_ctl"]["kind"] == "control"])):
        clr, log, bg = packs.readouts("frangieh2021", ids)
        for r in mel_all:
            scale[r][g] = V.scale(log[r], bg, clr[r])
    sc_of = lambda g: {r: scale[r][g] for r in mel_all}

    # --- the event graph: declared edges, and the edges admitted in run 43
    sites = sorted((set(un.tolist()) | {u for g in ("mel_co", "mel_ctl") for u in S[g]["unit"].tolist()}) - {"control"})
    declared, tab = V.declared(sites)
    r43 = json.loads((B.RUN.parent / "one_loop_20261125" / "results" / "retro.json").read_text())["edges"]
    on_record = [e for e in r43 if e["status"] != "declared"]       # what run 43 proposed: admitted, refused, or still proposed for want of a row
    admitted = [e for e in on_record if e["status"] == "admitted"]
    tab_used = {g: {**V.restrict(tab, B.USED["frangieh2021"]), **{(e["site"], e["readout"]): e["class"] for e in admitted if e["stratum"] == g}} for g in ("mel_co", "mel_ctl")}

    # --- the library under interferon, every readout
    YL, _log, _bg = packs.readouts("frangieh2021", lib_ids)
    ZL = V.unit_length(d["Z"].astype(np.float64))
    lib = V.library_law(ZL, YL, un, tab, sc_of("library"))
    flags = lib["law"].flags()
    proposed = EV.propose_edges(flags, known=declared + on_record)
    to_test = [e for e in on_record if e["status"] == "proposed"] + proposed
    print(f"library: {len(lib['law'].events)} events, {len(YL)} readouts; spread with no edge {lib['law'].sigma0:.3f} ({lib['law'].spread_from}); classes {dict((c, len(v)) for c, v in lib['classes'].items())}")
    print("  shares:", {f"{e}->{r}": round(v, 2) for (e, r), v in lib["shares"].items()})
    print("  flags:", [(f["event"], f["readout"], round(f["direct"], 2), f["flag"]) for f in flags])

    res = {"cells_file_sha256": sha(a.cells), "pack_a_sha256": sha(a.pack_a), "code": {"S1": "903532a", "S2": "d777070, runs/freeze_v5_20261126/scripts/_v5.py"},
           "scale": scale, "scale_measures": {r: {g: bool(v >= EV.MIN_SCALE) for g, v in x.items()} for r, x in scale.items()},
           "edges": {"declared": declared, "on_record_from_run_43": on_record, "proposed_here": proposed, "to_test_on_the_rows_of_stage_B": to_test},
           "library": {"events": len(lib["law"].events), "readouts": list(YL), "spread_of_a_direct_part_without_an_edge": lib["law"].sigma0, "spread_from": lib["law"].spread_from,
                       "share_of_pairs_with_a_part_of_their_own": lib["law"].own_share, "widening_of_the_control_prior": lib["law"].summary()["widening_of_the_control_prior"],
                       "classes": lib["classes"], "shares": [{"event": e, "readout": r, "share": v} for (e, r), v in lib["shares"].items()], "flags": flags},
           "strata": {}, "rows": []}
    bundle = {"lib_x": ZL, "lib_Y": np.column_stack([YL[r] for r in YL]), "lib_unit": un}
    meta = {"lib_readouts": list(YL), "edges": sorted([list(k) for k in tab]), "laws": []}

    def keep(key, study, readouts, x, Y, unit, xn, un_, placeholder=False):
        bundle.update({key + "_x": x, key + "_Y": np.column_stack([Y[r] for r in readouts]), key + "_unit": unit, key + "_xn": xn, key + "_un": un_})
        meta["laws"].append({"key": key, "study": study, "readouts": readouts, "placeholder": placeholder})

    def rows_of(g, group, law, rows, kinds):
        for r in rows:
            r.update(stratum=g, group=group, kind=kinds[r["event"]], law=f"{g}|{'used' if group == 'read before' else 'new'}")
        res["rows"].extend(rows)
        return {"recorded_events": [e for e in law.events if e != B.NOTHING], "spread_of_a_direct_part_without_an_edge": law.sigma0, "spread_from": law.spread_from,
                "widening_of_the_control_prior": law.summary()["widening_of_the_control_prior"], "pairs": len(rows)}

    # --- melanoma: the readouts read before, for the seven new units
    used_laws = {}
    for g in ("mel_co", "mel_ctl"):
        c = S[g]
        rec = c["kind"] != "new"
        Y = {r: np.array([used[r][i] for i in c["ids"][rec]]) for r in B.USED["frangieh2021"]}
        used_laws[g] = (V.stratum_law(c["Z"][rec], Y, c["unit"][rec], tab_used[g], sc_of(g), lib), Y, rec)
    twice = V.measured_twice(lib, [v[0] for v in used_laws.values()])
    res["pairs_measured_twice"] = {"pairs": len(twice), "spread_of_a_carried_share": EV.carried_spread(twice)}
    print("pairs measured twice:", len(twice), "| spread of a carried share:", res["pairs_measured_twice"]["spread_of_a_carried_share"])
    for g in ("mel_co", "mel_ctl"):
        c, ((law, ed, scl), Y, rec) = S[g], used_laws[g]
        new = c["kind"] == "new"
        kinds = dict(zip(c["unit"].tolist(), c["kind"].tolist()))
        rows = V.forecast(law, ed, scl, c["Z"][new], c["unit"][new], lib, twice)
        res["strata"].setdefault(g, {})["read before"] = rows_of(g, "read before", law, rows, kinds)
        keep(f"{g}|used", "frangieh2021", B.USED["frangieh2021"], c["Z"][rec], Y, c["unit"][rec], c["Z"][new], c["unit"][new])

    # --- melanoma: the readouts never read
    for g in ("mel_co", "mel_ctl"):
        c = S[g]
        rec, test = np.isin(c["kind"], ["control", "teacher"]), np.isin(c["kind"], ["read", "new"])
        assert packs.has("frangieh2021", c["ids"][rec]).all() and not packs.has("frangieh2021", c["ids"][test]).any(), "stage A holds a cell that is forecast, or lacks one the forecasts stand on"
        Y, _l, _b = packs.readouts("frangieh2021", c["ids"][rec], mel_new)
        x, unit, placeholder = c["Z"][rec], c["unit"][rec], not (c["kind"] == "teacher").any()
        if placeholder:
            x, Y, unit = B.with_placeholder(x, Y, unit)
            ed, scl = EV.by_scale(V.restrict(tab, mel_new), sc_of(g))
            law = AC.Law(x, Y, unit, edges=ed, types=V.PATHWAY, scale=scl).fit(sigma0=lib["law"].sigma0)
        else:
            law, ed, scl = V.stratum_law(x, Y, unit, tab, sc_of(g), lib)
        kinds = dict(zip(c["unit"].tolist(), c["kind"].tolist()))
        rows = V.forecast(law, ed, scl, c["Z"][test], c["unit"][test], lib, twice)
        res["strata"][g]["never read"] = rows_of(g, "never read", law, rows, kinds)
        keep(f"{g}|new", "frangieh2021", mel_new, x, Y, unit, c["Z"][test], c["unit"][test], placeholder)

    # --- THP-1: the three readouts never read; the spread with no edge from the THP-1 law on PD-L1
    c = S["thp1"]
    assert not (c["kind"] == "new").any()
    y = np.array([used["PDL1"][i] for i in c["ids"]])
    thp = AC.Law(c["Z"], {"PDL1": y}, c["unit"], types=V.PATHWAY).fit()
    thp_lib = {"law": thp, "classes": {}, "shares": {}}
    rec, test = np.isin(c["kind"], ["control", "teacher"]), c["kind"] == "read"
    assert packs.has("papalexi2021", c["ids"][rec]).all() and not packs.has("papalexi2021", c["ids"][test]).any()
    Y, _l, _b = packs.readouts("papalexi2021", c["ids"][rec], B.NEW["papalexi2021"])
    law, ed, scl = V.stratum_law(c["Z"][rec], Y, c["unit"][rec], {}, {}, thp_lib)
    kinds = dict(zip(c["unit"].tolist(), c["kind"].tolist()))
    rows = V.forecast(law, ed, scl, c["Z"][test], c["unit"][test], thp_lib, [])
    res["strata"]["thp1"] = {"never read": rows_of("thp1", "never read", law, rows, kinds),
                             "the_law_on_PD-L1_that_gives_the_spread": {"events": len(thp.events), "spread_of_a_direct_part_without_an_edge": thp.sigma0, "spread_from": thp.spread_from}}
    keep("thp1|new", "papalexi2021", B.NEW["papalexi2021"], c["Z"][rec], Y, c["unit"][rec], c["Z"][test], c["unit"][test])

    # --- S1, with the package of batch 2
    a.work.mkdir(parents=True, exist_ok=True)
    np.savez(a.work / "s1_bundle.npz", meta=json.dumps(meta), **bundle)
    subprocess.run([sys.executable, str(Path(__file__).resolve().parent / "s1.py"), "--package", str(a.old), "--bundle", str(a.work / "s1_bundle.npz"), "--out", str(a.work / "s1.json")], check=True)
    s1 = json.loads((a.work / "s1.json").read_text())
    by = {(r["law"], r["event"], r["readout"]): r for r in s1["rows"]}
    assert len(by) == len(s1["rows"]) == len(res["rows"])
    rows = []
    for r in res["rows"]:
        o = by[(r["law"], r["event"], r["readout"])]
        assert abs(o["read"] - r["read"]) < 1e-6, "the two packages read the control readout differently"
        untyped = r["type"] is None or r["type"].startswith("(")
        rows.append({"stratum": r["stratum"], "event": r["event"], "kind": r["kind"], "readout": r["readout"], "group": r["group"], "cells": r["cells"], "edge": r["edge"], "type": None if untyped else r["type"],
                     "taught_by": r["taught_by"], "scale": r["scale"] if r["edge"] and (r["event"], r["readout"]) in tab and not r["edge"].endswith(EV.AT_THE_BACKGROUND) else None,
                     "size_of_the_edge_from": r["size_of_the_edge_from"], "spread_of_a_carried_share": r["spread_of_a_carried_share"], "shift_times_noise": r["shift_times_noise"],
                     "S0": {"forecast": r["read"]},
                     "S1": {"forecast": o["forecast"], "sd": o["sd"], "interval": o["interval"], "from_the_recorded_events": o["from_the_recorded_events"], "by_its_edge": o["by_its_edge"], "no_forecast": o["no_forecast"]},
                     "S2": {"forecast": r["forecast"], "sd": r["sd"], "interval": r["interval"], "from_the_recorded_events": r["from_the_recorded_events"], "by_its_edge": r["by_its_edge"],
                            "no_forecast": r["no_forecast"], "parts_of_the_spread": r["parts_of_the_spread"]}})
    res["rows"], res["pairs"] = rows, len(rows)
    res["S1"] = {k: s1[k] for k in ("package", "library", "class_own_pooled", "laws")}
    key = lambda r: f"{r['stratum']}|{r['event']}|{r['readout']}"
    res["sets"] = {"the_systems_can_differ": [key(r) for r in rows if r["edge"] is not None or r["type"] == PATHWAY],
                   "H3": [key(r) for r in rows if r["edge"] == "own" and (r["size_of_the_edge_from"] or "").startswith("the pair")],
                   "H4": [key(r) for r in rows if r["type"] == PATHWAY and r["taught_by"] > 0 and r["edge"] is None and r["group"] == "never read"],
                   "H5": [key(r) for r in rows if r["type"] is None and r["edge"] is None and r["group"] == "never read"]}
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")

    print("pairs forecast:", res["pairs"], {g: sum(r["stratum"] == g for r in rows) for g in B.STRATA})
    for g in B.STRATA:
        for grp, v in res["strata"][g].items():
            if "pairs" in v:
                print(f"  {g:7s} {grp:11s} recorded {v['recorded_events']} | spread {v['spread_of_a_direct_part_without_an_edge']:.3f} ({v['spread_from']}) | S1 {s1['laws'][g + ('|used' if grp == 'read before' else '|new')]['spread_of_a_direct_part_without_an_edge']:.3f}")
    print("sets:", {k: len(v) for k, v in res["sets"].items()})
    print("edges to test on stage B:", [e["edge_id"] for e in to_test])
    for r in rows:
        if r["edge"] is not None:
            print(f"  {r['stratum']:7s} {r['event']:6s}->{r['readout']:5s} {r['edge']:24s} scale {r['scale'] if r['scale'] is None else round(r['scale'], 2)} | S0 {r['S0']['forecast']:+.2f} | S1 {r['S1']['forecast']:+.2f} "
                  f"| S2 {r['S2']['forecast']:+.2f} [{r['S2']['interval'][0]:+.2f}, {r['S2']['interval'][1]:+.2f}] {r['size_of_the_edge_from']}")


if __name__ == "__main__":
    main()
