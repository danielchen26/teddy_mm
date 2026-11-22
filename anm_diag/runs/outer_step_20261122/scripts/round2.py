#!/usr/bin/env python3
"""One round under package 2: MHC-I in melanoma under IFN-gamma, with the events admitted on the selection third.

  python round2.py --observe2 OBSERVE2 --work2 WORK2 --selection ../results/selection.json     (writes results/round2.json)

OBSERVE2 holds the grid mel_ifn_mhc_MHC1 (and the round's work folder is ANM_LOOP_WORK, for the co-culture cells whose
repair is the arm "other"). Environment as runs/teg_loop_20261120.

Cells: outside the selection third. Two events, so every reading is per event and resampled by cell: no statement about
a family. The change is the one of the round of 20 November: the output is replaced, for cells of this type, by a
linear readout of Z fitted on the type's control cells and updated once by the rows of its admitted events. Two
halves: every cell is read by a readout fitted without it. The type of a cell comes from the record.

Arms (protocol version 3, and the three arms of a coverage claim):
  rows          the rows of the admitted events (the own arm)
  controls      the same readout without rows: what a readout of this type's control cells reads
  rows_covered  rows of the five pathway genes of this type, which move MHC-I little: cells of a kind already covered
  other         the repair of co-culture applied to these cells: another map
  rows_other    the rows of co-culture's events on this type's readout
  cells         the readout fitted on every cell of the fitted half (the ceiling of this map)
Admission per event: A1 and A2 on the selection third (selection.json); A3 on the grid (the head loses what a readout
of the other cells of this type reads: the gap from the head's output to the `units` scope at Z above the step margin).
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bridge as B  # noqa: E402
from bridge import LP, R0, diff, val  # noqa: E402
from _common import forward_np, halves, interval, own_output_u, unit, weights  # noqa: E402
from build import OUTPUT  # noqa: E402

STEP, MARGIN, SEED = B.STEP, B.MARGIN, B.SEED
PATHWAY = ["IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1"]
ARMS = ["controls", "rows", "cells", "rows_covered", "rows_other", "other"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observe2", type=Path, required=True)
    ap.add_argument("--work2", type=Path, required=True)
    ap.add_argument("--selection", type=Path, required=True)
    a = ap.parse_args()
    out_file = B.RUN / "results" / "round2.json"
    if out_file.exists():
        raise SystemExit("round2.json exists: this reading was made once and is not made again")
    sel = json.loads(a.selection.read_text())["MHC1"]
    sites = sel["admitted"]
    P0, names, _j = weights()
    j = names.index(OUTPUT["MHC1"])
    head = lambda x: forward_np(P0, x)[1][:, j]
    d = np.load(a.work2 / "mel_ifn_mhc.npz")
    c = dict(ids=d["ids"].astype(str), x=unit(d["Z"]).astype(np.float64), A=d["A"].astype(int), tg=d["unit"].astype(str), y=d["y_MHC1"].astype(np.float64))
    c["o"] = head(c["x"])
    first = halves(np.where(c["A"] == 1, c["tg"], "control"), np.random.default_rng(SEED))
    co = B.load_cells(Path(os.environ["ANM_LOOP_WORK"]), "mel_co", "MHC1", head)
    whole_co = R0.Readout(co["x"], co["y"], co["o"], np.where(co["A"] == 0)[0], [(s, np.where((co["A"] == 1) & (co["tg"] == s))[0]) for s in PATHWAY])
    obs = a.observe2 / "mel_ifn_mhc_MHC1"
    by = {r["event"]["unit"]: r for r in json.loads((obs / "scope.json").read_text())["results"]}
    G, draws = {u: R0.grid_row(r) for u, r in by.items()}, np.load(obs / "scope_draws.npz")
    n = len(c["A"])
    out, inside, fits = {k: np.full(n, np.nan) for k in ARMS}, np.full(n, np.nan), []
    for fold in (0, 1):
        fit = first if fold == 0 else ~first
        rd = ~fit
        C = np.where(fit & (c["A"] == 0))[0]
        ev = lambda genes: [(s, np.where(fit & (c["A"] == 1) & (c["tg"] == s))[0]) for s in genes]
        R = R0.Readout(c["x"], c["y"], c["o"], C, ev(sites), (co["x"], co["y"], co["A"], co["tg"], PATHWAY))
        Rc = R0.Readout(c["x"], c["y"], c["o"], C, ev(PATHWAY))
        fits.append({"rows": R.log, "rows_covered": Rc.log["rows"]})
        assert not set(c["ids"][rd]) & set(c["ids"][fit])
        for k in ("controls", "rows", "cells", "rows_other"):
            out[k][rd] = R.out(k, c["x"][rd])
        out["rows_covered"][rd] = Rc.out("rows", c["x"][rd])
        out["other"][rd] = whole_co.out("rows", c["x"][rd])
        inside[fit] = R.out("rows", c["x"][fit])
    res = {"readout": "MHC1", "group": "mel_ifn", "type": {"cell_type": "melanoma", "stratum": "IFN-gamma"}, "events_admitted_on_the_selection_third": sites,
           "cells": {k: int(v) for k, v in zip(*np.unique(c["tg"], return_counts=True))}, "margin": MARGIN, "step_margin": STEP, "events": {}, "fits": fits,
           "note": "per event, resampled by cell; two events give no statement about a family"}
    ledger = []
    for s in sites:
        g = G[s]
        m = np.where(((c["A"] == 1) & (c["tg"] == s)) | (c["A"] == 0))[0]
        u0, s0 = own_output_u(c["y"][m], c["o"][m], c["A"][m], None)
        assert np.allclose(u0[0], draws[f"{s}|shift|system|frozen|system"][0], atol=2e-4), "the resamples do not pair"
        ref = draws[f"{s}|shift|Z|units|{g['scopes']['units']['class']}"].astype(np.float64)
        own = draws[f"{s}|shift|Z|own|{g['scopes']['own']['class']}"].astype(np.float64)
        gap = g["gaps"]["from_frozen"]["units"]
        e = {"response": g["response"], "n_cells": g["n_treated"], "n_control": g["n_control"], "selection_third": next(x for x in sel["top"] if x["event"] == s),
             "grid": {k: B.cell(v) for k, v in g["scopes"].items()},
             "admission": {"A1_A2": "on the selection third", "A3": bool(gap["ci"][0] > STEP), "gap_head_to_units": gap,
                           "attainable": "absolute" if g["scopes"]["units"]["verdict"] == "support" else "relative"},
             "reference": {"units": val(ref), "own": val(own), "own_beyond_units": diff(own, ref)},
             "as_trained": {**val(u0), "slope": s0, "beyond_reference": diff(u0, ref), "beyond_own": diff(u0, own)}, "arms": {}}
        U = {}
        for k in ARMS:
            p = out[k]
            u, slope = own_output_u(c["y"][m], p[m], c["A"][m], None)
            U[k] = u
            T_, C_ = c["A"][m] == 1, c["A"][m] == 0
            gain = np.abs(u0) - np.abs(u)
            e["arms"][k] = {**val(u), "slope": slope, "valid": bool(slope > 0 and np.sign(p[m][T_].mean() - p[m][C_].mean()) == np.sign(c["y"][m][T_].mean() - c["y"][m][C_].mean())),
                            "gain": float(gain[0]), "gain_ci": interval(gain[1:]), "beyond_reference": diff(u, ref), "beyond_own": diff(u, own)}
        dd = np.abs(U["controls"]) - np.abs(U["rows"])
        iv = interval(dd[1:])
        e["rows_minus_controls"] = {"gain_difference": float(dd[0]), "ci": iv, "reading": "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")}
        ui, _s = own_output_u(c["y"][m], inside[m], c["A"][m], None)
        g3, g0 = diff(ui, ref), e["as_trained"]["beyond_reference"]
        e["gate3"] = {"on_the_cells_it_was_fitted_on": val(ui), "against_the_reference": g3, "passes": bool(g3["ci"][1] < STEP),
                      "control_the_output_as_trained": {**g0, "passes": bool(g0["ci"][1] < STEP)},
                      "guard": {"fall": 0.0, "kept": True, "why": "the type comes from the record: a bone-marrow cell keeps the head's own output"}}
        e["closure"] = B.grade({"arms": e["arms"], "as_trained": e["as_trained"], "gate3": e["gate3"]}, e["admission"]["attainable"]) if e["admission"]["A3"] else \
            {"grade": "not a target for a repair of the head", "why": "A3: the head does not lose what a readout of the other cells of this type reads; what is lost is lost before the head"}
        ceil = g["ceiling"]["coverage"]
        at_least = round(ceil["removes"] - STEP, 3)
        H = LP.hypothesis_event(f"coverage:mel_ifn:MHC1:{s}:round2", {"event": f"{s}, melanoma under IFN-gamma", "readout": "MHC1", "target": "shift", "quantity": "fall of |u| of the system's own output"},
                                "coverage", {"rows": ceil["ci"]}, {"own": "rows", "at_least": at_least, "controls": {k: at_least for k in ("controls", "rows_covered", "other")}},
                                ["scope:mel_ifn_mhc_MHC1"], rivals=["function", "upstream", "unreadable"],
                                note="protocol version 3; the forecast is the grid's ceiling for this event; the bound is that ceiling minus the step margin")
        arms = {k: {"gain": v["gain"], "ci": v["gain_ci"], "runs": [f"round2:{s}"], "admitted": bool(e["gate3"]["passes"])} for k, v in e["arms"].items() if k in ("rows", "controls", "rows_covered", "other")}
        e["hypothesis"] = H
        if e["admission"]["A3"]:
            e["test"] = LP.test(H, arms)
            ledger += [H, e["test"]]
        else:                                                   # OUTER_LOOP.md, Section 10: where the reference shows nothing to remove, the outcome is "not testable", not "refuted"
            e["test"] = {"type": "test", "hypothesis": H["id"], "outcome": "not testable", "arms": {k: {"gain": v["gain"], "ci": v["ci"]} for k, v in arms.items()},
                         "why": f"nothing to remove: the grid's ceiling for a repair of the head is {ceil['removes']:.2f}; the reference stage reads no more than the head (A3)"}
            ledger += [H]
        res["events"][s] = e
    res["pool_of_the_two_events"] = {"grid": {k: B.cell(v) for k, v in G["pool"]["scopes"].items()}, "gap_head_to_units": G["pool"]["gaps"]["from_frozen"]["units"],
                                     "A3": bool(G["pool"]["gaps"]["from_frozen"]["units"]["ci"][0] > STEP), "note": "resampled over two genes: not an interval to rest anything on"}
    res["ledger_check"] = LP.ledger_check(ledger)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda v: f"{v:6.3f}"
    for s, e in res["events"].items():
        print(f"\n{s}: response {e['response']['value']:.2f} SD (snr {e['response']['snr']:.1f}), n {e['n_cells']} | grid " + " ".join(f"{k} {v['u']:.2f} {v['verdict']}" for k, v in e["grid"].items())
              + f" | A3 {e['admission']['A3']} (gap {e['admission']['gap_head_to_units']['value']:.2f} {[round(x, 2) for x in e['admission']['gap_head_to_units']['ci']]})")
        print(f"  as trained {f(e['as_trained']['u'])} | reference units {f(e['reference']['units']['u'])} {e['reference']['units']['verdict']} own {f(e['reference']['own']['u'])} {e['reference']['own']['verdict']}")
        for k, v in e["arms"].items():
            print(f"  {k:12s} {f(v['u'])} [{f(v['ci'][0])}, {f(v['ci'][1])}] {v['verdict']:10s} gain {f(v['gain'])} [{f(v['gain_ci'][0])}, {f(v['gain_ci'][1])}] beyond units {f(v['beyond_reference']['added'])} "
                  f"{v['beyond_reference']['reading']:13s} beyond own {f(v['beyond_own']['added'])} {v['beyond_own']['reading']}")
        print(f"  gate 3 {f(e['gate3']['against_the_reference']['added'])} passes {e['gate3']['passes']} | rows - controls {f(e['rows_minus_controls']['gain_difference'])} {e['rows_minus_controls']['reading']}"
              f" | test {e['test']['outcome']} ({e['test']['why']}) | closure {e['closure']['grade']}")
    print("\nledger:", [(x["code"], x["severity"]) for x in res["ledger_check"]])


if __name__ == "__main__":
    main()
