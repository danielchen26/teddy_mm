#!/usr/bin/env python3
"""The bridge reading: the cells and the model of the round of 20 November, read again under package 2 (DESIGN.md).
No new data and no change to the model.

  python bridge.py --observe OBSERVE1 --observe2 OBSERVE2 --work2 WORK2            (writes results/bridge.json)

OBSERVE1 holds the grids of the round (runs/teg_loop_20261120), OBSERVE2 the grids of this step (thp1a_PDL1,
mel_co_MHC1_z2), WORK2 the folder of build2.py. Environment as the round.

What differs from the round's reading (protocol version 3):
  3e  a cell is read by the readout of its recorded type; there is no reader of types in the path
  3a  gate 3 against the reference stage: on the cells the repair was fitted on, |u of the output| - |u of the reference|
      has its interval below the step margin; the output as trained is put through the same gate
  3b  the test has two controls: the repair of another group, and the readout fitted on control cells without rows;
      the difference between the rows arm and the controls-only arm is read on paired resamples
and the THP-1 pool is read three ways: as one pool of five (the schema of the round), as type A (four genes, their own
grid), and as type B (IRF1 with its own row; its reference is a reader of Z with its own cells).
The halves, the seeds and the fits are those of the round.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE.parent
RUNS = RUN.parent
sys.path.insert(0, str(RUNS.parent))
sys.path.insert(0, str(RUNS / "teg_training_20261119" / "scripts"))
sys.path.insert(0, str(RUNS / "teg_loop_20261120" / "scripts"))
import round as R0  # noqa: E402
from anm_diag import loop as LP  # noqa: E402
from _common import forward_np, halves, interval, own_output_u, unit, weights  # noqa: E402
from build import GROUPS, OUTPUT  # noqa: E402

STEP, MARGIN, SEED = R0.STEP, R0.MARGIN, R0.SEED
TYPE_A = ["IFNGR1", "IFNGR2", "JAK2", "STAT1"]
step_reading = lambda iv: "loss" if iv[0] > STEP else ("no added loss" if iv[1] < STEP else "unresolved")
cell = lambda v: {k: v[k] for k in ("u", "ci", "verdict", "class")}


def val(u):
    iv = interval(u[1:])
    return {"u": float(u[0]), "ci": iv, "verdict": R0.verdict(iv)}


def diff(a, b):
    """|a| - |b| on paired resamples."""
    d = np.abs(a) - np.abs(b)
    iv = interval(d[1:])
    return {"added": float(d[0]), "ci": iv, "reading": step_reading(iv)}


def load_cells(work, g, readout, head):
    d = np.load(work / f"{g}.npz")
    c = dict(ids=d["ids"].astype(str), x=unit(d["Z"]).astype(np.float64), A=d["A"].astype(int), tg=d["unit"].astype(str), y=d[f"y_{readout}"].astype(np.float64))
    c["o"] = head(c["x"])
    c["first"] = halves(np.where(c["A"] == 1, c["tg"], "control"), np.random.default_rng(SEED))     # the halves of the round
    return c


def read(c, sites, G, draws, ref_key, other=None, whole_other=None, by_event=True, name=""):
    """One target under package 2. c: the cells; sites: the events whose rows repair it; G, draws: its grid; ref_key: the draws
    of the reference stage; other: (cells of the other group, its admitted sites) for the arm rows_other; whole_other: the
    other group's repair fitted on all its cells, for the arm other."""
    n = len(c["A"])
    arms = ["controls", "rows", "cells"] + (["rows_other"] if other else []) + (["other"] if whole_other is not None else [])
    out, inside, fits = {a: np.full(n, np.nan) for a in arms}, np.full(n, np.nan), []
    for fold in (0, 1):
        fit = c["first"] if fold == 0 else ~c["first"]
        rd = ~fit
        C = np.where(fit & (c["A"] == 0))[0]
        ev = [(s, np.where(fit & (c["A"] == 1) & (c["tg"] == s))[0]) for s in sites]
        oth = (other[0]["x"], other[0]["y"], other[0]["A"], other[0]["tg"], other[1]) if other else None
        R = R0.Readout(c["x"], c["y"], c["o"], C, ev, oth)
        fits.append(R.log)
        assert not set(c["ids"][rd]) & set(c["ids"][fit])
        for a in ("controls", "rows", "cells", "rows_other"):
            if a in out and a in R.th:
                out[a][rd] = R.out(a, c["x"][rd])                    # 3e: the readout of the cell's recorded type
        if whole_other is not None:
            out["other"][rd] = whole_other.out("rows", c["x"][rd])
        inside[fit] = R.out("rows", c["x"][fit])
    cl = c["tg"] if len(sites) > 1 else None
    u0, s0 = own_output_u(c["y"], c["o"], c["A"], cl)
    assert np.allclose(u0[0], draws[ref_key.split("|")[0] + "|shift|system|frozen|system"][0], atol=2e-4), f"{name}: the resamples do not pair"
    ref = draws[ref_key].astype(np.float64)
    e = {"events": sites, "resampled_by": "gene" if cl is not None else "cell", "reference": {"draws": ref_key, **val(ref)},
         "as_trained": {**val(u0), "slope": s0, "beyond_reference": diff(u0, ref)}, "arms": {}, "fits": fits}
    U = {}
    for a in arms:
        p = out[a]
        assert np.isfinite(p).all(), (name, a)
        u, slope = own_output_u(c["y"], p, c["A"], cl)
        U[a] = u
        T_, C_ = c["A"] == 1, c["A"] == 0
        gain = np.abs(u0) - np.abs(u)
        e["arms"][a] = {**val(u), "slope": slope, "valid": bool(slope > 0 and np.sign(p[T_].mean() - p[C_].mean()) == np.sign(c["y"][T_].mean() - c["y"][C_].mean())),
                        "gain": float(gain[0]), "gain_ci": interval(gain[1:]), "beyond_reference": diff(u, ref), "events": {}}
        if by_event and len(sites) > 1:
            for s in sites:
                m = np.where(((c["A"] == 1) & (c["tg"] == s)) | (c["A"] == 0))[0]
                ue, _s = own_output_u(c["y"][m], p[m], c["A"][m], None)
                re = draws[f"{s}|shift|Z|units|{G[s]['scopes']['units']['class']}"].astype(np.float64)
                e["arms"][a]["events"][s] = {**val(ue), "beyond_reference": diff(ue, re)}
    d = np.abs(U["controls"]) - np.abs(U["rows"])                   # what the rows add over a readout of control cells
    iv = interval(d[1:])
    e["rows_minus_controls"] = {"gain_difference": float(d[0]), "ci": iv, "reading": "a step" if iv[0] > STEP else ("no step" if iv[1] < STEP else "unresolved")}
    ui, _s = own_output_u(c["y"], inside, c["A"], cl)
    g3, g0 = diff(ui, ref), e["as_trained"]["beyond_reference"]
    e["gate3"] = {"on_the_cells_it_was_fitted_on": val(ui), "absolute_as_in_version_2": bool(abs(ui[0]) < MARGIN), "against_the_reference": g3, "passes": bool(g3["ci"][1] < STEP),
                  "control_the_output_as_trained": {**g0, "passes": bool(g0["ci"][1] < STEP)}, "no_audit_cell_in_a_fit": True,
                  "guard": {"fall": 0.0, "kept": True, "why": "the type comes from the record: a bone-marrow cell keeps the head's own output"}}
    return e, U, out, ref


def grade(e, attainable):
    r = e["arms"]["rows"]
    base_loss = e["as_trained"]["beyond_reference"]["ci"][0] > STEP
    if not e["gate3"]["passes"]:
        g = "not tested (gate 3)"
    elif r["verdict"] == "support" and r["valid"] and e["as_trained"]["verdict"] == "reject":
        g = "absolute"
    elif r["beyond_reference"]["reading"] == "no added loss" and r["valid"] and base_loss:
        g = "relative"
    elif r["beyond_reference"]["reading"] == "unresolved":
        g = "loss no longer resolved"
    else:
        g = "not closed"
    order = ["not closed", "loss no longer resolved", "relative", "absolute"]
    return {"grade": g, "attainable": attainable, "reached_what_was_attainable": bool(g in order and order.index(g) >= order.index(attainable))}


def test(e, G, draws, hid, readout, basis):
    ceil = G["pool"]["ceiling"]["coverage"]
    at_least = round(ceil["removes"] - STEP, 3)
    fc = {"rows": ceil["ci"]}                                       # the grid has no forecast for a readout of control cells alone
    if "other" in e["arms"] and "other" in G["pool"]["scopes"]:
        fo = draws[f"pool|shift|Z|other|{G['pool']['scopes']['other']['class']}"].astype(np.float64)
        fc["other"] = interval((np.abs(draws["pool|shift|system|frozen|system"]) - np.abs(fo))[1:])
    ctl = {k: at_least for k in ("other", "controls") if k in e["arms"]}
    H = LP.hypothesis_event(hid, {"event": "the admitted knockouts, pooled", "readout": readout, "target": "shift", "quantity": "fall of |u| of the system's own output"},
                            "coverage", fc, {"own": "rows", "at_least": at_least, "controls": ctl}, basis, rivals=["function", "upstream"], retrospective=True,
                            note="protocol version 3, amendment 3b: two controls. Retrospective: the cells were read in the round")
    arms = {k: {"gain": v["gain"], "ci": v["gain_ci"], "runs": [hid + ":arms"], "admitted": bool(e["gate3"]["passes"])} for k, v in e["arms"].items() if k in ("rows", "other", "controls")}
    return H, LP.test(H, arms)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--observe", type=Path, required=True)
    ap.add_argument("--observe2", type=Path, required=True)
    ap.add_argument("--work2", type=Path, required=True)
    a = ap.parse_args()
    work = Path(os.environ["ANM_LOOP_WORK"])
    P0, names, _j = weights()
    res = {"note": "the round of 20 November read again under package 2; exploratory; no new data, no change to the model", "targets": {}}
    for readout, groups in (("PDL1", ["mel_ifn", "mel_co", "thp1"]), ("MHC1", ["mel_co"])):
        j = names.index(OUTPUT[readout])
        head = lambda x, j=j: forward_np(P0, x)[1][:, j]
        cells = {g: load_cells(work, g, readout, head) for g in GROUPS if readout in GROUPS[g]["readouts"]}
        grid, draws, genes = {}, {}, {}
        for g in cells:
            obs = a.observe2 / "mel_co_MHC1_z2" if (g, readout) == ("mel_co", "MHC1") else a.observe / f"{g}_{readout}"
            if not (obs / "scope.json").exists():
                continue
            d = json.loads((obs / "scope.json").read_text())
            by = {r["event"]["unit"]: r for r in d["results"]}
            grid[g], draws[g] = {u: R0.grid_row(r) for u, r in by.items()}, np.load(obs / "scope_draws.npz")
        adm = {g: R0.admit(g, readout, grid[g]) for g in grid}
        genes = {g: [x["site"] for x in adm[g]["events"] if x["admission"]["status"] == "admitted"] for g in grid}
        whole = {g: R0.Readout(c["x"], c["y"], c["o"], np.where(c["A"] == 0)[0], [(s, np.where((c["A"] == 1) & (c["tg"] == s))[0]) for s in genes[g]])
                 for g, c in cells.items() if genes.get(g)}
        for g in groups:
            c, G, og = cells[g], grid[g], GROUPS[g]["other"]
            other = (cells[og], genes[og]) if genes.get(og) else None
            key = f"pool|shift|Z|units|{G['pool']['scopes']['units']['class']}"
            e, U, out, ref = read(c, genes[g], G, draws[g], key, other, whole.get(og), name=f"{g} {readout}")
            e["schema"] = "one pool of the group's admitted events (as the round)"
            e["grid"] = {s: cell(v) for s, v in G["pool"]["scopes"].items()}
            e["grid_by_event"] = {s: {k: cell(G[s]["scopes"][k]) for k in ("units", "own")} for s in genes[g]}
            e["grid_by_class"] = {k: {cl: {q: x[q] for q in ("u", "ci")} for cl, x in v.get("by_class", {}).items()} for k, v in
                                  json.loads(((a.observe2 / "mel_co_MHC1_z2") if (g, readout) == ("mel_co", "MHC1") else (a.observe / f"{g}_{readout}")).joinpath("scope.json").read_text())
                                  ["results"][-1]["targets"]["shift"]["stages"]["Z"].items()}
            e["closure"] = grade(e, adm[g]["attainable"])
            e["hypothesis"], e["test"] = test(e, G, draws[g], f"coverage:{g}:{readout}:bridge", readout, [f"scope:{g}_{readout}"])
            res["targets"][f"{readout} {g}"] = e
            if (g, readout) == ("thp1", "PDL1"):
                # type A: the four genes, their own grid (the reference may not contain IRF1); the halves of the round
                d = json.loads((a.observe2 / "thp1a_PDL1" / "scope.json").read_text())
                GA = {r["event"]["unit"]: R0.grid_row(r) for r in d["results"]}
                dA = np.load(a.observe2 / "thp1a_PDL1" / "scope_draws.npz")
                m = np.isin(c["tg"], TYPE_A + ["control"])
                ka = np.load(a.work2 / "thp1a.npz")
                assert np.array_equal(ka["ids"].astype(str), c["ids"][m])
                cA = {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in c.items()}
                keyA = f"pool|shift|Z|units|{GA['pool']['scopes']['units']['class']}"
                eA, UA, _u, _r = read(cA, TYPE_A, GA, dA, keyA, other, whole.get(og), name="thp1 type A")
                eA["schema"] = "provisional split, type A: IFNGR1, IFNGR2, JAK2, STAT1"
                eA["grid"] = {s: cell(v) for s, v in GA["pool"]["scopes"].items()}
                eA["grid_by_event"] = {s: {k: cell(GA[s]["scopes"][k]) for k in ("units", "own")} for s in TYPE_A}
                eA["closure"] = grade(eA, "absolute" if GA["pool"]["scopes"]["units"]["verdict"] == "support" else "relative")
                eA["hypothesis"], eA["test"] = test(eA, GA, dA, "coverage:thp1a:PDL1:bridge", readout, ["scope:thp1_PDL1"])
                res["targets"]["PDL1 thp1 type A"] = eA
                # type B: IRF1 with its own row; the reference is a reader of Z with its own cells
                mb = np.isin(c["tg"], ["IRF1", "control"])
                cB = {k: (v[mb] if isinstance(v, np.ndarray) else v) for k, v in c.items()}
                keyB = f"IRF1|shift|Z|own|{G['IRF1']['scopes']['own']['class']}"
                eB, UB, _o, refB = read(cB, ["IRF1"], G, draws[g], keyB, None, None, by_event=False, name="thp1 type B")
                eB["schema"] = "provisional split, type B: IRF1"
                eB["grid"] = {k: cell(G["IRF1"]["scopes"][k]) for k in ("frozen", "train", "other", "units", "own") if k in G["IRF1"]["scopes"]}
                eB["closure"] = grade(eB, "relative")
                # the same event under the schema of the round: the pooled rows, read against the same reference
                mi = np.where(mb)[0]
                up, _s = own_output_u(c["y"][mi], out["rows"][mi], c["A"][mi], None)
                eB["with_the_pooled_rows"] = {**val(up), "beyond_reference": diff(up, refB), "against_its_own_row": diff(up, UB["rows"])}
                rna = a.observe / "thp1_PDL1_rna"
                if (rna / "scope.json").exists():
                    rr = {r["event"]["unit"]: r for r in json.loads((rna / "scope.json").read_text())["results"]}["IRF1"]["targets"]["shift"]["stages"]["RNA512"]
                    dr = np.load(rna / "scope_draws.npz")
                    zin = draws[g][keyB].astype(np.float64), dr[f"IRF1|shift|RNA512|own|{rr['own']['class']}"].astype(np.float64)
                    eB["input_view"] = {"reader": "512 components of RNA, ridge (the round)", "own": {k: rr["own"][k] for k in ("u", "ci", "verdict")},
                                        "units": {k: rr["units"][k] for k in ("u", "ci", "verdict")}, "from_the_input_view_to_Z_own_cells": diff(*zin)}
                res["targets"]["PDL1 thp1 type B"] = eB
    (RUN / "results").mkdir(exist_ok=True)
    (RUN / "results" / "bridge.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    f = lambda v: f"{v:6.3f}"
    for k, e in res["targets"].items():
        print(f"\n{k} | {e['schema']} | reference {f(e['reference']['u'])} {e['reference']['verdict']} | as trained {f(e['as_trained']['u'])}")
        for arm, v in e["arms"].items():
            print(f"  {arm:10s} {f(v['u'])} [{f(v['ci'][0])}, {f(v['ci'][1])}] {v['verdict']:10s} gain {f(v['gain'])} [{f(v['gain_ci'][0])}, {f(v['gain_ci'][1])}] beyond ref {f(v['beyond_reference']['added'])} "
                  f"[{f(v['beyond_reference']['ci'][0])}, {f(v['beyond_reference']['ci'][1])}] {v['beyond_reference']['reading']} | " + " ".join(f"{s}:{q['u']:.2f}" for s, q in v["events"].items()))
        g3 = e["gate3"]
        print(f"  gate 3: {f(g3['against_the_reference']['added'])} [{f(g3['against_the_reference']['ci'][0])}, {f(g3['against_the_reference']['ci'][1])}] passes {g3['passes']} | as trained passes {g3['control_the_output_as_trained']['passes']}"
              f" | rows - controls {f(e['rows_minus_controls']['gain_difference'])} {[round(x, 3) for x in e['rows_minus_controls']['ci']]} {e['rows_minus_controls']['reading']}"
              f" | closure {e['closure']['grade']}" + (f" | test {e['test']['outcome']} ({e['test']['why']})" if "test" in e else ""))
        if "input_view" in e:
            print("  input view:", json.dumps(e["input_view"], default=float)[:400])
        if "with_the_pooled_rows" in e:
            print("  with the pooled rows:", e["with_the_pooled_rows"])


if __name__ == "__main__":
    main()
