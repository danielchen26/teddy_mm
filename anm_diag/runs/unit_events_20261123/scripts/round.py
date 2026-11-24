#!/usr/bin/env python3
"""The round on eight donors: every protein with an admitted donor, by the general code (anm_diag.cycle, anm_diag.outer).

  python round.py --work WORK --observe OBSERVE [--part K --of N]     (writes results/round.json; refuses to run when it exists)

WORK: the folder build.py wrote. OBSERVE: the grid of `anm_diag scope` for WORK/donors.yaml (one run, `each_readout`).
With --part the proteins are split over N processes; each writes results/round_part_K.json and `--merge` joins them.
A donor is read against the other donors (`cycle.one_units`): its cells by cell, the other donors by unit. DESIGN.md, Section 6.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN, RUNS = HERE.parent, HERE.parents[1]
sys.path.insert(0, str(RUNS.parent))
from anm_diag import cycle, outer  # noqa: E402

TYPE = {"study": "hao2021", "cells": "PBMC", "event": "a donor against the other donors"}
WORDS = {"self event": "the event acts on the readout itself", "upstream": "a change of the encoder"}
unit_length = lambda z: z / (np.linalg.norm(z, axis=1, keepdims=True) + 1e-6)


def one_protein(name, x, d, sel, G):
    j, events = sel["index"], sel["admitted"]
    by, draws = cycle.load(G, readout=j)
    T = cycle.one_units(name, "hao2021", TYPE, x, d["Y"][:, j].astype(np.float64), d["frozen"][:, j].astype(np.float64), d["unit"].astype(str), events, by, draws, "Z",
                        inputs=["RNA512"], control="PERM", evidence=f"round.json:{name}", confounded=[("donor", "batch")], pool={"more_units_of_the_type": None})
    t = T["target"]
    ev = {}
    for s, e in T["reading"].items():
        p = T["placed"]["events"][s]
        e["cells"] = {"event": int((d["unit"].astype(str) == s).sum()), "other_donors": int((d["unit"].astype(str) != s).sum())}
        e["selection_third"] = next(q for q in sel["table"] if q["event"] == s)
        e["grid"] = {st: {sc: {k: v[k] for k in ("u", "ci", "verdict", "class")} for sc, v in by[s][st].items()} for st in by[s]}
        e["placed"] = {"input": p["input"], "model_stage": p["model_stage"], "from_the_input_to_the_model_stage": p["from_the_input_to_the_model_stage"], "control": p["control"]}
        e["status"] = outer.target_status(cycle.target_by_event(name, "hao2021", TYPE, {s: e}, T["placed"]))
        ev[s] = e
    return {"readout": name, "index": j, "events": ev, "fits": T["fits"], "status": outer.target_status(t), "target": t,
            "positive_control": T["placed"]["positive_control"], "control_view": T["placed"]["control_view"]}


def summary(res):
    P = [(n, s, e) for n, t in res["targets"].items() for s, e in t["events"].items()]
    share = lambda xs: (float(np.mean(xs)) if len(xs) else None)
    A3 = [x for x in P if x[2]["A3"]["holds"]]
    rej = [x for x in P if x[2]["reference"]["verdict"] == "reject"]
    it = lambda n, s, r: any(i["readout"] == n and i["event"] == s and i["reading"] == r for i in res["open_items"])
    return {"pairs": len(P), "proteins": len(res["targets"]),
            "as_trained": {"rejected": share([e["as_trained"]["verdict"] == "reject" for _n, _s, e in P]), "median_abs_u": float(np.median([abs(e["as_trained"]["u"]) for _n, _s, e in P]))},
            "reference": {v: share([e["reference"]["verdict"] == v for _n, _s, e in P]) for v in ("support", "unresolved", "reject")},
            "A3": {"holds": share([e["A3"]["holds"] for _n, _s, e in P]), "pairs": len(A3)},
            "among_A3": {"closed": share([e["closure"]["grade"] in ("relative", "absolute") for _n, _s, e in A3]), "absolute": share([e["closure"]["grade"] == "absolute" for _n, _s, e in A3]),
                         "grades": {g: int(sum(e["closure"]["grade"] == g for _n, _s, e in A3)) for g in sorted({e["closure"]["grade"] for _n, _s, e in A3})},
                         "median_fall_to_rows": float(np.median([e["arms"]["rows"]["gain"] for _n, _s, e in A3])) if A3 else None,
                         "own_row_needed": share([e["rows_minus_unseen"]["reading"] == "a step" for _n, _s, e in A3]),
                         "own_row_not_needed": share([e["rows_minus_unseen"]["reading"] == "no step" for _n, _s, e in A3]),
                         "median_unseen_minus_rows": float(np.median([e["rows_minus_unseen"]["gain_difference"] for _n, _s, e in A3])) if A3 else None,
                         "unseen_no_added_loss": share([e["arms"]["unseen"]["beyond_reference"]["reading"] == "no added loss" for _n, _s, e in A3]),
                         "rows_add_a_step_over_controls": share([e["rows_minus_controls"]["reading"] == "a step" for _n, _s, e in A3]),
                         "median_abs_u": {a: float(np.median([abs(e["arms"][a]["u"]) for _n, _s, e in A3])) for a in ("controls", "rows", "cells", "unseen")} if A3 else None},
            "among_reference_rejected": {"pairs": len(rej), "unreadable": share([it(n, s, "unreadable") for n, s, _e in rej]), "upstream": share([it(n, s, "upstream") for n, s, _e in rej])},
            "control_view_reads_nothing_for_every_pair": bool(all(abs(e["placed"]["control"]["u"]) >= 0.9 for _n, _s, e in P)),
            "least_abs_u_of_the_control_view": float(min(abs(e["placed"]["control"]["u"]) for _n, _s, e in P)),
            "positive_control_holds_for_some_protein": bool(any(t["positive_control"]["holds"] for t in res["targets"].values())),
            "statuses": {s: int(sum(t["status"] == s for t in res["targets"].values())) for s in sorted({t["status"] for t in res["targets"].values()})},
            "statuses_by_pair": {s: int(sum(e["status"] == s for _n, _s, e in P)) for s in sorted({e["status"] for _n, _s, e in P})},
            "pairs_by_donor": {u: int(sum(s == u for _n, s, _e in P)) for u in sorted({s for _n, s, _e in P})}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", type=Path)
    ap.add_argument("--observe", type=Path)
    ap.add_argument("--part", type=int, default=0)
    ap.add_argument("--of", type=int, default=1)
    ap.add_argument("--merge", action="store_true")
    a = ap.parse_args()
    out_file = RUN / "results" / "round.json"
    if out_file.exists():
        raise SystemExit("round.json exists: this reading was made once and is not made again")
    sel = json.loads((RUN / "results" / "selection.json").read_text())
    todo = [n for n, x in sel["proteins"].items() if x["admitted"]]
    part = lambda k: RUN / "results" / f"round_part_{k}.json"
    if not a.merge:
        if part(a.part).exists():
            raise SystemExit(f"{part(a.part).name} exists: this part was read once")
        d = np.load(a.work / "donors.npz")
        d = {k: d[k] for k in ("unit", "Z", "Y", "frozen")}
        x = unit_length(d["Z"].astype(np.float64))                # the head reads the state at unit length
        G = {"json": json.loads((a.observe / "scope.json").read_text()), "draws": np.load(a.observe / "scope_draws.npz")}
        got = {}
        for k, name in enumerate(todo):
            if k % a.of != a.part:
                continue
            got[name] = one_protein(name, x, d, sel["proteins"][name], G)
            print(name, got[name]["status"], {s: (round(e["as_trained"]["u"], 2), e["closure"]["grade"], e["rows_minus_unseen"]["reading"]) for s, e in got[name]["events"].items()}, flush=True)
            part(a.part).with_suffix(".tmp").write_text(json.dumps(got, default=float))
        part(a.part).with_suffix(".tmp").rename(part(a.part))
        return
    got = {}
    for k in range(a.of):
        got.update(json.loads(part(k).read_text()))
    assert sorted(got) == sorted(todo), "a protein with an admitted donor was not read"
    res = {"margin": cycle.MARGIN, "step_margin": cycle.STEP,
           "note": "exploratory; development donors; every pair of a donor and a protein on its own; the other donors count as units (seven: every interval is approximate)",
           "targets": {n: got[n] for n in todo}}
    res["not_targets"] = [{"readout": n, "why": "no donor passes A1 on the selection third", "largest": {k: x["table"][0][k] for k in ("event", "cells", "change_in_control_sd", "snr")}}
                          for n, x in sel["proteins"].items() if not x["admitted"]]
    T = [t["target"] for t in res["targets"].values()]
    res["open_items"] = outer.open_items(T, {"typed_right": None, "selection_on_audited_cells": None, "final_units": 133, "undeclared_cases": [], "words": WORDS})
    res["summary"] = summary(res)
    out_file.write_text(json.dumps(res, indent=1, default=float) + "\n")
    for k in range(a.of):
        part(k).unlink()
    print(json.dumps(res["summary"], indent=1))
    print("open items by reading:", {r: sum(i["reading"] == r for i in res["open_items"]) for r in sorted({i["reading"] for i in res["open_items"]})})


if __name__ == "__main__":
    main()
