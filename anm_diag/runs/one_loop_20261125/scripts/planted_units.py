#!/usr/bin/env python3
"""The account on a planted study of units in which no cell is a control cell. results/planted_units.json.

  python planted_units.py [seeds]

The system is `anm_diag.toy.planted_units` (run 32): ordinary units, and nine units that each differ from the others in
a known place. Every unit is read against the other units at the input view `view_a` and the state `z`, the other units'
rows as `others`. Two sizes (40 ordinary units of 300 cells; 16 of 150), twelve seeds each.

The rule, written before the run. For every special unit, what the account must say:

  h1 h2 h3     a factor that varies inside units and is in the state      no part open
  early enc    a factor the model drops before the state                  open: lost before the state, and nothing else
  modb post    a factor outside `view_a`; the readout itself              open: not in the input, and nothing else; the other
                                                                         units do not forecast it; the plan's lever is G
  solo         the readout itself, with a marker that tells the unit      the state holds it as the identity of the unit; the
               apart completely                                           plan says so; the other units do not forecast it
  weak         barely differs                                             no change to account for

right: exactly that. undecided: nothing false is said, and what should be open is unresolved (or, where nothing should
be open, a part is unresolved; or whether the others forecast is unresolved). wrong: anything else. The two controls of
the measure (a permuted stage carries nothing: |share| below 0.25; a made-up unit has no change) are counted apart.
"""
import json
import sys
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RUN.parents[1]))
from anm_diag import account as AC  # noqa: E402
from anm_diag import cycle, toy  # noqa: E402

SIZES = {"large": dict(n_ordinary=40, n_cell=300), "small": dict(n_ordinary=16, n_cell=150)}
EXPECT = {"h1": None, "h2": None, "h3": None, "early": "lost_before_the_state", "enc": "lost_before_the_state", "modb": "not_in_the_input", "post": "not_in_the_input"}


def score(e, a):
    ch, z = a["chain"], a["stages"]["z"]
    op, un = [p["part"] for p in ch["open"]], [p["part"] for p in ch["unresolved"]]
    if e == "weak":
        return "right" if ch["says"].startswith("no change to account for") else "wrong"
    if e == "solo":
        ident = z["carried_by"] == "identity" and a["plan"] and a["plan"][0]["reading"] == "held_as_identity"
        ae = (ch.get("across_events") or {}).get("reading")
        return "wrong" if (not ident or ae == "the other events forecast it") else ("right" if ae == "the other events do not forecast it" else "undecided")
    want = EXPECT[e]
    if want is None:
        return "wrong" if op else ("right" if not un else "undecided")
    if op == [want]:
        if want != "not_in_the_input":
            return "right"
        ae = ch["across_events"]["reading"]
        lever = a["plan"][0]["lever"]
        return "wrong" if ae == "the other events forecast it" else ("right" if ae == "the other events do not forecast it" and lever == "G: the graph" else "undecided")
    return "undecided" if (not op and want in un) else "wrong"


def main():
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    res = {"seeds": seeds, "sizes": SIZES, "by_size": {}, "rows": []}
    for size, kw in SIZES.items():
        tally = {e: {"right": 0, "undecided": 0, "wrong": 0} for e in list(EXPECT) + ["solo", "weak"]}
        controls = {"checks": 0, "hold": 0}
        for seed in range(seeds):
            d = toy.planted_units(seed, **kw)
            un = d["unit"]
            first = cycle.halves(un, seed)
            units = sorted(set(un.tolist()))
            for e in tally:
                a = AC.account({"view_a": d["view_a"], "z": d["z"]}, d["y"], un, first, e, others=[u for u in units if u != e], inputs=["view_a"], state="z", o=d["o"],
                               control=None, cluster=un, draws=200, seed=seed, controls=e in ("post", "early", "h1"))
                s = score(e, a)
                tally[e][s] += 1
                if "must_fail" in a:
                    mf = a["must_fail"]
                    controls["checks"] += 2
                    controls["hold"] += int(abs(mf["permuted_stage"]["carried"]["share"]) < 0.25) + int(mf["made_up_event"]["no_change_to_account_for"])
                res["rows"].append({"size": size, "seed": seed, "unit": e, "score": s, "says": a["chain"]["says"], "carried_by": a["stages"]["z"]["carried_by"],
                                    "parts": {k: v["share"] for k, v in a["chain"]["parts"].items()}})
            print(size, "seed", seed, {e: next(k for k, v in t.items() if v) if sum(t.values()) == 1 else dict(t) for e, t in tally.items()} if seed == 0 else "", flush=True)
        tot = {k: sum(t[k] for t in tally.values()) for k in ("right", "undecided", "wrong")}
        res["by_size"][size] = {"by_unit": tally, "total": tot, "controls_of_the_measure": controls}
        print(size, tot, controls, flush=True)
    res["total"] = {k: sum(v["total"][k] for v in res["by_size"].values()) for k in ("right", "undecided", "wrong")}
    (RUN / "results" / "planted_units.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("total", res["total"])


if __name__ == "__main__":
    main()
